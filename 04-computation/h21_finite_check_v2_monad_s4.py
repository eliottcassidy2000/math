#!/usr/bin/env python3
"""h21_finite_check_v2_monad_s4.py — optimized exhaustive closure of the H=21
finite window (HYP-2193 reduction). See h21_finite_check_monad_s4.py for the
full reduction; this is the faster engine.

Two improvements over v1:
  (1) DFS-pruned vertex extension. When we add a new vertex by deciding which of
      the existing vertices it beats, we decide those bits IN ORDER (0,1,...,k-1)
      via DFS, accumulating the new 3-cycles {i,j,new} as soon as both endpoints
      i<j are decided, and PRUNE the whole subtree the instant partial c3 > CAP.
      Most of the 2^k orientations blow past CAP early, so this is a big win.
  (2) Incremental per-level scan + flush: right after generating level k we scan
      its strong classes for H, print, and flush — so m=11 / m=12 results land
      (and can be checkpointed) as soon as each level is done, instead of only at
      the very end.

Validity of the c3<=CAP pruning: c3 is monotone non-decreasing under vertex
addition, so every induced subtournament of a c3<=CAP tournament is c3<=CAP;
building vertex-by-vertex and pruning loses nothing. Completeness of the
isomorph-free generation is validated against A000568 (no-cap) before the run.

Session: monad-compute-2026-06-04-S4.
"""
import sys, time
sys.stdout.reconfigure(line_buffering=True)
from itertools import permutations

CAP = 10

# ---- core invariants -------------------------------------------------------

def c3_count(beats, m):
    cnt = 0
    for i in range(m):
        bi = beats[i]
        for j in range(i+1, m):
            for k in range(j+1, m):
                di = ((bi >> j) & 1) + ((bi >> k) & 1)
                dj = ((beats[j] >> i) & 1) + ((beats[j] >> k) & 1)
                dk = ((beats[k] >> i) & 1) + ((beats[k] >> j) & 1)
                if di == 1 and dj == 1 and dk == 1:
                    cnt += 1
    return cnt

def is_strong(beats, m):
    full = (1 << m) - 1
    def reach(fwd):
        seen = 1; fr = [0]
        while fr:
            x = fr.pop()
            if fwd:
                nb = beats[x] & ~seen
            else:
                nb = 0
                for y in range(m):
                    if (not (seen >> y & 1)) and ((beats[y] >> x) & 1):
                        nb |= 1 << y
            while nb:
                y = (nb & -nb).bit_length() - 1
                seen |= 1 << y; fr.append(y); nb &= nb - 1
        return seen == full
    return reach(True) and reach(False)

def H_count(beats, m):
    size = 1 << m
    dp = [0] * (size * m)
    for v in range(m):
        dp[(1 << v) * m + v] = 1
    for mask in range(size):
        base = mask * m
        for v in range(m):
            c = dp[base + v]
            if c:
                av = beats[v] & ~mask
                while av:
                    w = (av & -av).bit_length() - 1
                    dp[((mask | (1 << w)) * m) + w] += c
                    av &= av - 1
    last = (size - 1) * m
    return sum(dp[last:last + m])

# ---- canonical form (colour refinement + within-class min) -----------------

def _refine_colors(beats, m):
    col = [bin(beats[v]).count("1") for v in range(m)]
    inn = [0] * m
    for v in range(m):
        bv = beats[v]
        while bv:
            w = (bv & -bv).bit_length() - 1
            inn[w] |= 1 << v
            bv &= bv - 1
    for _ in range(m):
        sig = []
        for v in range(m):
            outc = sorted(col[w] for w in range(m) if (beats[v] >> w) & 1)
            inc = sorted(col[w] for w in range(m) if (inn[v] >> w) & 1)
            sig.append((col[v], tuple(outc), tuple(inc)))
        order = {s: i for i, s in enumerate(sorted(set(sig)))}
        newcol = [order[s] for s in sig]
        if newcol == col:
            break
        col = newcol
    return col

def _encode(perm, beats, m):
    code = 0; bit = 0
    for a in range(m):
        oa = perm[a]
        ba = beats[oa]
        for b in range(a + 1, m):
            if (ba >> perm[b]) & 1:
                code |= 1 << bit
            bit += 1
    return code

def canon(beats, m):
    col = _refine_colors(beats, m)
    classes = {}
    for v in range(m):
        classes.setdefault(col[v], []).append(v)
    class_lists = [classes[c] for c in sorted(classes)]
    if all(len(cl) == 1 for cl in class_lists):
        return _encode([cl[0] for cl in class_lists], beats, m)
    best = [None]
    def rec(idx, acc):
        if idx == len(class_lists):
            code = _encode(acc, beats, m)
            if best[0] is None or code < best[0]:
                best[0] = code
            return
        for p in permutations(class_lists[idx]):
            rec(idx + 1, acc + list(p))
    rec(0, [])
    return best[0]

def beats_from_canon(code, m):
    beats = [0] * m; bit = 0
    for a in range(m):
        for b in range(a + 1, m):
            if (code >> bit) & 1:
                beats[a] |= 1 << b
            else:
                beats[b] |= 1 << a
            bit += 1
    return tuple(beats)

# ---- DFS-pruned extension --------------------------------------------------

def extend(parent_beats, k, cap, out_set):
    """parent_beats: tuple length k (beats among 0..k-1). Add new vertex index k.
    Decide bits 0..k-1 (does new beat i?) in order via DFS, prune partial c3>cap.
    For each completed orientation with c3<=cap, add canon(child) to out_set."""
    pb = parent_beats
    # precompute existing orientation between i,j (i<j<k): does i beat j?
    base_c3 = c3_count(pb, k)
    new = k  # new vertex index
    # store the new vertex's beat-mask over existing verts as we build
    # decided[i] in {0,1}: 1 => new beats i ; 0 => i beats new
    decided = [0] * k
    def dfs(i, c3, beat_mask):
        if c3 > cap:
            return
        if i == k:
            # build child beats
            child = list(pb) + [beat_mask]
            for v in range(k):
                if not ((beat_mask >> v) & 1):
                    child[v] |= 1 << new
            out_set.add(canon(tuple(child), k + 1))
            return
        # decide bit i: case new beats i (1), case i beats new (0)
        for bit in (1, 0):
            add = 0
            # triangles {j, i, new} for j < i, with both j and i now decided
            # triple vertices: j, i, new. orientation:
            #   j-i: existing
            #   new-i: bit (new beats i iff bit)
            #   new-j: decided[j]
            for j in range(i):
                ji = (pb[j] >> i) & 1  # j beats i
                ij = 1 - ji
                ni = bit              # new beats i
                nj = decided[j]       # new beats j
                # out-degrees within the triple {j, i, new}
                dj = ji + (1 - nj)    # j beats i? + j beats new?
                di = ij + (1 - ni)    # i beats j? + i beats new?
                dn = ni + nj          # new beats i? + new beats j?
                if dj == 1 and di == 1 and dn == 1:
                    add += 1
            nc = c3 + add
            if nc <= cap:
                decided[i] = bit
                dfs(i + 1, nc, beat_mask | (bit << i))
    dfs(0, base_c3, 0)

# ---- driver ----------------------------------------------------------------

def generate_and_scan(maxm, cap=CAP):
    R = {1: {0}}
    results = {}
    for k in range(1, maxm):
        cur = R[k]
        nxt = set()
        t0 = time.time()
        for code in cur:
            pb = beats_from_canon(code, k)
            extend(pb, k, cap, nxt)
        R[k + 1] = nxt
        gen_t = time.time() - t0
        m = k + 1
        # incremental scan
        if m >= 3:
            t1 = time.time()
            strong = 0; h21 = 0; minH = None; Hdist = {}
            for ccode in nxt:
                beats = beats_from_canon(ccode, m)
                if is_strong(beats, m):
                    strong += 1
                    H = H_count(beats, m)
                    Hdist[H] = Hdist.get(H, 0) + 1
                    if minH is None or H < minH:
                        minH = H
                    if H == 21:
                        h21 += 1
            lo = sorted(Hdist)[:6]
            results[m] = (len(nxt), strong, h21, minH, lo)
            print(f"  m={m:>2}: classes={len(nxt):>10,}  strong={strong:>9,}  "
                  f"H=21:{h21}  minH={minH}  lowestH={lo}   "
                  f"(gen {gen_t:.1f}s, scan {time.time()-t1:.1f}s)")
        else:
            print(f"  m={m:>2}: classes={len(nxt):>10,}   (gen {gen_t:.1f}s)")
        # free the previous level
        del R[k]
    return results

def validate():
    print("\n  [validation] no-cap iso-class counts vs A000568:")
    A = {1:1,2:1,3:2,4:4,5:12,6:56,7:456}
    R = {1: {0}}
    ok = True
    for k in range(1, 7):
        nxt = set()
        for code in R[k]:
            extend(beats_from_canon(code, k), k, 10**9, nxt)
        R[k+1] = nxt; del R[k]
        exp = A.get(k+1)
        mark = "OK" if exp == len(nxt) else "MISMATCH"
        if exp != len(nxt): ok = False
        print(f"      n={k+1}: {len(nxt):>5}  expected {exp}  [{mark}]")
    print("  validation:", "PASS" if ok else "FAIL")
    return ok

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxm", type=int, default=12)
    ap.add_argument("--cap", type=int, default=CAP)
    ap.add_argument("--no-validate", dest="validate", action="store_false")
    ap.set_defaults(validate=True)
    args = ap.parse_args()
    print("=" * 74)
    print(f"  H=21 FINITE WINDOW v2: exhaustive strong c_3<={args.cap}, m<= {args.maxm}")
    print("=" * 74)
    if args.validate:
        if not validate():
            print("ABORT: validation failed"); sys.exit(1)
    print(f"\n  [main] generate + incremental scan up to m={args.maxm}:")
    res = generate_and_scan(args.maxm, args.cap)
    print("\n  SUMMARY (strong, c_3<=%d):" % args.cap)
    any21 = False
    for m in sorted(res):
        cl, st, h21, mh, lo = res[m]
        if h21: any21 = True
        print(f"    m={m:>2}: strong={st:>9,}  H=21:{h21}  minH={mh}")
    print("\n  H=21 FOUND ANYWHERE:" , any21)
    print("  DONE.")
