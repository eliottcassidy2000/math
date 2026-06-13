#!/usr/bin/env python3
"""h21_finite_check_monad_s4.py — close the FINITE remaining window of the H=21 proof.

BACKGROUND (HYP-2193 / S617 reduction, corrected from HYP-2187):
  H(T)=21 => (H multiplicative over strong components; 21=3*7; 7 is NOT a strong
  H-value, THM-029) => some STRONG component S has H(S)=21.
  H(S)=21 = I(Omega,2) = 1 + 2*alpha_1 + 4*alpha_2 + ... => alpha_1 <= 10.
  3-cycles are odd cycles => c_3(S) <= alpha_1 <= 10.
  Moon (1968): a strong tournament on m vertices has c_3 >= m-2.
  => m-2 <= 10 => m <= 12.
  THM-079 Part G proved H=21 impossible EXHAUSTIVELY for m<=8 (2^28 at n=8).
  REMAINING: STRONG tournaments on m in {9,10,11,12} with c_3 <= 10. FINITE.

THIS SCRIPT: isomorph-free exhaustive enumeration of ALL tournaments with
c_3 <= CAP (default 10) up to m vertices, via canonical augmentation (orderly
generation). For each STRONG one we compute H by Held-Karp and check H != 21.
Equivalently we report min(alpha_1) / min(H): alpha_1 >= 11 (i.e. H >= 23) on
the whole window finishes H=21 (since H=21 needs alpha_1<=10).

Pruning is valid because c_3 is MONOTONE non-decreasing under vertex addition
(adding a vertex only creates new triangles), so every induced subtournament of
a c_3<=10 tournament also has c_3<=10. Hence building vertex-by-vertex and
pruning partial c_3 > CAP loses nothing.

No C compiler / numba / nauty on this node: pure Python. Canonicalization is
cheap here because low-c_3 tournaments are near-transitive (scores nearly all
distinct => refinement gives singleton color classes => unique labeling).

Session: monad-compute-2026-06-04-S4.
"""
import sys, time
sys.stdout.reconfigure(line_buffering=True)

CAP = 10

# ---------------------------------------------------------------------------
# Tournament representation: `beats` = tuple of out-neighbor bitmasks, length m.
# beats[a] has bit b set  <=>  a -> b.  (a beats b)
# ---------------------------------------------------------------------------

def c3_count(beats, m):
    cnt = 0
    for i in range(m):
        bi = beats[i]
        for j in range(i+1, m):
            for k in range(j+1, m):
                # triple {i,j,k}: cyclic iff each has out-degree 1 in the triple
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
    """Number of Hamiltonian paths = Redei H(T), via Held-Karp DP."""
    size = 1 << m
    dp = [0] * (size * m)  # dp[mask*m + v]
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

def alpha1_count(beats, m):
    """Total number of directed odd cycles (alpha_1). Brute DFS over simple
    directed cycles, counting only odd lengths. Used as a cheaper early reject
    (alpha_1 >= 11 => H >= 23). Counts each directed cycle once."""
    cnt = 0
    # enumerate simple cycles with smallest vertex as the start to avoid dup
    for start in range(m):
        # DFS paths start->...->start, vertices > start only as interior
        stack = [(start, 1 << start, 1)]  # (cur, visited, length)
        while stack:
            cur, vis, ln = stack.pop()
            bc = beats[cur]
            # try to close
            if ln >= 3 and (ln & 1) and ((bc >> start) & 1):
                cnt += 1
            nb = bc & ~vis
            # only allow interior vertices > start
            nb &= ~((1 << (start + 1)) - 1)
            while nb:
                w = (nb & -nb).bit_length() - 1
                stack.append((w, vis | (1 << w), ln + 1))
                nb &= nb - 1
    return cnt

# ---------------------------------------------------------------------------
# Canonical form via iterated color refinement + within-class minimization.
# ---------------------------------------------------------------------------

def _refine_colors(beats, m):
    # initial color = out-degree
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
        # relabel signatures to small ints, preserving order
        order = {s: i for i, s in enumerate(sorted(set(sig)))}
        newcol = [order[s] for s in sig]
        if newcol == col:
            break
        col = newcol
    return col

def _encode(perm, beats, m):
    # perm: new_position -> old_vertex.  Encode upper triangle.
    pos = [0] * m
    for newp, old in enumerate(perm):
        pos[old] = newp
    code = 0
    bit = 0
    for a in range(m):
        for b in range(a + 1, m):
            # does new-vertex a beat new-vertex b?
            oa, ob = perm[a], perm[b]
            if (beats[oa] >> ob) & 1:
                code |= 1 << bit
            bit += 1
    return code

def canon(beats, m):
    col = _refine_colors(beats, m)
    # group vertices by color, order classes by color value
    from itertools import permutations
    classes = {}
    for v in range(m):
        classes.setdefault(col[v], []).append(v)
    ordered_colors = sorted(classes)
    class_lists = [classes[c] for c in ordered_colors]
    # candidate perms: product of permutations within each class, concatenated
    # (fast path when all singletons)
    if all(len(cl) == 1 for cl in class_lists):
        perm = [cl[0] for cl in class_lists]
        return _encode(perm, beats, m)
    best = None
    # iterate product of within-class perms
    def rec(idx, acc):
        nonlocal best
        if idx == len(class_lists):
            code = _encode(acc, beats, m)
            if best is None or code < best:
                best = code
            return
        for p in permutations(class_lists[idx]):
            rec(idx + 1, acc + list(p))
    rec(0, [])
    return best

def beats_from_canon(code, m):
    beats = [0] * m
    bit = 0
    for a in range(m):
        for b in range(a + 1, m):
            if (code >> bit) & 1:
                beats[a] |= 1 << b
            else:
                beats[b] |= 1 << a
            bit += 1
    return tuple(beats)

# ---------------------------------------------------------------------------
# Canonical-augmentation generation: level by level.
# R[k] = set of canonical codes of tournaments on k vertices with c3<=CAP.
# ---------------------------------------------------------------------------

def generate(maxm, cap=CAP, verbose=True):
    # start: m=1, single vertex, code 0
    R = {1: {0}}
    for k in range(1, maxm):
        cur = R[k]
        nxt = set()
        t0 = time.time()
        for code in cur:
            beats = list(beats_from_canon(code, k)) + [0]
            base = tuple(beats[:k])
            # extend: new vertex = index k; choose subset S of {0..k-1} it beats
            for mask in range(1 << k):
                nb = list(base) + [mask]
                # set incoming arcs to new vertex
                for i in range(k):
                    if not ((mask >> i) & 1):
                        nb[i] |= 1 << k
                nbt = tuple(nb)
                if c3_count(nbt, k + 1) <= cap:
                    nxt.add(canon(nbt, k + 1))
                # reset (rebuild each loop anyway)
        R[k + 1] = nxt
        if verbose:
            print(f"    m={k+1}: iso-classes c3<={cap}: {len(nxt):>12,}   "
                  f"(gen {time.time()-t0:.1f}s)")
    return R

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxm", type=int, default=12)
    ap.add_argument("--cap", type=int, default=CAP)
    ap.add_argument("--no-validate", dest="validate", action="store_false",
                    help="skip the A000568 self-validation")
    ap.set_defaults(validate=True)
    args = ap.parse_args()

    print("=" * 70)
    print("  H=21 FINITE WINDOW: exhaustive strong c_3<=%d tournaments" % args.cap)
    print("=" * 70)

    if args.validate:
        print("\n  [validation] iso-class counts with NO c3 cap vs A000568:")
        A000568 = {1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880}
        Rv = generate(min(args.maxm, 7), cap=10**9, verbose=False)
        ok = True
        for k in sorted(Rv):
            exp = A000568.get(k, "?")
            mark = "OK" if exp == len(Rv[k]) else "MISMATCH"
            if exp != len(Rv[k]) and exp != "?":
                ok = False
            print(f"      n={k}: generated {len(Rv[k]):>6}  expected {exp}  [{mark}]")
        print("  validation:", "PASS" if ok else "FAIL")

    print(f"\n  [main] generating c3<={args.cap} up to m={args.maxm}:")
    R = generate(args.maxm, cap=args.cap)

    print(f"\n  [scan] strong + H over c3<={args.cap} classes:")
    for k in sorted(R):
        if k < 3:
            continue
        strong = 0; h21 = 0; minH = None; minA1 = None
        Hdist = {}
        for code in R[k]:
            beats = beats_from_canon(code, k)
            if is_strong(beats, k):
                strong += 1
                H = H_count(beats, k)
                Hdist[H] = Hdist.get(H, 0) + 1
                if minH is None or H < minH:
                    minH = H
                if H == 21:
                    h21 += 1
        lo = sorted(Hdist)[:6] if Hdist else []
        print(f"      m={k}: classes={len(R[k]):>10,}  strong={strong:>9,}  "
              f"H=21:{h21}  minH={minH}  lowestH={lo}")
    print("\n  DONE.")
