#!/usr/bin/env python3
"""death-star-2026-07-16-S23 (HYP-7046): (1) THE DILATE RATE LAW referee; (2) THE
DIAGONAL-POLYNOMIAL n=5 METAGRAPH TEST (2607.14068 Lemma 2.3).

PART 1 — THEOREM (dilate rate law): for E = c*E0, the movie is the base movie traversed c
times: merged events n = c*n0, states V = V0, so
    rate = 1 - (V0-1)/(c*n0),  k = (c-1)*n0 + k0,  d = 2 for c >= 2,
    palette code = (c-1)*n0 parallel-edge doublets (+) lifted base code.
Referee: c-ladders on two bases; near-dilate (consecutive) V(c) measured for contrast.

PART 2 — for W subset of F2^6 (tilings of an n=5 iso class or invariant level set), let
d(W) = min degree with the F2 evaluation matrix (monomials deg <= d on W) having rank |W|
(delta functions realizable). Lemma 2.3 applied to the induced subgraph Q6[W] (labels =
tilings themselves: distinct, Hamming-1 adjacency automatic) gives
    avg internal wiggly degree of W <= 2*d(W).
Census: all 12 iso classes, H-levels, c3-parity sets. Compare bound vs exact.
"""
from fractions import Fraction as Fr
from itertools import combinations, product as iproduct
import sys, time

# ---------------- Part 1: movies (S22 machinery, condensed) ----------------

def movie_stats(E):
    mov = [f for f in E if f > 0]
    evs = []
    for i, f in enumerate(mov):
        for w in range(7 * f):
            evs.append((Fr(w, 7 * f), i, w))
    evs.sort()
    pos = []
    i = 0
    while i < len(evs):
        p = evs[i][0]
        while i < len(evs) and evs[i][0] == p:
            i += 1
        pos.append(p)
    n = len(pos)
    states = set()
    seq = []
    for j in range(n):
        a = pos[j]; b = pos[j + 1] if j + 1 < n else Fr(1)
        mid = (a + b) / 2
        st = tuple(int((f * mid % 1) * 7) for f in mov)
        seq.append(st); states.add(st)
    V = len(states)
    # parallel-edge count for d=2 check
    stid = {s: i for i, s in enumerate(states)}
    epairs = {}
    par = 0
    for j in range(n):
        u, v = stid[seq[j]], stid[seq[(j + 1) % n]]
        key = (min(u, v), max(u, v))
        par += epairs.get(key, 0) and 1
        epairs[key] = epairs.get(key, 0) + 1
    npar = sum(1 for kv in epairs.values() if kv >= 2)
    return n, V, n - V + 1, npar

def part1():
    print("PART 1: THE DILATE RATE LAW — referee")
    for base, name in [([2, 3, 5, 7, 11, 13], "B1=[2,3,5,7,11,13]"),
                       ([3, 4, 7, 9, 10, 11], "B2=[3,4,7,9,10,11]")]:
        n0, V0, k0, _ = movie_stats([0] + base)
        print(f"  base {name}: n0={n0} V0={V0} k0={k0}")
        for c in [2, 3, 4, 5]:
            E = [0] + [c * f for f in base]
            n, V, k, npar = movie_stats(E)
            pred_n, pred_V, pred_k = c * n0, V0, (c - 1) * n0 + k0
            ok = (n == pred_n and V == pred_V and k == pred_k)
            print(f"    c={c}: [n,V,k] = [{n},{V},{k}] vs law [{pred_n},{pred_V},{pred_k}] "
                  f"{'OK' if ok else '*** MISMATCH'}  parallel-pairs={npar} (d=2: {npar>0})")
        sys.stdout.flush()
    print("  contrast — NEAR-dilate (consecutive) V(c): law does NOT apply, measured:")
    for c in [10, 20, 30, 50]:
        n, V, k, _ = movie_stats([0] + list(range(c, c + 6)))
        print(f"    [{c}..{c+5}]: n={n} V={V} rate={k/n:.3f}  (1-(V-1)/n consistency: {1-(V-1)/n:.3f})")

# ---------------- Part 2: n=5 tilings, classes, diagonal degree ----------------

def tournament_from_tiling(bits):
    """n=5: base path 5->4->3->2->1; tiles (a,b), a>=b+2, order: klein/canonical
    enumeration (x from 5 down, per y) — any fixed order is fine for our purposes."""
    n = 5
    T = [[0] * (n + 1) for _ in range(n + 1)]
    for a in range(n, 1, -1):
        T[a][a - 1] = 1
    tiles = []
    for y in range(1, n - 1):
        for x in range(n, y + 1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    assert len(tiles) == 6
    for (a, b), bit in zip(tiles, bits):
        if bit == 0:
            T[a][b] = 1
        else:
            T[b][a] = 1
    return T, tiles

def canon(T):
    """canonical form of a 5-tournament under relabeling: min adjacency bitstring."""
    from itertools import permutations
    n = 5
    best = None
    for perm in permutations(range(1, n + 1)):
        s = 0
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                if i != j and T[perm[i - 1]][perm[j - 1]]:
                    s = s * 2 + 1
                else:
                    s = s * 2
        if best is None or s < best:
            best = s
    return best

def c3_count(T):
    n = 5
    c = 0
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            for k in range(1, n + 1):
                if i < j and T[i][j] and T[j][k] and T[k][i] and k != i and k != j:
                    pass
    # simpler: c3 = C(n,3) - sum C(outdeg,2)... use score formula
    outs = [sum(T[i][j] for j in range(1, n + 1)) for i in range(1, n + 1)]
    from math import comb
    return comb(n, 3) - sum(comb(o, 2) for o in outs)

def diag_degree(W, m=6):
    """min d such that deg<=d monomial evaluations on W have F2-rank |W|."""
    rows = []
    for d in range(0, m + 1):
        # build evaluation matrix incrementally: monomials of degree exactly d
        for sup in combinations(range(m), d):
            row = 0
            for i, w in enumerate(W):
                val = 1
                for b in sup:
                    val &= (w >> b) & 1
                row |= val << i
            rows.append(row)
        # rank over F2 of rows (as bit-ints of length |W|)
        rk = 0
        basis = []
        for r in rows:
            x = r
            for b in basis:
                x = min(x, x ^ b)
            if x:
                basis.append(x)
                basis.sort(reverse=True)
        # proper elimination
        basis = []
        for r in rows:
            x = r
            for b in basis:
                hb = b.bit_length() - 1
                if (x >> hb) & 1:
                    x ^= b
            if x:
                basis.append(x)
                basis.sort(key=lambda z: -z.bit_length())
        if len(basis) == len(W):
            return d
    return m

def part2():
    print("\nPART 2: THE DIAGONAL-POLYNOMIAL n=5 METAGRAPH TEST (Lemma 2.3)")
    classes = {}
    info = {}
    for bits_int in range(64):
        bits = [(bits_int >> i) & 1 for i in range(6)]
        T, tiles = tournament_from_tiling(bits)
        key = canon(T)
        classes.setdefault(key, []).append(bits_int)
        if key not in info:
            info[key] = c3_count(T)
    print(f"  #classes = {len(classes)} (A000568(5) = 12 expected)")
    # internal wiggly degrees per class: edges within W in Q6
    print(f"  {'class(c3)':>10} {'|W|':>4} {'int-edges':>9} {'avg-deg':>8} {'d(W)':>5} {'2d':>3} {'tight?':>6}")
    rows = []
    for key, W in sorted(classes.items(), key=lambda kv: (info[kv[0]], len(kv[1]))):
        Wset = set(W)
        ie = sum(1 for w in W for b in range(6) if (w ^ (1 << b)) in Wset) // 2
        avg = 2 * ie / len(W)
        dW = diag_degree(W)
        rows.append((info[key], len(W), ie, avg, dW))
        print(f"  {info[key]:>10} {len(W):>4} {ie:>9} {avg:>8.3f} {dW:>5} {2*dW:>3} "
              f"{'TIGHT' if abs(avg - 2*dW) < 1e-9 else ''}")
    # H-levels: c3 determines H at n=5? (H = ... use c3 as the invariant); c3-parity sets
    for par in [0, 1]:
        W = [w for key, ws in classes.items() for w in ws if info[key] % 2 == par]
        Wset = set(W)
        ie = sum(1 for w in W for b in range(6) if (w ^ (1 << b)) in Wset) // 2
        avg = 2 * ie / len(W)
        dW = diag_degree(W)
        print(f"  c3-parity {par}: |W|={len(W)} int-edges={ie} avg-deg={avg:.3f} "
              f"d(W)={dW} bound 2d={2*dW}")
    sys.stdout.flush()

if __name__ == "__main__":
    t0 = time.time()
    part1()
    part2()
    print(f"[total {time.time()-t0:.1f}s]")
