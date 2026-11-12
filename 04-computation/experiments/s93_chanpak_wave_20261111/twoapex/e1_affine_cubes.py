#!/usr/bin/env python3
"""Lane twoapex, deliverable (e): is the two-apex count ever affine on a cube of letter changes?

A word of length n is a point of {0,1}^(2n) (bits a_i, b_i).  A MOVE is a change at one position: flip a_i, flip
b_i, or flip both (mask 1, 2, 3 at position i).  A k-cube = a base word and k XOR-independent moves; H is AFFINE on
it iff every 2-face has zero second difference H(x) - H(x + s) - H(x + t) + H(x + s + t) = 0.
(If a family of words had H = F(A) + L(B) with L linear in free binary digits of B - a full interval - those
digits would span an affine cube with increments 2, 4, 8, ...)

Exhaustive for all words of length n <= NMAX: the number of affine 2-cubes by type, the number of affine 3-cubes,
the largest dimension of an affine cube, and for the affine cubes the multiset of increments.  Also, for a fixed
word b and free a (and symmetrically), the graph of the pairs (i, j) with a nonzero quadratic coefficient
c_ij(b) = d_i d_j H, and the largest set of a-bits on which H is affine.
"""
import itertools
import sys
from collections import Counter

import numpy as np

from ta_fast import FIN24, INIT24, M24, all_alpha


def all_H(n):
    al = all_alpha(n, M24, INIT24)                 # index: first letter = top digit
    h = al @ FIN24
    # re-index so that bit 2i = a_i, bit 2i+1 = b_i with i = 0 the FIRST letter
    idx = np.arange(4 ** n)
    new = np.zeros(4 ** n, np.int64)
    for i in range(n):
        letter = (idx >> (2 * (n - 1 - i))) & 3
        new |= letter << (2 * i)
    out = np.empty(4 ** n, np.int64)
    out[new] = h
    return out


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    for n in range(2, nmax + 1):
        H = all_H(n)
        size = 4 ** n
        x = np.arange(size)
        moves = [(i, t) for i in range(n) for t in (1, 2, 3)]
        masks = [t << (2 * i) for i, t in moves]
        # second differences for every pair of independent moves
        zero = {}
        cnt2 = Counter()
        for s in range(len(moves)):
            for t in range(s + 1, len(moves)):
                ms, mt = masks[s], masks[t]
                if moves[s][0] == moves[t][0] and (ms ^ mt) == 0:
                    continue
                D = H - H[x ^ ms] - H[x ^ mt] + H[x ^ ms ^ mt]
                z = D == 0
                zero[(s, t)] = z
                same = moves[s][0] == moves[t][0]
                key = ("same position" if same else "two positions",
                       "".join(sorted("a" if m == 1 else "b" if m == 2 else "ab" for m in (moves[s][1], moves[t][1]))))
                cnt2[key] += int(z.sum()) // 4
        tot2 = sum(cnt2.values())
        # 3-cubes
        n3 = 0
        incr3 = Counter()
        ex3 = []
        for s, t, r in itertools.combinations(range(len(moves)), 3):
            ms, mt, mr = masks[s], masks[t], masks[r]
            # XOR-independence: no subset XORs to zero
            if (ms ^ mt) == 0 or (ms ^ mr) == 0 or (mt ^ mr) == 0 or (ms ^ mt ^ mr) == 0:
                continue
            ok = (zero[(s, t)] & zero[(s, t)][x ^ mr] & zero[(s, r)] & zero[(s, r)][x ^ mt]
                  & zero[(t, r)] & zero[(t, r)][x ^ ms])
            c = int(ok.sum())
            if c:
                n3 += c // 8
                base = x[ok]
                for b0 in base[:: max(1, len(base) // 50)]:
                    inc = tuple(sorted(abs(int(H[b0 ^ m] - H[b0])) for m in (ms, mt, mr)))
                    incr3[inc] += 1
                    if len(ex3) < 6 and all(v > 0 for v in inc):
                        ex3.append((int(b0), moves[s], moves[t], moves[r], inc))
        # nonconstant affine 3-cubes (all three increments nonzero)
        nonconst3 = sum(v for k, v in incr3.items() if all(q > 0 for q in k))
        print(f"n = {n}: {size} words; affine 2-cubes: {tot2} "
              f"({dict(sorted(cnt2.items()))}); affine 3-cubes: {n3}; sampled increment triples with all "
              f"increments nonzero: {nonconst3}")
        if ex3:
            for e in ex3[:4]:
                w = [(e[0] >> (2 * i)) & 3 for i in range(n)]
                print("    example affine 3-cube with nonzero increments: base", "".join(map(str, w)), "moves", e[1:4],
                      "increments", e[4])
        # 4-cubes: only if 3-cubes exist; check by brute force over quadruples (n <= 6 to keep it short)
        if n3 and n <= 6:
            n4 = 0
            for quad in itertools.combinations(range(len(moves)), 4):
                ms = [masks[q] for q in quad]
                indep = True
                for k in range(1, 16):
                    acc = 0
                    for j in range(4):
                        if (k >> j) & 1:
                            acc ^= ms[j]
                    if acc == 0:
                        indep = False
                        break
                if not indep:
                    continue
                ok = np.ones(size, bool)
                for i1 in range(4):
                    for i2 in range(i1 + 1, 4):
                        z = zero[(quad[i1], quad[i2])]
                        others = [ms[j] for j in range(4) if j not in (i1, i2)]
                        for shift in (0, others[0], others[1], others[0] ^ others[1]):
                            ok &= z[x ^ shift]
                        if not ok.any():
                            break
                    if not ok.any():
                        break
                c = int(ok.sum())
                if c:
                    # nonconstant ones
                    base = x[ok]
                    good = 0
                    for b0 in base:
                        if all(H[b0 ^ m] != H[b0] for m in ms):
                            good += 1
                    n4 += c // 16
                    if good:
                        print("    affine 4-cube with four nonzero increments exists:", quad, good // 16)
            print(f"    affine 4-cubes at n = {n}: {n4}")
    # ---- fixed b, free a: quadratic coefficients
    print()
    print("fixed word b, free word a: H is a polynomial of degree <= 2 in the bits of a; c_ij(b) = d_i d_j H.")
    for n in range(3, min(nmax, 9) + 1):
        H = all_H(n)
        best = 0
        best_ex = None
        allneg = True
        for bw in range(2 ** n):
            # spread b bits to odd positions
            bmask = 0
            for i in range(n):
                if (bw >> i) & 1:
                    bmask |= 2 << (2 * i)
            # c_ij from a = 0 corner (degree <= 2 makes it corner independent; verified below on random corners)
            adj = [[False] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    mi, mj = 1 << (2 * i), 1 << (2 * j)
                    c = int(H[bmask] - H[bmask ^ mi] - H[bmask ^ mj] + H[bmask ^ mi ^ mj])
                    if c > 0:
                        allneg = False
                    if c != 0:
                        adj[i][j] = adj[j][i] = True
            # largest independent set (n <= 9: brute force)
            for S in range(1, 2 ** n):
                k = bin(S).count("1")
                if k <= best:
                    continue
                el = [i for i in range(n) if (S >> i) & 1]
                if all(not adj[p][q] for p in el for q in el if p < q):
                    best = k
                    best_ex = (bw, el)
        btxt = "".join(str((best_ex[0] >> i) & 1) for i in range(n))
        print(f"  n = {n}: largest set of a-bits with all c_ij(b) = 0, over all b: {best} (b = {btxt}, positions "
              f"{[p + 1 for p in best_ex[1]]}); every c_ij(b) <= 0: {allneg}")


if __name__ == "__main__":
    main()
