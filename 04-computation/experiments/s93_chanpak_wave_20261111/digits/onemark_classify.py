#!/usr/bin/env python3
"""Lane "digits": exact affine classes of the one-mark block system (state (h, e, B, g), readout h).

For every block op N of order <= KMAX (onemark_matrices_KMAX.json) and every integer eigenvalue a >= 1 of N:
the space C(N; a) of states whose orbit obeys h -> a h + b exactly (b = lambda s, lambda = L N - a L).
For every PAIR of ops and every pair of integer eigenvalues: C(N1, N2; a1, a2), and the functional
    kappa = (a2 - 1) lambda_1 - (a1 - 1) lambda_2,
which vanishes at s iff the two affine maps x -> a1 x + b1, x -> a2 x + b2 of the class of s commute
(same fixed point).  A digit system needs two NON-commuting maps, i.e. a state s in C with kappa s != 0.
Categories of a pair: (i) C = 0; (ii) C != 0 and kappa = 0 on all of C; (iii) kappa != 0 somewhere on C.
Actual states: every marked tournament (T, z) with T of order <= 6.

Usage: onemark_classify.py KMAX
"""
import json
import sys
import time
from collections import Counter, defaultdict

from affine import class_space, dot, int_eigs
from onemark import state
from tour import iso_classes


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    t0 = time.time()
    ops = json.load(open(f"onemark_matrices_{kmax}.json"))
    mats = [o["N"] for o in ops]
    L = [1, 0, 0, 0]
    classes = iso_classes(6)
    states = {}
    for n in range(1, 7):
        for key, adj in classes[n].items():
            for z in range(n):
                s = state(n, adj, z)
                states.setdefault(s, (n, [int(x) for x in adj], z))
    S = sorted(states)
    print(f"{len(mats)} block matrices (blocks of order <= {kmax}); {len(S)} distinct actual states (h, e, B, g) from marked tournaments of order <= 6")
    eigs = [[a for a in int_eigs(G) if a >= 1] for G in mats]
    cnt = Counter()
    for e in eigs:
        for a in e:
            cnt[a] += 1
    print("integer eigenvalues >= 1 of the block matrices (eigenvalue: number of matrices):", dict(sorted(cnt.items())))

    def desc(i):
        o = ops[i]
        return f"#{i}(k={o['k']},X={o['adjX']},W={o['W']:0{o['k']}b},m={o['m']})"

    # ---------------- single operations
    print()
    print("SINGLE OPERATIONS with an exact class h -> a h + b on which h is not fixed (a = integer eigenvalue):")
    print("  [census states in C with a h + b != h; offsets b seen with multiplicities]")
    summ = defaultdict(list)
    for i, G in enumerate(mats):
        for a in eigs[i]:
            V, lam = class_space(L, [G], [a])
            dimC = 4 - V.dim()
            if dimC == 0:
                continue
            inC = [s for s in S if all(dot(r, s) == 0 for r in V.rows)]
            mov = [s for s in inC if a * s[0] + dot(lam[0], s) != s[0]]
            bs = Counter(dot(lam[0], s) for s in mov)
            mult = all(x == 0 for x in lam[0])
            kind = "multiplicative on every state" if mult else ("b = 0 on the class" if all(b == 0 for b in bs) else "AFFINE, b != 0")
            if mov:
                summ[(a, kind)].append(i)
                bl = sorted(bs.items())
                bshow = bl if len(bl) <= 10 else bl[:5] + ["..."] + bl[-4:]
                print(f"  {desc(i)}: a = {a}, dim C = {dimC}, moving states {len(mov)}, b: {bshow}  [{kind}]")
    print("summary (multiplier, kind): number of ops")
    for k2, v in sorted(summ.items()):
        print(f"   {k2}: {len(v)}")
    # ---------------- pairs
    print()
    cat = Counter()
    free_pairs = []
    free_with_states = []
    for i in range(len(mats)):
        for j in range(i + 1, len(mats)):
            for a1 in eigs[i]:
                for a2 in eigs[j]:
                    V, lam = class_space(L, [mats[i], mats[j]], [a1, a2])
                    dimC = 4 - V.dim()
                    if dimC == 0:
                        cat["(i) C = 0"] += 1
                        continue
                    kappa = [(a2 - 1) * x - (a1 - 1) * y for x, y in zip(lam[0], lam[1])]
                    if V.contains(kappa):
                        cat["(ii) C != 0, the two maps commute on all of C"] += 1
                        continue
                    cat["(iii) C != 0, kappa != 0 somewhere on C"] += 1
                    inC = [s for s in S if all(dot(r, s) == 0 for r in V.rows)]
                    good = [(s, dot(lam[0], s), dot(lam[1], s)) for s in inC if dot(kappa, s) != 0]
                    free_pairs.append((i, j, a1, a2, dimC, len(inC), len(good), V.annihilator()))
                    if good:
                        free_with_states.append((i, j, a1, a2, dimC, good))
    print("PAIRS OF OPERATIONS (op1 < op2, integer eigenvalues a1, a2 >= 1):")
    for k2, v in sorted(cat.items()):
        print(f"   {k2}: {v}")
    print(f"category (iii) pairs: {len(free_pairs)}; of these with a census state s in C and kappa s != 0: {len(free_with_states)}")
    for (i, j, a1, a2, dimC, ninC, ngood, basis) in free_pairs[:60]:
        print(f"   {desc(i)} & {desc(j)}: a = ({a1}, {a2}), dim C = {dimC}, census states in C: {ninC}, with kappa != 0: {ngood}, C basis {basis}")
    for (i, j, a1, a2, dimC, good) in free_with_states[:40]:
        print(f"   NON-COMMUTING: {desc(i)} & {desc(j)}: a = ({a1}, {a2}); states (h,e,B,g), b1, b2: {good[:6]}")
    json.dump({"categories": dict(cat),
               "free_pairs": [[i, j, a1, a2, dimC, ninC, ngood, basis] for (i, j, a1, a2, dimC, ninC, ngood, basis) in free_pairs],
               "free_with_states": [[i, j, a1, a2, dimC, [[list(s), b1, b2] for (s, b1, b2) in good[:50]]] for (i, j, a1, a2, dimC, good) in free_with_states]},
              open(f"onemark_pairs_{kmax}.json", "w"))
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
