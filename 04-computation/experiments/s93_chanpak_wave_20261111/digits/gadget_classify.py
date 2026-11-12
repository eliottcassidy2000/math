#!/usr/bin/env python3
"""Lane "digits", part (2): exact affine classes of the gadget system of THM-5130 (A).

State s = (r, M, h, E) of a pair (B, O); gadget (X, P) acts by G(X, P); readout h = H(B) (third coordinate), so
h+ = H_X r + Pi_X M is the path count of the closed tournament B => X => v.
For every gadget matrix with an integer eigenvalue a >= 1: the class space C(G; a).  For every pair of gadget
matrices and integer eigenvalues: C(G1, G2; a1, a2) and whether the two affine maps h -> a_i h + b_i commute on
all of C (kappa = (a2 - 1) lambda_1 - (a1 - 1) lambda_2 vanishes on C).
Census of actual states: every pair (B, O), B of order <= 6 (one representative per isomorphism class, every
nonempty O), and the empty base (1, 0, 1, 0).

Usage: gadget_classify.py NMAX [WORDS]
  WORDS = 1: also products G2 G1 of two gadgets of order <= WORDMAX (default 4) as single operations.
"""
import json
import sys
import time
from collections import Counter, defaultdict

from affine import class_space, dot, int_eigs, matmul
from gadgets import pair_state
from tour import iso_classes


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    t0 = time.time()
    gad = json.load(open(f"gadget_matrices_{nmax}.json"))
    mats = [g["G"] for g in gad]
    L = [0, 0, 1, 0]
    classes = iso_classes(6)
    states = {(1, 0, 1, 0): (0, [], 0)}
    for n in range(1, 7):
        for key, adj in classes[n].items():
            for O in range(1, 1 << n):
                s = pair_state(n, adj, O)
                states.setdefault(s, (n, [int(x) for x in adj], O))
    S = sorted(states)
    print(f"{len(mats)} gadget matrices (X of order <= {nmax}); {len(S)} distinct actual states (r, M, h, E) of pairs of order <= 6")
    eigs = [[a for a in int_eigs(G) if a >= 1] for G in mats]
    cnt = Counter()
    for e in eigs:
        for a in e:
            cnt[a] += 1
    print("integer eigenvalues >= 1 (eigenvalue: number of matrices):", dict(sorted(cnt.items())))
    withe = [i for i in range(len(mats)) if eigs[i]]
    print(f"matrices with an integer eigenvalue >= 1: {len(withe)} of {len(mats)}")

    def desc(i):
        g = gad[i]
        return f"#{i}(n={g['n']},X={g['adj']},P={g['P']:0{g['n']}b})"

    print()
    print("SINGLE GADGETS with an exact class h -> a h + b (census states in C on which h moves):")
    summ = defaultdict(list)
    for i in withe:
        G = mats[i]
        for a in eigs[i]:
            V, lam = class_space(L, [G], [a])
            dimC = 4 - V.dim()
            if dimC == 0:
                continue
            inC = [s for s in S if all(dot(r, s) == 0 for r in V.rows)]
            mov = [s for s in inC if a * s[2] + dot(lam[0], s) != s[2]]
            bs = Counter(dot(lam[0], s) for s in mov)
            kind = "b = 0" if all(b == 0 for b in bs) else "AFFINE b != 0"
            summ[(a, dimC, "no census state" if not mov else kind)].append(i)
            if mov and len(summ[(a, dimC, kind)]) <= 6:
                bl = sorted(bs.items())
                bshow = bl if len(bl) <= 10 else bl[:5] + ["..."] + bl[-4:]
                print(f"  {desc(i)}: a = {a}, dim C = {dimC}, C basis {V.annihilator()}, moving census states {len(mov)}, b: {bshow} [{kind}] e.g. {mov[:3]}")
    print("summary (multiplier a, dim C, kind): number of gadget matrices")
    for k2, v in sorted(summ.items()):
        print(f"   {k2}: {len(v)}")
    print()
    cat = Counter()
    free_pairs = []
    for x in range(len(withe)):
        i = withe[x]
        for y in range(x + 1, len(withe)):
            j = withe[y]
            for a1 in eigs[i]:
                for a2 in eigs[j]:
                    V, lam = class_space(L, [mats[i], mats[j]], [a1, a2])
                    dimC = 4 - V.dim()
                    if dimC == 0:
                        cat["(i) C = 0"] += 1
                        continue
                    kappa = [(a2 - 1) * p - (a1 - 1) * q for p, q in zip(lam[0], lam[1])]
                    if V.contains(kappa):
                        cat["(ii) C != 0, the two maps commute on all of C"] += 1
                        continue
                    cat["(iii) C != 0, kappa != 0 somewhere on C"] += 1
                    inC = [s for s in S if all(dot(r, s) == 0 for r in V.rows)]
                    good = [(s, dot(lam[0], s), dot(lam[1], s)) for s in inC if dot(kappa, s) != 0]
                    free_pairs.append((i, j, a1, a2, dimC, len(inC), good, V.annihilator()))
    print("PAIRS OF GADGETS (both with an integer eigenvalue >= 1):")
    for k2, v in sorted(cat.items()):
        print(f"   {k2}: {v}")
    nst = sum(1 for fp in free_pairs if fp[6])
    print(f"category (iii) pairs: {len(free_pairs)}; with a census state s in C, kappa s != 0: {nst}")
    for (i, j, a1, a2, dimC, ninC, good, basis) in free_pairs[:80]:
        print(f"   {desc(i)} & {desc(j)}: a = ({a1}, {a2}), dim C = {dimC}, C basis {basis}, census states in C: {ninC}, non-commuting: {good[:4]}")
    json.dump({"categories": dict(cat),
               "free_pairs": [[i, j, a1, a2, dimC, ninC, [[list(s), b1, b2] for (s, b1, b2) in good[:20]], basis] for (i, j, a1, a2, dimC, ninC, good, basis) in free_pairs]},
              open(f"gadget_pairs_{nmax}.json", "w"))
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
