#!/usr/bin/env python3
"""Verifier (lane digits), check V2: class spaces of the one-mark block system, own exact linear algebra.

(a) integer eigenvalues >= 1 of the 1,922 block matrices (rational root theorem on the exact characteristic polynomial);
(b) single operations, CENSUS-FREE: for every (N, a): is the offset functional lambda identically 0 on C(N; a)?
    (the lane decides "b = 0 on the class" on census states of order <= 6 only);
(c) pairs: categories (i) C = 0, (ii) kappa = 0 on C, (iii) kappa != 0 somewhere on C;
(d) pairs of translations (a1 = a2 = 1) commute trivially, so (ii) says nothing there: the offsets are examined directly.
"""
import json
import sys
import time
from collections import Counter, defaultdict

from vlin import class_space, dot, int_eigs, rowmat


def ann_basis(V, d):
    """basis of { s : rho s = 0 for all rows } by brute elimination over Q (small d)."""
    from fractions import Fraction
    rows = [[Fraction(x) for x in r] for r in V.rows.values()]
    piv = []
    R = []
    for r in rows:
        for rr, p in zip(R, piv):
            if r[p] != 0:
                f = r[p]
                r = [x - f * y for x, y in zip(r, rr)]
        p = next((i for i, x in enumerate(r) if x != 0), None)
        if p is None:
            continue
        f = r[p]
        r = [x / f for x in r]
        for i in range(len(R)):
            if R[i][p] != 0:
                g = R[i][p]
                R[i] = [x - g * y for x, y in zip(R[i], r)]
        R.append(r)
        piv.append(p)
    out = []
    for f in range(d):
        if f in piv:
            continue
        v = [Fraction(0)] * d
        v[f] = Fraction(1)
        for rr, p in zip(R, piv):
            v[p] = -rr[f]
        den = 1
        for x in v:
            den = den * x.denominator
        out.append([int(x * den) for x in v])
    return out


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    t0 = time.time()
    here = sys.path[0]
    ops = json.load(open(f"{here}/../onemark_matrices_{kmax}.json"))
    mats = [o["N"] for o in ops]
    L = [1, 0, 0, 0]
    eigs = [[a for a in int_eigs(G) if a >= 1] for G in mats]
    cnt = Counter(a for e in eigs for a in e)
    print(f"{len(mats)} matrices; integer eigenvalues >= 1 (exact): {dict(sorted(cnt.items()))}")
    print(f"   (matrix, eigenvalue) incidences: {sum(cnt.values())}; matrices with 1 / 2 / 3+ such eigenvalues: "
          f"{sum(1 for e in eigs if len(e) == 1)} / {sum(1 for e in eigs if len(e) == 2)} / {sum(1 for e in eigs if len(e) > 2)}")

    def desc(i):
        o = ops[i]
        return f"#{i}(k={o['k']},X={o['adjX']},W={o['W']:0{o['k']}b},m={o['m']})"

    # ---------------- (b) single ops, census-free
    kinds = Counter()
    affine_capable = []
    for i, G in enumerate(mats):
        for a in eigs[i]:
            V, lam = class_space(L, [G], [a])
            dimC = 4 - V.dim()
            if dimC == 0:
                kinds[(a if a == 1 else "a>1", "C = 0")] += 1
                continue
            if not any(lam[0]):
                kinds[(a if a == 1 else "a>1", "lambda = 0 identically")] += 1
            elif V.contains(lam[0]):
                kinds[(a if a == 1 else "a>1", "lambda = 0 on C (b = 0 for EVERY state of the class)")] += 1
            else:
                kinds[(a if a == 1 else "a>1", "lambda != 0 somewhere on C")] += 1
                affine_capable.append((i, a, dimC, ann_basis(V, 4), lam[0]))
    print("(b) single operations, census-free classification of (matrix, a):")
    for k2, v in sorted(kinds.items(), key=str):
        print(f"      {k2}: {v}")
    print(f"    (matrix, a) with lambda not identically 0 on C: {len(affine_capable)}")
    # refine: does H move on C?  (L N - L) s != 0 for some s in C;  can C meet the cone of actual states?
    #   valid for every marked tournament: g >= 1, h >= 1, 0 <= e <= g, e <= h, B >= h - e + g   (cut a path at z)
    from scipy.optimize import linprog
    import numpy as np
    summary = Counter()
    survivors = []
    for (i, a, dimC, basis, lam) in affine_capable:
        G = mats[i]
        V, _ = class_space(L, [G], [a])
        mu = [x - y for x, y in zip(rowmat(L, G), L)]          # H' - H
        moving = not V.contains(mu)
        o = ops[i]
        full = o["W"] == (1 << o["k"]) - 1
        trans = sorted(o["adjX"]) == sorted((1 << j) - 1 for j in range(o["k"]))
        ray = full and trans and o["m"] == -1 and a == 2 ** o["k"]
        g0 = all(b[3] == 0 for b in basis)
        # LP: s = sum t_j basis_j in the cone;  maximise / minimise lambda s and mu s
        Bm = np.array(basis, dtype=float).T          # 4 x dimC
        # constraints A s <= c:  -g <= -1, -h <= -1, -e <= 0, e - g <= 0, e - h <= 0, h - e + g - B <= 0
        A = np.array([[0, 0, 0, -1], [-1, 0, 0, 0], [0, -1, 0, 0], [0, 1, 0, -1], [-1, 1, 0, 0], [1, -1, -1, 1]], dtype=float)
        c = np.array([-1, -1, 0, 0, 0, 0], dtype=float)
        feas = linprog(np.zeros(dimC), A_ub=A @ Bm, b_ub=c, bounds=[(None, None)] * dimC, method="highs")
        feasible = feas.status == 0
        lam_nonzero = False
        mu_nonzero = False
        if feasible:
            for fun, name in ((lam, "lam"), (mu, "mu")):
                obj = np.array(fun, dtype=float) @ Bm
                for sgn in (1, -1):
                    r = linprog(sgn * obj, A_ub=A @ Bm, b_ub=c, bounds=[(None, None)] * dimC, method="highs")
                    if r.status == 3 or (r.status == 0 and abs(r.fun) > 1e-7):
                        if name == "lam":
                            lam_nonzero = True
                        else:
                            mu_nonzero = True
        key = ("ray a^k" if ray else "other", "a=1" if a == 1 else "a>1", "H moves on C" if moving else "H fixed on C",
               "C in {g=0}" if g0 else "C not in {g=0}", "cone feasible" if feasible else "cone EMPTY",
               "b != 0 possible in cone" if lam_nonzero else "b = 0 on cone", "H moves in cone" if mu_nonzero else "H fixed in cone")
        summary[key] += 1
        if (not ray) and feasible and lam_nonzero and mu_nonzero:
            survivors.append((i, a, dimC, basis, lam))
    print("    refined classification of the (matrix, a) with lambda not identically 0 on C  [LP over the cone of actual states]:")
    for k2, v in sorted(summary.items()):
        print(f"      {k2}: {v}")
    print(f"    NON-RAY operations whose class space meets the cone with b != 0 and H moving (LP): {len(survivors)}")
    for (i, a, dimC, basis, lam) in survivors[:30]:
        print(f"       {desc(i)} a = {a}, dim C = {dimC}, C basis (h,e,B,g) {basis}, lambda {lam}")
    # ---------------- (c) pairs
    cat = Counter()
    free = []
    trans_pairs = 0
    for i in range(len(mats)):
        if not eigs[i]:
            continue
        for j in range(i + 1, len(mats)):
            for a1 in eigs[i]:
                for a2 in eigs[j]:
                    V, lam = class_space(L, [mats[i], mats[j]], [a1, a2])
                    if V.dim() == 4:
                        cat["(i) C = 0"] += 1
                        continue
                    kappa = [(a2 - 1) * x - (a1 - 1) * y for x, y in zip(lam[0], lam[1])]
                    if V.contains(kappa):
                        cat["(ii) kappa = 0 on C"] += 1
                        if a1 == 1 and a2 == 1:
                            trans_pairs += 1
                            # two translations: is either offset nonzero somewhere on C?
                            if not (V.contains(lam[0]) and V.contains(lam[1])):
                                cat["   of (ii): both a = 1 and some offset not identically 0 on C"] += 1
                    else:
                        cat["(iii) kappa != 0 somewhere on C"] += 1
                        free.append((i, j, a1, a2))
    print("(c) pairs (i < j, integer eigenvalues a1, a2 >= 1):")
    for k2, v in sorted(cat.items()):
        print(f"      {k2}: {v}")
    print(f"    pairs with a1 = a2 = 1 inside (ii): {trans_pairs}")
    print(f"    category (iii): {len(free)} {[(desc(i), desc(j), a1, a2) for (i, j, a1, a2) in free[:10]]}")
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
