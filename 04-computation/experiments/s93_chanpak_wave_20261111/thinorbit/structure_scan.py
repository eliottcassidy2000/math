#!/usr/bin/env python3
"""Lane thinorbit (S93).  Search for 2 x 2 (continuant-like) structure inside the base-2 port-chain block library.

For a library of block matrices N (rows (h', E', O') in (h, E, O)), exact arithmetic:
 (1) parity lemma check: N = [[1,0,0],[*,x,y],[*,y,x]] (mod 2)                      (PROOFS.txt, Lemma P)
 (2) singular blocks of rank 2: characteristic polynomial t (t^2 - a t + b); the induced 2 x 2 action has
     determinant b; check a odd, b even; list the (a, b) that occur and the minimum of |b|
 (3) rank-2 families with a common kernel line: induced 2 x 2 matrices on Z^3 / kernel
 (4) unimodular blocks: characteristic polynomials, rational eigenvalues, invariant lines and planes;
     pairs of unimodular blocks sharing an invariant line or plane (a reducible unimodular sub-semigroup)
 (5) blocks whose seam matrix has e-column (0, +-1, 0)^T (the (h, B) part would evolve autonomously and unimodularly)
Usage: structure_scan.py MATFILE
"""
import sys
from collections import Counter, defaultdict
from math import gcd

import sympy as sp


def primitive(v):
    v = [sp.Rational(x) for x in v]
    den = 1
    for x in v:
        den = sp.ilcm(den, x.q)
    v = [int(x * den) for x in v]
    g = 0
    for x in v:
        g = gcd(g, abs(x))
    if g:
        v = [x // g for x in v]
    for x in v:
        if x != 0:
            if x < 0:
                v = [-y for y in v]
            break
    return tuple(v)


def main():
    path = sys.argv[1]
    mats = []
    reps = []
    for line in open(path + ".reps"):
        a, b = line.split("|")
        mats.append(sp.Matrix(3, 3, list(map(int, a.split()))))
        reps.append(b.strip())
    n = len(mats)
    t = sp.symbols("t")
    print(f"{n} block matrices from {path}")
    # (1) parity lemma
    pat = Counter()
    ok = True
    for M in mats:
        m = [int(x) % 2 for x in M]
        if not (m[0] == 1 and m[1] == 0 and m[2] == 0 and m[4] == m[8] and m[5] == m[7]):
            ok = False
        pat[(m[4], m[5])] += 1
    print("(1) parity lemma N = [[1,0,0],[*,x,y],[*,y,x]] mod 2 holds for all:", ok, " (x, y) counts:", dict(pat))
    dets = [int(M.det()) for M in mats]
    ranks = [M.rank() for M in mats]
    print("    determinants:", sorted(Counter(dets).items()))
    print("    ranks:", sorted(Counter(ranks).items()))
    print("    det parity = x + y for all:", all((dets[i] - (int(mats[i][1, 1]) + int(mats[i][1, 2]))) % 2 == 0 for i in range(n)))
    print("    trace odd for all:", all(int(M.trace()) % 2 == 1 for M in mats))
    # (2) rank 2
    ab = Counter()
    for i, M in enumerate(mats):
        if ranks[i] == 2:
            a = int(M.trace())
            b = int(sum(M.minor(j, j) for j in range(3)))
            ab[(a, b)] += 1
    print("(2) rank-2 blocks: t (t^2 - a t + b); (a, b): count")
    for (a, b), c in sorted(ab.items()):
        print(f"      a = {a:4d}  b = {b:5d}  x{c}")
    print("    all a odd:", all(a % 2 == 1 for (a, b) in ab), "; all b even:", all(b % 2 == 0 for (a, b) in ab),
          "; min |b| over b != 0:", min(abs(b) for (a, b) in ab if b != 0), "; b = 0 cases:", sum(c for (a, b), c in ab.items() if b == 0))
    # (3) common kernel families
    fam = defaultdict(list)
    for i, M in enumerate(mats):
        if ranks[i] == 2:
            fam[primitive(list(M.nullspace()[0]))].append(i)
        elif ranks[i] == 1:
            pass
    print("(3) rank-2 blocks by kernel line: kernel -> number of blocks; induced 2 x 2 matrices (distinct) and their determinants")
    for k, idx in sorted(fam.items(), key=lambda kv: -len(kv[1])):
        # quotient basis: complete k to a unimodular basis: find integer matrix U with last column k, det +-1
        kv = sp.Matrix(3, 1, list(k))
        # functionals vanishing on k: integer basis of the lattice {f : f.k = 0}
        F = sp.Matrix([list(k)]).nullspace()
        Fi = [primitive(list(f)) for f in F]
        Fm = sp.Matrix([list(Fi[0]), list(Fi[1])])
        # saturate: use Smith/Hermite via sympy
        from sympy.matrices.normalforms import smith_normal_form
        snf = smith_normal_form(Fm, domain=sp.ZZ)
        sat = (abs(snf[0, 0]) == 1 and abs(snf[1, 1]) == 1)
        induced = Counter()
        for i in idx:
            M = mats[i]
            # coordinates y = Fm s on the quotient (valid if Fm is saturated); M induces A with Fm M = A Fm
            A = (Fm * M) * Fm.pinv()
            A = sp.Matrix(2, 2, [sp.nsimplify(x) for x in A])
            assert Fm * M == A * Fm
            induced[tuple(A)] += 1
        print(f"    kernel {k}: {len(idx)} blocks; quotient coordinates {Fi} (saturated: {sat}); {len(induced)} distinct induced matrices")
        for A, c in sorted(induced.items(), key=lambda kv: (abs(kv[0][0] * kv[0][3] - kv[0][1] * kv[0][2]), kv[0])):
            d = A[0] * A[3] - A[1] * A[2]
            print(f"        {list(A)}  det {d}  x{c}")
    # (4) unimodular blocks
    uni = [i for i in range(n) if abs(dets[i]) == 1]
    print(f"(4) unimodular blocks: {len(uni)}")
    lines = defaultdict(list)
    planes = defaultdict(list)
    for i in uni:
        M = mats[i]
        cp = sp.factor(M.charpoly(t).as_expr())
        info = []
        for lam in (1, -1):
            ns = (M - lam * sp.eye(3)).nullspace()
            for v in ns:
                pv = primitive(list(v))
                lines[pv].append((i, lam))
                info.append(f"eigenvalue {lam}: line {pv}")
            ns = (M.T - lam * sp.eye(3)).nullspace()
            for v in ns:
                pv = primitive(list(v))
                planes[pv].append((i, lam))
                info.append(f"plane normal {pv} (factor {lam})")
        print(f"    {i:3d} {list(M)} det {dets[i]} charpoly {cp} [{reps[i]}] {'; '.join(info)}")
    print("    invariant lines shared by >= 2 unimodular blocks:", {k: v for k, v in lines.items() if len(v) >= 2})
    print("    invariant planes shared by >= 2 unimodular blocks:", {k: v for k, v in planes.items() if len(v) >= 2})
    # all invertible blocks: rational eigen-lines and planes (any rational eigenvalue)
    lines = defaultdict(list)
    planes = defaultdict(list)
    for i in range(n):
        if dets[i] == 0:
            continue
        M = mats[i]
        for lam, mult in sp.roots(M.charpoly(t).as_expr(), t, filter="Q").items():
            for v in (M - lam * sp.eye(3)).nullspace():
                lines[primitive(list(v))].append((i, lam))
            for v in (M.T - lam * sp.eye(3)).nullspace():
                planes[primitive(list(v))].append((i, lam))
    print("    invertible blocks: invariant lines shared by >= 2:", {k: v for k, v in lines.items() if len(v) >= 2})
    print("    invertible blocks: invariant planes shared by >= 2:", {k: v for k, v in planes.items() if len(v) >= 2})
    # (5) seam e-column
    P = sp.Matrix([[1, 0, 0], [0, -1, 1], [0, 1, 1]])
    Pinv = P.inv()
    hits = []
    zero_e = 0
    for i, M in enumerate(mats):
        S = P * M * Pinv
        col = (S[0, 1], S[1, 1], S[2, 1])
        if col[0] == 0 and col[2] == 0:
            if col[1] == 0:
                zero_e += 1
            if abs(col[1]) == 1:
                hits.append((i, list(S), dets[i]))
    print(f"(5) seam matrices with e-column (0, c, 0): c = 0 for {zero_e} blocks (the (h, B)-closed singular blocks); c = +-1:")
    for h in hits:
        print("     ", h, reps[h[0]])
    if not hits:
        print("      none")


if __name__ == "__main__":
    main()
