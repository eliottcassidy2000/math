#!/usr/bin/env python3
"""Lane twoapex, deliverable (a), part 2: the minimal linear representation of the two-apex count over Z.

1. forward lattice  L_f = Z-span{ e0 M(w) } in Z^36 (rank d_f), invariant under the letter matrices;
2. backward lattice L_b = Z-span{ F(s) g } in the coordinates of L_f (rank r);
3. minimal representation: y(w) in Z^r, y(w l) = y(w) A_l, H(w) = y(w) . rho, with integer A_l;
4. Hankel rank from the counts alone (exact rational rank of [H(p s)], |p|, |s| <= 3): equals r, so r is minimal;
5. the bilinear (seam) law H(w1 w2) = <alpha(w1), beta(w2)>, alpha(w1) = y0 A(w1), beta(w2) = A(w2) rho, and its rank.
Writes minrep17.json (matrices, initial vector, readout).
"""
import itertools
import json
import random
import sys

from ta_core import FIN36, H36, INIT36, M36, dot, matvec, vecmat, words
from ta_lattice import Lattice, det_int, rank_frac


def main():
    # ---- 1. forward lattice
    Lf = Lattice(36)
    frontier = [list(INIT36)]
    Lf.add(INIT36)
    depth = 0
    stable = 0
    while stable < 3:
        nxt = []
        changed = False
        for v in frontier:
            for l in range(4):
                w = vecmat(v, M36[l])
                nxt.append(w)
                if Lf.add(w):
                    changed = True
        # keep the frontier small: only the vectors of the current basis matter for closure
        frontier = Lf.basis()
        depth += 1
        stable = 0 if changed else stable + 1
        if depth > 40:
            break
    Bf = Lf.basis()
    df = len(Bf)
    inv_ok = all(Lf.contains(vecmat(r, M36[l])) for r in Bf for l in range(4))
    print(f"forward lattice: rank {df}, stable after depth {depth}, invariant under the four letters: {inv_ok}")
    F = [[Lf.coords(vecmat(r, M36[l])) for r in Bf] for l in range(4)]       # row i: Bf[i] M_l in the basis
    c0 = Lf.coords(INIT36)
    g = [dot(r, FIN36) for r in Bf]

    # ---- 2. backward lattice (column vectors in Z^df)
    Lb = Lattice(df)
    Lb.add(g)
    stable = 0
    depth = 0
    while stable < 3:
        changed = False
        for v in Lb.basis():
            for l in range(4):
                if Lb.add(matvec(F[l], v)):
                    changed = True
        depth += 1
        stable = 0 if changed else stable + 1
        if depth > 40:
            break
    Bb = Lb.basis()
    r = len(Bb)
    inv_ok = all(Lb.contains(matvec(F[l], b)) for b in Bb for l in range(4))
    print(f"backward lattice inside the forward space: rank {r}, invariant: {inv_ok}")

    # ---- 3. minimal representation
    # F_l b_j = sum_k b_k A_l[k][j]
    A = []
    for l in range(4):
        cols = [Lb.coords(matvec(F[l], b)) for b in Bb]          # cols[j][k] = A_l[k][j]
        A.append([[cols[j][k] for j in range(r)] for k in range(r)])
    y0 = [dot(c0, b) for b in Bb]
    rho = Lb.coords(g)

    def Hmin(word):
        y = y0
        for l in word:
            y = vecmat(y, A[l])
        return dot(y, rho)

    bad = 0
    cnt = 0
    for n in range(0, 8):
        for w in words(n):
            cnt += 1
            if Hmin(w) != H36(w):
                bad += 1
    rng = random.Random(7)
    for _ in range(2000):
        w = [rng.randrange(4) for _ in range(rng.randrange(8, 60))]
        cnt += 1
        if Hmin(w) != H36(w):
            bad += 1
    print(f"minimal representation (dimension {r}) against the 36-state transfer: {cnt} words (all of length <= 7, "
          f"2000 random of length 8..59), mismatches {bad}")

    # ---- 4. Hankel rank from the counts alone
    for L in (1, 2, 3):
        ws = [w for n in range(L + 1) for w in words(n)]
        rows = [[H36(p + s) for s in ws] for p in ws]
        print(f"Hankel matrix [H(p s)], |p|, |s| <= {L} ({len(ws)} x {len(ws)}): exact rank {rank_frac(rows)}")

    # ---- 5. seam law
    bad = 0
    for _ in range(3000):
        w1 = [rng.randrange(4) for _ in range(rng.randrange(0, 12))]
        w2 = [rng.randrange(4) for _ in range(rng.randrange(0, 12))]
        al = y0
        for l in w1:
            al = vecmat(al, A[l])
        be = rho
        for l in reversed(w2):
            be = matvec(A[l], be)
        if dot(al, be) != H36(w1 + w2):
            bad += 1
    print(f"seam law H(w1 w2) = <alpha(w1), beta(w2)> on 3000 random pairs: violations {bad}")
    # the prefix vectors span Q^r and the suffix vectors span Q^r (so the pairing has rank exactly r)
    ws = [w for n in range(5) for w in words(n)]
    al_rows = []
    be_rows = []
    for w in ws:
        al = y0
        for l in w:
            al = vecmat(al, A[l])
        al_rows.append(al)
        be = rho
        for l in reversed(w):
            be = matvec(A[l], be)
        be_rows.append(be)
    print(f"rank of the prefix vectors alpha(w), |w| <= 4: {rank_frac(al_rows)}; of the suffix vectors beta(w): "
          f"{rank_frac(be_rows)}")
    # lattice indices: do the prefix vectors generate Z^r?  do the suffix vectors?
    La = Lattice(r)
    for v in al_rows:
        La.add(v)
    Lbe = Lattice(r)
    for v in be_rows:
        Lbe.add(v)
    ia = 1
    for row, c in zip(La.basis(), sorted(La.rows)):
        ia *= row[c]
    ib = 1
    for row, c in zip(Lbe.basis(), sorted(Lbe.rows)):
        ib *= row[c]
    print(f"index in Z^{r} of the lattice generated by the prefix vectors: {ia}; by the suffix vectors: {ib}")

    out = {"dimension": r, "A": A, "y0": y0, "rho": rho,
           "convention": "row vectors: y(w l) = y(w) A[l]; H(w) = y(w) . rho; y(empty) = y0; letter l = a + 2 b"}
    with open("minrep17_raw.json", "w") as f:
        json.dump(out, f)
    print("entries of the raw matrices: min", min(min(min(row) for row in M) for M in A), "max",
          max(max(max(row) for row in M) for M in A))
    print("y0 =", y0)
    print("rho =", rho)


if __name__ == "__main__":
    main()
