#!/usr/bin/env python3
"""Auditor proof-5135, script a7: two supporting checks.

  [L] Lemma behind THM-5135 (3)(d).  Claim to be supplied: an integral quaternary quadratic form with the class
      counts of the seam a1 a2 + b1 b2 modulo p (already: with its number of zeros modulo p) is non-degenerate
      modulo p of hyperbolic type; so a form with the seam counts at every prime is even unimodular of rank 4,
      hence indefinite.  Odd p: every form is diagonalisable over F_p; zero counts of all diag(a1..a4), a_i in
      F_p (zeros allowed), against p^3 + p^2 - p.  (p = 2, 3: all forms were enumerated in a2_squares.py [Q].)
  [O] Orders.  The norm forms on the orders spanned by the four units e_a (integer coordinates):
        #S = 3: x0^2 + x1^2 + x2^2 + x3^2 (Lipschitz order)            -> r_4
        #S = 1: x0^2 + x1^2 - x2^2 - x3^2 (index 4 in M_2(Z))          vs the determinant a1 a2 + b1 b2 on M_2(Z)
        #S = 2: ((a+d)^2 + (b-c)^2)((a-d)^2 + (b+c)^2) on Z[1,i,j,ij]  vs |z1|^2 |z2|^2 on Z[i]^2 (index 4)
      Class counts modulo 2, 4, 8 and odd p of the #S = 1 diagonal form against the seam; unit groups.
"""
import sys
import itertools

import numpy as np


def table(m):
    cnt = [0] * m
    for a in range(m):
        for b in range(m):
            cnt[(a * b) % m] += 1
    return cnt


def cyc_conv(a, b):
    m = len(a)
    out = [0] * m
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[(i + j) % m] += ai * bj
    return out


def main():
    out = []
    P = out.append

    P("[L] zero counts of diagonal quaternary forms modulo odd p against the seam's p^3 + p^2 - p")
    for p in (3, 5, 7, 11, 13):
        target = p ** 3 + p ** 2 - p
        sq = [0] * p
        for a in range(p):
            sq[(a * a) % p] += 1
        qr = set((a * a) % p for a in range(1, p))
        n_match = 0
        n_match_good = 0
        n_good = 0
        for coeffs in itertools.product(range(p), repeat=4):
            dist = None
            for c in coeffs:
                d = [0] * p
                for t in range(p):
                    d[(c * t) % p] += sq[t]
                dist = d if dist is None else cyc_conv(dist, d)
            nondeg = all(c % p for c in coeffs)
            prod = 1
            for c in coeffs:
                prod = (prod * c) % p
            good = nondeg and (prod in qr)
            n_good += good
            if dist[0] == target:
                n_match += 1
                n_match_good += good
        P("  p = %2d: %d diagonal forms; non-degenerate with square discriminant: %d; with exactly p^3 + p^2 - p zeros: %d, all of them non-degenerate with square discriminant: %s"
          % (p, p ** 4, n_good, n_match, n_match == n_match_good == n_good))

    P("")
    P("[O] norm forms on the orders spanned by the units")
    for m in (2, 4, 8, 16, 3, 9, 5, 7):
        sq = [0] * m
        nsq = [0] * m
        for a in range(m):
            sq[(a * a) % m] += 1
            nsq[(-a * a) % m] += 1
        split_basis = cyc_conv(cyc_conv(sq, sq), cyc_conv(nsq, nsq))
        seam = cyc_conv(table(m), table(m))
        four = cyc_conv(cyc_conv(sq, sq), cyc_conv(sq, sq))
        P("  modulo %2d: x0^2+x1^2-x2^2-x3^2 equals the seam a1a2+b1b2: %s; equals four squares: %s%s"
          % (m, split_basis == seam, split_basis == four,
             ("; counts %s vs seam %s" % (split_basis, seam)) if m in (2, 4) else ""))
    # units
    R = 2
    units_basis = 0
    units_max = 0
    for a, b, c, d in itertools.product(range(-R, R + 1), repeat=4):
        if ((a + d) ** 2 + (b - c) ** 2) * ((a - d) ** 2 + (b + c) ** 2) == 1:
            units_basis += 1
    for a, b, c, d in itertools.product(range(-R, R + 1), repeat=4):
        if (a * a + b * b) * (c * c + d * d) == 1:
            units_max += 1
    P("  bicomplex units: order Z[1,i,j,ij]: %d (the extension group C_4 x C_2); maximal order Z[i]^2: %d (C_4 x C_4)" % (units_basis, units_max))
    P("  index of Z[1,i,j,ij] in Z[i]^2: |det [[1,0,0,1],[0,1,-1,0],[1,0,0,-1],[0,1,1,0]]| = %d"
      % abs(round(np.linalg.det(np.array([[1, 0, 0, 1], [0, 1, -1, 0], [1, 0, 0, -1], [0, 1, 1, 0]], dtype=float)))))

    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
