#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9665 (A). Direct test of the ingredient called "missing":
# joint equidistribution, on cells of B_n (i top digits and i bottom digits prescribed, r = n - 2i free digits),
# of (q(h) mod 2^(j-1) m, h mod m).  For every cell and every pair of residues the count should be
# 2^r / (2^(j-1) m^2) up to O((2 rho)^r), rho = max(rho_m, 2^(-1/2)).
# Reported: D(n, i) = max over cells and residue pairs of |count - 2^r/(2^(j-1) m^2)|, and D / 2^r.
import numpy as np, sys
from f1core import build_psi


def main():
    nmax = 22
    psi = build_psi(1 << nmax)
    h = np.arange(1 << nmax, dtype=np.int64)
    q = h - psi
    print("# mixed_cells.py (auditor hypotheses; own code)")
    for (j, m) in ((2, 3), (3, 3), (2, 5), (2, 9), (2, 7), (2, 15)):
        Mq = (1 << (j - 1)) * m
        ncl = Mq * m
        rows = []
        for n in (14, 18, 22):
            i = 3
            r = n - 2 * i
            lo, hi = 1 << (n - 1), 1 << n
            hh = h[lo:hi]
            qq = q[lo:hi]
            top = hh >> (n - i)             # i top digits (leading 1 included)
            bot = hh & ((1 << i) - 1)       # i bottom digits
            cell = (top - (1 << (i - 1))) * (1 << i) + bot
            ncell = (1 << (i - 1)) * (1 << i)
            key = cell * ncl + (qq % Mq) * m + (hh % m)
            cnt = np.bincount(key, minlength=ncell * ncl).reshape(ncell, ncl)
            assert cnt.sum(axis=1).min() == cnt.sum(axis=1).max() == (1 << r)
            D = np.abs(cnt - (1 << r) / ncl).max()
            rows.append((n, r, D, D / (1 << r) * ncl))
        dec = (rows[2][3] / rows[0][3]) ** (1.0 / (rows[2][1] - rows[0][1]))
        print("modulus 2^%d * %d (q mod %d, h mod %d): " % (j, m, Mq, m) + "; ".join(
            "n=%d r=%d: max dev %.1f, relative to the mean cell count %.4f" % row for row in rows) + ";  decay per free digit %.4f" % dec)


if __name__ == "__main__":
    main()
