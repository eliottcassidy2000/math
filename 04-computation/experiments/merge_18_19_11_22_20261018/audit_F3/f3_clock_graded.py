#!/usr/bin/env python3
"""Audit F3, item (e) extension: independent clock-graded ideal ladder (written from THM-4690's preimage rules and the
band conventions of clock_duality's 3c; no session code imported).

Depth-d node t = r (mod 3^(d+1)) in band b: M_d/2^(b+1) < t/N <= M_d/2^b, M_d = 4 (3/2)^d, weight 2^-(b+d+2) per class.
Clock-graded coefficient vector V_d(r, b)[tau] = sum over depth-d backward words of clock tau of the signed source values:
  V_0(r, b) = {0: h(r, e = -b)}                  (source band (2^(1-b), 2^(2-b)], representative 3*2^-b)
  V_d(r, b) = sum_j (-1)^(a_j) shift_(a_j + 2j) V_(d-1)(r_j, b + 1 - a_j - 2j)   (terms with negative band vanish),
with rho_j = B^j(r) mod 3^(d+1), a_j = 2 / 1 / none for rho_j = 1 / 2 / 0 (mod 3), r_j = (2^(a_j) rho_j - 1)/3 mod 3^d.
Depth-d meet sums (pairs of depth-d words of one node with different first-level children) by clock lag D:
  sum_b w_b sum_r [ A(V_d(r,b))[D] - sum_j A(V_(d-1)(r_j, b_j))[D] ]   (A = autocorrelation; D = 0 halved),
which equals (E_d - E_(d-1))/2 in total.  Reported in THM-4667's units (twice the pair sums), as in the note's table.
Usage: python3 f3_clock_graded.py DMAX BMAX
"""
import sys, os
import numpy as np
from fractions import Fraction as Fr
from functools import lru_cache

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from f3_band_exact import h                    # own ideal source coefficient (audit code)

DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
BMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 80
TLEN = BMAX + DMAX + 4                          # clock range (a depth-d word in band b has clock <= b + d)


def children(d, r, b):
    """(sign, clock, child residue mod 3^d, child band) for a depth-d node r mod 3^(d+1) in band b."""
    mod = 3 ** (d + 1)
    out = []
    rho, j = r % mod, 0
    while b + 1 - 1 - 2 * j >= 0:              # a >= 1, so the child band b + 1 - a - 2j < 0 for all further j
        if rho % 3 == 1:
            a = 2
        elif rho % 3 == 2:
            a = 1
        else:
            a = None
        if a is not None:
            bc = b + 1 - a - 2 * j
            if bc >= 0:
                v = (2 ** a * rho - 1) % mod
                out.append(((-1) ** a, a + 2 * j, (v // 3) % (mod // 3), bc))
        rho = (4 * rho + 1) % mod
        j += 1
    return out


@lru_cache(maxsize=None)
def V(d, r, b):
    vec = np.zeros(TLEN, dtype=np.int64)
    if d == 0:
        vec[0] = h(r % 3, -b)
        return vec
    for sg, tau, rc, bc in children(d, r, b):
        sub = V(d - 1, rc, bc)
        vec[tau:] += sg * sub[:TLEN - tau]
    return vec


@lru_cache(maxsize=None)
def AC(d, r, b):
    v = V(d, r, b)
    nz = np.flatnonzero(v)
    if len(nz) == 0:
        return np.zeros(TLEN, dtype=np.int64)
    lo, hi = nz[0], nz[-1] + 1
    w = v[lo:hi]
    c = np.correlate(w, w, 'full')[len(w) - 1:]
    out = np.zeros(TLEN, dtype=np.int64)
    out[:len(c)] = c
    return out


def main():
    out = []
    E = {}
    for d in range(0, DMAX + 1):
        mod = 3 ** (d + 1)
        Enum = 0
        lagnum = np.zeros(TLEN, dtype=object)
        for b in range(BMAX):
            sc = 2 ** (BMAX - b)                # weight 2^-(b+d+2) = sc / 2^(BMAX+d+2)
            tot_b = 0
            lag_b = np.zeros(TLEN, dtype=np.int64)
            for r in range(mod):
                v = V(d, r, b)
                s = int(v.sum())
                tot_b += s * s
                if d >= 1:
                    a = AC(d, r, b).copy()
                    for sg, tau, rc, bc in children(d, r, b):
                        a -= AC(d - 1, rc, bc)
                    lag_b += a
            Enum += tot_b * sc
            if d >= 1:
                lagnum += lag_b.astype(object) * sc
        den = 2 ** (BMAX + d + 2)
        E[d] = Fr(Enum, den)
        line = 'd = %d: E_d = %.12f' % (d, float(E[d]))
        if d >= 1:
            lags = [Fr(int(x), den) for x in lagnum]
            lags[0] = lags[0] / 2
            cross = 2 * sum(lags)
            line += ' ; cross term E_d - E_(d-1) = %.6f (direct %.6f); twice the pair sums by lag: ' % (
                float(cross), float(E[d] - E[d - 1]))
            line += ', '.join('lag %d %+.6f' % (D, float(2 * lags[D])) for D in range(0, 8))
        out.append(line)
        print(line, flush=True)
    open(os.path.join(HERE, 'f3_clock_graded.out'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
