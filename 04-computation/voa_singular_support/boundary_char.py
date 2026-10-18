#!/usr/bin/env python3
"""Test the boundary-level product formula ch L_k(g) = Dhat(q^u, z)/Dhat(q, z) (unflavoured: total dims per weight)
for g = sl_n at k = n/u - n, against direct engine dims.  Dhat(q,z) = prod_n (1-q^n)^r prod_{a>0} (1 - e^a q^n)(1 - e^-a q^(n-1)).
We evaluate the z-dependence exactly as Laurent polynomials in the root lattice."""
import sys, itertools
from collections import defaultdict

def roots_sl(n):
    pos = []
    for i in range(n):
        for j in range(i + 1, n):
            v = [0] * n; v[i] += 1; v[j] -= 1
            pos.append(tuple(v))
    return pos

def series_Dhat(n, u, Dmax, invert):
    r = n - 1
    pos = roots_sl(n)
    zero = tuple([0] * n)
    ser = {(0, zero): 1}
    def mulfac(ser, qp, vec, sign):
        # multiply by (1 - x) with x = q^qp e^vec, or by 1/(1-x)
        out = defaultdict(int)
        for (d, z), c in ser.items():
            if not invert:
                out[(d, z)] += c
                if d + qp <= Dmax:
                    out[(d + qp, tuple(a + b * sign for a, b in zip(z, vec)))] -= c
            else:
                kk = 0
                while d + kk * qp <= Dmax:
                    out[(d + kk * qp, tuple(a + kk * b * sign for a, b in zip(z, vec)))] += c
                    kk += 1
        return {kk: v for kk, v in out.items() if v}
    for m in range(1, Dmax // u + 2):
        qp = u * m
        if qp <= Dmax:
            for _ in range(r):
                ser = mulfac(ser, qp, zero, 1)
        for a in pos:
            if qp <= Dmax:
                ser = mulfac(ser, qp, a, 1)          # (1 - e^a q^{um})
            qq = u * (m - 1)
            if qq <= Dmax and qq > 0:
                ser = mulfac(ser, qq, a, -1)          # (1 - e^-a q^{u(m-1)})
    return ser

if __name__ == '__main__':
    n, u, Dmax = map(int, sys.argv[1:4])
    # numerator: Dhat(q^u) without the q^0 factor prod_{a>0}(1 - e^-a) (cancels with denominator's)
    num = series_Dhat(n, u, Dmax, invert=False)
    den_inv = series_Dhat(n, 1, Dmax, invert=True)
    prod = defaultdict(int)
    for (d1, z1), c1 in num.items():
        for (d2, z2), c2 in den_inv.items():
            if d1 + d2 <= Dmax:
                prod[(d1 + d2, tuple(a + b for a, b in zip(z1, z2)))] += c1 * c2
    tot = [0] * (Dmax + 1)
    for (d, z), c in prod.items():
        tot[d] += c
    print(f"sl{n}, u={u}, k={n}/{u}-{n}: dims", tot, " min coeff", min(prod.values()))
