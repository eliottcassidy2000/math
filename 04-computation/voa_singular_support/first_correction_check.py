#!/usr/bin/env python3
"""Exact check of the first quantum correction for Vir_{p,q}:
u = L_{-2}L_{-1}v - N L_{-3}v (v = singular vector, coeff(L_{-2}^N)=1) has degree-N symbol
  A*b_N + B*d^3(x0^N),  b_N = x0^{N-2}x1x2 + (2N-7)/12 x0^{N-3}x1^3,
  A = N(2N+1)(c - c_{2,2N+1})/8,  B = -(6N+c-10)/24."""
import sys
from fractions import Fraction
from vir_ss import singular_vector, apply_neg, central_charge

def apply_vec_neg(m, vec):
    out = {}
    for mono, co in vec.items():
        for mo, cc in apply_neg(m, mono):
            out[mo] = out.get(mo, 0) + co * cc
    return {k: v for k, v in out.items() if v}

def check(p, q):
    N = (p - 1) * (q - 1) // 2
    c = central_charge(p, q)
    v = singular_vector(p, q)
    u1 = apply_vec_neg(2, apply_vec_neg(1, v))
    u2 = apply_vec_neg(3, v)
    u = dict(u1)
    for k, val in u2.items():
        u[k] = u.get(k, 0) - N * val
    u = {k: val for k, val in u.items() if val}
    top = max(len(k) for k in u)
    degN = {k: val for k, val in u.items() if len(k) == N}
    m1 = (5,) + (2,) * (N - 1)
    m2 = (4, 3) + (2,) * (N - 2)
    m3 = (3, 3, 3) + (2,) * (N - 3)
    extra = {k: val for k, val in degN.items() if k not in (m1, m2, m3)}
    s1, s2, s3 = degN.get(m1, 0), degN.get(m2, 0), degN.get(m3, 0)
    A_pred = Fraction(N * (2 * N + 1), 8) * (c - (1 - Fraction(3 * (2 * N - 1) ** 2, 2 * N + 1)))
    B_pred = -(6 * N + c - 10) / 24
    # predicted coefficients
    p1 = B_pred * 6 * N
    p2 = A_pred + B_pred * 6 * N * (N - 1)
    p3 = A_pred * Fraction(2 * N - 7, 12) + B_pred * N * (N - 1) * (N - 2)
    ok = (top == N) and not extra and (s1, s2, s3) == (p1, p2, p3)
    print(f"({p},{q}) N={N} c={c}: top degree of u = {top}; deg-N coeffs = ({s1}, {s2}, {s3}); "
          f"predicted ({p1}, {p2}, {p3}); A={A_pred}; match={ok}")
    return ok

if __name__ == '__main__':
    cases = [(3,4),(3,5),(4,5),(3,7),(3,8),(4,7),(5,6),(2,7),(2,9),(2,11),(5,7),(3,10)]
    print(all(check(p, q) for p, q in cases))
