#!/usr/bin/env python3
"""Audit C: an independent (non-PARI) construction of the weight-3 newform g in S_3^new(Gamma_0(32), (-4/.)).

Construction
  f1(z) = eta(4z)^2 eta(8z)^2 in S_2(Gamma_0(32))           (Euler pentagonal products, exact integers)
  theta(z)^2 = sum r_2(n) q^n in M_1(Gamma_0(4), chi_-4)     (lattice-point count, exact)
  F_d = f1(z) * theta(dz)^2  (d = 1, 2, 4, 8)  in S_3(Gamma_0(32), chi_-4)  (cusp form x holomorphic form)
  plus eta(4z)^6, eta(8z)^6 (the CM oldforms from level 16).
  The span is checked to have rank 4 = dim S_3(Gamma_0(32), chi_-4) (Sturm bound 12; we compare 200 coefficients).
  T_3 is applied on q-expansions (a(3n) + chi(3) 3^2 a(n/3)), its matrix on the span is computed exactly, and g is
  the (normalised) eigenvector with eigenvalue 4i.  g = A + iB with A, B rational forms.
Output: b(p) = (Re, Im) for all primes p <= N, written to audit_newform_bp.json.
Checks: Hecke multiplicativity b(mn) = b(m) b(n) (coprime), b(p^2) = b(p)^2 - chi(p) p^2 for p^2 <= N,
        b(p^3) = b(p) b(p^2) - chi(p) p^2 b(p), and |b(p)| <= 2p (Deligne) for every p <= N.
Also: a_p(f1) = p + 1 - #E(F_p) for E: y^2 = x^3 - x, p < 2000 (the claim that f1 is the congruent-number curve).

Usage: python3 -u audit_newform.py N > audit_newform.out
"""
import sys
import json
import math
import time
from fractions import Fraction
import numpy as np
import sympy

OUT = '/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/audit_C_20261014/'


def euler_product(N, d, power):
    """prod_{n>=1} (1 - q^{d n})^power as int64 array of length N+1 (power >= 0)."""
    base = np.zeros(N + 1, dtype=np.int64)
    k = 0
    # pentagonal: prod (1 - x^n) = sum_k (-1)^k x^{k(3k-1)/2}, k in Z
    while True:
        e1 = k * (3 * k - 1) // 2
        e2 = k * (3 * k + 1) // 2
        if d * e1 > N and d * e2 > N:
            break
        s = -1 if k % 2 else 1
        if d * e1 <= N:
            base[d * e1] += s
        if k > 0 and d * e2 <= N:
            base[d * e2] += s
        k += 1
    out = np.zeros(N + 1, dtype=np.int64)
    out[0] = 1
    nz = np.nonzero(base)[0]
    for _ in range(power):
        new = np.zeros(N + 1, dtype=np.int64)
        for e in nz:
            new[e:] += base[e] * out[:N + 1 - e]
        out = new
    return out


def theta_sq(N, d):
    """theta(dz)^2 = sum_{x,y in Z} q^{d(x^2+y^2)} up to q^N."""
    out = np.zeros(N + 1, dtype=np.int64)
    m = int(np.sqrt(N / d)) + 1
    xs = np.arange(-m, m + 1)
    for x in xs:
        v = d * (x * x + xs * xs)
        v = v[v <= N]
        np.add.at(out, v, 1)
    return out


def mul_exact_fft(a, b, N):
    """exact product of integer series (check rounding)."""
    L = 1
    while L < 2 * (N + 1):
        L *= 2
    fa = np.fft.rfft(a.astype(np.float64), L)
    fb = np.fft.rfft(b.astype(np.float64), L)
    c = np.fft.irfft(fa * fb, L)[:N + 1]
    r = np.rint(c)
    dev = float(np.max(np.abs(c - r)))
    assert dev < 0.05, dev
    return r.astype(np.int64), dev


def shift(series, s, N):
    out = np.zeros(N + 1, dtype=np.int64)
    out[s:] = series[:N + 1 - s]
    return out


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
    t0 = time.time()
    # f1 = q * prod(1-q^{4n})^2 (1-q^{8n})^2
    f1 = shift(mul_exact_fft(euler_product(N, 4, 2), euler_product(N, 8, 2), N)[0], 1, N)
    # E: y^2 = x^3 - x point counts for p < 2000
    badE = []
    for p in sympy.primerange(3, 2000):
        cnt = 1 + sum(1 + (0 if (x ** 3 - x) % p == 0 else (1 if pow((x ** 3 - x) % p, (p - 1) // 2, p) == 1 else -1))
                      for x in range(p))
        if f1[p] != p + 1 - cnt:
            badE.append(p)
    print(f"f1 = eta(4z)^2 eta(8z)^2: first coefficients {f1[:20].tolist()}; a_p(f1) = p+1-#E(F_p) for y^2 = x^3 - x "
          f"at all odd p < 2000: failures {badE}", flush=True)
    eta4_6 = shift(euler_product(N, 4, 6), 1, N)
    eta8_6 = shift(euler_product(N, 8, 6), 2, N)
    basis = {}
    devs = []
    for d in (1, 2, 4, 8):
        prod, dev = mul_exact_fft(f1, theta_sq(N, d), N)
        basis[f'F{d}'] = prod
        devs.append(dev)
    basis['eta4^6'] = eta4_6
    basis['eta8^6'] = eta8_6
    names = list(basis)
    L = 200
    Mt = sympy.Matrix([[int(basis[nm][n]) for n in range(1, L + 1)] for nm in names])
    print(f"rank of span{{F1, F2, F4, F8, eta(4z)^6, eta(8z)^6}} on 200 coefficients: {Mt.rank()} "
          f"(dim S_3(32, chi_-4) = 4); FFT rounding max gap {max(devs):.2e}", flush=True)
    # choose a basis of 4 among F1, F2, F4, F8 if they are independent
    B4 = ['F1', 'F2', 'F4', 'F8']
    MB = sympy.Matrix([[int(basis[nm][n]) for n in range(1, L + 1)] for nm in B4])
    assert MB.rank() == 4
    # every other generator must lie in the span (consistency with dim 4)
    for nm in ('eta4^6', 'eta8^6'):
        aug = MB.col_join(sympy.Matrix([[int(basis[nm][n]) for n in range(1, L + 1)]]))
        assert aug.rank() == 4, nm
    # T_3 on q-expansions, compare on n = 1..L (needs coefficients to 3L)
    def T3(series):
        out = []
        for n in range(1, L + 1):
            v = int(series[3 * n])
            if n % 3 == 0:
                v += -1 * 9 * int(series[n // 3])     # chi_-4(3) = -1, p^(k-1) = 9
            out.append(v)
        return sympy.Matrix([out])
    # solve T3(F_j) = sum_i c_ij F_i using the first L coefficients
    MBt = MB.T                                          # L x 4
    Tmat = sympy.zeros(4, 4)
    for j, nm in enumerate(B4):
        rhs = T3(basis[nm]).T
        sol, params = MBt.gauss_jordan_solve(rhs)
        assert len(params) == 0
        assert MBt * sol == rhs, nm                       # exact fit on all L coefficients
        Tmat[:, j] = sol
    print("T_3 matrix on (F1, F2, F4, F8) (column j = coordinates of T_3 F_j):", Tmat.tolist())
    print("characteristic polynomial of T_3:", sympy.factor(Tmat.charpoly().as_expr()))
    I = sympy.I
    ev = (Tmat - 4 * I * sympy.eye(4)).nullspace()
    assert len(ev) == 1
    v = ev[0]
    gser = [sympy.nsimplify(sum(v[i] * int(basis[B4[i]][n]) for i in range(4))) for n in range(0, 40)]
    lead = gser[1]
    v = v / lead
    coeffs = [sympy.simplify(c) for c in v]
    print("g = sum c_i F_i with c =", coeffs)
    # real and imaginary rational parts
    cre = [sympy.re(c) for c in coeffs]
    cim = [sympy.im(c) for c in coeffs]
    den = sympy.ilcm(*[sympy.fraction(sympy.Rational(x))[1] for x in cre + cim])
    cre_i = [int(sympy.Rational(x) * den) for x in cre]
    cim_i = [int(sympy.Rational(x) * den) for x in cim]
    Aser = sum(cre_i[i] * basis[B4[i]].astype(object) for i in range(4))
    Bser = sum(cim_i[i] * basis[B4[i]].astype(object) for i in range(4))
    assert all(int(Aser[n]) % den == 0 and int(Bser[n]) % den == 0 for n in range(N + 1))
    Aser = [int(x) // den for x in Aser]
    Bser = [int(x) // den for x in Bser]
    b = lambda n: complex(Aser[n], Bser[n])
    print("g = " + " + ".join(f"({Aser[n]}{Bser[n]:+d}i)q^{n}" for n in range(1, 24) if Aser[n] or Bser[n]), flush=True)
    # Hecke checks
    chi4 = lambda n: 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
    nbad = 0
    nchk = 0
    rng = np.random.default_rng(20261014)
    for _ in range(200000):
        m = int(rng.integers(2, 400))
        n = int(rng.integers(2, N // m + 1)) if N // m >= 2 else 2
        if m * n > N or math.gcd(m, n) != 1:
            continue
        lhs = (Aser[m * n], Bser[m * n])
        bm, bn = (Aser[m], Bser[m]), (Aser[n], Bser[n])
        rhs = (bm[0] * bn[0] - bm[1] * bn[1], bm[0] * bn[1] + bm[1] * bn[0])
        nchk += 1
        nbad += lhs != rhs
    primes = list(sympy.primerange(2, N + 1))
    bad2 = []
    for p in primes:
        if p * p > N:
            break
        bp = (Aser[p], Bser[p])
        bp2 = (bp[0] ** 2 - bp[1] ** 2 - chi4(p) * p * p, 2 * bp[0] * bp[1])
        if (Aser[p * p], Bser[p * p]) != bp2:
            bad2.append(p)
        if p ** 3 <= N:
            bp3 = (bp[0] * bp2[0] - bp[1] * bp2[1] - chi4(p) * p * p * bp[0],
                   bp[0] * bp2[1] + bp[1] * bp2[0] - chi4(p) * p * p * bp[1])
            if (Aser[p ** 3], Bser[p ** 3]) != bp3:
                bad2.append(('cube', p))
    deligne = [p for p in primes if Aser[p] ** 2 + Bser[p] ** 2 > 4 * p * p]
    realimag = [p for p in primes if p > 2 and ((p % 4 == 1 and Bser[p] != 0) or (p % 4 == 3 and Aser[p] != 0))]
    even = [n for n in range(2, N + 1, 2) if Aser[n] or Bser[n]]
    print(f"Hecke multiplicativity on {nchk} random coprime pairs mn <= {N}: failures {nbad}; "
          f"p^2, p^3 relations: failures {bad2}; Deligne |b(p)| <= 2p violations: {deligne}; "
          f"b(p) real for p = 1 mod 4 and imaginary for p = 3 mod 4: violations {realimag[:5]}; "
          f"nonzero even-index coefficients: {len(even)}", flush=True)
    json.dump({str(p): [Aser[p], Bser[p]] for p in primes}, open(OUT + 'audit_newform_bp.json', 'w'))
    print(f"wrote b(p) for {len(primes)} primes p <= {N}  [{time.time() - t0:.1f}s]")


if __name__ == '__main__':
    main()
