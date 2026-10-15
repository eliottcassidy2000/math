#!/usr/bin/env python3
"""Wheels are power sums of Jacobi sums; W_4 and W_6 are level-8 modular forms.

PROVED (telescoping x_{i+1} = t_i x_i, then Fourier analysis on the cyclic group F_q^x):
    T_n := Tr((DS)^n) = sum_{x in F_q^n} prod_i chi(x_i) chi(x_{i+1} - x_i) = sum_{psi} (chi(-1) J(psi, chi))^n,
where psi runs over all q - 1 multiplicative characters and J(psi, chi) = sum_t psi(t) chi(1 - t).  U(W_n) = q T_n.
So the nonzero eigenvalues of DS are exactly the Jacobi sums chi(-1) J(psi, chi).
FINITE-EXACT here (q = 3 mod 4, all primes 7 <= q <= QMAX):
    T_4 = -(q - 1)(q + a_q),               a = eta(2z)^4 eta(4z)^4                         (weight 4, level 8)
    T_6 = (q - 1)(q^2 - q a_q - b_q),       b = the newform of S_6(Gamma_0(8)), b_3 = 20     (weight 6, level 8)
b is computed as the T_3-eigenvector with eigenvalue 20 in the span of a(z) E_2-type weight-2 products.
Session opus-2026-10-14-S35.
"""
import cmath, math, sys
from sympy import primerange, primitive_root
import sympy as sp

QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
M = QMAX + 2


def jacobi_traces(q, ns):
    g = primitive_root(q)
    dlog = [0] * q
    x = 1
    for k in range(q - 1):
        dlog[x] = k; x = x * g % q
    chi = [0] + [1 if pow(t, (q - 1) // 2, q) == 1 else -1 for t in range(1, q)]
    sgn = chi[q - 1]
    out = {n: 0 for n in ns}
    for j in range(q - 1):  # psi_j(g^k) = exp(2 pi i j k / (q - 1))
        J = 0j
        for t in range(2, q):
            if chi[(1 - t) % q]:
                J += cmath.exp(2j * math.pi * j * dlog[t] / (q - 1)) * chi[(1 - t) % q]
        J += chi[0] if False else 0  # t = 1 contributes chi(0) = 0; t = 0 contributes psi(0) = 0
        v = sgn * J
        for n in ns:
            out[n] += v ** n
    return {n: round(out[n].real) for n in ns}


def series_mul(A, B):
    out = [0] * M
    for i, a in enumerate(A):
        if a:
            for j in range(M - i):
                if B[j]:
                    out[i + j] += a * B[j]
    return out


def forms():
    # a = eta(2z)^4 eta(4z)^4 = q prod (1 - q^{2n})^4 (1 - q^{4n})^4
    P = [0] * M; P[0] = 1
    for n in range(1, M):
        for step in (2 * n, 4 * n):
            if step >= M:
                continue
            for _ in range(4):
                for i in range(M - 1, step - 1, -1):
                    P[i] -= P[i - step]
    a = [0] + P[:M - 1]
    sig = lambda n: sum(d for d in range(1, n + 1) if n % d == 0)
    E2 = [1] + [-24 * sig(n) for n in range(1, M)]
    sc = lambda F, d: [F[n // d] if n % d == 0 else 0 for n in range(M)]
    w2 = [[x - 2 * y for x, y in zip(E2, sc(E2, 2))], [x - 4 * y for x, y in zip(E2, sc(E2, 4))],
          [x - 8 * y for x, y in zip(E2, sc(E2, 8))]]
    gens = [series_mul(a, g) for g in w2]
    # T_3 on span(gens) (dimension 3 = S_6(Gamma_0(8))); eigenvector for eigenvalue 20
    B = 40
    Bm = sp.Matrix([g[1:B + 1] for g in gens])
    T = []
    for g in gens:
        v = sp.Matrix([[g[3 * n] + (243 * g[n // 3] if n % 3 == 0 else 0) for n in range(1, B + 1)]])
        sol, params = Bm.T.gauss_jordan_solve(v.T)
        assert params.shape[0] == 0
        T.append(list(sol))
    T3 = sp.Matrix(T)
    ns = (T3.T - 20 * sp.eye(3)).nullspace()
    assert len(ns) == 1
    v = ns[0]
    b = [sum(v[i] * gens[i][n] for i in range(3)) for n in range(M)]
    b = [sp.Rational(c) / b[1] for c in b]
    assert all(c.q == 1 for c in b)
    return a, [int(c) for c in b], sp.factor(T3.charpoly().as_expr())


if __name__ == '__main__':
    a, b, cp = forms()
    print('S_6(Gamma_0(8)) via a(z) x weight-2 forms: T_3 charpoly', cp)
    print('b = q + ...:', b[1:16])
    ok4 = ok6 = True; n = 0
    for q in primerange(7, QMAX + 1):
        if q % 4 != 3:
            continue
        T = jacobi_traces(q, (4, 6))
        e4 = T[4] == -(q - 1) * (q + a[q])
        e6 = T[6] == (q - 1) * (q * q - q * a[q] - b[q])
        ok4 &= e4; ok6 &= e6; n += 1
        if not (e4 and e6):
            print('FAIL at', q, T, -(q - 1) * (q + a[q]), (q - 1) * (q * q - q * a[q] - b[q]))
    print(f'primes q = 3 mod 4 in [7, {QMAX}]: {n};  T_4 law holds at all: {ok4};  T_6 law holds at all: {ok6}')
