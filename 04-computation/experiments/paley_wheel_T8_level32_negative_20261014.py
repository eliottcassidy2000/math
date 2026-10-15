#!/usr/bin/env python3
"""T_8 = Tr((DS)^8) (the wheel W_8; Greene's 8F7(1) up to normalisation) is NOT in the span of trivial-character
modular forms of level dividing 32 and weight <= 8 with Tate twists (FINITE-EXACT negative, overdetermined test).

Data: T_8 at all 39 primes q = 3 (mod 4) in [7, 400], from the Jacobi-sum formula T_n = sum_psi (chi(-1)J(psi,chi))^n
(paley_wheel_jacobi_20261014.py).  Basis functions (21):
  q^3, q^2, q, 1;  q^2 a, q a, a  (a = eta(2z)^4 eta(4z)^4);  q b6, b6 (weight-6 level-8 newform);  q e6, e6
  (e6 = eta(2z)^12, weight 6 level 4);  and the 24 coefficient functions q -> B_i[q] of a basis of S_8(Gamma_0(32))
  (dimension 24; at odd primes these span all newform eigenvalues of level dividing 32, every Galois embedding).
S_8(Gamma_0(32)) is built in full from products: a(dz) x weight-4 forms, b6(dz) x weight-2 forms, and the four
weight-3 level-32 chi_-4 eta quotients x theta powers.  The linear system (39 equations, 35 unknowns) is solved
exactly; inconsistency = T_8/(q - 1) is not in the span.  (The S35 audit did the same at levels 2/4/8 with 15 functions.)
Session opus-2026-10-14-S35.
"""
import importlib.util, sys, os
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + '.py'))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
W = _load('paley_wheel_ladder_fit_20261014')
J = _load('paley_wheel_jacobi_20261014')

M = 410
W.M = M
mul, eigenforms, eta_q = W.mul, W.eigenforms, W.eta_q


def sigma(n, k):
    return sum(d ** k for d in range(1, n + 1) if n % d == 0)


def scale(G, d):
    return [G[n // d] if n % d == 0 else 0 for n in range(M)]


def pw(G, k):
    R = [1] + [0] * (M - 1)
    for _ in range(k):
        R = mul(R, G)
    return R


def main():
    qs = [q for q in sp.primerange(7, 401) if q % 4 == 3]
    T8 = {q: J.jacobi_traces(q, (8,))[8] for q in qs}
    print(f'T_8 at {len(qs)} primes (Jacobi-sum formula); e.g. T_8(7) = {T8[7]}, T_8(103) = {T8[103]}', flush=True)
    E4 = [1] + [240 * sigma(n, 3) for n in range(1, M)]
    E2 = [1] + [-24 * sigma(n, 1) for n in range(1, M)]
    th = [0] * M
    for n in range(-25, 26):
        if n * n < M:
            th[n * n] += 1
    a4 = [f for (_, _, ev, f) in eigenforms(8, 4)][0]
    b6 = [f for (_, _, ev, f) in eigenforms(8, 6) if ev == 20][0]
    e6 = [f for (_, _, ev, f) in eigenforms(4, 6)][0]
    w3 = [eta_q(rr) for rr in [{1: -4, 2: 10, 4: -2, 8: 2}, {1: -2, 2: 5, 4: -2, 8: 7, 16: -2},
                                {2: -2, 4: 7, 8: -2, 16: 5, 32: -2}, {2: -2, 4: 7, 16: -1, 32: 2}]]
    th10, th2 = pw(th, 10), pw(th, 2)
    w4 = [scale(E4, d) for d in (1, 2, 4, 8, 16, 32)] + [scale(pw(th, 8), d) for d in (1, 2, 4)] + \
         [mul(pw(th, 4), scale(pw(th, 4), d)) for d in (2, 4, 8)]
    w2 = [[x - d * y for x, y in zip(E2, scale(E2, d))] for d in (2, 4, 8, 16, 32)] + [scale(pw(th, 4), d) for d in (1, 2, 4, 8)]
    gens = []
    for f in (a4, scale(a4, 2), scale(a4, 4)):
        gens += [mul(f, g) for g in w4]
    for f in (b6, scale(b6, 2), scale(b6, 4)):
        gens += [mul(f, g) for g in w2]
    for f in w3:
        gens += [mul(f, th10), mul(f, scale(th10, 2)), mul(mul(f, th2), pw(th, 8))]
    Bcols = 120
    cur = sp.zeros(0, Bcols); B = []
    for G in gens:
        if G[0] != 0:
            continue
        test = cur.col_join(sp.Matrix([G[:Bcols]]))
        if test.rank() > cur.rows:
            cur = test; B.append(G)
        if len(B) == 24:
            break
    print('span dimension (dim S_8(Gamma_0(32)) = 24):', len(B), flush=True)
    # Every form in S_8(Gamma_0(32)) is a rational combination of the basis B; at odd primes q its coefficient is a
    # combination of the eigenvalues a_q(f) of all newforms f of level dividing 32 (every Galois embedding).  So the
    # most general test uses the 24 basis coefficient functions q -> B_i[q] themselves (exact integers).
    cand = {'q^3': [q ** 3 for q in qs], 'q^2': [q ** 2 for q in qs], 'q': list(qs), '1': [1] * len(qs),
            'q^2 a': [q * q * a4[q] for q in qs], 'q a': [q * a4[q] for q in qs], 'a': [a4[q] for q in qs],
            'q b6': [q * b6[q] for q in qs], 'b6': [b6[q] for q in qs], 'q e6': [q * e6[q] for q in qs], 'e6': [e6[q] for q in qs]}
    for i, G in enumerate(B):
        cand[f'B{i}'] = [G[q] for q in qs]
    names = list(cand)
    A = sp.Matrix([cand[n] for n in names]).T
    b = sp.Matrix([sp.Rational(T8[q], q - 1) for q in qs])
    print(f'system: {A.rows} equations (primes), {A.cols} unknowns; rank(A) = {A.rank()}, rank([A|b]) = {A.row_join(b).rank()}')
    try:
        x, params = A.gauss_jordan_solve(b)
        print('CONSISTENT (unexpected): ', dict(zip(names, list(x))))
    except ValueError:
        print(f'INCONSISTENT: T_8/(q-1) is not in the span of these {A.cols} functions (all of S_8(Gamma_0(32)) at odd primes, plus Tate-twisted lower weights).')


if __name__ == '__main__':
    main()
