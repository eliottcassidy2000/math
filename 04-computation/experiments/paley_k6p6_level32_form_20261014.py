#!/usr/bin/env python3
"""The K_6 - P_6 terminal is the symmetric square of a level-32 weight-3 newform.

Claim (FINITE-EXACT, all primes 3 <= q <= 691):
    Y_q := U(K_6 - P_6)/(q(q-1)) = |a_q(f)|^2 - q^2 + 2 chi_{-4}(q) q,
where f is the newform in S_3(Gamma_0(32), chi_{-4}) with Hecke field Q(i) (T_3-eigenvalue a_3 = 4i; its Galois
conjugate has a_3 = -4i).  Equivalently Y_q = chi_{-4}(q) (tr Sym^2 rho_f(Frob_q) + 2q), and
    Y_q - 1 = (|a_q| - (q - chi_{-4}(q))) (|a_q| + (q - chi_{-4}(q))).

Construction of f: the weight-3 eta quotients of level 32 with character chi_{-4} (Newman/Ligozat conditions,
cuspidal at all cusps, exponents |r| <= 10) span a 4-dimensional Hecke-stable space; T_3 has characteristic polynomial
x^2 (x^2 + 16).  The eigenvalue-0 part is old (eta(4z)^6 and eta(8z)^6, CM by Q(i)); the x^2 + 16 part is the newform
pair.  The eigenvector is computed exactly over Q(i).
Input: Y_q from paley_k6p6_fast_table_20261014.out.  Session opus-2026-10-14-S35.
"""
import itertools, re, sys
from math import gcd
import sympy as sp

M = 700
OUT = sys.argv[1] if len(sys.argv) > 1 else '05-knowledge/results/paley_k6p6_fast_table_20261014.out'


def pentagonal(delta, M):
    E = [0] * M
    k = 0
    while True:
        hit = False
        for kk in ((k, -k) if k else (0,)):
            e = kk * (3 * kk - 1) // 2
            if e * delta < M:
                E[e * delta] += (-1) ** abs(kk); hit = True
        if not hit and k > 2:
            break
        k += 1
    return E


def mul(A, B):
    out = [0] * M
    for i, a in enumerate(A):
        if a:
            for j in range(M - i):
                if B[j]:
                    out[i + j] += a * B[j]
    return out


def inv(A):
    out = [0] * M; out[0] = 1
    for n in range(1, M):
        out[n] = -sum(A[k] * out[n - k] for k in range(1, n + 1) if A[k])
    return out


def eta_quotient(r):
    lead = sum(d * x for d, x in r.items()) // 24
    S = [1] + [0] * (M - 1)
    for d, x in r.items():
        if x == 0:
            continue
        E = pentagonal(d, M)
        if x < 0:
            E = inv(E); x = -x
        for _ in range(x):
            S = mul(S, E)
    return [0] * lead + S[:M - lead]


def quotients(N=32, RMAX=10):
    divs = [d for d in range(1, N + 1) if N % d == 0]
    for r in itertools.product(range(-RMAX, RMAX + 1), repeat=len(divs)):
        if sum(r) != 6:
            continue
        rr = dict(zip(divs, r))
        if sum(d * x for d, x in rr.items()) % 24 or sum((N // d) * x for d, x in rr.items()) % 24:
            continue
        orders = [sum(gcd(c, d) ** 2 * x * N // (gcd(c, N // c) * c * d) for d, x in rr.items()) for c in divs]
        if min(orders) <= 0:
            continue
        e2 = sum(x * {1: 0, 2: 1, 4: 2, 8: 3, 16: 4, 32: 5}[d] for d, x in rr.items())
        if e2 % 2 == 0:  # character chi_{-4}
            yield rr


def main():
    Q = list(quotients())
    B = 80
    basis, cur = [], sp.zeros(0, B)
    for rr in Q:
        F = eta_quotient(rr)
        test = cur.col_join(sp.Matrix([F[:B]]))
        if test.rank() > cur.rows:
            cur = test; basis.append((rr, F))
    print(f'{len(Q)} eta quotients of level 32, weight 3, character chi_-4; span dimension {len(basis)}')
    for rr, F in basis:
        print('  basis element', {d: x for d, x in rr.items() if x}, ' q-expansion', F[:12])
    # T_3 on the span
    Bp = 60
    Bm = sp.Matrix([F[1:Bp + 1] for _, F in basis])
    T = []
    for _, F in basis:
        v = sp.Matrix([[F[3 * n] + (-9 * F[n // 3] if n % 3 == 0 else 0) for n in range(1, Bp + 1)]])
        sol, params = Bm.T.gauss_jordan_solve(v.T)
        assert params.shape[0] == 0
        T.append(list(sol))
    T3 = sp.Matrix(T)
    print('T_3 characteristic polynomial:', sp.factor(T3.charpoly().as_expr()))
    # eigenvector for 4i (row-vector convention: f = x . basis, x T3 = 4i x)
    ns = (T3.T - 4 * sp.I * sp.eye(T3.rows)).nullspace()
    assert len(ns) == 1
    x = ns[0]
    f = [sp.nsimplify(sp.expand(sum(x[k] * basis[k][1][n] for k in range(len(basis))))) for n in range(M)]
    f = [sp.simplify(c / f[1]) for c in f]
    print('newform f = q + ... coefficients a_1..a_20:', [sp.nsimplify(c) for c in f[1:21]])
    # eigenform sanity: a_9 = a_3^2 - chi(3) 3^2 a_1, a_15 = a_3 a_5
    print('sanity a_9 = a_3^2 + 9:', sp.simplify(f[9] - (f[3] ** 2 + 9)) == 0, ';  a_15 = a_3 a_5:', sp.simplify(f[15] - f[3] * f[5]) == 0)
    # compare with Y_q
    Y = {}
    for l in open(OUT):
        m = re.match(r'q =\s*(\d+)\s+q mod 16 =\s*\d+\s+Y_q =\s*(-?\d+)', l)
        if m:
            Y[int(m.group(1))] = int(m.group(2))
    bad = []
    for q, y in sorted(Y.items()):
        if q >= M:
            continue
        a = sp.expand(f[q])
        norm = sp.simplify(a * sp.conjugate(a))
        chi = 1 if q % 4 == 1 else -1
        pred = norm - q * q + 2 * chi * q
        if pred != y:
            bad.append((q, y, pred))
    print(f'Y_q = |a_q(f)|^2 - q^2 + 2 chi_-4(q) q at {len([q for q in Y if q < M])} primes 3..691; failures: {bad}')
    print('a_q real for q = 1 mod 4 and purely imaginary for q = 3 mod 4:',
          all((sp.im(f[q]) == 0) if q % 4 == 1 else (sp.re(f[q]) == 0) for q in Y if q < M))
    print('primes q < 700 with a_q(f) = 0:', [q for q in Y if q < M and f[q] == 0])
    print('primes with |a_q| = q + 1 (Y_q = 1):', [q for q in Y if q < M and q % 4 == 3 and sp.simplify(f[q] * sp.conjugate(f[q])) == (q + 1) ** 2])


if __name__ == '__main__':
    main()
