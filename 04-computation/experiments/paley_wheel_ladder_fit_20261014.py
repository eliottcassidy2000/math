#!/usr/bin/env python3
"""Exact identification of the wheel traces T_n = Tr((DS)^n), n = 4, 6, 8 (q = 3 mod 4), by modular forms of level 8.

Input: 05-knowledge/results/paley_wheel_ladder_20261014.out (primes 7..103, q = 3 mod 4: 13 data points).
Basis: q^j and Tate twists q^((n-w)/2) f_q of the rational Hecke eigenforms f of weight w <= n in the spans of
holomorphic cuspidal eta quotients of levels 2, 4, 8 (trivial character).  The system is solved exactly and the
solution must reproduce every data point (more equations than unknowns).
Results (FINITE-EXACT, 13 primes):
  T_4 = -(q - 1)(q + a_q),                 a = eta(2z)^4 eta(4z)^4 (weight 4, level 8)
  T_6 = (q - 1)(q^2 - q a_q - b_q),        b = the weight-6 level-8 newform q + 20q^3 - 74q^5 - 24q^7 + ...
  T_8 = see output.
Session opus-2026-10-14-S35.
"""
import itertools, re, sys
from math import gcd
import sympy as sp

M = 110
DATA = sys.argv[1] if len(sys.argv) > 1 else '05-knowledge/results/paley_wheel_ladder_20261014.out'


def pentagonal(delta):
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


def eta_q(r):
    lead = sum(d * x for d, x in r.items()) // 24
    S = [1] + [0] * (M - 1)
    for d, x in r.items():
        if x == 0:
            continue
        E = pentagonal(d)
        if x < 0:
            E = inv(E); x = -x
        for _ in range(x):
            S = mul(S, E)
    return [0] * lead + S[:M - lead]


def quotients(N, k, RMAX):
    divs = [d for d in range(1, N + 1) if N % d == 0]
    for r in itertools.product(range(-RMAX, RMAX + 1), repeat=len(divs)):
        if sum(r) != 2 * k:
            continue
        rr = dict(zip(divs, r))
        if sum(d * x for d, x in rr.items()) % 24 or sum((N // d) * x for d, x in rr.items()) % 24:
            continue
        orders = [sum(gcd(c, d) ** 2 * x * N // (gcd(c, N // c) * c * d) for d, x in rr.items()) for c in divs]
        if min(orders) <= 0:
            continue
        if sum(x * (d.bit_length() - 1) for d, x in rr.items()) % 2 == 0:
            yield rr


def eigenforms(N, k):
    cur = sp.zeros(0, M); B = []
    for rr in quotients(N, k, 12 if N <= 4 else 9):
        F = eta_q(rr)
        test = cur.col_join(sp.Matrix([F]))
        if test.rank() > cur.rows:
            cur = test; B.append(F)
    if not B:
        return []
    Bp = 30
    Bm = sp.Matrix([F[1:Bp + 1] for F in B])
    T = []
    for F in B:
        v = sp.Matrix([[F[3 * n] + (3 ** (k - 1) * F[n // 3] if n % 3 == 0 else 0) for n in range(1, Bp + 1)]])
        sol, params = Bm.T.gauss_jordan_solve(v.T)
        T.append(list(sol))
    T3 = sp.Matrix(T)
    out = []
    for ev, mult, vecs in T3.T.eigenvects():
        if not ev.is_rational:
            continue
        for v in vecs:
            f = [sum(v[i] * B[i][n] for i in range(len(B))) for n in range(M)]
            if f[1] == 0:
                continue
            f = [c / f[1] for c in f]
            if any(f[n] != 0 for n in range(2, M, 2)) and N > 2:
                # keep only forms with a_2 = 0 normalisation-independent check is not needed; keep all
                pass
            out.append((N, k, ev, f))
    return out


def main():
    T = {}
    for l in open(DATA):
        m = re.match(r'q =\s*(\d+) \((\d) mod 4\): (.*)$', l.strip())
        if m:
            T[int(m.group(1))] = dict((int(a), int(b)) for a, b in re.findall(r'T(\d) = (-?\d+)', m.group(3)))
    qs = [q for q in sorted(T) if q % 4 == 3]
    forms = []
    seen = set()
    for k in (4, 6, 8):
        for N in (2, 4, 8):
            for (N_, k_, ev, f) in eigenforms(N, k):
                key = (k_, tuple(f[p] for p in (3, 5, 7, 11)))
                if key in seen:
                    continue
                seen.add(key); forms.append((N_, k_, ev, f))
                print(f'weight {k_}, level {N_}: rational eigenform a_3 = {ev}; a_p (p = 3,5,7,11,13) = {[f[p] for p in (3, 5, 7, 11, 13)]}', flush=True)
    for n in (4, 6, 8):
        names, cols = [], []
        for j in range(n // 2 + 1):
            names.append(f'q^{j}'); cols.append([q ** j for q in qs])
        for (N_, k_, ev, f) in forms:
            if k_ <= n:
                tw = (n - k_) // 2
                for e in range(tw + 1):
                    names.append(f'q^{e}*f[w{k_},N{N_},a3={ev}]'); cols.append([q ** e * f[q] for q in qs])
        A = sp.Matrix(cols).T
        b = sp.Matrix([sp.Rational(T[q][n], q - 1) for q in qs])
        print(f'\nT_{n}/(q - 1): {len(qs)} primes, {len(names)} basis functions')
        try:
            x, params = A.gauss_jordan_solve(b)
            if params.shape[0]:
                x = x.subs({p: 0 for p in params})
                print('  (underdetermined: free parameters set to 0)')
            terms = [f'{c}*{nm}' for c, nm in zip(list(x), names) if c != 0]
            print('  exact fit:', ' + '.join(terms))
        except ValueError:
            print('  no exact fit in this basis')


if __name__ == '__main__':
    main()
