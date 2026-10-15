#!/usr/bin/env python3
"""Audit (c) of S35: wheels.  T_n = Tr((DS)^n), D = diag chi(x), S_xy = chi(y - x).

1. Proposition T_n = sum_psi (chi(-1) J(psi, chi))^n: spot check (complex floats) against exact integer traces, both
   classes mod 4, odd and even n; U(W_n) = q T_n by brute force at small q (hub 0 -> rim, rim i -> i+1).
2. Level-8 eigenforms built here from theta products theta(z)^a theta(2z)^b (a + b = 2k, b even): these span
   M_{2k}(Gamma_0(8)) if the rank equals dim S_{2k}(Gamma_0(8)) + 4 (4 cusps).  T_3-eigenvectors give the systems
   weight 4: -4; weight 6: -12 (old, level 4), 20 (new); weight 8: 12 (old, level 2), 44, -84 (new).
3. T_4 = -(q-1)(q + a_q) and T_6 = (q-1)(q^2 - q a_q - b6_q) at all primes q = 3 mod 4 in [7, 400] (exact traces).
4. T_8: is T_8/(q-1) in the Z-span of q^0..q^4 and the Tate twists of ALL level-2/4/8 eigen-systems of weights <= 8
   (15 functions, including the a_3 = -84 newform that the session's fit basis lacked), at 39 primes?
"""
import cmath, math, time
import numpy as np
import sympy as sp

QMAX = 400


def chi_table(q):
    c = np.full(q, -1, dtype=np.int64)
    c[0] = 0
    for x in range(1, q):
        c[(x * x) % q] = 1
    return c


def DS(q):
    chi = chi_table(q)
    x = np.arange(q)
    S = chi[(x[None, :] - x[:, None]) % q]
    return chi[:, None] * S


def traces(q):
    A = DS(q)
    M2 = A @ A
    M4 = M2 @ M2
    T = {}
    T[2] = int(np.sum(A * A.T))
    T[3] = int(np.sum(M2 * A.T))
    T[4] = int(np.sum(M2 * M2.T))
    T[5] = int(np.sum(M4 * A.T))
    T[6] = int(np.sum(M4 * M2.T))
    T[8] = int(np.sum(M4 * M4.T))
    return T


def jacobi_powers(q, ns):
    g = sp.primitive_root(q)
    dlog = {}
    x = 1
    for k in range(q - 1):
        dlog[x] = k
        x = x * g % q
    chi = chi_table(q)
    cm1 = int(chi[q - 1])
    out = {n: 0j for n in ns}
    for j in range(q - 1):
        J = 0j
        for t in range(1, q):
            c = int(chi[(1 - t) % q])
            if c:
                J += cmath.exp(2j * math.pi * j * dlog[t] / (q - 1)) * c
        for n in ns:
            out[n] += (cm1 * J) ** n
    return {n: out[n] for n in ns}


def U_wheel_brute(q, n):
    chi = chi_table(q)
    tot = 0
    for xs in np.ndindex(*([q] * (n + 1))):
        h = xs[0]
        p = 1
        for i in range(1, n + 1):
            p *= chi[(xs[i] - h) % q] * chi[(xs[i % n + 1] - xs[i]) % q]
            if p == 0:
                break
        tot += p
    return tot


def theta_terms(d, N):
    terms, m = [], 0
    while d * m * m < N:
        terms.append((d * m * m, 1 if m == 0 else 2))
        m += 1
    return terms


def sparse_mul(A, terms, N):
    out = np.zeros(N, dtype=object)
    for e, c in terms:
        if e < N:
            out[e:] += c * A[:N - e]
    return out


def level8_systems(k2, N):
    """eigen-systems (T_3 eigenvalue, a_n list) of S_{k2}(Gamma_0(8)) from theta products; returns dict ev -> coeffs."""
    forms = []
    for b in range(0, k2 * 2 + 1, 2):
        a = 2 * k2 - b
        A = np.zeros(N, dtype=object)
        A[0] = 1
        for _ in range(a):
            A = sparse_mul(A, theta_terms(1, N), N)
        for _ in range(b):
            A = sparse_mul(A, theta_terms(2, N), N)
        forms.append([int(v) for v in A])
    B = 120
    rows = sp.Matrix([f[:B] for f in forms])
    rank = rows.rank()
    T = []
    for f in forms:
        v = [f[3 * n] + (3 ** (k2 - 1) * f[n // 3] if n % 3 == 0 else 0) for n in range(B)]
        sol, params = rows.T.gauss_jordan_solve(sp.Matrix(v))
        assert params.shape[0] == 0
        T.append(list(sol))
    T3 = sp.Matrix(T)
    lam = sp.symbols('lam')
    cp = sp.factor(T3.charpoly(lam).as_expr())
    systems = {}
    for ev, mult, vecs in T3.T.eigenvects():
        if ev == 1 + 3 ** (k2 - 1):
            continue                      # Eisenstein
        for vv in vecs:
            f = [sum(vv[i] * forms[i][n] for i in range(len(forms))) for n in range(N)]
            if f[1] != 0:
                systems[int(ev)] = [sp.Rational(c) / f[1] for c in f]
                break
    return rank, cp, systems


def main():
    t0 = time.time()
    # 1. proposition spot checks
    for q in (7, 11, 13, 17, 19):
        T = traces(q)
        J = jacobi_powers(q, (2, 3, 4, 5, 6, 8))
        ok = all(abs(J[n] - T[n]) < 1e-6 * max(1, abs(T[n])) for n in (2, 3, 4, 5, 6, 8))
        print(f'q = {q:2d}: exact traces T_2..T_8 = {[T[n] for n in (2, 3, 4, 5, 6, 8)]}; Jacobi power sums agree: {ok}', flush=True)
    for q, n in ((5, 4), (7, 4), (5, 5), (7, 6), (11, 4)):
        Ub = U_wheel_brute(q, n)
        T = traces(q)
        Tn = T[n] if n in T else None
        print(f'U(W_{n}) brute force at q = {q}: {Ub};  q * T_{n} = {q * Tn}; equal: {Ub == q * Tn}', flush=True)
    # 2. eigenforms
    N = QMAX + 2
    sysw = {}
    for k2 in (4, 6, 8):
        rank, cp, systems = level8_systems(k2, N)
        dimS = {4: 1, 6: 3, 8: 5}[k2]
        print(f'weight {k2}: theta-product span rank {rank} (dim M_{k2}(Gamma_0(8)) = {dimS + 4}); T_3 charpoly {cp}; '
              f'cusp systems a_3 -> a_5, a_7: {[(ev, int(f[5]), int(f[7])) for ev, f in systems.items()]}', flush=True)
        sysw[k2] = systems
    a4 = sysw[4][-4]
    b6 = sysw[6][20]
    print('b6 = q + ...:', [int(c) for c in b6[1:16]])
    # 3. T_4, T_6 laws and T_8 data
    primes = [q for q in sp.primerange(7, QMAX + 1) if q % 4 == 3]
    data = {}
    ok4 = ok6 = True
    for q in primes:
        T = traces(q)
        data[q] = T
        ok4 &= T[4] == -(q - 1) * (q + int(a4[q]))
        ok6 &= T[6] == (q - 1) * (q * q - q * int(a4[q]) - int(b6[q]))
    print(f'{len(primes)} primes q = 3 mod 4 in [7, {QMAX}]: T_4 law: {ok4}; T_6 law: {ok6}  [{time.time() - t0:.0f} s]', flush=True)
    # 4. T_8 span test
    names, cols = [], []
    for j in range(5):
        names.append(f'q^{j}')
        cols.append([q ** j for q in primes])
    for k2, systems in sysw.items():
        for ev, f in systems.items():
            for e in range((8 - k2) // 2 + 1):
                names.append(f'q^{e} f[w{k2}, a3={ev}]')
                cols.append([q ** e * f[q] for q in primes])
    A = sp.Matrix(cols).T
    bvec = sp.Matrix([sp.Rational(data[q][8], q - 1) for q in primes])
    rA = A.rank()
    rAb = A.row_join(bvec).rank()
    print(f'T_8/(q-1) at {len(primes)} primes vs {len(names)} basis functions {names}: rank {rA}, augmented rank {rAb} '
          f'-> {"in the span" if rA == rAb else "NOT in the span (no exact fit)"}')
    # also the session's 13-prime, 14-function setting: underdetermined?
    sub = [i for i, q in enumerate(primes) if q <= 103]
    names14 = [nm for nm in names if 'a3=-84' not in nm]
    cols14 = [c for nm, c in zip(names, cols) if 'a3=-84' not in nm]
    A14 = sp.Matrix([[c[i] for c in cols14] for i in sub])
    b14 = sp.Matrix([bvec[i] for i in sub])
    print(f'session setting (primes <= 103: {len(sub)} equations; {len(names14)} functions without a3 = -84): rank {A14.rank()}, '
          f'augmented rank {A14.row_join(b14).rank()} -> an exact fit {"exists (no information)" if A14.rank() == A14.row_join(b14).rank() else "does not exist"}')


if __name__ == '__main__':
    main()
