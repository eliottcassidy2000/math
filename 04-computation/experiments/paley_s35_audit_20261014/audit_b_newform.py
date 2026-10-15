#!/usr/bin/env python3
"""Audit (b1) of opus-2026-10-14-S35: the weight-3 level-32 newform b with character chi_{-4}, rebuilt independently.

Construction (different from the session's eta quotients): theta products
    theta(z)^a theta(2z)^b theta(4z)^c theta(8z)^d,  a+b+c+d = 6,  b+d even,
i.e. theta series of diagonal senary forms of level dividing 32 with character (-prod d_i / .) = chi_{-4}.
Cohen-Oesterle gives dim M_3(Gamma_0(32), chi_{-4}) = 8 (dim S_3 = 4, Eisenstein 4); if the theta span has rank 8 it
is the whole space, hence Hecke-stable.  T_3 (a(3n) - 9 a(n/3)) is computed exactly; the 4i-eigenvector is the newform.
Then: the eta-theta closed form of audit C (HYP-9276) is checked against it far beyond the Sturm bound (12), and used
to tabulate b(n) for n <= 524,288 by exact sparse convolutions (int64).  Checks on the big table: Hecke relations,
Ramanujan bound, real/imaginary pattern.
Outputs b(q) for all primes q <= 524,288 to audit_b_newform_bq.json (q -> [Re, Im]).
Audit of S35, 2026-10-14.  Independent code.
"""
import json, time, itertools
import numpy as np
import sympy as sp

NH = 1500          # coefficients for the Hecke computation
NBIG = 524289      # big table


def theta(d, N):
    t = np.zeros(N, dtype=np.int64)
    m = 0
    while d * m * m < N:
        t[d * m * m] += 1 if m == 0 else 2
        m += 1
    return t


def sparse_mul(A, sparse_terms, N):
    """A (dense int64, length N) times sum c x^e over sparse_terms [(e, c)], truncated to N."""
    out = np.zeros(N, dtype=np.int64)
    for e, c in sparse_terms:
        if e >= N:
            continue
        out[e:] += c * A[:N - e]
    return out


def theta_terms(d, N):
    terms = []
    m = 0
    while d * m * m < N:
        terms.append((d * m * m, 1 if m == 0 else 2))
        m += 1
    return terms


def theta_product(exps, N):
    A = np.zeros(N, dtype=np.int64)
    A[0] = 1
    for d, e in zip((1, 2, 4, 8), exps):
        for _ in range(e):
            A = sparse_mul(A, theta_terms(d, N), N)
    return A


def cohen_oesterle_dim_M3():
    # k = 3, N = 32 = 2^5, chi conductor 4 = 2^2: (k-1)/12 * index = 2/12 * 48 = 8; lambda(r=5, s=2, p=2) = 2 * 2^2 = 8
    # dim S_3 - dim M_{-1} = 8 - 8/2 + 0 + 0 = 4 (no solutions of x^2+1 or x^2+x+1 mod 32);
    # Eisenstein part (k >= 3) = prod lambda = 8.  So dim M_3 = 12.
    return 4, 8


def pentagonal_terms(scale, N):
    """prod_{m>=1} (1 - x^{scale m}) as sparse terms up to N."""
    terms = []
    k = 0
    while True:
        added = False
        for kk in ((0,) if k == 0 else (k, -k)):
            e = kk * (3 * kk - 1) // 2 * scale
            if e < N:
                terms.append((e, -1 if kk % 2 else 1))
                added = True
        if not added and k > 0:
            break
        k += 1
    return terms


def closed_form_auditC(N):
    """g = eta(4z)^2 eta(8z)^2 (theta(4z)^2 + i (theta(2z)^2 - theta(4z)^2)); returns (Re, Im) int64 arrays."""
    G = np.zeros(N, dtype=np.int64)
    G[1] = 1                                   # q^{(2*4 + 2*8)/24} = q
    for scale in (4, 4, 8, 8):
        G = sparse_mul(G, pentagonal_terms(scale, N), N)
    t4 = theta_terms(4, N)
    t2 = theta_terms(2, N)
    Re = sparse_mul(sparse_mul(G, t4, N), t4, N)
    Im = sparse_mul(sparse_mul(G, t2, N), t2, N) - Re
    return Re, Im


def main():
    t0 = time.time()
    dS, dE = cohen_oesterle_dim_M3()
    combos = [e for e in itertools.product(range(7), repeat=4) if sum(e) == 6 and (e[1] + e[3]) % 2 == 0]
    forms = {e: theta_product(e, NH) for e in combos}
    # greedy basis
    basis, rows = [], []
    for e in combos:
        test = sp.Matrix(rows + [list(forms[e][:200])])
        if test.rank() > len(rows):
            rows.append(list(forms[e][:200]))
            basis.append(e)
    print(f'{len(combos)} theta products (a,b,c,d) with b+d even; span rank {len(basis)}; Cohen-Oesterle dim M_3 = {dS + dE} '
          f'(S_3: {dS}, Eisenstein: {dE});  basis exponents: {basis}', flush=True)
    assert len(basis) == dS + dE
    nb = len(basis)
    B = 400
    Bm = sp.Matrix([[int(forms[e][n]) for n in range(0, B)] for e in basis])     # include constant term
    T = []
    for e in basis:
        f = forms[e]
        v = [int(f[3 * n]) - (9 * int(f[n // 3]) if n % 3 == 0 else 0) for n in range(0, B)]
        sol, params = Bm.T.gauss_jordan_solve(sp.Matrix(v))
        assert params.shape[0] == 0
        T.append(list(sol))
    T3 = sp.Matrix(T)       # row convention: T_3 f_i = sum_j T3[i, j] f_j
    lam = sp.symbols('lam')
    print('T_3 characteristic polynomial on M_3(Gamma_0(32), chi_-4):', sp.factor(T3.charpoly(lam).as_expr(), extension=sp.I), flush=True)
    ns = (T3.T - 4 * sp.I * sp.eye(nb)).nullspace()
    assert len(ns) == 1
    x = ns[0]
    coeffs = [sp.expand(sum(x[j] * int(forms[basis[j]][n]) for j in range(nb))) for n in range(B)]
    c1 = coeffs[1]
    bq = [sp.nsimplify(sp.simplify(c / c1)) for c in coeffs]
    assert bq[0] == 0
    print('newform (4i-eigenvector, a_1 = 1): b(1..20) =', bq[1:21], flush=True)
    # closed form of audit C against the eigenvector, n < B (Sturm bound 12)
    Re, Im = closed_form_auditC(B)
    ok = all(sp.simplify(bq[n] - (int(Re[n]) + sp.I * int(Im[n]))) == 0 for n in range(B))
    print(f'eta-theta closed form (audit C) = theta-product eigenvector for all n < {B}: {ok}  (Sturm bound 12)', flush=True)
    assert ok
    # big table
    Re, Im = closed_form_auditC(NBIG)
    print(f'big table to n = {NBIG - 1} built [{time.time() - t0:.0f} s]', flush=True)
    primes = list(sp.primerange(3, NBIG))
    # pattern and Ramanujan
    bad_pattern = [p for p in primes if (p % 4 == 1 and Im[p] != 0) or (p % 4 == 3 and Re[p] != 0)]
    bad_ram = [p for p in primes if int(Re[p]) ** 2 + int(Im[p]) ** 2 > 4 * p * p]
    print(f'real for p = 1 mod 4 / imaginary for p = 3 mod 4 fails at: {bad_pattern[:5]}; Ramanujan |b(p)| <= 2p fails at: {bad_ram[:5]}', flush=True)
    # Hecke relations: b(p^2) = b(p)^2 - chi(p) p^2 ; b(mn) = b(m) b(n) for coprime m, n (sample)

    def bc(n):
        return complex(int(Re[n]), int(Im[n]))

    bad_h = []
    for p in primes:
        if p * p >= NBIG:
            break
        chi = 1 if p % 4 == 1 else -1
        lhs = (int(Re[p * p]), int(Im[p * p]))
        br, bi = int(Re[p]), int(Im[p])
        rhs = (br * br - bi * bi - chi * p * p, 2 * br * bi)
        if lhs != rhs:
            bad_h.append(p)
    rng = np.random.default_rng(1)
    bad_m = 0
    for _ in range(20000):
        m = int(rng.integers(3, 2000)) | 1
        n = int(rng.integers(3, 260)) | 1
        if sp.gcd(m, n) != 1 or m * n >= NBIG:
            continue
        a1, b1, a2, b2 = int(Re[m]), int(Im[m]), int(Re[n]), int(Im[n])
        if (int(Re[m * n]), int(Im[m * n])) != (a1 * a2 - b1 * b2, a1 * b2 + a2 * b1):
            bad_m += 1
    print(f'Hecke b(p^2) = b(p)^2 - chi(p) p^2 failures (p < 724): {bad_h}; multiplicativity failures in 20000 random coprime pairs: {bad_m}', flush=True)
    # even coefficients: level 32 => b(2) = 0
    print('b(2), b(4), b(6) =', (int(Re[2]), int(Im[2])), (int(Re[4]), int(Im[4])), (int(Re[6]), int(Im[6])))
    json.dump({str(p): [int(Re[p]), int(Im[p])] for p in primes}, open('audit_b_newform_bq.json', 'w'))
    print(f'wrote b(p) for {len(primes)} odd primes < {NBIG}  [{time.time() - t0:.0f} s]')


if __name__ == '__main__':
    main()
