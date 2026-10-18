#!/usr/bin/env python3
"""Audit check (b'): THM-4710's recursion over the rational function field Q(N) (stable range N >= 9).

For N >= 9 every space Sym^N(V_2)_L, L <= 9, has the N-independent basis mu (partitions of L), and d, D_n are given by
exact formulas with coefficients in Z[N] (only d involves N, through the factor N - l(mu)).  So the whole recursion
can be run over K = Q(N).  This script
  1. computes S_L over K (L = 1..9, generator pattern {3,5,7,9}, conditions D_1 g = 0 and D_m g in M_(L-m),
     2 <= m <= L), and the generic dimensions;
  2. converts the generic G_L(N) (L = 5, 7, 9) to Cayley coordinates and compares each coefficient, as a rational
     function, with the closed forms of THM-4710 / witt_rigidity_interpolation_9_48.json;
  3. computes, for each L, an explicit nonzero maximal minor Delta_L(N) of the condition matrix (rows/columns chosen
     at a random large integer point), whose nonvanishing at an integer N0 >= 9 certifies dim S_L(N0) <= 1 there,
     and lists the integer roots N0 >= 1 of Delta_L (and of the denominators met) for exact rechecking.
Output: witt_symbolic.out / witt_symbolic.json.
"""
import json, os, sys, time
from fractions import Fraction
from math import comb
import sympy as sp
from sympy import QQ, ZZ
from sympy.polys.matrices import DomainMatrix

HERE = os.path.dirname(os.path.abspath(__file__))
Nsym = sp.Symbol('N')
K = QQ.frac_field(Nsym)
PATTERN = (3, 5, 7, 9)


def partitions(n, maxpart=None, minpart=1):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), minpart - 1, -1):
        for rest in partitions(n - k, k, minpart):
            yield (k,) + rest


def basis(L):
    return list(partitions(L)) if L > 0 else [()]


KN = K.from_sympy(Nsym)
K0 = K.zero
K1 = K.one


def kint(x):
    return K.from_sympy(sp.Rational(x))


def d_op(vec):
    out = {}
    for mu, c in vec.items():
        new = tuple(sorted(mu + (1,), reverse=True))
        out[new] = out.get(new, K0) + c * (KN - kint(len(mu)))
        for i, j in enumerate(mu):
            new = tuple(sorted(mu[:i] + (j + 1,) + mu[i + 1:], reverse=True))
            out[new] = out.get(new, K0) + c * kint(j + 1)
    return {k: v for k, v in out.items() if v != K0}


def D_op(vec, n):
    out = {}
    for mu, c in vec.items():
        for i, j in enumerate(mu):
            if j >= n:
                rest = mu[:i] + ((j - n,) if j - n > 0 else ()) + mu[i + 1:]
                new = tuple(sorted(rest, reverse=True))
                out[new] = out.get(new, K0) + c * kint(j + n + 2)
    return {k: v for k, v in out.items() if v != K0}


def to_row(d, idx):
    r = [K0] * len(idx)
    for k, v in d.items():
        r[idx[k]] = r[idx[k]] + v
    return r


def condition_matrix(L, M, conds):
    cols = basis(L)
    n = len(cols)
    blocks = []
    nb = 0
    for m in conds:
        blocks.append((m, basis(L - m), M[L - m], nb))
        nb += len(M[L - m])
    rows = []
    t1 = basis(L - 1)
    t1i = {m: i for i, m in enumerate(t1)}
    img1 = [to_row(D_op({mu: K1}, 1), t1i) for mu in cols]
    for r in range(len(t1)):
        rows.append([img1[j][r] for j in range(n)] + [K0] * nb)
    for (m, tb, Mv, off) in blocks:
        ti = {x: i for i, x in enumerate(tb)}
        img = [to_row(D_op({mu: K1}, m), ti) for mu in cols]
        Mvec = [to_row(v, ti) for v in Mv]
        for r in range(len(tb)):
            row = [img[j][r] for j in range(n)] + [K0] * nb
            for k, mv in enumerate(Mvec):
                row[n + off + k] = -mv[r]
            rows.append(row)
    return rows, n, nb


def span_basis(dicts, cols):
    if not dicts:
        return []
    idx = {m: i for i, m in enumerate(cols)}
    A = DomainMatrix([to_row(d, idx) for d in dicts], (len(dicts), len(cols)), K)
    R, piv = A.rref()
    Rl = R.to_Matrix().tolist() if False else R.rep.to_ddm() if False else None
    out = []
    Rr = R.to_list()
    for i in range(len(piv)):
        out.append({cols[j]: Rr[i][j] for j in range(len(cols)) if Rr[i][j] != K0})
    return out


def closure_symbolic(Lmax=9, conds_all=True):
    M = {0: [{(): K1}]}
    S = {0: [{(): K1}]}
    cond_mats = {}
    for L in range(1, Lmax + 1):
        conds = list(range(2, L + 1)) if conds_all else [2]
        rows, n, nb = condition_matrix(L, M, conds)
        A = DomainMatrix(rows, (len(rows), n + nb), K)
        NS = A.nullspace()           # rows span the nullspace
        nsl = NS.to_list()
        proj = [r[:n] for r in nsl]
        cols = basis(L)
        if proj:
            P_ = DomainMatrix(proj, (len(proj), n), K)
            R, piv = P_.rref()
            Rr = R.to_list()
            SL = [{cols[j]: Rr[i][j] for j in range(n) if Rr[i][j] != K0} for i in range(len(piv))]
        else:
            SL = []
        S[L] = SL
        cond_mats[L] = (rows, n, nb)
        dM = [v for v in (d_op(v) for v in M[L - 1]) if v]
        M[L] = span_basis(dM + (SL if L in PATTERN else []), cols)
    return S, M, cond_mats


def cayley_A(j):
    f = {}
    for i in range(0, j + 1):
        coef = Fraction(comb(j + 3, i)) * Fraction(-1, 4) ** i
        parts = [1] * i + ([j - i] if j - i > 0 else [])
        m = tuple(sorted(parts, reverse=True))
        f[m] = f.get(m, 0) + coef
    return {k: v for k, v in f.items() if v}


def poly_mul(f, g):
    out = {}
    for a, ca in f.items():
        for b, cb in g.items():
            m = tuple(sorted(a + b, reverse=True))
            out[m] = out.get(m, 0) + ca * cb
    return {k: v for k, v in out.items() if v}


def cayley_coords_symbolic(G, L):
    nus = list(partitions(L, minpart=2))
    polys = []
    for nu in nus:
        f = {(): Fraction(1)}
        for j in nu:
            f = poly_mul(f, cayley_A(j))
        polys.append(f)
    mus = basis(L)
    # augmented system over K
    top = G[(L,)]
    rows = [[kint(poly.get(mu, 0)) for poly in polys] + [G.get(mu, K0) / top] for mu in mus]
    A = DomainMatrix(rows, (len(mus), len(nus) + 1), K)
    R, piv = A.rref()
    assert len(nus) not in piv
    Rr = R.to_list()
    sol = {nus[pc]: Rr[i][-1] for i, pc in enumerate(piv)}
    return {nu: K.to_sympy(sol[nu]) for nu in nus}


def main():
    t0 = time.time()
    S, M, cond_mats = closure_symbolic()
    dims = {L: len(S[L]) for L in S}
    print(f'generic dims over Q(N): {dims}  ({time.time()-t0:.1f}s)', flush=True)
    fits = json.load(open(os.path.join(HERE, '..', '..', 'quantum_transfer_s44', 'witt_rigidity_interpolation_9_48.json')))['fits']
    formulas = {}
    for key, ent in fits.items():
        L = int(key.split()[0].split('_')[1])
        nu = tuple(int(x) for x in key.split(' ', 1)[1].strip('()').split(',') if x.strip())
        formulas[(L, nu)] = sp.sympify(ent['expr'], locals={'N': Nsym})
    out = {'generic_dims': dims, 'cayley': {}, 'match_closed_forms': {}}
    allmatch = True
    for L in (3, 5, 7, 9):
        cay = cayley_coords_symbolic(S[L][0], L)
        out['cayley'][L] = {str(k): str(sp.factor(v)) for k, v in cay.items()}
        for nu, v in cay.items():
            if nu == (L,):
                assert sp.simplify(v - 1) == 0
                continue
            f = formulas.get((L, nu))
            ok = f is not None and sp.simplify(sp.together(v - f)) == 0
            out['match_closed_forms'][f'G_{L} {nu}'] = ok
            allmatch &= ok
            print(f'G_{L} {list(nu)}: {sp.factor(v)}   matches THM-4710 closed form: {ok}', flush=True)
    out['all_closed_forms_match_identically'] = allmatch
    print(f'ALL CLOSED FORMS MATCH AS RATIONAL FUNCTIONS: {allmatch}', flush=True)
    # explicit maximal minors certifying generic uniqueness
    N0 = 1000003
    cert = {}
    for L in (3, 5, 7, 9):
        rows, n, nb = cond_mats[L]
        nvar = n + nb
        # numeric matrix at N0 to choose rows/columns of a nonsingular maximal minor
        num = [[sp.Rational(K.to_sympy(x).subs(Nsym, N0)) for x in r] for r in rows]
        Mnum = sp.Matrix(num)
        rnk = Mnum.rank()
        assert rnk == nvar - 1, (L, rnk, nvar)
        # choose pivot rows (via rref of transpose) and pivot columns
        _, pc = Mnum.rref()
        _, pr = Mnum.T.rref()
        sub = [[rows[i][j] for j in pc] for i in pr]
        # clear denominators row by row, determinant over ZZ[N]
        subsym = sp.Matrix([[K.to_sympy(x) for x in r] for r in sub])
        lcmden = [sp.lcm([sp.denom(sp.together(x)) for x in subsym.row(i)]) for i in range(subsym.rows)]
        polym = sp.Matrix([[sp.cancel(subsym[i, j] * lcmden[i]) for j in range(subsym.cols)] for i in range(subsym.rows)])
        Dm = DomainMatrix.from_Matrix(polym).convert_to(ZZ[Nsym])
        det = Dm.det()
        detexpr = sp.factor(ZZ[Nsym].to_sympy(det))
        roots = sorted(set(int(r) for r in sp.Poly(sp.expand(ZZ[Nsym].to_sympy(det)), Nsym).ground_roots().keys()
                           if r.is_integer))
        dens = sorted(set(int(r) for den in lcmden for r in sp.Poly(den, Nsym).ground_roots().keys() if r.is_integer)) \
            if any(sp.Poly(den, Nsym).degree() > 0 for den in lcmden) else []
        cert[L] = {'nvar': nvar, 'rank_at_N0': rnk, 'minor_size': len(pc), 'det_factored': str(detexpr),
                   'integer_roots_of_det': roots, 'integer_roots_of_row_denominators': dens}
        print(f'L={L}: maximal minor {len(pc)}x{len(pc)}; det = {detexpr}; integer roots {roots}; '
              f'row-denominator integer roots {dens}', flush=True)
    out['uniqueness_certificates'] = cert
    out['time_s'] = time.time() - t0
    json.dump(out, open(os.path.join(HERE, 'witt_symbolic.json'), 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
