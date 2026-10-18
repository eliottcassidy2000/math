#!/usr/bin/env python3
"""Degree-N generators of the singular support as sl_2 lowest-weight vectors (Cayley semi-invariants), S44.

Setting.  gr_F V^c = C[x_0, x_1, ...], x_j <-> L_(-j-2).  The global conformal sl_2 = <L_(-1), L_0, L_1> acts on the
Li-graded algebra by derivations preserving the polynomial degree: d x_j = (j+1) x_(j+1) and D_1 x_j = (j+3) x_(j-1)
(D_1 x_0 = 0), [D_1, d] = 2 D_0 with D_0 x_j = (j+2) x_j.  Since I is L_1-stable and L_1 preserves PBW length, grI_N
(the degree-N part of grI) is an sl_2-submodule of Sym^N(V_2) (V_2 = lowest weight 2 Verma module), hence a direct sum
of lowest-weight Verma modules; its lowest weight vectors are the degree-N minimal generators (HYP-9315(b)).

Cayley coordinates.  In the stable range (N >= L) a degree-N monomial of Li degree L is x_0^(N - l(mu)) x_mu with mu a
partition of L (parts >= 1).  The semi-invariants (ker D_1) are spanned by x_0^(N - |nu|) prod A_(nu_i) over partitions
nu of L into parts >= 2, where (Tschirnhaus normalisation, x_1 killed)
    A_j = sum_(i=0..j) C(j+3, i) (-1/4)^i x_0^(j-1-i) x_1^i x_(j-i).
This script computes, for each model and each Li degree L = 2k+1 <= 2p-3, the lowest weight vector of grI_(2N+L, N),
normalised by its x_0^(N-1) x_L coefficient, in Cayley coordinates, exactly (CRT over three primes + rational
reconstruction), and compares models.
Usage: python semiinvariant_generators_20261018.py
"""
import os, sys, json, itertools
from fractions import Fraction
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'voa_singular_support'))
import flint
from vir_ss import SingularSupport, block, block_index, central_charge, ratrec

PRIMES = (2147483647, 2147483629, 2147483587)


def partitions(n, maxpart=None, minpart=1):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), minpart - 1, -1):
        for rest in partitions(n - k, k, minpart):
            yield (k,) + rest


def D1_block(vec):
    """D_1 on a dict over partitions lam (parts >= 2, part k <-> x_(k-2)): x_j -> (j+3) x_(j-1), i.e. part k -> (k+1) * part (k-1), k >= 3."""
    out = {}
    for lam, co in vec.items():
        lst = list(lam)
        for i, k in enumerate(lst):
            if k >= 3:
                new = tuple(sorted(lst[:i] + [k - 1] + lst[i + 1:], reverse=True))
                out[new] = out.get(new, 0) + co * (k + 1)
    return out


def to_mu(lam):
    return tuple(sorted((k - 2 for k in lam if k > 2), reverse=True))


# ---- Cayley semi-invariants in mu-coordinates (x_0 -> 1)
def poly_mul(f, g):
    out = {}
    for a, ca in f.items():
        for b, cb in g.items():
            m = tuple(sorted(a + b, reverse=True))
            out[m] = out.get(m, 0) + ca * cb
    return {k: v for k, v in out.items() if v}


def cayley_A(j):
    f = {}
    for i in range(0, j + 1):
        coef = Fraction(comb(j + 3, i)) * Fraction(-1, 4) ** i
        # monomial x_1^i x_(j-i) (x_0 powers dropped; x_0 itself contributes nothing)
        parts = [1] * i + ([j - i] if j - i > 0 else [])
        m = tuple(sorted(parts, reverse=True))
        f[m] = f.get(m, 0) + coef
    return {k: v for k, v in f.items() if v}


def cayley_basis(L):
    nus = [nu for nu in partitions(L, minpart=2)]
    polys = []
    for nu in nus:
        f = {(): Fraction(1)}
        for j in nu:
            f = poly_mul(f, cayley_A(j))
        polys.append(f)
    return nus, polys


def check_cayley_semiinvariant(L):
    """D_1 A_nu = 0 in mu-coordinates: D_1 x_j = (j+3) x_(j-1), x_1 -> 4 x_0 = 4."""
    nus, polys = cayley_basis(L)
    for f in polys:
        out = {}
        for mu, co in f.items():
            lst = list(mu)
            for i, j in enumerate(lst):
                new = lst[:i] + ([j - 1] if j > 1 else []) + lst[i + 1:]
                m = tuple(sorted(new, reverse=True))
                out[m] = out.get(m, 0) + co * (j + 3)
        assert all(v == 0 for v in out.values()), (L, out)
    return True


def lowest_weight_vector(ss, D, N, P):
    """Return list of nullspace vectors (dicts over block(D,N)) of D_1 restricted to grI_(D,N)."""
    g = ss.grI.get((D, N), [])
    if not g:
        return []
    tgt = block(D - 1, N)
    tidx = block_index(D - 1, N)
    cols = []
    for v in g:
        w = D1_block(v)
        col = [0] * len(tgt)
        for lam, co in w.items():
            col[tidx[lam]] = co % P
        cols.append(col)
    r = len(g)
    M = flint.nmod_mat(len(tgt), r, [cols[j][i] for i in range(len(tgt)) for j in range(r)], P)
    ns, nullity = M.nullspace()
    out = []
    for t in range(nullity):
        a = [int(ns[j, t]) for j in range(r)]
        vec = {}
        for aj, v in zip(a, g):
            if aj:
                for lam, co in v.items():
                    vec[lam] = (vec.get(lam, 0) + aj * co) % P
        out.append({k: x for k, x in vec.items() if x})
    return out


def solve_cayley_modp(vmu, L, P):
    """Coordinates of vmu (dict mu -> int mod P) in the Cayley basis; checks consistency."""
    nus, polys = cayley_basis(L)
    mus = list(partitions(L))
    midx = {m: i for i, m in enumerate(mus)}
    def red(x):
        return x.numerator % P * pow(x.denominator % P, -1, P) % P
    rows = len(mus)
    colsn = len(nus)
    A = [[0] * (colsn + 1) for _ in range(rows)]
    for j, f in enumerate(polys):
        for m, co in f.items():
            A[midx[m]][j] = red(co)
    for m, co in vmu.items():
        A[midx[m]][colsn] = co % P
    M = flint.nmod_mat(rows, colsn + 1, [x for row in A for x in row], P)
    R, rk = M.rref()
    Mc = flint.nmod_mat(rows, colsn, [A[i][j] for i in range(rows) for j in range(colsn)], P)
    assert Mc.rank() == colsn
    assert rk == colsn, 'vector is not in the Cayley span (not a semi-invariant?)'
    sol = [0] * colsn
    Rl = R.tolist()
    for i in range(rk):
        row = [int(x) for x in Rl[i]]
        piv = next(j for j in range(colsn + 1) if row[j])
        sol[piv] = row[colsn]
    return nus, sol


def crt(residues, primes):
    x, m = 0, 1
    for r, p in zip(residues, primes):
        t = (r - x) * pow(m, -1, p) % p
        x += m * t
        m *= p
    return x, m


def ratrec_big(a, m):
    """Rational reconstruction with |num|, den <= sqrt(m/2)."""
    import math
    bound = math.isqrt(m // 2)
    r0, r1 = m, a % m
    s0, s1 = 0, 1
    while r1 > bound:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)


def model_generators(p, q):
    N = (p - 1) * (q - 1) // 2
    Lmax = 2 * p - 3
    per_prime = []
    for P in PRIMES:
        ss = SingularSupport(p, q, 2 * N + Lmax, P)
        ss.run(verbose=False)
        res = {}
        for L in range(0, Lmax + 1):
            lw = lowest_weight_vector(ss, 2 * N + L, N, P)
            entry = {'nullity': len(lw)}
            if len(lw) == 1 and L >= 2:
                v = lw[0]
                vmu = {}
                for lam, co in v.items():
                    vmu[to_mu(lam)] = (vmu.get(to_mu(lam), 0) + co) % P
                top = vmu.get((L,), 0)
                if top:
                    inv = pow(top, -1, P)
                    vmu = {k: x * inv % P for k, x in vmu.items()}
                    nus, sol = solve_cayley_modp(vmu, L, P)
                    entry['cayley'] = dict(zip(nus, sol))
                else:
                    entry['top_zero'] = True
            res[L] = entry
        per_prime.append(res)
    out = {}
    for L in range(0, Lmax + 1):
        nul = [r[L]['nullity'] for r in per_prime]
        assert len(set(nul)) == 1
        e = {'nullity': nul[0]}
        if 'cayley' in per_prime[0][L]:
            co = {}
            for nu in per_prime[0][L]['cayley']:
                x, m = crt([r[L]['cayley'][nu] for r in per_prime], PRIMES)
                co[nu] = ratrec_big(x, m)
            e['cayley'] = co
        out[L] = e
    return N, out


if __name__ == '__main__':
    for L in range(2, 10):
        check_cayley_semiinvariant(L)
    print('Cayley A_j products are D_1-invariant in mu-coordinates for L = 2..9: OK')
    models = [(4, 5), (4, 7), (4, 9), (4, 11), (5, 6), (5, 7), (5, 8), (5, 9), (6, 7), (3, 4), (3, 5), (3, 7), (3, 8)]
    table = {}
    for (p, q) in models:
        N, out = model_generators(p, q)
        c = central_charge(p, q)
        table[f'{p},{q}'] = {'N': N, 'c': str(c),
                             'L': {str(L): {'nullity': e['nullity'],
                                            'cayley': {str(k): str(v) for k, v in e.get('cayley', {}).items()}}
                                   for L, e in out.items()}}
        nz = {L: e['nullity'] for L, e in out.items() if e['nullity']}
        print(f'Vir({p},{q}) N={N} c={c}: lowest weight vectors in grI_N at Li degrees {nz}')
        for L, e in out.items():
            if 'cayley' in e:
                print(f'    L={L}: ' + ', '.join(f'A{list(k)}: {v}' for k, v in e['cayley'].items()))
    json.dump(table, open(os.path.join(HERE, 'semiinvariant_generators.json'), 'w'), indent=1)
