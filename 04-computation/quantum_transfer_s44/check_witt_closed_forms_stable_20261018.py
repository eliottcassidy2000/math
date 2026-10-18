#!/usr/bin/env python3
"""Stable-range evaluation of the pre-registered tests P9, P10 of THM-4710 (PREREGISTRATION_witt_closed_forms_20261018.md).

Why a second implementation.  The pre-registered run (witt_rigidity_interpolation / closure in the PBW basis
block(2N+L, N)) died with MemoryError at N = 50: vir_ss.block enumerates every partition of 2N+L into parts >= 2
(about 10^8 at N = 50) before filtering by length.  For N >= L the degree-N part of Li degree L has the faithful
coordinates mu (partitions of L, monomial x_0^(N - l(mu)) x_mu), in which every operator is an exact formula in N:
    d(mu)   = (N - l(mu)) (mu + {1}) + sum_i (mu_i + 1) (mu_i -> mu_i + 1)
    D_n(mu) = sum_(mu_i >= n) (mu_i + n + 2) (mu_i -> mu_i - n, zero parts dropped)          (n >= 1).
The recursion (generators at Li 3, 5, 7, 9; S_L = admissible lowest-weight vectors) is identical.  This implementation
is VALIDATED first: for N = 9..49 it must reproduce dims {3,5,7,9: 1} and the closed forms that the PBW-basis runs
reproduced exactly on 9..48 (and N = 49, the first pre-registered value, PASSED in the PBW basis before the crash).
Then P9 (N = 49..64) and P10 (Vir_(4,13), PBW basis, weight <= 43) are evaluated.
Usage: python check_witt_closed_forms_stable_20261018.py
"""
import os, sys, json
from fractions import Fraction
import sympy as sp
import flint
HERE = os.path.dirname(os.path.abspath(__file__))
from semiinvariant_generators_20261018 import partitions, solve_cayley_modp, crt, ratrec_big, model_generators

PRIMES = (2147483647, 2147483629, 2147483587, 2147483579, 2147483563, 2147483549)
Nsym = sp.Symbol('N')


def d_mu(vec, N, P):
    out = {}
    for mu, c in vec.items():
        new = tuple(sorted(mu + (1,), reverse=True))
        out[new] = (out.get(new, 0) + c * (N - len(mu))) % P
        lst = list(mu)
        for i, j in enumerate(lst):
            new = tuple(sorted(lst[:i] + [j + 1] + lst[i + 1:], reverse=True))
            out[new] = (out.get(new, 0) + c * (j + 1)) % P
    return {k: v for k, v in out.items() if v}


def Dn_mu(vec, n, P):
    out = {}
    for mu, c in vec.items():
        lst = list(mu)
        for i, j in enumerate(lst):
            if j >= n:
                rest = lst[:i] + ([j - n] if j - n > 0 else []) + lst[i + 1:]
                new = tuple(sorted(rest, reverse=True))
                out[new] = (out.get(new, 0) + c * (j + n + 2)) % P
    return {k: v for k, v in out.items() if v}


def basis(L):
    return list(partitions(L)) if L > 0 else [()]


def to_row(vec, idx):
    row = [0] * len(idx)
    for k, v in vec.items():
        row[idx[k]] = v
    return row


def rref(rows, n, P):
    if not rows:
        return []
    M = flint.nmod_mat(len(rows), n, [x % P for r in rows for x in r], P)
    R, rk = M.rref()
    return [[int(x) for x in r] for r in R.tolist()[:rk]]


def closure_mu(N, P, Lmax=9, pattern=(3, 5, 7, 9)):
    M = {0: [{(): 1}]}
    S = {0: [{(): 1}]}
    for L in range(1, Lmax + 1):
        cols = basis(L)
        idx = {m: i for i, m in enumerate(cols)}
        n = len(cols)
        dM = [v for v in (d_mu(v, N, P) for v in M[L - 1]) if v]
        # lowest-weight space: nullspace of D_1: mu-coords of L -> L-1
        tgt = basis(L - 1)
        tidx = {m: i for i, m in enumerate(tgt)}
        img = [to_row(Dn_mu({m: 1}, 1, P), tidx) for m in cols]
        A = flint.nmod_mat(len(tgt), n, [img[j][i] for i in range(len(tgt)) for j in range(n)], P)
        ns, nul = A.nullspace()
        lw = [[int(ns[i, t]) for i in range(n)] for t in range(nul)]
        cond = []
        for m in range(2, L + 1):
            cc = basis(L - m)
            cidx = {x: i for i, x in enumerate(cc)}
            Mrows = rref([to_row(v, cidx) for v in M[L - m]], len(cc), P)
            piv = [next(j for j, x in enumerate(r) if x) for r in Mrows]
            ims = []
            for vec in lw:
                g = {cols[i]: vec[i] for i in range(n) if vec[i]}
                row = to_row(Dn_mu(g, m, P), cidx)
                for r, pj in zip(Mrows, piv):
                    if row[pj]:
                        f = row[pj]
                        row = [(a - f * b) % P for a, b in zip(row, r)]
                ims.append(row)
            for c in range(len(cc)):
                cond.append([ims[t][c] for t in range(len(lw))])
        if lw and cond:
            C = flint.nmod_mat(len(cond), len(lw), [x for r in cond for x in r], P)
            ns2, nul2 = C.nullspace()
            sols = [[int(ns2[t, s]) for t in range(len(lw))] for s in range(nul2)]
        elif lw:
            sols = [[1 if t == s else 0 for t in range(len(lw))] for s in range(len(lw))]
        else:
            sols = []
        Sv = []
        for a in sols:
            v = [0] * n
            for t, at in enumerate(a):
                for i in range(n):
                    v[i] = (v[i] + at * lw[t][i]) % P
            Sv.append({cols[i]: v[i] for i in range(n) if v[i]})
        S[L] = Sv
        M[L] = dM + (Sv if L in pattern else [])
        rows = rref([to_row(v, idx) for v in M[L]], n, P)
        M[L] = [{cols[i]: r[i] for i in range(n) if r[i]} for r in rows]
    return S


def forced(N):
    per = [closure_mu(N, P) for P in PRIMES]
    dims = {L: len(per[0][L]) for L in per[0]}
    assert all({L: len(s[L]) for L in s} == dims for s in per)
    out = {}
    for L in (5, 7, 9):
        if dims.get(L) == 1:
            co = []
            for s, P in zip(per, PRIMES):
                v = s[L][0]
                top = v.get((L,), 0)
                inv = pow(top, -1, P)
                vmu = {k: x * inv % P for k, x in v.items()}
                nus, sol = solve_cayley_modp(vmu, L, P)
                co.append(dict(zip(nus, sol)))
            out[L] = {nu: ratrec_big(*crt([c[nu] for c in co], PRIMES)) for nu in co[0]}
    return dims, out


def compare(N, g, formulas):
    mism = []
    for (L, nu), f in formulas.items():
        val = g.get(L, {}).get(nu)
        pred = sp.Rational(f.subs(Nsym, N))
        if val is None or sp.Rational(val.numerator, val.denominator) != pred:
            mism.append((L, nu, str(val), str(pred)))
    return mism


if __name__ == '__main__':
    fits = json.load(open(os.path.join(HERE, 'witt_rigidity_interpolation_9_48.json')))['fits']
    formulas = {}
    for key, ent in fits.items():
        L = int(key.split()[0].split('_')[1])
        nu = tuple(int(x) for x in key.split(' ', 1)[1].strip('()').split(',') if x.strip())
        formulas[(L, nu)] = sp.sympify(ent['expr'], locals={'N': Nsym})
    res = {'validation_9_49': [], 'P9': [], 'P10': None}
    # validation against the PBW-basis results (exact closed forms on 9..48, and N = 49 PASSED in the PBW basis)
    vok = True
    for N in range(9, 50):
        dims, g = forced(N)
        okd = all(dims.get(L) == 1 for L in (3, 5, 7, 9))
        mism = compare(N, g, formulas)
        vok &= okd and not mism
        res['validation_9_49'].append({'N': N, 'dims_ok': okd, 'mismatches': mism})
    print(f'validation N = 9..49 (stable-range implementation vs PBW-basis results): {"OK" if vok else "FAILED"}', flush=True)
    # P9
    ok9 = True
    for N in range(49, 65):
        dims, g = forced(N)
        okd = all(dims.get(L) == 1 for L in (3, 5, 7, 9))
        mism = compare(N, g, formulas)
        ok9 &= okd and not mism
        res['P9'].append({'N': N, 'dims_ok': okd, 'mismatches': mism})
        print(f'N={N}: dims_ok={okd}, coefficient mismatches={len(mism)}', flush=True)
    print(f'P9 {"PASSED" if ok9 else "FAILED"}', flush=True)
    # P10 (PBW basis; small weights, no memory issue)
    N, r = model_generators(4, 13)
    gens = {L: e['nullity'] for L, e in r.items() if e['nullity']}
    cay5 = r[5].get('cayley', {})
    want = {(5,): Fraction(1), (3, 2): Fraction(-89, 73)}
    ok10 = (N == 18 and gens == {0: 1, 3: 1, 5: 1} and dict(cay5) == want)
    print(f'P10: Vir(4,13) N={N}, lowest-weight generators at Li {gens}, G_5 = {cay5}: {"PASSED" if ok10 else "FAILED"}', flush=True)
    from vir_ss import SingularSupport
    from semiinvariant_generators_20261018 import lowest_weight_vector
    Pp = 2147483647
    ss = SingularSupport(4, 13, 2 * 18 + 7, Pp)
    ss.run(verbose=False)
    pos7 = {L: len(lowest_weight_vector(ss, 2 * 18 + L, 18, Pp)) for L in range(0, 8)}
    pos7 = {L: n for L, n in pos7.items() if n}
    ok10b = pos7 == {0: 1, 3: 1, 5: 1}
    print(f'P10 (Li-7 clause): lowest-weight vectors of grI_18 through Li 7 at {pos7}: {"PASSED" if ok10b else "FAILED"}', flush=True)
    res['P10'] = {'N': N, 'positions': gens, 'G5': {str(k): str(v) for k, v in cay5.items()}, 'passed': ok10,
                  'positions_through_Li7': pos7, 'passed_Li7_clause': ok10b}
    json.dump(res, open(os.path.join(HERE, 'check_witt_closed_forms_stable.json'), 'w'), indent=1, default=str)
