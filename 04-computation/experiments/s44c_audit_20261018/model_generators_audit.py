#!/usr/bin/env python3
"""Audit check (c) of THM-4710: degree-N generators of grI for actual minimal models, with the independent engine.

For Vir_(p,q), N = (p-1)(q-1)/2:
  1. the singular vector v of weight 2N (kernel of L_1, L_2 on the PBW basis of V^c); checks it is unique and that its
     top length is N with symbol x_0^N;
  2. for each Li degree L <= Lmax: I_(2N+L) = span{L_(-lam_1) ... L_(-lam_s) v : lam partition of L}; the subspace of
     elements with no PBW component of length > N, and their length-N components = grI_N(L) (mu-coordinates);
  3. the lowest-weight vectors (kernel of D_1 on grI_N(L)), their number per L (the generator positions), and, for a
     single one with N >= L, its Cayley coordinates normalized by the A_L coefficient;
  4. comparison with THM-4710's closed forms G_5(N), G_7(N), G_9(N).
Exact over Q (P = 0) or modulo primes (P > 0; several primes give CRT + rational reconstruction for the Cayley
coefficients, and the predicted rationals are also compared directly modulo each prime).
Usage: python model_generators_audit.py p q Lmax P1 [P2 ...]      (P1 = 0 for exact rational arithmetic)
"""
import sys, os, json, time
from fractions import Fraction
from math import comb
import flint
from vir_pbw import VirVacuum, partitions_ge2, partitions, central_charge
from witt_closure_exact import cayley_basis, rref as rref_Q, nullspace as nullspace_Q

HERE = os.path.dirname(os.path.abspath(__file__))


def closed_forms(N):
    """THM-4710 closed forms (as written in the theorem file / JSON) evaluated at N, as Fractions."""
    import sympy as sp
    Ns = sp.Symbol('N')
    fits = json.load(open(os.path.join(HERE, '..', '..', 'quantum_transfer_s44', 'witt_rigidity_interpolation_9_48.json')))['fits']
    out = {}
    for key, ent in fits.items():
        L = int(key.split()[0].split('_')[1])
        nu = tuple(int(x) for x in key.split(' ', 1)[1].strip('()').split(',') if x.strip())
        val = sp.Rational(sp.sympify(ent['expr'], locals={'N': Ns}).subs(Ns, N))
        out.setdefault(L, {})[nu] = Fraction(int(val.p), int(val.q))
    for L in out:
        out[L][(L,)] = Fraction(1)
    out[3] = {(3,): Fraction(1)}
    return out


def red(x, P):
    if P == 0:
        return x
    x = Fraction(x)
    return x.numerator % P * pow(x.denominator % P, -1, P) % P


def nullspace_modp(rows, ncols, P):
    if not rows:
        return [[1 if i == j else 0 for i in range(ncols)] for j in range(ncols)]
    M = flint.nmod_mat(len(rows), ncols, [x % P for r in rows for x in r], P)
    X, nul = M.nullspace()
    return [[int(X[i, t]) for i in range(ncols)] for t in range(nul)]


def rref_modp(rows, ncols, P):
    if not rows:
        return [], []
    M = flint.nmod_mat(len(rows), ncols, [x % P for r in rows for x in r], P)
    R, rk = M.rref()
    Rl = [[int(x) for x in r] for r in R.tolist()[:rk]]
    piv = [next(j for j, x in enumerate(r) if x) for r in Rl]
    return Rl, piv


def nullspace_any(rows, ncols, P):
    if P == 0:
        return nullspace_Q(rows, ncols) if rows else [[Fraction(int(i == j)) for i in range(ncols)] for j in range(ncols)]
    return nullspace_modp(rows, ncols, P)


def rref_any(rows, ncols, P):
    if P == 0:
        return rref_Q(rows, ncols) if rows else ([], [])
    return rref_modp(rows, ncols, P)


def singular_vector(V, N, P):
    W = 2 * N
    monos = list(partitions_ge2(W))
    idx = {m: i for i, m in enumerate(monos)}
    t1 = list(partitions_ge2(W - 1)); t1i = {m: i for i, m in enumerate(t1)}
    t2 = list(partitions_ge2(W - 2)); t2i = {m: i for i, m in enumerate(t2)}
    nr = len(t1) + len(t2)
    if P == 0:
        rows = [[Fraction(0)] * len(monos) for _ in range(nr)]
        for j, m in enumerate(monos):
            for m2, a in V.act(1, m).items():
                rows[t1i[m2]][j] += a
            for m2, a in V.act(2, m).items():
                rows[len(t1) + t2i[m2]][j] += a
        ns = nullspace_any(rows, len(monos), P)
    else:
        Mt = flint.nmod_mat(nr, len(monos), P)      # zero matrix, filled sparsely
        for j, m in enumerate(monos):
            for m2, a in V.act(1, m).items():
                Mt[t1i[m2], j] = (int(Mt[t1i[m2], j]) + a) % P
            for m2, a in V.act(2, m).items():
                Mt[len(t1) + t2i[m2], j] = (int(Mt[len(t1) + t2i[m2], j]) + a) % P
        X, nul = Mt.nullspace()
        ns = [[int(X[i, t]) for i in range(len(monos))] for t in range(nul)]
    assert len(ns) == 1, f'singular vector space has dim {len(ns)}'
    v = ns[0]
    top = (2,) * N
    t = v[idx[top]]
    assert t != 0, 'L_-2^N coefficient vanishes'
    if P == 0:
        vec = {m: v[i] / t for i, m in enumerate(monos) if v[i]}
    else:
        inv = pow(int(t), -1, P)
        vec = {m: int(v[i]) * inv % P for i, m in enumerate(monos) if v[i]}
    return vec, len(monos)


def D1_mu(vec):
    out = {}
    for mu, c in vec.items():
        for i, j in enumerate(mu):
            rest = mu[:i] + ((j - 1,) if j - 1 > 0 else ()) + mu[i + 1:]
            new = tuple(sorted(rest, reverse=True))
            out[new] = out.get(new, 0) + c * (j + 3)
    return out


def to_mu(mono):
    return tuple(k - 2 for k in mono if k > 2)


def main(p, q, Lmax, primes):
    N = (p - 1) * (q - 1) // 2
    c = central_charge(p, q)
    pred = closed_forms(N) if N >= 9 else None
    results = {}
    for P in primes:
        t0 = time.time()
        V = VirVacuum(c, P if P else None)
        v, dimW = singular_vector(V, N, P)
        top_lengths = sorted(set(len(m) for m in v), reverse=True)
        print(f'[P={P}] Vir({p},{q}) N={N} c={c}: singular vector unique (dim V_2N = {dimW}), top length {top_lengths[0]}'
              f' ({time.time()-t0:.1f}s)', flush=True)
        assert top_lengths[0] == N
        w = {(): v}
        res = {}
        for L in range(1, Lmax + 1):
            lams = list(partitions(L))
            for lam in lams:
                w[lam] = V.apply(-lam[0], w[lam[1:]])
            # conditions: components of length > N vanish
            highmonos = sorted(set(m for lam in lams for m in w[lam] if len(m) > N))
            hidx = {m: i for i, m in enumerate(highmonos)}
            zero = Fraction(0) if P == 0 else 0
            rows = [[zero] * len(lams) for _ in highmonos]
            for j, lam in enumerate(lams):
                for m, a in w[lam].items():
                    if len(m) > N:
                        rows[hidx[m]][j] += a
            ns = nullspace_any(rows, len(lams), P)
            # length-N components in mu-coordinates
            mus = [m for m in partitions(L) if len(m) <= N]
            midx = {m: i for i, m in enumerate(mus)}
            gvecs = []
            for a in ns:
                g = [zero] * len(mus)
                for j, lam in enumerate(lams):
                    if a[j]:
                        for m, x in w[lam].items():
                            if len(m) == N:
                                g[midx[to_mu(m)]] += a[j] * x
                if P:
                    g = [x % P for x in g]
                gvecs.append(g)
            R, piv = rref_any(gvecs, len(mus), P)
            dimgr = len(R)
            # lowest-weight vectors: combinations of R killed by D_1
            t1 = [m for m in partitions(L - 1) if len(m) <= N] if L > 1 else [()]
            t1i = {m: i for i, m in enumerate(t1)}
            imgs = []
            for r in R:
                vv = {mus[i]: x for i, x in enumerate(r) if x}
                dv = D1_mu(vv)
                row = [zero] * len(t1)
                for m, x in dv.items():
                    row[t1i[m]] += x
                imgs.append(row)
            mat = [[imgs[j][i] for j in range(len(R))] for i in range(len(t1))]
            kn = nullspace_any(mat, len(R), P) if R else []
            entry = {'dim_grI_N': dimgr, 'n_lowest_weight': len(kn)}
            if len(kn) == 1 and N >= L and L >= 3:
                coeffs = kn[0]
                g = [zero] * len(mus)
                for j, r in enumerate(R):
                    for i, x in enumerate(r):
                        g[i] += coeffs[j] * x
                if P:
                    g = [x % P for x in g]
                gdict = {mus[i]: x for i, x in enumerate(g) if x}
                nus, polys = cayley_basis(L)
                topc = gdict[(L,)]
                if P == 0:
                    gn = {k: x / topc for k, x in gdict.items()}
                    rowsC = [[poly.get(mu, Fraction(0)) for poly in polys] + [gn.get(mu, Fraction(0))] for mu in mus]
                    RC, pivC = rref_Q(rowsC, len(nus) + 1)
                    assert len(nus) not in pivC
                    cay = {nus[pc]: RC[i][-1] for i, pc in enumerate(pivC)}
                else:
                    inv = pow(int(topc), -1, P)
                    gn = {k: x * inv % P for k, x in gdict.items()}
                    rowsC = [[red(poly.get(mu, 0), P) for poly in polys] + [gn.get(mu, 0)] for mu in mus]
                    RC, pivC = rref_modp(rowsC, len(nus) + 1, P)
                    assert len(nus) not in pivC
                    cay = {nus[pc]: RC[i][-1] for i, pc in enumerate(pivC)}
                entry['cayley'] = cay
                if pred is not None and L in pred:
                    entry['matches_closed_form'] = all(cay.get(nu, 0) == red(pred[L][nu], P) for nu in nus)
            res[L] = entry
            print(f'   L={L}: dim grI_N = {dimgr}, lowest-weight vectors {len(kn)}'
                  + (f", matches G_{L}(N) closed form: {entry.get('matches_closed_form')}" if 'matches_closed_form' in entry else '')
                  + (f"  cayley {dict((k, str(v)) for k, v in entry['cayley'].items())}" if (P == 0 and 'cayley' in entry) else '')
                  + f'  ({time.time()-t0:.1f}s)', flush=True)
        results[P] = res
        V.cache.clear()
    positions = {L: e['n_lowest_weight'] for L, e in results[primes[0]].items() if e['n_lowest_weight']}
    positions[0] = 1
    summary = {'p': p, 'q': q, 'N': N, 'c': str(c), 'Lmax': Lmax, 'primes': primes,
               'generator_positions': dict(sorted(positions.items())),
               'per_prime': {str(P): {str(L): {k: (str(v) if k != 'cayley' else {str(a): str(b) for a, b in v.items()})
                                               for k, v in e.items()} for L, e in r.items()} for P, r in results.items()}}
    json.dump(summary, open(os.path.join(HERE, f'model_generators_{p}_{q}_L{Lmax}.json'), 'w'), indent=1)
    print(f'Vir({p},{q}): generator positions through Li {Lmax}: {summary["generator_positions"]}')


if __name__ == '__main__':
    a = [int(x) for x in sys.argv[1:]]
    main(a[0], a[1], a[2], a[3:] or [0])
