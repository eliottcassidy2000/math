#!/usr/bin/env python3
"""First homotopy of the derived arc space of a quasi-homogeneous plane curve: the Cramer syzygy of a contact locus
(S44; explains kind-pasteur's cusp "curiosity", singular_support_arc_spaces_classical_freeness_20261018.md 5a).

For f(x, y) quasi-homogeneous of weight d (wt x = a, wt y = b), the derived arc space is the Koszul complex of the
jets f_n = d^n f (weight d + n) on C[x_j, y_j] (weights a + j, b + j).  Its Euler characteristic is the product
    E(t) = prod_{n>=0}(1 - t^(d+n)) / prod_{j>=0}((1 - t^(a+j))(1 - t^(b+j))),
and the arc-space Hilbert series H_0(t) equals E(t) below the first weight w_1 with H_1 != 0.

Mechanism (predicted w_1).  For a monomial contact locus C(i, k) = {x_0..x_(i-1) = y_0..y_(k-1) = 0} the jets
f_0, ..., f_(m-1) vanish on C, m = ord_C f = min over monomials x^u y^v of f of (i u + k v).  If m > i + k (surplus),
the K + 1 = i + k + 1 jets f_0..f_K lie in the K-generated ideal (x_0..x_(i-1), y_0..y_(k-1)), so Cramer's rule gives
the determinantal syzygy sum_l (-1)^l det(A_l) f_l = 0 of weight
    w(C) = sum_{l=0}^{i+k} (d + l) - sum_{j<i}(a + j) - sum_{j<k}(b + j).
Prediction: w_1 = min over surplus contact loci of w(C) (the minimum is attained near the lct-computing locus).
This script computes H_0 by exact rank (mod P) weight by weight and reports the first weight where H_0 != E, against
the prediction.  Controls: the node xy (log canonical: no surplus, H_0 = E forever).
PRE-REGISTERED (committed before the long runs; predictor = minimal nonzero Cramer weight over monomial ideals):
    E6 y^3 - x^4 (wt 3,4): first H_0 != E at weight 29;  A4 y^2 - x^5 (wt 2,5): first H_0 != E at weight 31.
    (Already observed before registration: D4 10, x^2y+y^3 10, A3 16, x^4+y^4 11; cusp 34 by kind-pasteur.)
Usage: python derived_arc_curves_first_syzygy_20261018.py [Wmax] [curve-filter]
Session opus-2026-10-18-S44.
"""
import sys, os, json, itertools, functools
import flint
P = 2147483647


def predicted_w1(a, b, d, monos, imax=12, kmax=12):
    """Contact-locus Cramer weight (monomial contact loci only)."""
    best = None
    arg = None
    for i in range(0, imax + 1):
        for k in range(0, kmax + 1):
            if i + k == 0:
                continue
            m = min(i * u + k * v for (u, v) in monos)
            if m > i + k:
                w = sum(d + l for l in range(i + k + 1)) - sum(a + j for j in range(i)) - sum(b + j for j in range(k))
                if best is None or w < best:
                    best, arg = w, (i, k, m)
    return best, arg


def predicted_w1_ideals(a, b, d, fs, jmax=4, emax=4):
    """Cramer weight minimized over monomial ideals Q = (x_j^(e_j), y_j^(g_j)) (at most one generator per jet
    variable, j <= jmax, exponents <= emax) containing the first K+1 jets f_0..f_K, K = number of generators."""
    import itertools
    WT = {'x': a, 'y': b}
    def in_Q(poly, gens):
        for mo in poly:
            cnt = {}
            for v in mo:
                cnt[v] = cnt.get(v, 0) + 1
            if not any(cnt.get(var, 0) >= e for var, e in gens):
                return False
        return True
    best, arg = None, None
    choices = []
    for var in [(v, j) for v in ('x', 'y') for j in range(jmax + 1)]:
        choices.append([None] + [(var, e) for e in range(1, emax + 1)])
    # restrict to "downward closed in j" supports to keep the search small: generators on x_0..x_(i-1), y_0..y_(k-1)
    for i in range(0, jmax + 2):
        for k in range(0, jmax + 2):
            K = i + k
            if K == 0:
                continue
            xs = [('x', j) for j in range(i)]
            ys = [('y', j) for j in range(k)]
            for ex in itertools.product(range(1, emax + 1), repeat=K):
                gens = list(zip(xs + ys, ex))
                if len(fs) < K + 1:
                    continue
                if all(in_Q(fs[l], gens) for l in range(K + 1)):
                    w = sum(d + l for l in range(K + 1)) - sum(e * (WT[v[0]] + v[1]) for v, e in gens)
                    if (best is None or w < best) and full_rank(fs, gens, K):
                        best, arg = w, tuple(('%s%d^%d' % (v[0], v[1], e)) for v, e in gens)
    return best, arg


def full_rank(fs, gens, K, trials=2):
    """Is the (K+1) x K coefficient matrix A (f_l = sum_j A_lj g_j, each monomial assigned to the first generator
    dividing it) of rank K over the jet field?  Random evaluation mod P.  Rank K makes the Cramer syzygy nonzero."""
    import random
    for _ in range(trials):
        val = {}
        def ev(cnt):
            r = 1
            for var, c in cnt.items():
                if c:
                    if var not in val:
                        val[var] = random.randrange(1, P)
                    r = r * pow(val[var], c, P) % P
            return r
        rows = []
        for l in range(K + 1):
            row = [0] * K
            for mo, co in fs[l].items():
                cnt = {}
                for v in mo:
                    cnt[v] = cnt.get(v, 0) + 1
                for j, (var, e) in enumerate(gens):
                    if cnt.get(var, 0) >= e:
                        rest = dict(cnt)
                        rest[var] -= e
                        row[j] = (row[j] + co * ev(rest)) % P
                        break
            rows.append(row)
        M = flint.nmod_mat(K + 1, K, [x for r in rows for x in r], P)
        if M.rank() == K:
            return True
    return False


def run_curve(name, a, b, f_terms, Wmax):
    """f_terms: list of (coeff, u, v) for monomials x^u y^v (all of weight d)."""
    d = a * f_terms[0][1] + b * f_terms[0][2]
    assert all(a * u + b * v == d for _, u, v in f_terms)
    WT = {'x': a, 'y': b}

    @functools.lru_cache(maxsize=None)
    def monos(w):
        vars_ = [(v, j) for v in ('x', 'y') for j in range(0, w) if WT[v] + j <= w]
        vars_.sort(key=lambda t: (WT[t[0]] + t[1], t[0]))
        out = []
        def rec(i, rem, cur):
            if rem == 0:
                out.append(tuple(sorted(cur)))
                return
            for kk in range(i, len(vars_)):
                v = vars_[kk]
                wv = WT[v[0]] + v[1]
                if wv > rem:
                    continue
                cur.append(v)
                rec(kk, rem - wv, cur)
                cur.pop()
        rec(0, w, [])
        return tuple(out)

    def mul(x, y):
        return tuple(sorted(x + y))

    def deriv(poly):
        out = {}
        for mo, co in poly.items():
            for i, (v, j) in enumerate(mo):
                new = tuple(sorted(mo[:i] + ((v, j + 1),) + mo[i + 1:]))
                out[new] = (out.get(new, 0) + co * (j + 1)) % P
        return {k: v for k, v in out.items() if v}

    f0 = {}
    for co, u, v in f_terms:
        mo = tuple(sorted([('x', 0)] * u + [('y', 0)] * v))
        f0[mo] = (f0.get(mo, 0) + co) % P
    fs = [f0]
    while d + len(fs) - 1 < Wmax:
        fs.append(deriv(fs[-1]))
    # Euler product
    E = [0] * (Wmax + 1)
    E[0] = 1
    def mulser(s, k, sign):
        # multiply series s by (1 - t^k)^(sign) truncated
        if sign > 0:
            for i in range(Wmax, k - 1, -1):
                s[i] -= s[i - k]
        else:
            for i in range(k, Wmax + 1):
                s[i] += s[i - k]
    for j in range(0, Wmax + 1):
        if a + j <= Wmax:
            mulser(E, a + j, -1)
        if b + j <= Wmax:
            mulser(E, b + j, -1)
    for n in range(0, Wmax + 1):
        if d + n <= Wmax:
            mulser(E, d + n, +1)
    H0 = []
    first = None
    for w in range(0, Wmax + 1):
        M = monos(w)
        idx = {m: i for i, m in enumerate(M)}
        rows = []
        for n, fn in enumerate(fs):
            rest = w - d - n
            if rest < 0:
                break
            for m in monos(rest):
                row = {}
                for mo, co in fn.items():
                    kk = idx[mul(mo, m)]
                    row[kk] = (row.get(kk, 0) + co) % P
                rows.append(row)
        if rows and M:
            mat = flint.nmod_mat(len(rows), len(M), [r.get(j, 0) for r in rows for j in range(len(M))], P)
            rk = mat.rank()
        else:
            rk = 0
        h = len(M) - rk
        H0.append(h)
        if first is None and h != E[w]:
            first = (w, h, E[w])
    pred, arg = predicted_w1(a, b, d, [(u, v) for _, u, v in f_terms])
    fs_long = [f0]
    while len(fs_long) < 14:
        fs_long.append(deriv(fs_long[-1]))
    pred2, arg2 = predicted_w1_ideals(a, b, d, fs_long)
    print('%-22s wt(x,y)=(%d,%d) d=%2d | contact-locus w = %s %s | monomial-ideal w = %s %s | observed first H0 != E: %s'
          % (name, a, b, d, pred, arg, pred2, arg2, first if first else 'none through %d' % Wmax), flush=True)
    return dict(name=name, a=a, b=b, d=d, predicted_contact=pred, locus=arg, predicted_ideal=pred2, ideal=arg2,
                observed=first, H0=H0, E=E)


def main():
    Wmax = int(sys.argv[1]) if len(sys.argv) > 1 else 22
    curves = [
        ('node xy', 1, 1, [(1, 1, 1)]),
        ('D4: x^3 + y^3', 1, 1, [(1, 3, 0), (1, 0, 3)]),
        ('x^2 y + y^3 (D4)', 1, 1, [(1, 2, 1), (1, 0, 3)]),
        ('A3: y^2 - x^4', 1, 2, [(1, 0, 2), (-1, 4, 0)]),
        ('x^4 + y^4', 1, 1, [(1, 4, 0), (1, 0, 4)]),
        ('E6: y^3 - x^4', 3, 4, [(1, 0, 3), (-1, 4, 0)]),
        ('A4: y^2 - x^5', 2, 5, [(1, 0, 2), (-1, 5, 0)]),
        ('cusp y^2 - x^3', 2, 3, [(1, 0, 2), (-1, 3, 0)]),
    ]
    filt = sys.argv[2] if len(sys.argv) > 2 else None
    res = []
    for name, a, b, terms in curves:
        if filt and filt not in name:
            continue
        res.append(run_curve(name, a, b, terms, Wmax))
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'derived_arc_curves_first_syzygy_%s_W%d.json' % ((filt or 'all').split(':')[0].replace(' ', ''), Wmax)), 'w'), indent=1)


if __name__ == '__main__':
    main()
