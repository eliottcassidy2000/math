#!/usr/bin/env python3
"""Independent arc-space Hilbert series (audit2).  Pure commutative algebra, no vertex-algebra code, no author code.

R = C[x,y]/(f_1..f_r), wt x = 2, wt y = 3.  Arc algebra C[J_infty X] = C[x_j, y_j : j >= 0]/(F_{i,n}), where
f_i(sum_j x_j t^j, sum_j y_j t^j) = sum_n F_{i,n} t^n.  Bigraded by weight W (x_j: 2+j, y_j: 3+j) and Li degree
L (= j), F_{i,n} has bidegree (wt f_i + n, n).  HS(W) = sum_L (#monomials(W,L) - rank span{m F_{i,n}}).
Ranks mod P with python-flint; polynomials given as dict {(a,b): coeff (Fraction or int)} meaning coeff x^a y^b.
"""
import sys
from fractions import Fraction
from functools import lru_cache
import flint

WX, WY = 2, 3


def Fp(x, P):
    x = Fraction(x)
    return x.numerator % P * pow(x.denominator % P, P - 2, P) % P


class Arc:
    def __init__(self, rels, Wmax, P=2147483647):
        self.P = P
        self.Wmax = Wmax
        self.rels = []
        for f in rels:
            ws = {WX * a + WY * b for (a, b) in f}
            assert len(ws) == 1, ('not weighted homogeneous', f)
            self.rels.append((ws.pop(), {k: Fp(v, P) for k, v in f.items() if Fraction(v) != 0}))
        self.Fs = [(w, self.series(f, Wmax - w)) for (w, f) in self.rels]
        self._monos = {}

    # monomials: sorted tuples of ('x', j) / ('y', j)
    def series(self, f, N):
        """coefficients F_0..F_N of f(x(t), y(t))"""
        P = self.P
        xs = [{(('x', j),): 1} for j in range(N + 1)]
        ys = [{(('y', j),): 1} for j in range(N + 1)]

        def mulpoly(A, B):
            out = {}
            for ma, ca in A.items():
                for mb, cb in B.items():
                    m = tuple(sorted(ma + mb))
                    out[m] = (out.get(m, 0) + ca * cb) % P
            return {m: c for m, c in out.items() if c}

        def mulser(S, T):
            out = [dict() for _ in range(N + 1)]
            for i in range(N + 1):
                if not S[i]:
                    continue
                for j in range(N + 1 - i):
                    if not T[j]:
                        continue
                    pr = mulpoly(S[i], T[j])
                    o = out[i + j]
                    for m, c in pr.items():
                        o[m] = (o.get(m, 0) + c) % P
            return [{m: c for m, c in o.items() if c} for o in out]
        one = [{(): 1}] + [dict() for _ in range(N)]
        powx = [one]
        powy = [one]
        maxa = max(a for (a, b) in f)
        maxb = max(b for (a, b) in f)
        for _ in range(maxa):
            powx.append(mulser(powx[-1], xs))
        for _ in range(maxb):
            powy.append(mulser(powy[-1], ys))
        tot = [dict() for _ in range(N + 1)]
        for (a, b), co in f.items():
            pr = mulser(powx[a], powy[b])
            for n in range(N + 1):
                for m, c in pr[n].items():
                    tot[n][m] = (tot[n].get(m, 0) + c * co) % P
        return [{m: c for m, c in o.items() if c} for o in tot]

    def monos(self, W, L):
        key = (W, L)
        if key in self._monos:
            return self._monos[key]
        # variables with weight <= W and Li <= L
        vars_ = []
        for j in range(0, L + 1):
            if WX + j <= W:
                vars_.append(('x', j))
            if WY + j <= W:
                vars_.append(('y', j))
        vars_.sort()
        out = []

        def rec(i, w, l, cur):
            if w == 0 and l == 0:
                out.append(tuple(cur))
                return
            if w <= 0:
                return
            for k in range(i, len(vars_)):
                v = vars_[k]
                wv = (WX if v[0] == 'x' else WY) + v[1]
                lv = v[1]
                if wv > w or lv > l:
                    continue
                cur.append(v)
                rec(k, w - wv, l - lv, cur)
                cur.pop()
        rec(0, W, L, [])
        self._monos[key] = out
        return out

    def block_rank(self, W, L):
        P = self.P
        cols = self.monos(W, L)
        if not cols:
            return 0, 0
        ci = {m: i for i, m in enumerate(cols)}
        rows = []
        for (wf, F) in self.Fs:
            for n in range(0, L + 1):
                Wm, Lm = W - wf - n, L - n
                if Wm < 0 or n >= len(F):
                    continue
                for m in self.monos(Wm, Lm):
                    row = {}
                    for mo, co in F[n].items():
                        k = ci[tuple(sorted(mo + m))]
                        row[k] = (row.get(k, 0) + co) % P
                    rows.append(row)
        if not rows:
            return len(cols), 0
        nr, nc = len(rows), len(cols)
        flat = [0] * (nr * nc)
        for r, row in enumerate(rows):
            base = r * nc
            for k, v in row.items():
                flat[base + k] = v
        M = flint.nmod_mat(nr, nc, flat, P)
        return nc, M.rank()

    def hs(self, W):
        tot = 0
        for L in range(0, W + 1):
            nc, rk = self.block_rank(W, L)
            tot += nc - rk
        return tot

    def hs_list(self, Wmax=None, verbose=False):
        Wmax = Wmax or self.Wmax
        out = []
        for W in range(Wmax + 1):
            out.append(self.hs(W))
            if verbose:
                print(W, out[-1], flush=True)
        return out


def a2_product(q, W):
    a = [0] * (W + 1)
    for n in range(1, W + 1):
        r = n % q
        a[n] = (0 if r in (0, 1, q - 1) else 1) + (0 if r in (0, 1, q - 1, 2, q - 2) else 1)
    return pe_from_exponents(a, W)


def pe_from_exponents(a, W):
    ser = [1] + [0] * W
    for n in range(1, W + 1):
        for _ in range(a[n]):
            for i in range(n, W + 1):
                ser[i] += ser[i - n]
    return ser


def cusp_ci(W):
    # prod_{n>=2} (1-t^n)^-1 / ((1-t^3)(1-t^4)(1-t^5))
    a = [0] * (W + 1)
    for n in range(2, W + 1):
        a[n] += 1
    for n in (3, 4, 5):
        if n <= W:
            a[n] += 1
    return pe_from_exponents(a, W)


if __name__ == '__main__':
    pass
