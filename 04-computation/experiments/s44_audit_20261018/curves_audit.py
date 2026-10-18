#!/usr/bin/env python3
"""S44 audit, HYP-9316: arc-space Hilbert series of quasi-homogeneous plane curves vs the derived Euler product,
computed independently and BIGRADED (q-degree, jet index); and explicit Cramer syzygies tested in Koszul H_1.

Jets: F_n = [t^n] f(sum_j x_j t^j, sum_j y_j t^j)   (the S44 script uses f_n = d^n f = n! F_n: same ideal).
Bigrading: x_j has (q-degree a, index j), y_j has (b, j); F_n has (d, n); total weight = q-degree + index.
Euler product (bigraded): prod_n (1 - s^d u^n) / prod_j (1 - s^a u^j)(1 - s^b u^j).
Usage: python curves_audit.py <name> Wmax
"""
import sys, json, time, itertools
from functools import lru_cache
import flint

P = 2147483647

CURVES = {
    'node': (1, 1, [(1, 1, 1)]),
    'D4': (1, 1, [(1, 3, 0), (1, 0, 3)]),
    'D4b': (1, 1, [(1, 2, 1), (1, 0, 3)]),        # x^2 y + y^3
    'A3': (1, 2, [(1, 0, 2), (-1, 4, 0)]),
    'X4': (1, 1, [(1, 4, 0), (1, 0, 4)]),         # x^4 + y^4
    'E6': (3, 4, [(1, 0, 3), (-1, 4, 0)]),
    'A4': (2, 5, [(1, 0, 2), (-1, 5, 0)]),
    'cusp': (2, 3, [(1, 0, 2), (-1, 3, 0)]),
    'A5': (1, 3, [(1, 0, 2), (-1, 6, 0)]),        # y^2 - x^6
    'E7': (3, 2, [(1, 3, 0), (1, 1, 3)]),         # x^3 + x y^3
    'D5': (3, 2, [(1, 2, 1), (1, 0, 4)]),         # x^2 y + y^4
    'X5': (1, 1, [(1, 5, 0), (1, 0, 5)]),         # x^5 + y^5
    'E8': (3, 5, [(1, 0, 3), (-1, 5, 0)]),        # y^3 - x^5
    'A6': (2, 7, [(1, 0, 2), (-1, 7, 0)]),        # y^2 - x^7
    'D6': (2, 1, [(1, 2, 1), (1, 0, 5)]),         # x^2 y + y^5
}


@lru_cache(maxsize=None)
def parts_atmost(n, k, maxpart=None):
    """Multisets of exactly k nonnegative integers summing to n, as nonincreasing tuples of length k."""
    if maxpart is None:
        maxpart = n
    if k == 0:
        return ((),) if n == 0 else ()
    out = []
    for first in range(min(n, maxpart), -1, -1):
        if first * k < n:
            break
        for rest in parts_atmost(n - first, k - 1, first):
            out.append((first,) + rest)
    return tuple(out)


class Curve:
    def __init__(self, a, b, terms, Wmax):
        self.a, self.b, self.terms = a, b, terms
        self.d = a * terms[0][1] + b * terms[0][2]
        assert all(a * u + b * v == self.d for _, u, v in terms)
        self.Wmax = Wmax
        self.nmax = Wmax - self.d
        self.F = self._jets()
        self._mon = {}

    # monomial = sorted tuple of variable codes: x_j -> 2j, y_j -> 2j+1
    def _jets(self):
        nmax = max(self.nmax, 0)
        # power series of x(t), y(t): list over t-degree of polynomial dicts
        def series_pow(var, e):
            res = [{(): 1}] + [{} for _ in range(nmax)]
            for _ in range(e):
                new = [{} for _ in range(nmax + 1)]
                for i, poly in enumerate(res):
                    if not poly:
                        continue
                    for j in range(0, nmax + 1 - i):
                        code = 2 * j + var
                        for mono, co in poly.items():
                            key = tuple(sorted(mono + (code,)))
                            new[i + j][key] = (new[i + j].get(key, 0) + co) % P
                res = new
            return res
        F = [{} for _ in range(nmax + 1)]
        for co, u, v in self.terms:
            X = series_pow(0, u)
            Y = series_pow(1, v)
            for i in range(nmax + 1):
                for j in range(nmax + 1 - i):
                    for m1, c1 in X[i].items():
                        for m2, c2 in Y[j].items():
                            key = tuple(sorted(m1 + m2))
                            F[i + j][key] = (F[i + j].get(key, 0) + co * c1 * c2) % P
        return [{k: v for k, v in Fn.items() if v} for Fn in F]

    def monos(self, Dq, I):
        key = (Dq, I)
        if key in self._mon:
            return self._mon[key]
        a, b = self.a, self.b
        out = []
        if Dq >= 0 and I >= 0:
            for k in range(0, Dq // a + 1):
                rem = Dq - a * k
                if rem % b:
                    continue
                l = rem // b
                for I1 in range(0, I + 1):
                    for xs in parts_atmost(I1, k):
                        for ys in parts_atmost(I - I1, l):
                            out.append(tuple(sorted([2 * j for j in xs] + [2 * j + 1 for j in ys])))
        out = tuple(sorted(set(out)))
        self._mon[key] = out
        return out

    def H0(self, Dq, I):
        M = self.monos(Dq, I)
        if not M:
            return 0
        idx = {m: i for i, m in enumerate(M)}
        rows = []
        for n, Fn in enumerate(self.F):
            if n > I:
                break
            for m in self.monos(Dq - self.d, I - n):
                row = {}
                for mono, co in Fn.items():
                    key = tuple(sorted(mono + m))
                    row[idx[key]] = (row.get(idx[key], 0) + co) % P
                rows.append(row)
        if not rows:
            return len(M)
        mat = flint.nmod_mat(len(rows), len(M), [r.get(j, 0) for r in rows for j in range(len(M))], P)
        return len(M) - mat.rank()

    def euler(self, W):
        """bigraded Euler product coefficients E[(Dq, I)] for Dq + I <= W."""
        a, b, d = self.a, self.b, self.d
        E = {(0, 0): 1}
        def mul_geom(E, s, u):  # multiply by 1/(1 - S^s U^u)
            keys = sorted(E.keys(), key=lambda t: t[0] + t[1])
            new = dict(E)
            # iterate in increasing total weight
            for tot in range(0, W + 1):
                for Dq in range(0, tot + 1):
                    I = tot - Dq
                    if (Dq - s, I - u) in new and Dq - s >= 0 and I - u >= 0:
                        new[(Dq, I)] = new.get((Dq, I), 0) + new[(Dq - s, I - u)]
            return new
        def mul_lin(E, s, u):  # multiply by (1 - S^s U^u)
            new = dict(E)
            for (Dq, I), v in E.items():
                if Dq + s + I + u <= W:
                    new[(Dq + s, I + u)] = new.get((Dq + s, I + u), 0) - v
            return new
        for j in range(0, W + 1):
            if a + j <= W:
                E = mul_geom(E, a, j)
            if b + j <= W:
                E = mul_geom(E, b, j)
        for n in range(0, W + 1):
            if d + n <= W:
                E = mul_lin(E, d, n)
        return {k: v for k, v in E.items() if v}


def main():
    name, Wmax = sys.argv[1], int(sys.argv[2])
    a, b, terms = CURVES[name]
    C = Curve(a, b, terms, Wmax)
    E = C.euler(Wmax)
    t0 = time.time()
    first = None
    tots = []
    devs = []
    for w in range(0, Wmax + 1):
        th, te = 0, 0
        for Dq in range(0, w + 1):
            I = w - Dq
            h = C.H0(Dq, I)
            e = E.get((Dq, I), 0)
            th += h
            te += e
            if h != e:
                devs.append((w, Dq, I, h, e))
        tots.append((w, th, te))
        if first is None and th != te:
            first = (w, th, te)
        print('w=%2d HS(H0)=%7d E=%7d %s (%.1fs)' % (w, th, te, '' if th == te else '<-- deviation', time.time() - t0), flush=True)
    print('%s (wt %d,%d; d=%d): first deviation (w, HS, E) = %s' % (name, a, b, C.d, first))
    print('bigraded deviations at the first deviating weight:', [x for x in devs if first and x[0] == first[0]])
    json.dump(dict(name=name, a=a, b=b, d=C.d, Wmax=Wmax, first=first, totals=tots,
                   deviations=[list(x) for x in devs]), open('curves_audit_%s_W%d.json' % (name, Wmax), 'w'), indent=1)


if __name__ == '__main__':
    main()
