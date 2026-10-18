#!/usr/bin/env python3
"""Hilbert series of the arc algebra C[J_infty X], X = Spec C[x,y]/(f_1..f_r), wt x = 2, wt y = 3
(jets x_j, y_j of weight 2+j, 3+j; d x_j = (j+1) x_(j+1)).  Pure commutative algebra, independent of any VOA code.
Usage: python arc_hs.py q Wmax   (uses the C_2-relations of W_3(3,q) found by the engine, rescaled)."""
import sys
from fractions import Fraction
import flint
P = 2147483647
WT = {'x': 2, 'y': 3}

def monos(w, maxvar=None):
    # monomials of weight w: sorted tuples of (var, j); generate by partition into parts labelled x/y
    out = []
    vars_ = [(v, j) for v in ('x', 'y') for j in range(0, w)]
    vars_ = [v for v in vars_ if WT[v[0]] + v[1] <= w]
    vars_.sort(key=lambda t: (WT[t[0]] + t[1], t[0]))
    def rec(i, rem, cur):
        if rem == 0:
            out.append(tuple(sorted(cur))); return
        for k in range(i, len(vars_)):
            v = vars_[k]; wv = WT[v[0]] + v[1]
            if wv > rem: continue
            cur.append(v); rec(k, rem - wv, cur); cur.pop()
    rec(0, w, [])
    return out

def mul(a, b): return tuple(sorted(a + b))

def deriv(poly):
    out = {}
    for mo, co in poly.items():
        for i, (v, j) in enumerate(mo):
            new = tuple(sorted(mo[:i] + ((v, j + 1),) + mo[i + 1:]))
            out[new] = (out.get(new, 0) + co * (j + 1)) % P
    return {k: v for k, v in out.items() if v}

def poly(terms):  # terms: list of (coeff Fraction, nx, ny) at jet level 0
    out = {}
    for co, nx, ny in terms:
        mo = tuple(sorted([('x', 0)] * nx + [('y', 0)] * ny))
        co = Fraction(co)
        out[mo] = co.numerator % P * pow(co.denominator % P, -1, P) % P
    return out

RELS = {  # C_2 relations of W_3(3,q) (engine, w3_c2_relations.py)
    7: [[(1, 1, 1)], [(1, 0, 2), (Fraction(-7, 39), 3, 0)]],
    8: [[(1, 0, 2), (Fraction(-16, 279), 3, 0)], [(1, 2, 1)]],
    10: [[(1, 1, 2), (Fraction(-2, 123), 4, 0)], [(1, 0, 3), (Fraction(-8, 41), 3, 1)]],
    11: [[(1, 0, 3), (Fraction(-11, 138), 3, 1)], [(1, 2, 2), (Fraction(-11, 1035), 5, 0)]],
    13: [[(1, 1, 3), (Fraction(-13, 448), 4, 1)], [(1, 0, 4), (Fraction(-13, 56), 3, 2), (Fraction(169, 150528), 6, 0)]],  # CRT over two primes (audit 2, MISTAKE-633)
}

def arc_hs(rels, W):
    gens = [poly(r) for r in rels]
    gw = [sum(WT[v] + j for (v, j) in next(iter(g))) for g in gens]
    # derivatives of generators up to needed order
    ders = []
    for g, w in zip(gens, gw):
        lst = [g]
        while w + len(lst) - 1 < W:
            lst.append(deriv(lst[-1]))
        ders.append((w, lst))
    hs = []
    for w in range(0, W + 1):
        M = monos(w); idx = {m: i for i, m in enumerate(M)}
        rows = []
        for (gw0, lst) in ders:
            for n, dg in enumerate(lst):
                rest = w - gw0 - n
                if rest < 0: break
                for m in monos(rest):
                    row = {}
                    for mo, co in dg.items():
                        k = idx[mul(mo, m)]; row[k] = (row.get(k, 0) + co) % P
                    rows.append(row)
        if rows and M:
            mat = flint.nmod_mat(len(rows), len(M), [r.get(j, 0) for r in rows for j in range(len(M))], P)
            rk = mat.rank()
        else:
            rk = 0
        hs.append(len(M) - rk)
    return hs

def a2_product(q, W):
    a = [0] * (W + 1)
    for n in range(1, W + 1):
        r = n % q
        a[n] = (0 if r in (0, 1, q - 1) else 1) + (0 if r in (0, 1, q - 1, 2, q - 2) else 1)
    ser = [1] + [0] * W
    for n in range(1, W + 1):
        for _ in range(a[n]):
            for i in range(n, W + 1): ser[i] += ser[i - n]
    return ser

if __name__ == '__main__':
    q, W = int(sys.argv[1]), int(sys.argv[2])
    hs = arc_hs(RELS[q], W)
    pr = a2_product(q, W)
    print(f"q={q}: arc HS   ", hs)
    print(f"      A2 product", pr)
    print("      EQUAL" if hs == pr else f"      DIFFER at {[i for i in range(W+1) if hs[i] != pr[i]][:5]}")
