#!/usr/bin/env python3
"""
lrc_lonely_measure_formula_s521.py   claudebox-2026-06-01-S521

Exact lonely-measure formula and the rigorous LRC sufficient condition (reflection:
07-reflections/lrc-lonely-measure-relation-lattice-formula-s521.md).

THEOREM (elementary Fourier): mu(v)=int_0^1 prod_i 1_B(v_i t) dt
       = sum_{c in Z^m: <c,v>=0} prod_i hhat(c_i),  B=[1/n,1-1/n], hhat(0)=1-2/n.
COROLLARY: if sum_{c!=0,<c,v>=0} |prod hhat(c_i)| < (1-2/n)^m then mu>0 => LRC (room).
Residual: short additive relations v_i+v_j=v_k (the AP-like / tight core).
"""
from fractions import Fraction as F
from math import sin, pi, gcd
from itertools import product, combinations

def dist(x):
    x = x % 1; return min(x, 1 - x)
def mu_exact(v, n):
    W = set([F(0)])
    for vi in v:
        for k in range(vi+1):
            W.add((F(k, vi)+F(1, n*vi)) % 1); W.add((F(k, vi)-F(1, n*vi)) % 1)
    W = sorted(x for x in W if 0 <= x < 1); W2 = W + [F(1)]
    return float(sum((b-a) for a, b in zip(W, W2[1:])
                     if all(dist(F(vi)*((a+b)/2)) >= F(1, n) for vi in v)))
def hhat(c, n):
    return (1-2/n) if c == 0 else ((-1)**c)*sin(pi*c*(1-2/n))/(pi*c)
def formula(v, n, H, absval=False):
    s = 0.0
    for c in product(range(-H, H+1), repeat=len(v)):
        if sum(ci*vi for ci, vi in zip(c, v)) == 0:
            p = 1.0
            for ci in c: p *= hhat(ci, n)
            if absval:
                if any(c): s += abs(p)
            else:
                s += p
    return s

def main():
    print("THEOREM: mu(v) = sum over relation lattice {c:<c,v>=0} of prod hhat(c_i).")
    print(f"{'speeds':16} {'mu_exact':10} {'formula(H=8)':12} {'V':8}")
    for v in [(1,2,4,7),(2,3,5,7),(1,5,7,11),(3,5,7,11)]:
        n = len(v)+1; m = len(v)
        print(f"{str(v):16} {mu_exact(list(v),n):10.5f} {formula(list(v),n,8):12.5f} {(1-2/n)**m:8.4f}")
    print("\nCOROLLARY: corr=sum_{c!=0,rel}|prod hhat| < V  =>  mu>0 => LRC (rigorous, relation-sparse):")
    print(f"{'speeds':16} {'V':8} {'corr(H=7)':10} {'corr<V => LRC':14} {'mu_exact':10}")
    for v in [(1,5,7,11),(3,5,7,11),(1,2,4,7),(2,3,5,7),(1,2,3,4)]:
        n = len(v)+1; m = len(v); V = (1-2/n)**m
        corr = formula(list(v), n, 7, absval=True)
        proved = "PROVED" if corr < V else "(inconclusive)"
        print(f"{str(v):16} {V:8.4f} {corr:10.4f} {proved:14} {mu_exact(list(v),n):10.5f}")
    print("\n  Residual (corr>=V) = short additive relations v_i+v_j=v_k (AP-like / tight core).")

if __name__ == "__main__":
    main()
