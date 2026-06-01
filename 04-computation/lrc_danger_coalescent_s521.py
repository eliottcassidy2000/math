#!/usr/bin/env python3
"""
lrc_danger_coalescent_s521.py   claudebox-2026-06-01-S521

The danger-graph coalescent model for LRC, made computational.
(Reflection: 07-reflections/lrc-danger-coalescent-formalization-s521.md)

n=m+1, threshold 1/n. Danger count N(t)=#{i: ||v_i t||<1/n}; observer lonely <=> N(t)=0.
N(0)=m, mean E[N]=2(n-1)/n.  Two views computed:
  (A) N(t) excursion over cells: does it reach 0 in an interior cell? (tight sets only at boundary)
  (B) discretized covering of Z/N_* (N_*=n*lcm): LRC <=> bad sets don't cover; report uncovered count.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def dist(x):
    x = x % 1; return min(x, 1 - x)

def N_cells(sp, n):
    ev = set()
    for v in sp:
        for k in range(v):
            ev.add(F(k*n - 1, n*v) % 1); ev.add(F(k*n + 1, n*v) % 1)
    pts = sorted(e for e in ev if 0 <= e < 1)
    aug = [F(0)] + pts + [F(1)]
    return [sum(1 for v in sp if dist(F(v)*((a+b)/2)) < F(1, n))
            for a, b in zip(aug, aug[1:]) if 0 < (a+b)/2 < 1]

def lcm(a, b): return a*b//gcd(a, b)
def covering(sp):
    n = len(sp)+1; L = reduce(lcm, sp); Ns = n*L
    covered = set()
    for v in sp:
        for j in range(Ns):
            r = (v*j) % Ns
            if r < L or r > Ns - L: covered.add(j)
    return n, Ns, len(set(range(Ns)) - covered)

def main():
    print("(A) N(t) danger-count excursion (LRC: N reaches 0):")
    print(f"{'speeds':22} {'minN':5} {'maxN':5} {'E[N]':6} {'cells N=0':10}")
    from collections import Counter
    for sp in [(1,2,3,4),(1,2,4,7),(1,3,4,5,9),(2,3,5,7,11),(1,5,6,11,16),(1,2,3,4,5,6)]:
        n = len(sp)+1; cells = N_cells(list(sp), n); c = Counter(cells)
        print(f"{str(sp):22} {min(cells):5} {max(cells):5} {float(2*(n-1)/n):6.2f} {c.get(0,0):10}")
    print("\n(B) discretized covering of Z/N_* (uncovered = lonely discrete times):")
    print(f"{'speeds':22} {'n':3} {'N_*':7} {'uncovered':10}")
    for sp in [(1,2,3),(1,2,3,4),(1,2,4,7),(1,3,4,5,9),(2,3,5,7)]:
        n, Ns, unc = covering(list(sp))
        print(f"{str(sp):22} {n:3} {Ns:7} {unc:10}")
    print("\n(C) fraction of sets lonely only at the boundary (tight/extremal), by n:")
    random.seed(0)
    for m in (5, 6, 7):
        n = m+1; only_bdry = tested = 0
        for _ in range(600):
            sp = sorted(random.sample(range(1, 2*n), m))
            g = 0
            for v in sp: g = gcd(g, v)
            if g != 1: continue
            tested += 1
            if min(N_cells(sp, n)) > 0: only_bdry += 1
        print(f"  n={n}: {only_bdry}/{tested} sets reach N=0 only at the boundary wall (extremal)")

if __name__ == "__main__":
    main()
