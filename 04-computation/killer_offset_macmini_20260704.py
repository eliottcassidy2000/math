#!/usr/bin/env python3
"""
The killer-offset mechanism + single-killer formula (mac-mini-2026-07-04-S43). THM-618.
Covering forces a killer a with (n-1)|a; at t=1/(n-1) (small-base optimum) the killer is at 0 => hiding offsets.
Formula (proved): M({1..n-2, X}) = X/((n-1)(X+1)) for (n-1)n | X, minimized at X=(n-1)n=182 => 14/183.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
def gcd_all(xs): return reduce(gcd, xs)
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def M_exact(sp):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)): Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0); bt = None
    for q in Q:
        if q < 2: continue
        for a in range(1, q):
            t = F(a, q); m = min(nd(v*t) for v in sp)
            if m > best: best, bt = m, t
    return best, bt
if __name__ == "__main__":
    n = 14; phi6 = n*n-n+1
    print(f"covering-min = 1/(n-1) - 1/((n-1)Phi6) = {F(1,n-1)-F(1,(n-1)*phi6)} = n/Phi6 = {F(n,phi6)}")
    print(f"single-killer formula M({{1..{n-2},X}}) = X/((n-1)(X+1)):")
    for X in [182, 364, 546, 728, 2366]:
        M, t = M_exact(list(range(1, n-1))+[X])
        pred = F(X, (n-1)*(X+1))
        print(f"  X={X:>5}: M={M} pred={pred} match={M==pred} t*={t} (offset 1/13-t*={F(1,n-1)-t})")
    print("=> monotone increasing in X => min at smallest killer X=(n-1)n=182 => 14/183 (deep well).")
