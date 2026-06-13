#!/usr/bin/env python3
"""
Signed-LRC OBSERVER-INCLUSIVE mutual gap Gstar_full: is inf_S = 2/(4n-5) stable
at high speed bound B?   monad-explorer-2026-06-06-S3.

W_full(A,B) = {v_i : all i} u {|v_i-v_j| same side} u {v_i+v_j across}.
Gstar_full(S) = max over cuts M(W_full).  (Observer pairs = bare speeds, sign-blind.)

Families run found inf_S Gstar_full = 2/7,2/11,2/15,2/19 for n=3..6 = 2/(4n-5).
This script tests whether that closed form survives larger B (set-level pruning).
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd
from functools import reduce
import sys

def norm(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def candidate_times(Ws):
    T = set()
    for w in Ws:
        for a in range(0, w): T.add(F(2*a+1, 2*w))
    for x, y in combinations(Ws, 2):
        for d in (x+y, x-y):
            if d == 0: continue
            d = abs(d)
            for a in range(1, d): T.add(F(a, d))
    return T

def maximin_gt(Ws, thresh):
    best = F(0)
    for t in candidate_times(Ws):
        ok=True; m=F(1,2)
        for w in Ws:
            nv = norm(w*t)
            if nv <= best: ok=False; break
            if nv < m: m = nv
        if ok and m > best:
            best = m
            if best > thresh: return best, True
    return best, False

def cutWs(V, side):
    r=len(V); W=list(V)
    for i,j in combinations(range(r),2):
        W.append(abs(V[i]-V[j]) if side[i]==side[j] else V[i]+V[j])
    return sorted(set(W))

def gstar_full_capped(V, inf):
    r=len(V); best=F(0)
    for tail in product([0,1],repeat=r-1):
        side=(0,)+tail
        m,reached = maximin_gt(cutWs(V,side), inf)
        if m>best: best=m
        if reached: return best, False
    return best, True

def enum_sets(r,B):
    return [c for c in combinations(range(1,B+1),r) if reduce(gcd,c)==1]

def search(n,B):
    r=n-1; sets=enum_sets(r,B); inf=F(1,2); args=[]
    for V in sets:
        g,cand=gstar_full_capped(V,inf)
        if cand and g<inf: inf=g; args=[V]
        elif cand and g==inf: args.append(V)
    return inf,args,len(sets)

def main():
    if len(sys.argv)>1:
        n=int(sys.argv[1]); B=int(sys.argv[2])
        g,a,ns=search(n,B)
        pred=F(2,4*n-5)
        print(f"n={n} B<={B} ({ns} sets): inf Gstar_full={g}={float(g):.6f}  "
              f"2/(4n-5)={pred}={float(pred):.6f}  {'MATCH' if g==pred else 'DIFFERS'}  n*inf={float(n*g):.4f}")
        for V in a[:6]: print(f"   min {V}")
        return
    print("inf_S Gstar_full vs 2/(4n-5)  (monad-S3)")
    for n,Bs in [(4,[12,16,20]),(5,[10,13,16]),(6,[9,11,13])]:
        pred=F(2,4*n-5)
        print(f"\n=== n={n}  predicted 2/(4n-5)={pred}={float(pred):.6f} ===")
        for B in Bs:
            g,a,ns=search(n,B)
            print(f"  B<={B:2d} ({ns:5d} sets): inf={g}={float(g):.6f}  {'MATCH' if g==pred else 'LOWER' if g<pred else 'HIGHER'}  "
                  f"n*inf={float(n*g):.4f}  min e.g. {a[:3]}")
            sys.stdout.flush()

if __name__=="__main__":
    main()
