#!/usr/bin/env python3
"""The recursive fractal concept of TRANSLATED APs. (a) literal translates
{c+1,..,c+n-1}: M(c) structure & multiples; (b) flip=translation: worry-set = AP with
elements translated to antipodes; (c) binary IFS {x->2x, x->2x+1} & the doubling map
x->2x mod n dynamics (order of 2). opus-2026-06-03-S580."""
from fractions import Fraction as F
from math import gcd
def dist(x): x%=1; return min(x,1-x)
def Mexact(V):
    cands=set()
    for i in range(len(V)):
        for j in range(len(V)):
            if i==j: continue
            for D in (V[i]+V[j],abs(V[i]-V[j])):
                if D:
                    for m in range(1,D): cands.add(F(m,D))
        cands.add(F(1,2*V[i]))
    best=F(0)
    for t in cands:
        mn=min(dist(v*t) for v in V)
        if mn>best: best=mn
    return best
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def main():
    print("(a) LITERAL translated APs {c+1,..,c+(n-1)} (k=n-1 consecutive): M vs delta=1/n")
    for n in [8,12,14]:
        delta=F(1,n); row=[]
        for c in range(0,2*n):
            V=tuple(range(c+1,c+n)); 
            if len(V)!=n-1: continue
            Vp=prim(V); M=Mexact(Vp)
            mult=any(v%n==0 for v in V)
            row.append((c,float(M),float(M-delta),mult))
        print(f"  n={n}: (c, M, M-delta, has-mult-of-n)")
        for c,M,mar,mult in row[:8]:
            print(f"     c={c:2d}: M={M:.4f} margin={mar:+.4f} mult={mult}  V={tuple(range(c+1,c+n))}")
    print()
    print("(c) binary IFS {x->2x, x->2x+1} from 1 = the integers (binary tree); AP = truncation")
    print("    doubling map x->2x mod n on {1..n-1}: orbit of 1, and order/degeneracy")
    for n in [7,13,14,16]:
        orb=[1]; x=1; seen={1}
        for _ in range(2*n):
            x=(2*x)%n
            if x==0 or x in seen: orb.append(x); break
            orb.append(x); seen.add(x)
        # multiplicative order of 2 mod n (if gcd(2,n)=1)
        if n%2==1:
            o=1; y=2%n
            while y!=1 and o<2*n: y=(2*y)%n; o+=1
            ordstr=f"ord_2(mod {n})={o}"
        else:
            ordstr=f"2 NOT invertible mod {n} (degenerate doubling)"
        print(f"   n={n:2d}: orbit of 1 under x->2x mod n: {orb}; {ordstr}")
if __name__=='__main__': main()
