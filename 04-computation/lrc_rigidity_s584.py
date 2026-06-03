#!/usr/bin/env python3
"""Rigidity in LRC. (1) symmetry group G(V)={u in (Z/n)* : uV=V mod n}; AP has full
(Z/n)*. (2) witness rigidity: is M achieved only on the clock {j/n} (measure-0, RIGID)
or on an interval (FLEXIBLE)? (3) duality: high symmetry => rigid clock witness => M=1/n
exactly; low symmetry => flexible => M>1/n. opus-2026-06-03-S584."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
def dist(x): x%=1; return min(x,1-x)
def units(n): return [u for u in range(1,n) if gcd(u,n)==1]
def symmetry(V,n):
    Vm=sorted(set(v%n for v in V))
    G=[u for u in units(n) if sorted(set((u*v)%n for v in Vm))==Vm]
    return G
def Mexact(V):
    cands=set()
    for i in range(len(V)):
        for j in range(len(V)):
            if i==j: continue
            for D in (V[i]+V[j],abs(V[i]-V[j])):
                if D:
                    for m in range(1,D): cands.add(F(m,D))
        cands.add(F(1,2*V[i]))
    best=F(0); arg=[]
    for t in cands:
        mn=min(dist(v*t) for v in V)
        if mn>best: best=mn; arg=[t]
        elif mn==best: arg.append(t)
    return best,arg
def on_clock(args,n):  # are ALL optimal times of form j/n?
    return all((t*n)==int(t*n) for t in args)
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def main():
    print("RIGIDITY SURVEY: symmetry group |G(V)| vs M vs witness-on-clock (rigid)")
    import random; rng=random.Random(3)
    for n in [6,8,10,14]:
        delta=F(1,n); m=n-1
        AP=tuple(range(1,n))
        G=symmetry(AP,n); M,args=Mexact(AP)
        print(f"  n={n}: AP {AP}: |sym|={len(G)} (phi(n)={len(units(n))}), full={len(G)==len(units(n))}; "
              f"M={float(M):.4f}(delta={float(delta):.3f}); witness-on-clock(rigid)={on_clock(args,n)}, #witnesses={len(args)}")
        # spectrum: random configs, correlate symmetry with M and rigidity
        buckets={}
        B=2*n+4
        for _ in range(2500):
            V=prim(tuple(sorted(rng.sample(range(1,B+1),m))))
            if len(V)!=m: continue
            g=len(symmetry(V,n)); M2,a2=Mexact(V)
            key=g
            buckets.setdefault(key,[]).append((float(M2-delta),on_clock(a2,n)))
        line=[]
        for g in sorted(buckets):
            vals=buckets[g]; mn=min(v[0] for v in vals); rig=sum(1 for v in vals if v[1])/len(vals)
            line.append(f"|sym|={g}: minMargin={mn:+.3f},clock-rigid%={100*rig:.0f}(n={len(vals)})")
        print("       "+"; ".join(line))
if __name__=='__main__': main()
