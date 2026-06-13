"""The RESIDUAL of LRC after the shell-partner lemma + clock witness: configs with a multiple of n
AND no shell-partner. Count, min M, is it still loose (M≥2/(2n-1)), and via what witness. n=5,6,7.
opus-2026-06-07-S700b."""
from itertools import combinations
from math import gcd
from fractions import Fraction as F
def Mexact(V):
    ds=set(V)
    for a,b in combinations(V,2): ds.add(a+b); ds.add(abs(a-b))
    ds.discard(0); best=F(-1)
    for d in ds:
        for m in range(d):
            t=F(m,d); v=min(min((x*t)%1,1-(x*t)%1) for x in V)
            if v>best: best=v
    return best
def forbidden_count(V,C):
    Fset=set([0])
    for v in V:
        for m in range(1,C):
            if (v*m)%C in (1,C-1): Fset.add(m)
    return len(Fset)
def has_shell(V,C): return any((V[i]+V[j])%C==0 for i in range(len(V)) for j in range(i+1,len(V)))
def main():
    print("RESIDUAL: multiple of n, NO shell-partner, gcd 1.  n | count | min M | ≥2/(2n-1)? | discrete-witness fails")
    for n in range(5,8):
        C=2*n-1; B=2*n; thr=F(2,C); cnt=0; minM=F(2); dfail=0; examples=[]
        for V in combinations(range(1,B+1),n-1):
            if not any(v%n==0 for v in V): continue
            g=0
            for v in V: g=gcd(g,v)
            if g!=1: continue
            if has_shell(V,C): continue
            cnt+=1; m=Mexact(V)
            if m<minM: minM=m; 
            if forbidden_count(V,C)>=C: dfail+=1; examples.append(V)
        print(f" {n} | {cnt} | minM={minM}={float(minM):.4f} | loose(≥2/(2n-1)): {minM>=thr} | discrete-fails={dfail} {examples[:2]}")
    print("\n=> If minM ≥ 2/(2n-1) on the residual, then ALL multiples configs are loose (C' holds):")
    print("   shell-partner ones by the LEMMA (discrete witness), no-shell ones by a coarse witness.")
    print("   The discrete-fails are exactly where the 2n-1 tick can't fire (need the coarse t=1/k).")
if __name__=='__main__': main()
