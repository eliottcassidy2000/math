"""Fast minimizer-structure check: do all M=2/(2n-1) C2b minimizers contain {n, n-1} (shell
partners, forbidden j's collide)? Float screen + exact confirm. opus-2026-06-03-S599l."""
from itertools import combinations
from fractions import Fraction as F
from math import gcd
def dist_f(x):
    x-=int(x); 
    if x<0:x+=1
    return min(x,1-x)
def M_float(V):
    ds=set(V)
    for a,b in combinations(V,2): ds.add(a+b); ds.add(abs(a-b))
    ds.discard(0); best=-1.0
    for d in ds:
        for mm in range(d):
            val=min(dist_f(v*mm/d) for v in V)
            if val>best: best=val
    return best
def M_exact(V):
    ds=set(V)
    for a,b in combinations(V,2): ds.add(a+b); ds.add(abs(a-b))
    ds.discard(0); best=F(-1)
    for d in ds:
        for mm in range(d):
            t=F(mm,d); val=min(min((v*t)%1,1-(v*t)%1) for v in V)
            if val>best: best=val
    return best
def main():
    print(" n | #minimizers (M=2/(2n-1)) | all contain {n,n-1}? | counterexamples")
    for n in range(3,8):
        m=2*n-1; B=3*n; thr=2.0/m; mins=[]
        for V in combinations(range(1,B+1),n-1):
            if not any(v%n==0 for v in V): continue
            g=0
            for v in V: g=gcd(g,v)
            if g!=1: continue
            if abs(M_float(V)-thr)<1e-9 and M_exact(V)==F(2,m): mins.append(V)
        bad=[V for V in mins if not (n in V and (n-1) in V)]
        print(f" {n:2d} | {len(mins):3d} | {len(bad)==0} | {bad[:4]}")
if __name__=='__main__': main()
