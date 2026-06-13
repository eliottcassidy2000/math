"""C2b ANGLE 1: the DISCRETE witness at resolution 2n-1. M(S)>=2/(2n-1) iff some dilate jS
avoids {0,±1} mod (2n-1). Forbidden j = union {±v_i^{-1}}; collisions happen at shell-partners
(v_i+v_k≡0, THM-401). Test: does the discrete witness exist for all C2b configs? Compare M_exact.
opus-2026-06-03-S599i."""
from itertools import combinations
from fractions import Fraction as F
from math import gcd
def nrm(a,m): a%=m; return min(a,m-a)
def discrete_best(V,m):
    # max over j of min_i ||v_i j / m||  (in units of 1/m)
    best=0
    for j in range(1,m):
        val=min(nrm(v*j,m) for v in V)
        if val>best: best=val
    return best   # integer units of 1/m; witness certifies M >= best/m
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
    print("C2b ANGLE 1 — discrete witness at 2n-1. 'good j' = dilate jS avoids {0,±1} (best>=2).")
    print(" n | #C2b | discrete fails(best<2) | min M_exact | =2/(2n-1)? | when discrete fails: M_exact")
    for n in range(3,9):
        m=2*n-1; B=2*n; cnt=0; failj=0; minM=F(2); failMs=set()
        for V in combinations(range(1,B+1),n-1):
            if not any(v%n==0 for v in V): continue
            g=0
            for v in V: g=gcd(g,v)
            if g!=1: continue
            cnt+=1
            db=discrete_best(V,m)
            if db<2: failj+=1; failMs.add(M_exact(V))
            mm=M_exact(V)
            if mm<minM: minM=mm
        print(f" {n:2d} | {cnt:5d} | {failj:4d} | {minM}={float(minM):.4f} | {minM==F(2,m)} | {sorted(map(str,failMs))[:6]}")
    print("\nReading: if discrete-witness 'fails' count is 0, then jS avoids {0,±1} for ALL C2b")
    print("configs => M>=2/(2n-1) PROVEN by the 2n-1 clock tick. Else: the failing configs need")
    print("a finer witness; their M_exact tells us the true floor.")
if __name__=='__main__': main()
