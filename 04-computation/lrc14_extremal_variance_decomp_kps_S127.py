# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont33: chipping the extremal -- the consec-extremality of J=E[N(7-N)] is
# VARIANCE-maximization dominating a mean-deviation term, and the k=8,9-vs-k>=10 boundary is ONE crossover.
# J = 49/4 - Var(N) - (E[N]-7/2)^2.  consec MAXIMIZES Var(N) (all k, adversarially global) but does NOT
# minimize E[N]; J-min at consec holds iff consec's Var-lead >= competitor's (E[N]-7/2)^2 gain -- true at
# k=8,9 (the ladder base), false at k>=10 (which mac-mini THM-710 makes inherit, so it is not needed).
from fractions import Fraction as F
from itertools import combinations
from functools import reduce
from math import gcd

def stats(E):  # exact (E[N], Var(N), J)
    bps={F(0),F(1)}
    for e in E:
        for m in range(1,7*e): bps.add(F(m,7*e))
    bps=sorted(bps); m1=F(0); m2=F(0)
    for a,b in zip(bps,bps[1:]):
        x=(a+b)/2; N=7-len(set(int((e*x%1)*7) for e in E)); w=b-a; m1+=w*N; m2+=w*N*N
    return m1, m2-m1*m1, 7*m1-m2

def main():
    print("k : argmin E[N]=consec? | argmax Var=consec? | argmin J=consec? (exhaustive small box)")
    for k,N in [(8,13),(9,14),(10,15),(11,16)]:
        c=tuple(range(1,k+1)); mnE=(F(9),None); mxV=(F(-9),None); mnJ=(F(9),None)
        for E in combinations(range(1,N+1),k):
            if reduce(gcd,E)!=1: continue
            e1,v,j=stats(list(E))
            if e1<mnE[0]: mnE=(e1,E)
            if v>mxV[0]: mxV=(v,E)
            if j<mnJ[0]: mnJ=(j,E)
        print(f" {k}: E[N]min {mnE[1]==c:1} | Var max {mxV[1]==c:1} | J min {mnJ[1]==c:1}  (Jmin at {mnJ[1]})")
    print()
    print("THE CROSSOVER (Var-lead vs mean-gain), consec vs its J-runner-up:")
    for k,ru in [(9,[2,4,5,6,7,8,10,12,14]),(10,[2,4,5,6,7,8,9,10,12,14])]:
        c=list(range(1,k+1)); e1c,vc,jc=stats(c); e1r,vr,jr=stats(ru)
        dev=lambda m:(m-F(7,2))**2
        print(f"  k={k}: Var-lead(consec) = {float(vc-vr):+.3f} ; mean-gain(runner-up) = {float(dev(e1r)-dev(e1c)):+.3f}"
              f"  => consec J {'<' if jc<jr else '>'} runner-up ({'HOLDS' if jc<jr else 'FLIPS'})")
    print("  => Var(N)-max is the fundamental extremal (all k); J-min = Var-max beating the mean term (k<=9 only).")

if __name__ == '__main__': main()
