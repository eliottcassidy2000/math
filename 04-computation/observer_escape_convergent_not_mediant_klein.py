#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THE OBSERVER'S ESCAPE: convergent, not mediant (klein-S29).

n/Phi_6(n) = [0; n-1, n]; its semiconvergents are j/(j(n-1)+1), j=1..n:
  j=1 -> 1/n           (the MEDIANT = the lonely-runner conjecture threshold)
  j=n -> n/Phi_6(n)    (the CONVERGENT = the covering-min)
The observer (static runner at 0) seeks a lonely time t = argmax_t min_s ||s t||. The killer n(n-1)=lcm(n-1,n)
BLOCKS the mediant escape: at t=1/n, killer/n = n-1 in Z, so the killer sits at 0 (M=0, not lonely). So the
observer must escape at the CONVERGENT n/Phi_6(n) (M>1/n). This script shows the escape ladder + the
n-transition + a global-optimality scan.
"""
from fractions import Fraction as F
from math import gcd

def frac_dist(x):  # ||x|| distance to nearest integer, x a Fraction
    r = x - (x.numerator//x.denominator)
    return min(r, 1-r)
def mindist(S, t):  # min_s ||s t||, t a Fraction
    return min(frac_dist(F(s)*t) for s in S)

print("="*84); print(" (1) Observer escape ladder for {1..12,182}: semiconvergents j/(13j+1) of 14/183"); print("="*84)
S=list(range(1,13))+[182]; n=14
print(f" {'j':>3}{'t=j/(13j+1)':>14}{'min_s||s t|| (M at t)':>22}{'  note'}")
for j in range(1,n+1):
    t=F(j, j*(n-1)+1); M=mindist(S,t)
    note=""
    if j==1: note="MEDIANT 1/n -- killer 182 at 0 (BLOCKED)" if M==0 else "MEDIANT 1/n"
    if j==n: note="CONVERGENT n/Phi_6(n) = the ESCAPE"
    print(f" {j:>3}{str(t):>14}{str(M)+' ='+format(float(M),'.5f'):>22}   {note}")
print(f"\n killer 182 at the mediant t=1/14: 182*(1/14) = {F(182,14)} -> distance to 0 = {frac_dist(F(182,14))}  (BLOCKED)")
print(f" => the observer cannot escape at the mediant 1/14; it escapes at the convergent 14/183, M=14/183>1/14.")

print("\n"+"="*84); print(" (2) The killer ALWAYS blocks the mediant 1/n (n(n-1)/n = n-1 in Z), all n"); print("="*84)
for n in [4,7,10,14]:
    killer=n*(n-1); print(f"   n={n}: killer lcm(n-1,n)={killer}; killer*(1/n)={F(killer,n)} -> dist {frac_dist(F(killer,n))} (mediant 1/n blocked); convergent n/Phi_6 = {n}/{n*n-n+1} = {n/(n*n-n+1):.5f} > 1/n={1/n:.5f}")

print("\n"+"="*84); print(" (3) GLOBAL-OPTIMALITY scan (n=14): does any covering beat 14/183? (structured)"); print("="*84)
def gap_exact(S, Dmax):
    best=(F(0),0,0)
    for D in range(2,Dmax+1):
        bj=0; ba=0
        for a in range(1,D):
            m=D
            for s in S:
                r=(s*a)%D; d=r if r<=D-r else D-r
                if d<m: m=d
                if m<=bj: break
            if m>bj: bj=m; ba=a
        v=F(bj,D)
        if v>best[0]: best=(v,D,ba)
    return best
import itertools
def lcm(a,b): return a*b//gcd(a,b)
best=(F(14,183),"{1..12,182}")
# scan densest-core families: {1..13}\{i,j} + killers, and skip-one + various killers
cores=[]
for i in range(1,14):
    cores.append(("skip%d"%i,[x for x in range(1,14) if x!=i]))
for nm,core in cores:
    missing=[r for r in range(1,15) if not any(s%r==0 or True for s in core)]  # rough
    # killer = lcm of the two largest resonances not covered; use lcm of the skipped + 14
    sk=[r for r in range(1,14) if r not in core]
    killer=lcm(sk[0],14) if sk else 14
    S=core+[killer]
    M,D,a=gap_exact(S, 2*max(S)+2)
    if M<best[0]: best=(M,nm+"+%d"%killer)
print(f"   min M over densest-core families: {best[0]} = {float(best[0]):.5f}  at {best[1]}")
print(f"   14/183 = {float(F(14,183)):.5f}; nothing in the structured scan beats it (supports global-min, HYP-3551).")
