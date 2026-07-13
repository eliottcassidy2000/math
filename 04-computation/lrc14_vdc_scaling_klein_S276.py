#!/usr/bin/env python3
"""
lrc14_vdc_scaling_klein_S276.py
===============================
klein-2026-07-13-S276 (owner: work on the van der Corput O(k) bound).

Claim (S275): at NON-RESONANT w, Error*w = |Sum_endpoints eps_p G_{s(p)}(wp)| = O(k), independent of
Sigma-e'. DECISIVE scaling test (fixed PRIME grid Ng>>w -- NO Ng∝w aliasing):
 (1) fixed k=7, GROW diameter (=> grow Sigma-e'), clean w~20*diam: does err*w grow with Sigma-e'?
     -> distinguishes O(1)/O(sqrt(sum))/O(sum).
 (2) fixed diameter, GROW k (# offsets 3..11): err*w vs k -> is it O(k) (linear) or bounded?
 (3) DIOPHANTINE: fixed cluster, w from clean(prime) to resonant(near lcm-mult): err*w vs
     delta=min_e' ||w/e'||.
"""
import math
from math import gcd
from functools import reduce
NG=999983  # fixed large PRIME >> all w; independent of w -> no aliasing
def gcdall(xs): return reduce(gcd,[x for x in xs if x],0)
def lcm(xs):
    r=1
    for x in xs:
        if x: r=r*x//gcd(r,x)
    return r
def isprime(n):
    if n<2: return False
    for p in range(2,int(n**0.5)+1):
        if n%p==0: return False
    return True
def nextprime(n):
    n=max(2,n)
    while not isprime(n): n+=1
    return n
def occ_of(E,x):
    o=0
    for e in E:o|=1<<(int((e*x%1.0)*7.0)%7)
    return o
def moments3(E):
    s1=s2=s3=0
    for k in range(1,NG):
        N=7-bin(occ_of(E,k/NG)).count("1"); s1+=N;s2+=N*(N-1);s3+=N*(N-1)*(N-2)
    n=NG-1;return s1/n,s2/n,s3/n
def Phi(E):m1,m2,m3=moments3(E);return 1-(2/3)*m1+(47/252)*m2-(5/252)*m3
def Phi_inf(C):m1,m2,m3=moments3(C);return 1-(2/3)*(6/7)*m1+(47/252)*(5/7)*m2-(5/252)*(4/7)*m3
def errw(C,w): return abs(Phi(C+[w])-Phi_inf(C))*w
def dioph(C,w):  # min_e' ||w/e'|| over nonzero e'
    d=1.0
    for e in C:
        if e:
            f=(w/e)%1.0; d=min(d, f, 1-f)
    return d

print("="*72)
print("(1) fixed k=7, grow diameter D (=> grow Sigma-e'), clean w=nextprime(20D)")
print("="*72)
print(f"  {'D':>4} {'sum':>5} {'w':>7} {'delta':>7} {'err*w':>8} {'/sqrt(sum)':>11} {'/sum':>8}")
for D in [6,12,24,48,96,192,384]:
    C=[0,1,2,3,4,5,D] if D>=6 else [0,1,2,3,4,5,6]
    w=nextprime(20*D); s=sum(C); ew=errw(C,w)
    print(f"  {D:4d} {s:5d} {w:7d} {dioph(C,w):7.3f} {ew:8.3f} {ew/math.sqrt(s):11.3f} {ew/s:8.4f}")

print()
print("="*72)
print("(2) fixed diameter D=120, GROW k (offsets spread over [0,120]), clean w=2411")
print("="*72)
w=2411  # prime, ~20*120
print(f"  {'k':>3} {'sum':>5} {'delta':>7} {'err*w':>8} {'err*w/k':>8}")
for k in [3,4,5,6,7,8,9,10,11]:
    C=sorted(set(round(i*120/(k-1)) for i in range(k)))
    if len(C)<k: continue
    ew=errw(C,w); s=sum(C)
    print(f"  {k:3d} {s:5d} {dioph(C,w):7.3f} {ew:8.3f} {ew/k:8.3f}")

print()
print("="*72)
print("(3) DIOPHANTINE: fixed C={0,1,2,28,29,30,15} (lcm=12180), w clean->resonant; err*w vs delta")
print("="*72)
C=[0,1,2,28,29,30,15]; L=lcm(C)
print(f"  cluster lcm={L}")
print(f"  {'w':>7} {'delta=min||w/e||':>16} {'err*w':>8}")
cands=[1009,2003,4001,  L//4, L//2, L-1, L, 2*L]
for w in sorted(set(x for x in cands if 0<x<NG//30)):
    print(f"  {w:7d} {dioph(C,w):16.4f} {errw(C,w):8.3f}")
print("\ndone.")
