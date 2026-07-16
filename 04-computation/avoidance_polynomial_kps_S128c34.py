#!/usr/bin/env python3
"""avoidance_polynomial_kps_S128c34.py -- kind-pasteur S128 cont.34.
The avoidance polynomial of the tight scale-one system: speeds {1..13}, lambda = 1/14.
p(t) = sum_{k=0}^{13} (-1)^k S_k t^k, S_k = sum over k-subsets of mu(intersection of bad sets).
Exact integer arithmetic on the grid D = 14*lcm(1..13) = 14*360360.
Verdicts: p(1) = 0 identity; Bonferroni partial sums (death level); Newton table; roots."""
import sys
from fractions import Fraction as F
from itertools import combinations
from math import comb
sys.stdout.reconfigure(line_buffering=True)
L=360360; D=14*L
speeds=list(range(1,14))
# bad arcs of speed v: ((14j-1)/(14v), (14j+1)/(14v)) scaled by D -> ((14j-1)*L/v, (14j+1)*L/v)
arcs={}
for v in speeds:
    a=[]
    w=L//v
    for j in range(v):
        lo=(14*j-1)*w; hi=(14*j+1)*w
        if lo<0: a.append((0,hi)); a.append(((14*v-1)*w,14*v*w//14*14//14*0+D))  # wrap: split
        else: a.append((lo,hi))
    # fix wrap cleanly
    a=[]
    for j in range(v):
        lo=(14*j-1)*w; hi=(14*j+1)*w
        if lo<0:
            a.append((0,hi)); a.append((D+lo,D))
        else:
            a.append((lo,hi))
    arcs[v]=a
# sanity: mu(B_v) = 2/14 = 1/7 -> total length D/7
for v in speeds:
    tot=sum(h-l for l,h in arcs[v])
    assert tot==D//7,(v,tot,D//7)
S=[0]*14  # S[k] * ? : store integer measure on grid D (mu = S/D)
S[0]=D
def inter_measure(A):
    # measure (in grid units) of intersection of the unions arcs[v], v in A
    k=len(A)
    ev=[]
    for v in A:
        for l,h in arcs[v]:
            ev.append((l,1)); ev.append((h,-1))
    ev.sort()
    depth=0; last=0; tot=0
    for x,d in ev:
        if depth==k: tot+=x-last
        depth+=d
        if depth==k: last=x
    return tot
for k in range(1,14):
    s=0
    for A in combinations(speeds,k):
        s+=inter_measure(A)
    S[k]=s
    print("S_%d = %d/%d = %.6f   (equid: %.6f)"%(k,s,D,s/D,comb(13,k)/7**k),flush=True)
# p(1) identity
p1=sum((-1)**k*S[k] for k in range(14))
print("p(1)*D =", p1, "-> p(1) =", F(p1,D), "(EXPECT 0: tight system)")
# Bonferroni partial sums (in units of D)
print("Bonferroni partial sums B_m = sum_{k<=m} (-1)^k S_k  (death level = first m with B_m <= 0 ... sign scan):")
run=0
for m in range(14):
    run+= (-1)**m * S[m]
    print("  B_%-2d = %12d /D = %+.6f"%(m,run,run/D))
# Newton / log-concavity on a_k = S_k / C(13,k)
print("Newton table a_k = S_k/C(13,k), test a_k^2 >= a_{k-1} a_{k+1}:")
a=[F(S[k],comb(13,k)) for k in range(14)]
for k in range(1,13):
    lhs=a[k]*a[k]; rhs=a[k-1]*a[k+1]
    print("  k=%2d: a_k^2 %s a_{k-1}a_{k+1}  (ratio %.4f)"%(k,">=" if lhs>=rhs else "< ",float(lhs/rhs) if rhs>0 else float('inf')))
# roots of p(t) (numeric, Durand-Kerner on exact-rational coeffs)
co=[(-1)**k*S[k]/D for k in range(14)]
n=13
ws=[complex(0.4,0.9)**i*2 for i in range(n)]
for _ in range(2000):
    new=[]
    for i,wi in enumerate(ws):
        num=sum(co[j]*wi**j for j in range(n+1)); den=co[-1]
        for j,wj in enumerate(ws):
            if j!=i: den*=(wi-wj)
        new.append(wi-num/den if den!=0 else wi)
    ws=new
print("roots of p(t):")
for w in sorted(ws,key=lambda z:(abs(z.imag)>1e-7,z.real)):
    print("   %.6f%+.6fi  |t|=%.4f"%(w.real,w.imag,abs(w)))
nreal=sum(1 for w in ws if abs(w.imag)<1e-7)
print("real roots: %d / 13 ; REAL-ROOTED: %s"%(nreal,nreal==13))
print("DONE")
