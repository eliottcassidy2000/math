#!/usr/bin/env python3
"""relation_expansion_referee_kps_S128c36.py -- referee: (a) the relation expansion
mu(cap_A B) = sum_{h.v=0} prod c(h_i) (box-truncated) vs exact sweep; (b) the E_s lemma:
B5 = B5_eq + sum_{s=2}^{5} (-1)^s E_s * (relation mass at support s), E = [24/343, 24/49, -2/7, 1]."""
import sys, math
from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd, sin, pi
sys.stdout.reconfigure(line_buffering=True)
LAM=1/14
def c(h): return 2*LAM if h==0 else sin(2*pi*h*LAM)/(pi*h)
def sweep_mu(speeds):
    ev=[]
    for v in speeds:
        for j in range(v):
            lo=F(14*j-1,14*v); hi=F(14*j+1,14*v)
            if lo<0: ev.append((F(0),1)); ev.append((hi,-1)); ev.append((lo+1,1)); ev.append((F(1),-1))
            else: ev.append((lo,1)); ev.append((hi,-1))
    ev.sort()
    k=len(speeds); depth=0; last=F(0); tot=F(0)
    for x,d in ev:
        if depth==k: tot+=x-last
        depth+=d
        if depth==k: last=x
    return tot
print("== (a) expansion referee: subsets of the tight packet ==")
for A in [(1,2,3),(2,3,5),(1,2,3,4)]:
    exact=float(sweep_mu(list(A)))
    s=len(A); H=400
    # enumerate h with |h_i|<=H, h.v=0 by meet in middle on last coordinate
    tot=0.0
    if s==3:
        for h1 in range(-H,H+1):
            for h2 in range(-H,H+1):
                r=h1*A[0]+h2*A[1]
                if r% A[2]==0:
                    h3=-r//A[2]
                    if abs(h3)<=H: tot+=c(h1)*c(h2)*c(h3)
    else:
        for h1 in range(-H,H+1):
            for h2 in range(-H,H+1):
                for h3 in range(-H,H+1):
                    r=h1*A[0]+h2*A[1]+h3*A[2]
                    if r%A[3]==0:
                        h4=-r//A[3]
                        if abs(h4)<=H: tot+=c(h1)*c(h2)*c(h3)*c(h4)
    print("  A=%s: sweep %.8f  expansion(H=%d) %.8f  err %.1e"%(A,exact,H if s==3 else 60,tot,abs(exact-tot)))
print("== (b) E_s exact ==")
E={}
for s in range(2,6):
    Es=sum((-1)**j*comb(13-s,j)*F(1,7)**j for j in range(0,6-s))
    E[s]=Es
    print("  E_%d = %s = %.6f"%(s,Es,float(Es)))
B5eq=sum((-1)**k*comb(13,k)*F(1,7)**k for k in range(6))
print("  B5_equid = %s = %.6f  (ERRATUM: THM-930 quoted 0.0821; correct kill threshold = %s = %.3e; deep well 1/91 is %.0fx over)"%(
    B5eq,float(B5eq),B5eq/792,float(B5eq/792),(1/91)/float(B5eq/792)))
print("== (b2) E_s lemma referee on a small full system: speeds {1,2,3} at lambda=1/14, all levels ==")
# full B_m for {1,2,3}: check B_m = B_m^eq + sum over relations grouped by support with E_s^(m) weights
# (n=3 version: E_s^(m) = sum_{j<=m-s} (-1)^j C(3-s,j)(2lam)^j)
V=[1,2,3]
def depth_mu(speeds):
    ev=[]
    for v in speeds:
        for j in range(v):
            lo=F(14*j-1,14*v); hi=F(14*j+1,14*v)
            if lo<0: ev.append((F(0),1)); ev.append((hi,-1)); ev.append((lo+1,1)); ev.append((F(1),-1))
            else: ev.append((lo,1)); ev.append((hi,-1))
    ev.sort()
    mu=[F(0)]*(len(speeds)+1); depth=0; last=F(0)
    for x,d in ev:
        if x>last: mu[depth]+=x-last; last=x
        depth+=d
    if F(1)>last: mu[depth]+=F(1)-last
    return mu
mu=depth_mu(V)
S=[sum(comb(d,k)*mu[d] for d in range(4)) for k in range(4)]
H=3000
# relation masses by support (numeric)
mass={2:0.0,3:0.0}
for (i,j) in combinations(range(3),2):
    a,b=V[i],V[j]; g=gcd(a,b); pa,pb=b//g,a//g
    m=1
    while max(pa*m,pb*m)<=H:
        mass[2]+=2*c(pa*m)*c(pb*m)  # +-m
        m+=1
for h1 in range(-H,H+1):
    if h1==0: continue
    for h2 in range(-H,H+1):
        if h2==0: continue
        r=h1*1+h2*2
        if r%3==0:
            h3=-r//3
            if h3!=0 and abs(h3)<=H: mass[3]+=c(h1)*c(h2)*c(h3)
for m_ in range(1,4):
    Bm=float(sum((-1)**k*S[k] for k in range(m_+1)))
    Bmeq=float(sum((-1)**k*comb(3,k)*(1/7)**k for k in range(m_+1)))
    pred=Bmeq
    for s in (2,3):
        if s<=m_:
            Es=sum((-1)**j*comb(3-s,j)*(1/7)**j for j in range(0,m_-s+1))
            pred+=(+1)*Es*mass[s]   # sign: (-1)^s * (-1)^s from c-products? test empirically
    print("  m=%d: B_m exact %.8f  eq %.8f  eq+relmass %.8f"%(m_,Bm,Bmeq,pred))
print("DONE")
