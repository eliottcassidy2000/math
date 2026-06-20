#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quick per-subset means (no full brute), flushed line-by-line, to cross-check
the claimed E[c5] and E[alpha2-tri] closed forms fast."""
import sys
from itertools import combinations, product, permutations
from fractions import Fraction
from math import comb
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

def prob_dcycle(subset):
    sub = tuple(sorted(subset))
    pairs = list(combinations(sorted(sub, reverse=True), 2))  # (u,v) u>v
    free = [p for p in pairs if p[0]-p[1] != 1]
    forced = [p for p in pairs if p[0]-p[1] == 1]
    k = len(sub); s0 = sub[0]; rest = sub[1:]
    total = Fraction(0)
    for bits in product((0,1), repeat=len(free)):
        adj = {}
        for p in forced:
            adj[(p[0],p[1])] = 1; adj[(p[1],p[0])] = 0
        for idx,p in enumerate(free):
            if bits[idx]==0:
                adj[(p[0],p[1])]=1; adj[(p[1],p[0])]=0
            else:
                adj[(p[0],p[1])]=0; adj[(p[1],p[0])]=1
        c = 0
        for perm in permutations(rest):
            seq=(s0,)+perm; ok=True
            for i in range(k):
                if adj.get((seq[i],seq[(i+1)%k]),0)==0: ok=False;break
            if ok: c+=1
        total += c
    return total/(1<<len(free))

def Ec5(n):
    return sum(prob_dcycle(s) for s in combinations(range(1,n+1),5))

def Ec3(n):
    return sum(prob_dcycle(s) for s in combinations(range(1,n+1),3))

def Ea2tri(n):
    trips=list(combinations(range(1,n+1),3))
    pc={t:prob_dcycle(t) for t in trips}
    tot=Fraction(0)
    for i in range(len(trips)):
        for j in range(i+1,len(trips)):
            if set(trips[i]).isdisjoint(trips[j]):
                tot+=pc[trips[i]]*pc[trips[j]]
    return tot

def Ec5_formula(n):
    n=Fraction(n)
    return (Fraction(1,160)*n**5 - Fraction(1,16)*n**4 + Fraction(9,32)*n**3
            - Fraction(7,8)*n**2 + Fraction(147,80)*n - Fraction(7,4))

CLAIM_Ec5={5:Fraction(19,16),6:Fraction(49,8),7:Fraction(315,16),8:Fraction(199,4),
           9:Fraction(1727,16),10:Fraction(1683,8),11:Fraction(6055,16),12:Fraction(1279,2)}
CLAIM_Ea2={3:Fraction(0),4:Fraction(0),5:Fraction(0),6:Fraction(15,16),
           7:Fraction(93,16),8:Fraction(173,8),9:Fraction(495,8),10:Fraction(2395,16)}

print("n | Ec3_persub  Ec3_formula | Ec5_persub  formula  claim | match", flush=True)
for n in range(3,13):
    ec3=Ec3(n); ec3f=Fraction(comb(n,3)+(n-2),4)
    if n>=5:
        ec5=Ec5(n); ec5f=Ec5_formula(n); cl=CLAIM_Ec5[n]
        m = (ec5==ec5f==cl)
        print(f"{n:2d} | {str(ec3):>10} {str(ec3f):>10} | {str(ec5):>10} {str(ec5f):>10} {str(cl):>10} | {'OK' if (m and ec3==ec3f) else 'MISMATCH'}", flush=True)
    else:
        print(f"{n:2d} | {str(ec3):>10} {str(ec3f):>10} |  (c5 n/a)                          | {'OK' if ec3==ec3f else 'MISMATCH'}", flush=True)
print(flush=True)
print("n | Ea2tri_persub  claim | match", flush=True)
for n in range(3,11):
    a=Ea2tri(n); cl=CLAIM_Ea2[n]
    print(f"{n:2d} | {str(a):>12} {str(cl):>12} | {'OK' if a==cl else 'MISMATCH'}", flush=True)
