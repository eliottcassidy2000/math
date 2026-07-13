#!/usr/bin/env python3
"""mac-mini-S75: ATTACK the runner-1 positional bound on |core|=1 covering families.
For |core|=1: coreCover = meas(G' n D_1)/meas(G'), D_1=[0,1/14)u(13/14,1] (runner 1's arc),
G'=safe set of the 12 smooth (30030-smooth) non-core runners. coreCover<1 <=> G' has a point in
the MIDDLE [1/14,13/14]. Find the TIGHTEST |core|=1 covering family (max coreCover), characterize
what keeps it <1, and test the hypothesis: coreCover -> 1 iff -> AP (the non-covering wall)."""
from math import gcd
LV=1.0/14; P=[2,3,5,7,11,13]
def cop(v): return all(v%p!=0 for p in P)
def is_cov(S,n=14): return all(any(v%q==0 for v in S) for q in range(2,n+1))
def prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1
def nn(x): x%=1.0; return min(x,1-x)
def core1(S): 
    c=[v for v in S if cop(v)]
    return len(c)==1 and c==[1]
def coreCover(S,res=300000):
    nc=[v for v in S if not cop(v)]
    inGp=0; cd=0; mid=0
    for j in range(res):
        t=(j+0.5)/res
        if all(nn(w*t)>=LV for w in nc):
            inGp+=1
            if nn(t)<LV: cd+=1     # runner 1 danger
            else: mid+=1
    if inGp==0: return None
    return cd/inGp, inGp/res, mid/res

def dist_to_AP13(S):
    # distance to {1..13} (the AP) as a covering-family perturbation proxy: how many elements differ
    A=set(range(1,14)); return len(set(S)^A)//2

print(f"1/14={LV:.5f}. TIGHTEST |core|=1 covering families (max coreCover):\n")
# enumerate |core|=1 covering families: {1} + 12 smooth runners covering 2..14
# smooth pool up to a cap; require the family covering + primitive + |core|=1
import itertools
smoothpool=[v for v in range(2,100) if not cop(v)]  # 30030-smooth, 2..99
cands=[]
# structured: {1..k} + smooth outliers (the extremal shapes)
base_families=[
 [*range(1,13),182],[*range(1,12),13,84],[*range(1,12),13,168],[*range(1,11),22,13,84],
 [*range(1,13),364],[*range(1,12),26,84],[*range(1,10),20,11,13,84],[*range(1,12),13,84*2],
 [*range(1,11),11,12,13,84],[*range(1,13),182*2],[*range(1,12),13,182],
]
results=[]
for S in base_families:
    S=sorted(set(S))
    if len(S)!=13 or not(prim(S) and is_cov(S) and core1(S)): continue
    r=coreCover(S)
    if r: results.append((r[0],S,r[1],dist_to_AP13(S)))
# random |core|=1 covering search
import random; random.seed(3); tried=0
while tried<3000:
    tried+=1
    S=sorted(set([1]+random.sample(smoothpool,12)))
    if len(S)!=13 or not(prim(S) and is_cov(S) and core1(S)): continue
    r=coreCover(S)
    if r and r[0]>0.7: results.append((r[0],S,r[1],dist_to_AP13(S)))
results.sort(reverse=True)
print(f"{'coreCover':>9} | meas(G') | dist-to-AP | family")
for cc,S,mg,d in results[:12]:
    print(f"{cc:9.4f} | {mg:.5f} | {d:10d} | {S}")
if results:
    print(f"\nMAX coreCover among |core|=1 covering = {results[0][0]:.4f} at {results[0][1]}")
    print(f"  (dist-to-AP {results[0][3]}, meas(G') {results[0][2]:.4f})")
print("\nATTACK SUMMARY: coreCover<1 for all (safe pt in middle exists). The max is bounded away")
print("from 1 by the covering constraint (>=1 mult of 14 forces a far element => G' not fully in")
print("D_1). But NO uniform-margin proof: coreCover->1 as the family ->AP (blanket), and covering")
print("only bars the EXACT AP (mult-of-14), not near-AP => the margin is the covering-min gap itself.")
