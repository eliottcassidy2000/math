"""
CUSP-DEPTH LANDSCAPE (M dilation-invariant => function on PROJECTIVE speed-space).
Hunt: which integer 13-sets have the DEEPEST cusps (M near 1/14)? Is the AP dilation-class unique?
Any ABNORMAL deep cusps?  Plus: general APs, gcd structure, the '2 and 7' factor.
"""
import numpy as np, random
from math import gcd
from functools import reduce
def M_grid(S,Q=200000,n=14):
    t=np.arange(1,Q)/Q; f=np.ones(Q-1)
    for v in S: f=np.minimum(f,np.abs(((v*t+0.5)%1)-0.5))
    return f.max()
def g(S): return reduce(gcd,S)
print("deepest known: AP {1..13} & dilates, M=1/14=%.6f\n"%(1/14))
print("(A) general APs {a, a+d, ..., a+12d}: M and cusp depth:")
for a,d in [(1,1),(2,1),(1,2),(0,1),(1,3),(3,2),(5,3),(7,4),(1,7),(2,7),(1,14)]:
    S=[a+d*i for i in range(13)]
    if 0 in S or len(set(S))<13: continue
    print(f"   a={a},d={d}: {S[:4]}...{S[-1]}  M={M_grid(S):.6f}  gcd={g(S)}")
print("\n(B) random primitive 13-sets: any deep cusps (M<0.075)? min M found:")
random.seed(1); best=(1.0,None); deep=[]
for _ in range(300):
    S=random.sample(range(1,40),13)
    if g(S)!=1: continue
    M=M_grid(S,120000)
    if M<best[0]: best=(M,sorted(S))
    if M<0.075: deep.append((round(M,5),sorted(S)))
print(f"   min M over 300 random primitive sets: {best[0]:.6f} at {best[1]}")
print(f"   deep cusps (M<0.075): {len(deep)} found")
for M,S in sorted(deep)[:6]: print(f"      M={M}: {S}")
print("\n(C) near-AP perturbations: replace ONE speed of {1..13}, scan M (the cusp neighborhood):")
ap=list(range(1,14))
print(f"   drop k, add w -- deepest results:")
res=[]
for drop in range(1,14):
    for w in range(14,40):
        if w in ap: continue
        S=[x for x in ap if x!=drop]+[w]
        if len(set(S))<13: continue
        res.append((M_grid(S,120000),drop,w))
for M,drop,w in sorted(res)[:8]:
    print(f"      drop {drop:>2}, add {w:>2}: M={M:.6f}")
print("\n(D) the 2*7 structure: M of {1..13} vs scaling by 2, 7, 14:")
for lam in [1,2,7,14]:
    S=[lam*v for v in ap]; print(f"   {lam}*AP: M={M_grid(S):.6f}")
