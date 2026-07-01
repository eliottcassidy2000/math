"""
Chase Fourier-positivity + additive-energy lens.
 (A) danger sum D(t)=sum_v f(vt), f=1_{||.||<=1/n}; min_t D(t) = 1 iff tight (covering), <1 iff a hole.
 (B) additive energy A(S)=#{a+b=c+d in S^4}=sum_s r(s)^2. Does TIGHT correlate with HIGH A (=> Freiman => AP)?
"""
from fractions import Fraction
from itertools import combinations
def frac(x): x=x%1.0; return min(x,1-x)
def Dmin(S,n,G): return min(sum(1 for v in S if frac(v*i/G)<=1.0/n+1e-12) for i in range(G))
def Mf(S,Q):
    b=0.0
    for q in range(1,Q+1):
        for a in range(1,q):
            m=min(frac(s*a/q) for s in S)
            if m>b: b=m
    return b
def AE(S):
    from collections import Counter
    r=Counter()
    for a in S:
        for b in S: r[a+b]+=1
    return sum(c*c for c in r.values())
n=14; AP=list(range(1,n)); G=100*n
print(f"n={n}. (A) danger-sum min (=1 tight, <1 hole) and (B) additive energy A(S):")
sets={"AP {1..13}":AP, "GW {..11,13,24}":list(range(1,12))+[13,24],
      "swap 12->26 (non-tight)":[v for v in AP if v!=12]+[26],
      "swap 2->7 (n=5-style)":[v for v in AP if v!=2]+[7],
      "AP+1 {2..14}":list(range(2,n+1)), "random far":[1,3,7,8,12,15,17,19,20,23,25,27,30]}
Amax=AE(AP)
for name,S in sets.items():
    S=sorted(set(S))
    if len(S)!=n-1: 
        print(f"  {name}: |S|={len(S)} skip"); continue
    dm=Dmin(S,n,G); M=Mf(S,8*n); A=AE(S)
    print(f"  {name:>26}: Dmin={dm} (tight:{dm>=1}), M={M:.4f} (1/n={1/n:.4f}), A={A} (A/Amax={A/Amax:.3f})")
print()
print("=> if TIGHT sets all have A/Amax near 1 (high additive energy) and non-tight lower, Freiman gives the")
print("   lowness (high energy => AP-structure). Test correlation above.")
