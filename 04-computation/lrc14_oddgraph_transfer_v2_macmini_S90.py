#!/usr/bin/env python3
"""mac-mini-S90b: corrected build test. c3 = C(m,3)-sum C(s_i,2) (standard); H = Ham-path count (OCF
kin, odd by Redei). Q: does an ι-odd tournament invariant (c3 / H / score-shape) of the tight-time
difference tournament SEPARATE tight (AP,GW; L=0) from covering (deep well, near-AP; L>0), where the
ι-even Dedekind sum s(n,Phi6) is BLIND (HYP-3768)? And do any match f_14 (curve 14a) arithmetic?"""
from math import comb, gcd
from itertools import combinations
from functools import lru_cache
def diff_tour(S,num,den):
    m=len(S); adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i==j: continue
            fr=((S[i]-S[j])*num % den)/den
            adj[i][j]=1 if 0<fr<0.5 else 0
    return adj
def scores(adj): return [sum(r) for r in adj]
def c3cyc(adj):
    m=len(adj); s=scores(adj); return comb(m,3)-sum(comb(x,2) for x in s)
def Hpaths(adj):
    m=len(adj); full=(1<<m)-1
    @lru_cache(maxsize=None)
    def f(mask,last):
        if mask==full: return 1
        return sum(f(mask|(1<<nx),nx) for nx in range(m) if not (mask>>nx)&1 and adj[last][nx])
    r=sum(f(1<<s,s) for s in range(m)); f.cache_clear(); return r
fams=[("AP {1..13}          TIGHT(L=0)",list(range(1,14)),1,14),
      ("GW {1..11,13,24}    TIGHT(L=0)",[*range(1,12),13,24],1,14),
      ("deep well           COVER(L>0)",list(range(1,13))+[182],14,183),
      ("near-AP {1..11,13,84}COVER(L>0)",[*range(1,12),13,84],37,89),
      ("{1..14}\\{6} drop-6  COVER(L>0)",[v for v in range(1,15) if v!=6],1,14)]
print(f"{'family':34s} {'sortedscores':22s} {'c3cyc':>6} {'H(odd)':>10} {'scorevar':>9}")
res=[]
for name,S,num,den in fams:
    adj=diff_tour(S,num,den); sc=sorted(scores(adj)); m=len(S)
    cyc=c3cyc(adj); H=Hpaths(adj); var=sum((x-(m-1)/2)**2 for x in scores(adj))/m
    res.append((name,cyc,H,var)); print(f"{name:34s} {str(sc):22s} {cyc:>6} {H:>10} {var:>9.3f}")
print("\nSEPARATION verdict:")
tight=[r for r in res if 'TIGHT' in r[0]]; cover=[r for r in res if 'COVER' in r[0]]
print(f"  tight c3: {[r[1] for r in tight]}  vs  covering c3: {[r[1] for r in cover]}")
print(f"  tight H : {[r[2] for r in tight]}  vs  covering H : {[r[2] for r in cover]}")
mc=max(r[1] for r in tight); nc=min(r[1] for r in cover)
print(f"  c3 separates tight<covering cleanly? {mc < nc}  (max tight {mc} vs min cover {nc})")
# numerology vs f_14
print("\nf_14 / 14a numerology cross-check:")
print(f"  (Z/14)* order = {len([a for a in range(1,14) if gcd(a,14)==1])} = 6 = torsion(14a)=Z/6  [elementary phi(14)=6]")
print(f"  X_0(14) cusps = Klein-four V_4 = Atkin-Lehner <w_2,w_7> (14=2*7); tournament merged fold = Z_2 (complement) -- V_4 needs a 2nd fold (transpose)")
print(f"  H(AP diff-tour) = {tight[0][2]}; is it ~ a modular number? 14a: a_3=-2,a_13=-4; |disc(14a)|=14^? ; L(14a,1)!=0 (rank0)")
