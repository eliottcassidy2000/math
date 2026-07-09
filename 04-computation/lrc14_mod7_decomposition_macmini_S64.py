"""
mac-mini-2026-07-09-S64 -- CREATIVE ANGLE: the mod-7 decomposition of the 7-structured
dissociated good-period problem, connecting to THM-530's k<=7 pigeonhole.

STRUCTURAL FACT: at x = m/7, the phases {frac(e_i m/7)} = {(e_i*m mod 7)/7} all lie on the
7-point grid {0,1/7,...,6/7}. So maxgap(m/7):
  - if {e_i mod 7} MISSES a residue => an empty 1/7-slot => maxgap >= 2/7 > 1/7 (GOOD, big margin);
  - if the =0 mod 7 elements number |S7| >= k-6 => they COLLAPSE to 1 phase at m/7 => effective
    count k-|S7|+1 <= 7 => k<=7 pigeonhole => maxgap > 1/7 (THM-530 mechanism);
  - if {e_i mod 7} = Z/7 (all covered) => maxgap = 1/7 exactly (boundary).

SPLIT by gcd(7,Vmax):
  - 7 | Vmax: x=m/7 IS a grid point (j=m*Vmax/7); maxgap there governed by mod-7 coverage above.
  - gcd(7,Vmax)=1: grid decorrelates from the 1/7-resonances => samples the good bulk (mu high).

TEST: (1) over 7-structured dissociated sets, does {mod-7 misses a residue} OR {|S7|>=k-6} hold?
      (2) for 7|Vmax, is the grid point m*Vmax/7 good (mod-7 mechanism)?
      (3) for gcd(7,Vmax)=1, is existence robust?
"""
import numpy as np
from math import gcd
from functools import reduce
import random
random.seed(64)

def prim(E):
    E=sorted(E); return len(E)>=2 and reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
def longest_ap(E):
    S=set(E); best=2; E=sorted(E)
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i]; L=2; nx=E[j]+d
            while nx in S: L+=1; nx+=d
            bk=E[i]-d
            while bk in S: L+=1; bk-=d
            best=max(best,L)
    return best
def maxgap_int(E,j,V):
    ph=sorted({(e*j)%V for e in E}); m=len(ph)
    if m==1: return V
    mg=max((ph[(i+1)%m]-ph[i])%V for i in range(m-1)); return max(mg,ph[0]+V-ph[-1])
def good(E,j,V): return maxgap_int(E,j,V)*7>V
def has_good_period(E,V): return any(good(E,j,V) for j in range(1,V))
def make_7struct(k, s, min_s7):
    sevens=[x for x in range(0,s+1,7)]
    if len(sevens)<min_s7: return None
    S7=sorted(random.sample(sevens, random.randint(min_s7, min(len(sevens),k-1))))
    rest_pool=[x for x in range(1,s) if x%7!=0]
    need=k-len(S7)
    if need<0 or len(rest_pool)<need: return None
    E=sorted(set(S7+[0,s]+random.sample(rest_pool,max(0,need))))
    return E if len(E)==k else None

print("=== (1) structural dichotomy: 7-structured dissociated => misses-residue OR |S7|>=k-6? ===")
for k in (13,):
    n=0; dich_ok=0; both_fail=[]
    for _ in range(40000):
        s=random.randint(30,400); E=make_7struct(k,s,3)
        if E is None or not prim(E) or longest_ap(E)>k-6: continue
        res=set(e%7 for e in E); S7=sum(1 for e in E if e%7==0)
        n+=1
        misses = len(res)<7; collapse = S7>=k-6
        if misses or collapse: dich_ok+=1
        else: both_fail.append((tuple(E[:6]),s,len(res),S7))
    print(f"  k={k}: {n} 7-structured dissociated; dichotomy holds (misses-res OR |S7|>=k-6): "
          f"{dich_ok}/{n}; BOTH fail: {len(both_fail)}")
    if both_fail: print(f"    both-fail example (covers all 7 residues, |S7|<k-6): {both_fail[0]}")

print("\n=== (2) 7|Vmax: is grid point j=m*Vmax/7 good? (mod-7 mechanism) ===")
n=0; gp_ok=0
for _ in range(4000):
    s=random.randint(40,300); E=make_7struct(13,s,3)
    if E is None or not prim(E) or longest_ap(E)>7: continue
    V=7*((max(E)+7)//7 + random.randint(1,20))   # 7 | V, V>spread
    if V<=max(E): continue
    n+=1
    ok=any(good(E,(m*V//7)%V,V) for m in range(1,7) if (m*V//7)%V!=0)
    if ok: gp_ok+=1
print(f"  7|Vmax: {n} sets; grid point m*V/7 good for some m: {gp_ok}/{n}")

print("\n=== (3) gcd(7,Vmax)=1: existence robust? ===")
n=0; ex_ok=0
for _ in range(4000):
    s=random.randint(40,300); E=make_7struct(13,s,3)
    if E is None or not prim(E) or longest_ap(E)>7: continue
    V=max(E)+random.randint(1, max(1,max(E)//6))
    while V%7==0: V+=1
    if V<=max(E): continue
    n+=1
    if has_good_period(E,V): ex_ok+=1
print(f"  gcd(7,Vmax)=1: {n} sets; good period exists: {ex_ok}/{n}")
