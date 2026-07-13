#!/usr/bin/env python3
"""cont.56: the residues at a clean base q form a ROTATIONAL tournament. Explore the crux
families at their clean bases: AP {1..13} at q=14 (heptagon D7?), deep well {1..12,182} at
q=183. Does the tournament structure reveal WHY these are extremal (SC/resonant)?"""
from fractions import Fraction as F
from math import gcd
def residues(v,a,q): return [(x*a)%q for x in v]
def rot_tournament_info(res, q):
    # rotational tournament on Z/q (odd part): i beats j if (r_i - r_j) mod q in (0, q/2)
    # count cyclic triangles among the residue set
    n=len(res); wins={}
    for i in range(n):
        for j in range(n):
            if i!=j:
                d=(res[i]-res[j])%q
                wins[(i,j)] = 0 < d < q/2
    tri=0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                # cyclic if the 3 edges form a 3-cycle
                es=[wins.get((i,j)),wins.get((j,k)),wins.get((k,i))]
                if all(x is not None for x in es):
                    s=sum(1 for x in es if x)
                    if s in (0,3) and es[0]==es[1]==es[2]: tri+=1
    return tri
def M(v,Q=400):
    best=F(0); ba=None; bq=None
    for q in range(2,Q):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((x*a)%q, q-((x*a)%q)) for x in v)
            mm=F(m,q)
            if mm>best: best=mm; ba=a; bq=q
    return best,ba,bq
print("crux families: clean base q*, residues, rotational tournament structure:")
for nm,v in [("AP {1..13}",list(range(1,14))),
             ("V* {1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24]),
             ("deep well {1..12,182}",list(range(1,13))+[182])]:
    m,a,q=M(v, 200 if max(v)>50 else 30)
    res=residues(v,a,q)
    print(f"\n  {nm}: M={m}={float(m):.5f} at base q*={q}, a={a}")
    print(f"    residues mod {q}: {sorted(res)}")
    # the danger band is (0, q/14) and (13q/14, q); safe = residues avoid it
    band = q/F(14)
    print(f"    danger band radius q/14 = {float(band):.2f}; min |residue from 0| = {min(min(r,q-r) for r in res)}")
    # odd-part reduction (CRT 14=2*7 => mod 7)
    if q%7==0:
        res7=sorted(set(r%7 for r in res))
        print(f"    residues mod 7 (apex prime): {res7}  ({'ALL 7 classes' if len(res7)==7 else '%d classes'%len(res7)})")
print("\nOstrowski archimedean-finite bridge: M_k = k/(13k+1) = finite base q=13k+1 giving archimedean margin k/q:")
for k in [1,2,3,14]:
    q=13*k+1; print(f"  k={k}: base q={q}, margin k/q = {k}/{q} = {float(F(k,q)):.5f}  ({'AP/heptagon' if k==1 else 'deep well' if k==14 else ''})")
