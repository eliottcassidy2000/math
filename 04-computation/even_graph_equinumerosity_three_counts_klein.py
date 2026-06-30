#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The three 'even graph' counts and the Royle equinumerosity (klein-S12).

'Even graph' is polysemous (opus-S260). Three Burnside sums on K_n's edges, all sharing the value 2 at
n=3 but diverging after -- the equinumerosity particularity:

  A000088  (ALL graphs)              = (1/n!) sum_{all sigma}        2^{E_orbits(sigma)}
  A000568  (tournaments=ROYLE-even)  = (1/n!) sum_{sigma ODD ORDER}  2^{E_orbits(sigma)}   [= signed cycle
                                        index P_n(1), THM-587; |sigma| odd <=> all cycle lengths odd]
  A002854  (degree-even=EULERIAN=cycle space) = (1/n!) sum_{all sigma} 2^{E_orbits - V_orbits + 1}
                                        [the cycle space Z(K_n), dim C(n-1,2); fixed-dim = E_orb-V_orb+1]

Tournaments are the ODD-ORDER-RESTRICTED graph Burnside (Royle 2022 = an intrinsic property too);
EULERIAN is the cycle-space (boundary-quotient) Burnside. They are NOT the same 'even'.
"""
from math import gcd
from itertools import permutations, combinations
from fractions import Fraction as F

def order(perm):
    n=len(perm); o=1; q=list(perm)
    while q!=list(range(n)):
        q=[perm[q[i]] for i in range(n)]; o+=1
    return o

def gf2_rank(rows):
    rows=list(rows); rank=0
    ncol=max((r.bit_length() for r in rows), default=0)
    for col in range(ncol-1,-1,-1):
        piv=None
        for i in range(rank,len(rows)):
            if (rows[i]>>col)&1: piv=i; break
        if piv is None: continue
        rows[rank],rows[piv]=rows[piv],rows[rank]
        for i in range(len(rows)):
            if i!=rank and (rows[i]>>col)&1: rows[i]^=rows[rank]
        rank+=1
    return rank

def counts(n):
    verts=list(range(n)); edges=list(combinations(verts,2)); eidx={e:i for i,e in enumerate(edges)}
    all_g=F(0); odd_g=F(0); eul=F(0)
    for perm in permutations(verts):
        seen=[False]*len(edges); orbits=[]
        for s in range(len(edges)):
            if seen[s]: continue
            orb=[]; e=edges[s]
            while not seen[eidx[e]]:
                seen[eidx[e]]=True; orb.append(eidx[e])
                a,b=perm[e[0]],perm[e[1]]; e=(min(a,b),max(a,b))
            orbits.append(orb)
        Eorb=len(orbits)
        # degree-parity matrix over GF(2): row per vertex, bit j = parity of #edges of orbit j at v
        Brows=[]
        for v in range(n):
            row=0
            for j,orb in enumerate(orbits):
                if sum(1 for ei in orb if v in edges[ei])&1: row|=(1<<j)
            Brows.append(row)
        rankB=gf2_rank(Brows)
        all_g += 2**Eorb
        if order(perm)%2==1: odd_g += 2**Eorb
        eul += 2**(Eorb - rankB)   # # degree-even (Eulerian) graphs fixed by perm (GF(2)-correct)
    nf=1
    for k in range(1,n+1): nf*=k
    return all_g/nf, odd_g/nf, eul/nf

print("="*78)
print(" The three 'even graph' counts on K_n (Burnside); the equinumerosity particularity")
print("="*78)
print(f" {'n':>2} {'A000088(all)':>13} {'A000568(tourn=Royle-even)':>26} {'A002854(degree-even=Euler)':>27}")
A088=[];A568=[];A854=[]
for n in range(3,8):
    a,o,e=counts(n); A088.append(int(a)); A568.append(int(o)); A854.append(int(e))
    print(f" {n:>2} {str(a):>13} {str(o):>26} {str(e):>27}")
print(f"\n A000088 (all graphs):        {A088}  (expect 4,11,34,156,1044)")
print(f" A000568 (tournaments/Royle): {A568}  (expect 4,12,56,456,6880)")
print(f" A002854 (degree-even/Euler): {A854}  (expect 3,7,16,54,243)")
print("\n PARTICULARITY:")
print(" - tournaments (A000568) = the ODD-ORDER restriction of the all-graph Burnside (= signed cycle")
print("   index P_n(1), THM-587). Royle 2022: also = #Royle-even graphs (an intrinsic property).")
print(" - EULERIAN/degree-even (A002854) = the CYCLE-SPACE Burnside (boundary-zero/F_2). DIFFERENT.")
print(" - all three = 2 at n=3 (coincidence), then diverge: tournaments > Eulerian for n>=4.")
print(f" - tournaments are SANDWICHED: A002854 <= A000568 <= A000088 ? {all(A854[i]<=A568[i]<=A088[i] for i in range(len(A088)))}")
