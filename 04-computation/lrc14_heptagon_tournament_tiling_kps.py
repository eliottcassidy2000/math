#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE HEPTAGON BRIDGE: LRC atoms -> a p-gon -> a p-vertex tournament (tiling model) -> dihedral D_p.

kind-pasteur-2026-07-01-S9. The owner's cue: take the 6 unit atoms (Z/14)* at a/14, ADD the vertex
7/14, get 7 vertices; build a tournament on 7 vertices (vertex k <-> (2k+1)/14) and apply the tiling
model; consider tournaments and DIHEDRAL groups.  Findings tested here:

 (A) The 7 odd residues {1,3,5,7,9,11,13}/14 are EQUALLY SPACED by 1/7 = a REGULAR HEPTAGON; its
     symmetry is the DIHEDRAL group D_7 of ORDER 14 (= the LRC modulus).  vertex k <-> (2k+1)/14.
 (B) Two group actions on the 7 vertices: ADDITIVE heptagon rotation C_7 (k->k+1 = +1/7) with
     reflections = D_7 (geometry); MULTIPLICATIVE (Z/14)*=C_6 (fixes 7/14=1/2, permutes the 6 units).
     Inversion iota: t->-t = the reflection fixing the center 7/14 = TOURNAMENT COMPLEMENTATION T->T^op.
 (C) Tournaments on 7 vertices via the circle: ROTATIONAL R=C_7(1,2,3) (geometric forward-arc) and
     PALEY P=C_7(1,2,4) (QR mod 7, the H-MAXIMIZER).  TILING MODEL: fix base path 0-1-..-6; the 15
     tiles = gaps d>=2, forward iff d in the connection set.  Compute H, c3, scores, |Aut|, self-converse.
 (D) THE GENERAL BRIDGE (make the finish concrete): a full-orbit census class (Z/N)* with N=2p is the
     p-GON; the census infimum over full-orbit classes = an optimization over p-gon/D_p tournaments.
     N=6 (triangle p=3), N=10 (PENTAGON p=5 = the EXTREMIZER), N=14 (heptagon p=7 = TIGHT, meas 0).
"""
import sys, itertools
from fractions import Fraction as Fr
from math import gcd
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def units(N): return [a for a in range(1,N) if gcd(a,N)==1]

print("="*92); print(" (A)(B) THE 7 ODD POINTS = REGULAR HEPTAGON, symmetry D_7 (order 14)"); print("="*92)
pts=[Fr(2*k+1,14) for k in range(7)]
gaps=[pts[k+1]-pts[k] for k in range(6)]
print(f"  vertices k=0..6 -> (2k+1)/14 = {[str(p) for p in pts]}")
print(f"  consecutive gaps = {[str(g) for g in gaps]} all = 1/7? {all(g==Fr(1,7) for g in gaps)}  => REGULAR HEPTAGON")
U=units(14); center=7
print(f"  units (Z/14)* = {U} (6 units) + center 7 (=7/14=1/2, the added vertex, NON-unit)")
# multiplicative C_6 orbit of 1 under x3
orb=[1]; x=1
for _ in range(6):
    x=(x*3)%14; orb.append(x)
print(f"  MULT (Z/14)* = <x3> orbit: {orb[:6]} (C_6, fixes 7); ADD heptagon rotation = k->k+1 (+1/7).")
inv={a:next(b for b in U if (a*b)%14==1) for a in U}
print(f"  inversion iota (t->-t <=> a->14-a): 1<->13,3<->11,5<->9,7 fixed = reflection; a^-1 mod14: {inv}")

def circ(S7):  # circulant tournament C_7(S): A[i][j]=1 iff (j-i)%7 in S
    return [[1 if (j-i)%7 in S else 0 for j in range(7)] for i in range(7)]
def ham_paths(A,n=7):
    c=0
    for p in itertools.permutations(range(n)):
        if all(A[p[i]][p[i+1]] for i in range(n-1)): c+=1
    return c
def count3(A,n=7):
    c=0
    for i,j,k in itertools.combinations(range(n),3):
        # directed 3-cycle among i,j,k
        for (a,b,d) in [(i,j,k)]:
            pass
        # check both cyclic orientations
        if A[i][j] and A[j][k] and A[k][i]: c+=1
        if A[i][k] and A[k][j] and A[j][i]: c+=1
    return c
def scores(A,n=7): return sorted(sum(A[i]) for i in range(n))
def aut_order(A,n=7):
    c=0
    for p in itertools.permutations(range(n)):
        if all(A[i][j]==A[p[i]][p[j]] for i in range(n) for j in range(n)): c+=1
    return c
def self_converse(A,n=7):
    # exists permutation p with A[p[i]][p[j]] = A[j][i] (reverse)?
    for p in itertools.permutations(range(n)):
        if all(A[p[i]][p[j]]==A[j][i] for i in range(n) for j in range(n)): return True
    return False
def tiling_forward(S):  # base path 0-1-..-6 needs 1 in S; tiles = gaps d=2..6, forward iff d in S (mod7)
    return [d for d in range(2,7) if d in S], [d for d in range(2,7) if d not in S]

print("\n"+"="*92); print(" (C) TOURNAMENTS ON 7 VERTICES + TILING MODEL"); print("="*92)
for name,S in [("ROTATIONAL C_7(1,2,3)",{1,2,3}), ("PALEY C_7(1,2,4)=QR7",{1,2,4})]:
    A=circ(S); fwd,bwd=tiling_forward(S)
    print(f"  {name}: scores={scores(A)} (regular? {len(set(scores(A)))==1}); H(Ham paths)={ham_paths(A)}; "
          f"c3={count3(A)}; |Aut|={aut_order(A)}; self-converse(=iota-inv/complement)? {self_converse(A)}")
    print(f"     base path 0..6 needs 1 in S: {1 in S}; TILING: forward tiles gaps {fwd}, backward gaps {bwd} (15 tiles total across the staircase)")
print("  Paley_7 H=189 is the n=7 MAXIMUM (project canon). Both are SELF-CONVERSE => iota-invariant => the")
print("  dihedral reflection = tournament complementation T->T^op (the merged-metagraph Z_2, CLAUDE.md).")

print("\n"+"="*92); print(" (D) THE GENERAL BRIDGE: census class (Z/N)*, N=2p  <->  regular p-GON  <->  D_p  <->  p-vertex tournament"); print("="*92)
# measures from the census (HYP-3793 / S8, exact):
meas_by_N={6:Fr(1546,35035), 10:Fr(313,9702), 14:Fr(0,1)}
print(f"  {'N=2p':>5} {'p':>3} {'poly':>9} {'D_p order':>9} {'#units':>7} {'min meas(L_C)':>14} {'rot-tourn H(p)':>13}")
for N,p in [(6,3),(10,5),(14,7)]:
    poly={3:'triangle',5:'pentagon',7:'heptagon'}[p]
    # rotational tournament on p vertices C_p(1..(p-1)/2)
    Sp=set(range(1,(p-1)//2+1))
    Ap=[[1 if (j-i)%p in Sp else 0 for j in range(p)] for i in range(p)]
    Hp=0
    for perm in itertools.permutations(range(p)):
        if all(Ap[perm[i]][perm[i+1]] for i in range(p-1)): Hp+=1
    nu=len([a for a in units(N) if a!=N//2])  # units excluding center if present
    nu=len(units(N))
    m=meas_by_N[N]
    print(f"  {N:>5} {p:>3} {poly:>9} {2*p:>9} {nu:>7} {float(m):>14.5f} {Hp:>13}")
print("  => N=10 = PENTAGON (p=5) is the census EXTREMIZER (min meas 0.03226); N=14 = HEPTAGON (p=7) is")
print("     the TIGHT reference (meas 0, the AP {1..13}); N=6 = TRIANGLE (p=3) is looser (0.044).")
print("  The full-orbit census family = the REGULAR p-GONS; the two-clash family = the SPORADIC (non-polygon,")
print("  2-atom) sets = the Goddyn-Wong analog. Census infimum over full-orbit = optimization over p-gons,")
print("  attained at the PENTAGON. This is the tournament/dihedral image of the {AP, GW} dichotomy (THM-523).")
print("DONE.")
