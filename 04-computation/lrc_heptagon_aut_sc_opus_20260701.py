"""Aut groups, self-complementary (SC) structure, iso-class placement, and the LRC t=1/14 connection."""
import numpy as np
from itertools import permutations, combinations
V=list(range(7))
def rot_tour(conn):
    A=np.zeros((7,7),int)
    for i in V:
        for j in V:
            if i!=j and ((j-i)%7) in conn: A[i,j]=1
    return A
def perm_apply(A,p):
    B=np.zeros((7,7),int)
    for i in V:
        for j in V: B[p[i],p[j]]=A[i,j]
    return B
def aut_and_sc(A):
    Aop=np.where(np.eye(7,dtype=int)==1,0,1-A)  # converse/complement (reverse arcs)
    auts=[]; scs=[]
    for p in permutations(V):
        B=perm_apply(A,p)
        if np.array_equal(B,A): auts.append(p)
        if np.array_equal(B,Aop): scs.append(p)
    return auts,scs
for name,conn in [("ROTATIONAL {1,2,3}",{1,2,3}),("PALEY/QR {1,2,4}",{1,2,4})]:
    A=rot_tour(conn); auts,scs=aut_and_sc(A)
    print(f"{name}: |Aut|={len(auts)} (automorphisms), |SC-witnesses|={len(scs)} (converse-isos)")
    print(f"   self-complementary (SC = T ~ T^op)? {len(scs)>0}   full 'dihedral' group |Aut|+|SC| = {len(auts)+len(scs)}")
# the iota witness for rotational: k->6-k
iota=[6-k for k in V]; A=rot_tour({1,2,3}); Aop=np.where(np.eye(7,dtype=int)==1,0,1-A)
print(f"\n   iota=(k->6-k) is an SC-witness for rotational? {np.array_equal(perm_apply(A,iota),Aop)}")
print(f"   => Aut=C_7 (rotations preserve), reflections REVERSE (anti-aut) => D_7 acts, tournament is SC. On the repo's SC SPINE.")
# ---- LRC t=1/14 connection: the 13 runners split by parity ----
print("\n== LRC connection at t=1/14 (the AP's lonely time) ==")
runners=list(range(1,14))
odd=[v for v in runners if v%2==1]; even=[v for v in runners if v%2==0]
print(f"  13 runners v/14: ODD v={odd} -> odd/14 = the 7 heptagon vertices minus... units={[v for v in odd if np.gcd(v,14)==1]} (6 units) + v=7 -> 1/2 (antipode)")
print(f"  EVEN v={even} -> v/14 = {[f'{v//2}/7' for v in even]} = the non-trivial 7th roots of unity (the C_7 orbit)")
print(f"  So at t=1/14: ODD runners land on roots of z^7=-1 (the tournament vertices), EVEN runners on roots of z^7=+1 (7th roots).")
print(f"  M=1/14 is set by runners 1 & 13 = vertices 0 & 6 = the EXTREME iota-pair (nearest odd/14 points to origin 0).")
# gap structure: distances from origin 0 to the 7 vertices
d=sorted(min((2*k+1)/14,1-(2*k+1)/14) for k in V)
print(f"  distances origin->vertices: {[f'{x:.4f}' for x in d]}  min=1/14={1/14:.4f} (the gap M), max=1/2 (antipode)")
