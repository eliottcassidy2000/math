#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
MERGE: the minimal-flip metagraph DEPTH of the heptagon tournaments to transitive = their MIN FEEDBACK ARC SET.

kind-pasteur-2026-07-01-S10. opus-S14/klein-S70 (HYP-3802) found "R_7 = transitive + 6 flipped TILES = 6 units"
for the base labeling 0->1->...->6.  In my minimal-flip frame (HYP-3803) the TRUE metagraph geodesic from a
class to the transitive class is the MINIMUM FEEDBACK ARC SET (min arcs to reverse to make acyclic = min back-
arcs over ALL vertex orderings).  Is opus's 6 the true minimum, or does a cleverer ordering beat it?  And is
it phi(14)=6 (the LRC atom/unit count)?  Compute MFAS(R_7), MFAS(Paley_7), the base-ordering count, the
R_7<->Paley flip distance, and confirm the |Aut|/H/triangle invariants.
"""
import sys, itertools
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
n=7
def circ(conn):  # i beats j iff (j-i)%7 in conn
    A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i!=j and (j-i)%n in conn: A[i][j]=1
    return A
R7=circ({1,2,3}); Paley=circ({1,2,4})
def arcs(A): return [(i,j) for i in range(n) for j in range(n) if A[i][j]]
def backarcs(A, order):  # order = tuple, position of each vertex; count arcs a->b with pos[a]>pos[b]
    pos={v:k for k,v in enumerate(order)}
    return sum(1 for (a,b) in arcs(A) if pos[a]>pos[b])
def MFAS(A):
    best=10**9; bestord=None
    for order in itertools.permutations(range(n)):
        b=backarcs(A,order)
        if b<best: best=b; bestord=order
    return best,bestord
def three_cycles(A):
    c=0
    for i,j,k in itertools.combinations(range(n),3):
        # count cyclic triangles among the 3
        e=[(i,j),(j,i),(i,k),(k,i),(j,k),(k,j)]
        # cyclic if i->j->k->i or i->k->j->i
        if A[i][j] and A[j][k] and A[k][i]: c+=1
        elif A[i][k] and A[k][j] and A[j][i]: c+=1
    return c
def ham_paths(A):
    c=0
    for p in itertools.permutations(range(n)):
        if all(A[p[t]][p[t+1]] for t in range(n-1)): c+=1
    return c
def aut(A):
    c=0
    for p in itertools.permutations(range(n)):
        if all(A[i][j]==A[p[i]][p[j]] for i in range(n) for j in range(n)): c+=1
    return c
def flip_dist(A,B):  # min Hamming over relabelings of B
    arcsetA=set(arcs(A)); best=10**9
    for p in itertools.permutations(range(n)):
        Bp=[[B[p[i]][p[j]] for j in range(n)] for i in range(n)]
        d=sum(1 for i in range(n) for j in range(n) if i<j and (A[i][j]!=Bp[i][j]))
        if d<best: best=d
    return best

base=tuple(range(n))
print("="*94); print(" HEPTAGON TOURNAMENTS: minimal-flip depth to transitive = MIN FEEDBACK ARC SET"); print("="*94)
for name,A,connlabel in [("R_7 rotational {1,2,3}",R7,"{1,2,3}"),("Paley_7 QR {1,2,4}",Paley,"{1,2,4}")]:
    mf,ordr=MFAS(A); bc=backarcs(A,base); tc=three_cycles(A); hp=ham_paths(A); au=aut(A)
    print(f"\n {name}:")
    print(f"   #3-cycles={tc}  Ham-paths(H)={hp}  |Aut|={au}")
    print(f"   base-ordering 0->..->6 back-arcs (opus's flipped tiles) = {bc}")
    print(f"   TRUE MIN FEEDBACK ARC SET (min back-arcs over all 7! orders) = {mf}  (best order {ordr})")
    print(f"   => opus's 6 is the BASE count; the metagraph geodesic to transitive is MFAS={mf}"
          + (f" = phi(14)=6 (atom count)!" if mf==6 else f" (NOT 6; base count 6 is not minimal)"))
print("\n"+"="*94); print(" R_7 <-> Paley_7 minimal flip distance, and to transitive"); print("="*94)
dRP=flip_dist(R7,Paley)
print(f"   flip_dist(R_7, Paley_7) = {dRP}  (metagraph geodesic between the two heptagon classes)")
print(f"   MFAS(R_7)={MFAS(R7)[0]}, MFAS(Paley_7)={MFAS(Paley)[0]}  (depth of each below transitive)")
print("\n INTERPRET: MFAS = the minimal-flip metagraph depth; opus's '6 tiles' is the base-labeling Hamming")
print("   weight, an UPPER bound on the geodesic. The true depth is the feedback-arc-set invariant.")
print("DONE.")
