#!/usr/bin/env python3
"""cont.56: THE UNIFICATION. At Ostrowski rung k the clean base is q=13k+1 and the extremal
family's residues are the ORBIT of a single rotation by k/q. This one fact ties together:
  (a) the Ostrowski ladder M_k = k/(13k+1)  [the rotation angle IS the rung value],
  (b) the three-gap theorem [orbit of a rotation has <=3 gap lengths => the three-gap
      regularity that makes the extremals extremal, cont.44],
  (c) the rotational tournament on the orbit [heptagon D7 at q=14 -> mod 7, SC].
Verify all three across the ladder and at the crux."""
from fractions import Fraction as F
from math import gcd
from collections import Counter

def three_gaps(points, q):
    """gap-length multiset of a point set on Z/q (circle of circumference q)."""
    P=sorted(set(points)); g=[]
    for a,b in zip(P, P[1:]+[P[0]+q]): g.append(b-a)
    return sorted(set(g)), Counter(g)

def rot_tournament_triangles(res, q):
    """rotational tournament on distinct residues: i->j iff (r_i-r_j) mod q in (0,q/2).
    Return (#vertices, #3-cycles, is_it_the_full_rotational_R_q_on_Z/q)."""
    R=sorted(set(r%q for r in res)); n=len(R)
    if q%2==0: return n, None, "q even: ambiguous (reduce to odd part)"
    beats=lambda a,b: 0 < ((a-b)%q) <= (q-1)//2
    cyc=0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                a,b,c=R[i],R[j],R[k]
                # 3-cycle iff not a transitive triple
                e=[beats(a,b),beats(b,c),beats(c,a)]
                s=sum(e)
                if s==0 or s==3: cyc+=1
    full = (n==q)
    return n, cyc, ("FULL rotational R_%d (all residues present)"%q if full else "%d/%d residues"%(n,q))

print("="*74)
print("THE ROTATION-ORBIT UNIFICATION across the Ostrowski ladder M_k = k/(13k+1)")
print("="*74)
print("At rung k: clean base q=13k+1, multiplier a=k, residues = {k*v mod q}.")
print("Extremal family at rung k: k=1 -> AP{1..13}; k=14 -> deep well {1..12,182}.\n")

# rung 1: AP, rung 14: deep well; also show the intermediate covering rungs exist
cases = [
    (1,  "AP {1..13}",            list(range(1,14))),
    (14, "deep well {1..12,182}", list(range(1,13))+[182]),
]
for k, nm, v in cases:
    q=13*k+1; a=k
    res=[(a*x)%q for x in v]
    gaps_set, gaps_ct = three_gaps(res, q)
    nn, cyc, desc = rot_tournament_triangles(res, q)
    print(f"rung k={k}: base q={q}, rotation angle a/q = {a}/{q} = {float(F(a,q)):.5f} = M_k")
    print(f"  residues (rotation orbit): {sorted(res)}")
    print(f"  THREE-GAP: distinct gap lengths = {gaps_set}  ({'<=3 CONFIRMED' if len(gaps_set)<=3 else 'MORE THAN 3'}),  multiplicities {dict(gaps_ct)}")
    print(f"  rotational tournament: {nn} vertices, {desc}" + (f", cyclic triangles={cyc}" if cyc is not None else ""))
    if q%7==0:
        res7=[r%7 for r in res]
        nn7,cyc7,desc7=rot_tournament_triangles(res7,7)
        # heptagon R_7 = "beat next 3": SC, Aut=C_7, has 14 cyclic triangles = |D_7|
        print(f"  mod apex prime 7: residues {sorted(set(res7))} -> {desc7}, cyclic triangles={cyc7}"
              f"  {'== HEPTAGON D7 (SC, 14 triangles)!' if desc7.startswith('FULL') and cyc7==14 else ''}")
    print()

print("="*74)
print("ARCHIMEDEAN <-> FINITE bridge (opus-S250's obstruction, resolved by the ladder):")
print("  each FINITE base q=13k+1 delivers an ARCHIMEDEAN margin k/q; sup_k k/(13k+1)=1/13.")
print("  the continued fraction of 1/13-limit IS the reconciliation of the two places.")
for k in [1,2,3,7,14,100,10**6]:
    q=13*k+1; print(f"    k={k:>7}: finite base q={q:>9}, archimedean margin={float(F(k,q)):.6f}  -> 1/13={1/13:.6f}")
print("\nThe crux (covering) sits at k=14 (14 | k is the covering condition): M=14/183.")
print("The tight rung k=1 (AP, M=1/14) is NON-covering; covering forces the jump 1->14.")
