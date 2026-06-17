#!/usr/bin/env python3
"""
lrc14_perturbation_gap_structure_kps  (kind-pasteur, PROVE side)

GOAL: understand the STRUCTURE of the lonely set that opens when a tight config
(AP {1..13} or sporadic {1..11,13,24}) is perturbed by replacing speed e->w.

For a tight config C (L=0), the danger arcs of its 13 speeds cover [0,1) exactly
(every point is within 1/14 of an integer multiple of some speed). When we DROP
speed e, a set  G_e := [0,1) \ (union of danger arcs of C\{e})  opens up: the
points that ONLY speed e was protecting.  Then we ADD speed w, whose danger arcs
cover part of G_e.  The opened lonely measure is

    L(C\{e} ∪ {w}) = meas( G_e \ D_w )

where D_w is w's danger set. This is EXACT and is the structural object to study.

We compute G_e exactly (the "e-exclusive gap"), describe it as a union of arcs,
and then for each candidate w compute meas(G_e ∩ D_w) to see how much w re-covers.
The minimal opened L = meas(G_e) - max_w meas(G_e ∩ D_w).

This isolates the two competing pieces and explains 12->36, L=1/1260.
"""
from fractions import Fraction as F

# ---------- exact danger arcs of a single speed ----------
def danger(v):
    out=[]
    w=F(1,14*v)
    for k in range(v+1):
        lo=F(k,v)-w; hi=F(k,v)+w
        if lo<0: out += [(F(0),hi),(1+lo,F(1))]
        elif hi>1: out += [(lo,F(1)),(F(0),hi-1)]
        else: out.append((lo,hi))
    return [(a,b) for a,b in out if b>a]

def union(arcs):
    """merge list of (lo,hi) into sorted disjoint arcs"""
    arcs=sorted((a,b) for a,b in arcs if b>a)
    if not arcs: return []
    res=[]; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=max(ch,hi)
        else: res.append((cl,ch)); cl,ch=lo,hi
    res.append((cl,ch))
    return res

def total(arcs):
    return sum(b-a for a,b in union(arcs))

def complement(arcs):
    """complement of union(arcs) within [0,1)"""
    u=union(arcs); res=[]; prev=F(0)
    for lo,hi in u:
        if lo>prev: res.append((prev,lo))
        prev=max(prev,hi)
    if prev<1: res.append((prev,F(1)))
    return res

def intersect(A,B):
    """intersection of two unions-of-arcs (each must be sorted disjoint)"""
    A=union(A); B=union(B); res=[]; i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]); hi=min(A[i][1],B[j][1])
        if lo<hi: res.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return res

def subtract(A,B):
    """A minus B (both unions of arcs)"""
    return intersect(A, complement(B))

def L_exact(S):
    arcs=[]
    for v in set(S): arcs+=danger(v)
    return F(1)-total(arcs)

def gap_of_drop(C, e):
    """G_e = points covered ONLY by e among C (i.e. complement of union of others)"""
    others=[v for v in C if v!=e]
    cov_others=[]
    for v in others: cov_others+=danger(v)
    return complement(cov_others)   # union of arcs == the e-exclusive lonely gap

# ================================================================
TIGHT = {
    "AP {1..13}": list(range(1,14)),
    "sporadic {1..11,13,24}": list(range(1,12))+[13,24],
}

print("="*78)
for name,C in TIGHT.items():
    print(f"\n#### TIGHT CONFIG {name} : L = {L_exact(C)} (must be 0)")
    assert L_exact(C)==0, "not tight!"
    print(f"   For each dropped speed e, the e-exclusive gap G_e (arcs only e covered):")
    for e in C:
        Ge=gap_of_drop(C,e)
        mg=total(Ge)
        # number of arcs and their (representative) structure
        print(f"     drop e={e:3d}: meas(G_e)={str(mg):>10s} = {float(mg):.6e} ,  #arcs={len(Ge)}")
print("="*78)
