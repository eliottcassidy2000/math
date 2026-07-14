#!/usr/bin/env python3
"""mac-mini-S91: explore the V4/Atkin-Lehner thread on the METAGRAPH side. THM-584: complement =
antipodal on Q_{C(n,2)} = ONE involution (W_14/Fricke, ι-even bulk = R-even Perron). The Atkin-Lehner
V_4 = <W_2,W_7> needs a SECOND involution. Q: does the tournament ISO-CLASS metagraph have a second
natural involution (giving V_4), or only Z_2 (complement)? Test the candidate second generators:
grid-transpose (reversal relabel i->n+1-i) and any score/duality involution."""
from itertools import combinations, permutations
def canon(adj,n):
    """canonical form of a labeled tournament (min over relabelings) -- iso-class key."""
    best=None
    for p in permutations(range(n)):
        bits=0
        for i in range(n):
            for j in range(n):
                if i<j and adj[p[i]][p[j]]: bits|=1<<(i*n+j)
        if best is None or bits<best: best=bits
    return best
def all_classes(n):
    edges=list(combinations(range(n),2)); m=len(edges); seen={}
    for mask in range(1<<m):
        adj=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(edges):
            if (mask>>k)&1: adj[i][j]=1
            else: adj[j][i]=1
        key=canon(adj,n)
        if key not in seen: seen[key]=adj
    return seen
def complement(adj,n): return [[adj[j][i] for j in range(n)] for i in range(n)]
def reversal(adj,n):  # grid-transpose = relabel i -> n-1-i
    return [[adj[n-1-i][n-1-j] for j in range(n)] for i in range(n)]
for n in [3,4,5]:
    cls=all_classes(n); keys=set(cls); N=len(keys)
    # complement action
    comp_fix=sum(1 for k in keys if canon(complement(cls[k],n),n)==k)
    # reversal action on iso-classes
    rev_fix=sum(1 for k in keys if canon(reversal(cls[k],n),n)==k)
    rev_is_id=all(canon(reversal(cls[k],n),n)==k for k in keys)
    # does reversal EVER differ from identity or complement on iso-classes?
    rev_new=any(canon(reversal(cls[k],n),n) not in (k,canon(complement(cls[k],n),n)) for k in keys)
    print(f"n={n}: #iso-classes A000568={N}; SC (complement-fixed)={comp_fix}; "
          f"reversal-fixed={rev_fix} (all? {rev_is_id})")
    print(f"       reversal gives a NEW iso-class involution (!= id, != complement)? {rev_new}")
print()
print("Atkin-Lehner V_4 on X_0(14): W_14=Fricke, W_2=2-adic, W_7=apex-7 (klein-S10 dictionary).")
print("Fricke W_14 <-> complement (THM-584, PROVED metagraph involution).")
print("SECOND generator test above: is there a tournament-iso-class involution besides complement?")
print()
# The tiling-model grid symmetry (BLUE) -- is it an iso-class operation or a labeled/base-path artifact?
print("NOTE: grid-transpose = reversal relabel i->n-1-i is IN S_n => a relabeling => acts as IDENTITY")
print("on iso-classes. So the 'grid-symmetric/BLUE' structure is a TILING-MODEL (fixed base path)")
print("phenomenon, INVISIBLE at the iso-class level. => the metagraph iso-class group is Z_2 (complement),")
print("NOT V_4. The 2nd AL generator (W_2/W_7) has NO tournament-iso-class realization.")
