#!/usr/bin/env python3
"""
WHY does M3 (section majority vote at N=14) forbid the regular tournament
(score 2,2,2,2,2)? Investigate the mechanism.

M3 arc i->j iff  sum_{a in U14} sign( d(v_i a) - d(v_j a) ) > 0,
  where d(x)=min(x mod N, N - x mod N) is the "section depth" (distance to 0
  on the residue circle), U14=(Z/14)*={1,3,5,9,11,13}.

Hypothesis: the per-runner total "section-depth profile" induces a near-total
order because d(v a) summed over units is an INVARIANT-like score that is
hard to make all-equal -> majority vote rarely produces a balanced regular
tournament. Check:
  1. For each speed v, compute its DEPTH VECTOR D(v) = (d(v*a mod 14))_{a in U}.
  2. M3 compares depth vectors by a "wins more columns" majority rule (a
     Copeland-style / pairwise-majority order on vectors).
  3. The regular tournament needs a 3-cycle pattern AND perfect score balance.
     Show that the depth-vector multiset over residues mod 14 has limited
     distinct values -> Condorcet structure constrains achievable tournaments.
"""
from math import gcd
from itertools import combinations, permutations
from functools import reduce

N=14
U=[a for a in range(1,N) if gcd(a,N)==1]  # {1,3,5,9,11,13}
def depth(x): r=x%N; return min(r, N-r)

# depth vector depends only on v mod 14 (since a runs over residues)
def dvec(v): return tuple(depth((v*a)%N) for a in U)

print("Depth vectors D(v) for v mod 14 (v=1..13), U =", U)
for v in range(1,N):
    print(f"  v={v:2d}: sections {[ (v*a)%N for a in U]}  depths {dvec(v)}  sum={sum(dvec(v))}")
print("  v=0 (multiple of 14): sections all 0 -> depth vector all 0 (the PARKED runner).")

# The M3 pairwise rule: i beats j iff #{a: d_i(a)>d_j(a)} > #{a: d_i(a)<d_j(a)}.
def beats(vi, vj):
    Di, Dj = dvec(vi), dvec(vj)
    t=sum((Di[c]>Dj[c]) - (Di[c]<Dj[c]) for c in range(len(U)))
    return t  # >0 i beats j, <0 j beats i, 0 tie

print("\nPairwise majority margins between distinct residues 1..13:")
print("    ", " ".join(f"{v:3d}" for v in range(1,N)))
for vi in range(1,N):
    row=[]
    for vj in range(1,N):
        if vi==vj: row.append("  .")
        else: row.append(f"{beats(vi,vj):+3d}")
    print(f" {vi:2d}:", " ".join(row))

# KEY: the depth vector is the SAME for v and 14-v? and for the (Z/14)* orbit?
print("\nResidues grouped by their depth-MULTISET (sorted depth vector):")
groups={}
for v in range(1,N):
    key=tuple(sorted(dvec(v)))
    groups.setdefault(key,[]).append(v)
for key,vs in sorted(groups.items()):
    print(f"  depths(sorted)={key}: residues {vs}")

# Now: can a regular tournament on 5 distinct residues exist? Need 5 residues
# pairwise giving a balanced (2-2) Condorcet structure with a 5-cycle of beats.
# Exhaustively check ALL 5-subsets of residues 1..13 (distinct sections) and
# also allow repeated residues (speeds can share residue mod 14).
def m3_tournament(residues):
    k=len(residues)
    adj=[[False]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i==j: continue
            t=beats(residues[i],residues[j])
            if t>0: adj[i][j]=True
            elif t<0: adj[i][j]=False
            else:
                # tie broken by speed in original; here break by residue index
                adj[i][j]= residues[i]<residues[j]
    return adj
def score(adj,k): return tuple(sorted(sum(adj[i][j] for j in range(k) if j!=i) for i in range(k)))
def is_tour(adj,k):
    return all(adj[i][j]!=adj[j][i] for i in range(k) for j in range(i+1,k))

print("\n### Can ANY 5 residues (distinct, mod 14) give regular score (2,2,2,2,2)? ###")
found=False
for combo in combinations(range(1,N),5):
    adj=m3_tournament(list(combo))
    if is_tour(adj,5) and score(adj,5)==(2,2,2,2,2):
        print("  FOUND distinct-residue regular:", combo); found=True
if not found: print("  NONE among distinct-residue 5-subsets (C(13,5)=1287).")

print("\n### Allow REPEATED residues (speeds sharing a section): all multisets ###")
from itertools import combinations_with_replacement as cwr
found=False; cnt=0
for combo in cwr(range(1,N),5):
    cnt+=1
    adj=m3_tournament(list(combo))
    if is_tour(adj,5) and score(adj,5)==(2,2,2,2,2):
        found=True
        # only print first few
        print("  regular from residue-multiset:", combo)
        break
if not found: print(f"  NONE among all {cnt} residue-multisets of size 5.")

print("\n### The deep reason: pairwise margin is a WEAK ORDER? check transitivity of 'beats' ###")
# If 'beats' (majority margin sign) on residues is transitive (a weak order),
# then EVERY M3 tournament is transitive -> contradicts nontriviality.
# So it must have SOME 3-cycles but they are structurally limited.
cyc=0; tot=0
for a,b,c in combinations(range(1,N),3):
    tot+=1
    ab=beats(a,b); bc=beats(b,c); ca=beats(c,a)
    # treat 0 as tie (skip)
    if ab and bc and ca:
        if (ab>0 and bc>0 and ca>0) or (ab<0 and bc<0 and ca<0):
            cyc+=1
print(f"  3-cycles in the residue 'beats' relation: {cyc} cyclic / {tot} triples (strict, no ties)")

# Show the SCC / Condorcet ranking of residues 1..13 by Copeland score
copeland={v: sum(1 for w in range(1,N) if w!=v and beats(v,w)>0) for v in range(1,N)}
print("\nCopeland scores (how many residues each beats by majority):")
for v in sorted(range(1,N), key=lambda x:-copeland[x]):
    print(f"   residue {v:2d}: beats {copeland[v]} others, depthsum={sum(dvec(v))}")
