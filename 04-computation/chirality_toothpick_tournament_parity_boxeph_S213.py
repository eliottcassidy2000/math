#!/usr/bin/env python3
"""chirality_toothpick_tournament_parity_boxeph_S213.py -- boxeph-2026-07-21-S213

Chirality, the toothpick sequence A139250, and why tournament values are EVEN -- one law:

  a count graded by a Z/2 INVOLUTION satisfies  count = (#fixed) + 2*(#free pairs),  so
  count  ==  (#fixed points of the involution)   (mod 2)   [a Lefschetz/Euler parity].

Three realizations of the SAME reversal involution:
  * TOURNAMENTS: R: T -> T^op (converse = complement = antipodal map, THM-584). THM-587 (PROVED):
    P_n(1)=A000568 (tournament count), P_n(-1)=SC(n)=self-converse count = antipodal EULER/LEFSCHETZ number.
    A000568 is EVEN for n>=3 because SC(n) is even. count=SC+2*(chiral pairs).
  * TOOTHPICK A139250: the diagram has mirror (D4) symmetry; A139250(n) = (#axis toothpicks) + 2*(#pairs),
    and #axis toothpicks is ODD (the central seed), so A139250(n) is ODD for all n>=1.
  * LRC (boxeph S212): iota: t->1-t; chi(G_delta) == (#iota-fixed lonely points) (mod 2).

Same Lefschetz parity law; the toothpick DIAGRAM is its geometric picture (one achiral seed + chiral pairs).
"""
from itertools import combinations, permutations

def sep(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)

# ==========================================================================
sep("A  TOOTHPICK A139250: simulate the CA; verify OEIS, oddness, seed-parity, 2^k Jacobsthal formula")
# unit toothpicks on the even integer lattice (length 2): (x,y,orient) with orient in {'H','V'}.
def endpoints(tp):
    x,y,o=tp
    return ((x,y-1),(x,y+1)) if o=='V' else ((x-1,y),(x+1,y))
def toothpick_seq(gens):
    tps={(0,0,'V')}; totals=[1]
    for _ in range(gens):
        # endpoint multiplicities
        mult={}
        for tp in tps:
            for e in endpoints(tp): mult[e]=mult.get(e,0)+1
        # exposed endpoints (multiplicity 1) sprout a perpendicular toothpick centered there
        new=set()
        for tp in tps:
            x,y,o=tp
            for e in endpoints(tp):
                if mult[e]==1:
                    perp='H' if o=='V' else 'V'
                    new.add((e[0],e[1],perp))
        new-=tps
        if not new: totals.append(len(tps)); break
        tps|=new; totals.append(len(tps))
    return tps, totals
tps,totals=toothpick_seq(16)
oeis=[1,3,7,11,15,23,35,43,47,55,67,79,95,123,155,171,175]
print("  simulated A139250 (totals):", totals[:17])
print("  OEIS A139250            :", oeis)
print("  match OEIS?", totals[:16]==oeis[:16])
diffs=[totals[0]]+[totals[i]-totals[i-1] for i in range(1,len(totals))]
print("  first differences (added per gen):", diffs[:16], " -> #odd terms:", sum(1 for d in diffs if d%2==1))
print("  A139250(n) all ODD (n>=0 total>=1)?", all(t%2==1 for t in totals))
# 2^k formula A139250(2^k) = (2^(2k+1)+1)/3  (Jacobsthal A007583). totals[g]=A139250(g+1), so A139250(2^k)=totals[2^k-1].
print("  A139250(2^k) vs (2^(2k+1)+1)/3:", [(totals[2**k-1], (2**(2*k+1)+1)//3) for k in range(0,5) if 2**k-1<len(totals)])
# axis toothpicks (fixed by mirror x->-x): centered on x=0 axis  => the ODD 'fixed' part
axis=sum(1 for (x,y,o) in tps if x==0)
print(f"  final gen: total toothpicks={len(tps)} (odd? {len(tps)%2==1}); #on-axis(x=0) fixed-by-mirror={axis} (odd? {axis%2==1})")
print("  => total == #axis-toothpicks (mod 2): the ODD central seed forces A139250 odd. Lefschetz parity.")

# ==========================================================================
sep("B  TOURNAMENTS: A000568 = SC + 2*(chiral pairs); EVEN for n>=3 because SC (self-converse) is even")
def all_tournaments(n):
    pairs=list(combinations(range(n),2))
    for bits in range(2**len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for idx,(i,j) in enumerate(pairs):
            if (bits>>idx)&1: adj[i][j]=1
            else: adj[j][i]=1
        yield adj
def canon(adj,n,perms):
    best=None
    for p in perms:
        key=tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or key<best: best=key
    return best
def converse(adj,n):
    return [[adj[j][i] for j in range(n)] for i in range(n)]
for n in (3,4,5,6):
    perms=list(permutations(range(n)))
    classes={}
    for adj in all_tournaments(n):
        c=canon(adj,n,perms)
        if c not in classes: classes[c]=adj
    total=len(classes)
    sc=0
    for c,adj in classes.items():
        if canon(converse(adj,n),n,perms)==c: sc+=1
    pairs=(total-sc)//2
    print(f"  n={n}: A000568={total} (even? {total%2==0}) ; self-converse SC={sc} (even? {sc%2==0}) ; chiral pairs={pairs} ; check total==SC+2*pairs? {total==sc+2*pairs}")
print("  => THM-587: A000568=P_n(1), SC=P_n(-1) (antipodal Lefschetz). SC=2,2,8,12,88,176 (n=3..8), all even")
print("     => A000568 EVEN for n>=3. count parity = SC = self-converse (achiral) FIXED-POINT count (mod 2).")

# ==========================================================================
sep("C  THE ONE LAW + the toothpick correspondence")
print("""  count = (#involution FIXED) + 2*(#free pairs)  ==>  count == #fixed (mod 2)  [Lefschetz/Euler]:
    TOURNAMENTS  : A000568 == SC (mod 2); SC=P_n(-1) even (n>=3)  => count EVEN. fixed = self-converse(achiral).
    TOOTHPICK    : A139250 == #axis-toothpicks (mod 2); #axis odd (central seed) => total ODD.
    LRC (S212)   : chi(G_delta) == #iota-fixed lonely pts (mod 2); covering => 0 fixed => EVEN (mirror pairs).

  TOOTHPICK DIAGRAM <-> CHIRALITY STRUCTURE:
    central seed toothpick (on both mirror axes, D4-fixed)   <->  the achiral / SELF-CONVERSE tournaments (R-fixed)
    each mirror-symmetric toothpick PAIR (off axis)          <->  each CHIRAL pair {T, T^op}
    the D4 mirror symmetry of the whole diagram              <->  the converse involution R = complement = antipode
    self-similar doubling (A139250(2^k)=(2^(2k+1)+1)/3)      <->  the level-graded metagraph recursion (THM-587)
  So the toothpick fractal is the PICTURE of 'one achiral seed generating chiral mirror-pairs' -- the same
  Z/2 antipodal Euler class (P_n(+1) count / P_n(-1) fixed) that THM-587 proves and S212 uses for LRC.""")
