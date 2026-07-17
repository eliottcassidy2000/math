#!/usr/bin/env python3
"""wild_eight_analysis_klein_S317.py — what sees the 8 wild pairs at n=9?
Reconstructs the 8 exact-verified wild pairs from (class8, pattern) coordinates,
confirms non-isomorphism, and tests: deck (vertex-deleted cpA multiset), c3/c4,
|Aut|, reversal relations among the pairs (op-orbits), and tau vectors unsorted."""
import itertools, sys
sys.path.insert(0,'04-computation')
from n9_wild_hunt_klein_S317 import census7, panels_for_stack, canon, exact_panel, charpoly_frac
import numpy as np, hashlib

# rebuild reps8 exactly as in the hunt (deterministic)
from n9_wild_hunt_klein_S317 import build8
reps8,P8=build8()
reps8_A=np.zeros((6880,9,9),dtype=np.int64)
for ci,enc in enumerate(reps8):
    bit=0
    for i in range(8):
        for j in range(i+1,8):
            if (enc>>bit)&1: reps8_A[ci,i,j]=1
            else: reps8_A[ci,j,i]=1
            bit+=1

WILD=[(1168683,1491485),(1415199,1585693),(1422619,1422889),(1514556,1579291),
      (1518607,1563919),(1536271,1581854),(1700451,1708141)]
# find the 8th from the out file
import re
txt=open('05-knowledge/results/n9_wild_hunt_klein_S317.out').read()
gids=set()
for m in re.finditer(r"members gid=\[(\d+), (\d+)\]",txt):
    gids.add((int(m.group(1)),int(m.group(2))))
WILD=sorted(gids)
print(f"wild pairs recovered from census: {len(WILD)}")

def tourn(gid):
    ci,p=divmod(gid,256)
    M=reps8_A[ci].copy()
    for u in range(8):
        if (p>>u)&1: M[u,8]=1
        else: M[8,u]=1
    return [[int(M[i,j]) for j in range(9)] for i in range(9)]

def canon9(A):
    rows=tuple(int(sum((1<<v) if A[u][v] else 0 for v in range(9))) for u in range(9))
    return canon(rows,9)

def deck(A,n):
    cards=[]
    for r in range(n):
        B=[[A[u][v] for v in range(n) if v!=r] for u in range(n) if u!=r]
        cards.append(charpoly_frac(B))
    return tuple(sorted(cards))

def c3c4(A,n):
    A_=np.array(A)
    A2=A_@A_; A3=A2@A_; A4=A3@A_
    c3=int(np.trace(A3))//3
    return c3,int(np.trace(A4))

def autsize(A,n):
    rows=[tuple(A[u]) for u in range(n)]
    cnt=0
    sc=[sum(A[u]) for u in range(n)]
    from itertools import permutations
    for g in permutations(range(n)):
        if any(sc[v]!=sc[g[v]] for v in range(n)): continue
        if all(A[u][v]==A[g[u]][g[v]] for u in range(n) for v in range(n) if u!=v): cnt+=1
    return cnt

canons={}
print()
print("pair | non-iso | deck splits | c3 | tr(A^4) equal | |Aut| | op-partner")
pair_canons=[]
for g1,g2 in WILD:
    A1,A2=tourn(g1),tourn(g2)
    cf1,cf2=canon9(A1),canon9(A2)
    noniso=cf1!=cf2
    d1,d2=deck(A1,9),deck(A2,9)
    c31,t41=c3c4(A1,9); c32,t42=c3c4(A2,9)
    a1,a2=autsize(A1,9),autsize(A2,9)
    pair_canons.append(frozenset((cf1,cf2)))
    print(f"({g1},{g2}) | {noniso} | {d1!=d2} | {c31}={c32} | {t41==t42} | {a1},{a2}")

# reversal orbits among the pairs
def opA(A,n): return [[A[v][u] for v in range(n)] for u in range(n)]
print()
for i,(g1,g2) in enumerate(WILD):
    A1,A2=tourn(g1),tourn(g2)
    op_pair=frozenset((canon9(opA(A1,9)),canon9(opA(A2,9))))
    partner=[j for j,pc in enumerate(pair_canons) if pc==op_pair]
    print(f"pair {i}: op-image is pair {partner} (self-paired = op-closed pair)"
          if partner else f"pair {i}: op-image NOT among the 8 (?)")

# ---- deeper probes: full-panel deck, c5, SC status
def full_deck(A,n):
    cards=[]
    for r in range(n):
        B=[[A[u][v] for v in range(n) if v!=r] for u in range(n) if u!=r]
        cards.append(exact_panel(B,n-1))
    return tuple(sorted(cards))

def c5(A,n):
    A_=np.array(A); A2=A_@A_; A3=A2@A_; A5=A2@A3
    # closed 5-walks = tr A^5; 5-cycles = (tr A^5 - contributions)/5 for digraph:
    # for tournaments (no 2-cycles, no loops) closed 5-walks decompose over 3-cycles
    # with a tail? safest: count directed 5-cycles by brute force over 5-subsets
    from itertools import combinations, permutations
    cnt=0
    for S in combinations(range(n),5):
        for per in permutations(S[1:]):
            cyc=(S[0],)+per
            if all(A[cyc[i]][cyc[(i+1)%5]] for i in range(5)): cnt+=1
    return cnt

print()
print("pair | full-panel deck splits | c5 | member self-complementary?")
for g1,g2 in WILD:
    A1,A2=tourn(g1),tourn(g2)
    fd_split=full_deck(A1,9)!=full_deck(A2,9)
    c51,c52=c5(A1,9),c5(A2,9)
    sc1=canon9(opA(A1,9))==canon9(A1); sc2=canon9(opA(A2,9))==canon9(A2)
    print(f"({g1},{g2}) | {fd_split} | {c51} vs {c52} | SC: {sc1},{sc2}")

# ---- THE CHIRALITY MECHANISM: T2 = T1^op; panel-tie <=> cpL self-dual & tau balanced
print()
print("CHIRALITY VERIFICATION: member2 iso to member1^op; cpL_in == cpL_out; tau_in == tau_out")
ok=True
for g1,g2 in WILD:
    A1,A2=tourn(g1),tourn(g2)
    swap=canon9(opA(A1,9))==canon9(A2)
    p1=exact_panel(A1,9)
    # cpL_in = charpoly of L_in = D_in - A^T = charpoly of L_out(T^op)
    p1op=exact_panel(opA(A1,9),9)
    Lselfdual=(p1[1]==p1op[1])       # cpL_out(T) == cpL_out(T^op)
    taubal=(p1[3]==p1[4])            # tau_in == tau_out (sorted)
    print(f"({g1},{g2}): T2 ~= T1^op: {swap} | cpL self-dual: {Lselfdual} | tau balanced: {taubal}")
    ok&=swap and Lselfdual and taubal
print(f"ALL EIGHT are chirality pairs with self-dual L-spectrum and balanced tree vectors: {ok}")
