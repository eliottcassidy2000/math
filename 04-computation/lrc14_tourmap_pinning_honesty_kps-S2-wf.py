#!/usr/bin/env python3
"""
HONEST reconciliation. The residue-level exhaustion (with parked residue 0 and
free tie-breaks) CAN reach regular T_5. Yet over actual primitive speed sets
vmax<=25 with the SPEED tie-break, M3 never realized it. Resolve which is true
and state the forbiddance precisely.

Questions:
 1. Over primitive 5-speed sets, push vmax much higher: does regular T_5 EVER
    appear under M3 (N=14, speed tie-break)? If not even at large vmax, this is
    a robust EMPIRICAL forbiddance for the real LRC family.
 2. Find the residue+tiebreak witness that produces regular T_5 (from the proof
    script) and see WHY no primitive speed set realizes it (does it need parked
    residue 0 = a multiple of 14? does it need a specific tie-break that the
    speed order can't supply?).
 3. Same for M5's two forbidden classes.
"""
from math import gcd
from itertools import combinations, permutations, combinations_with_replacement as cwr, product
from functools import reduce
N=14
U=[a for a in range(1,N) if gcd(a,N)==1]
def depth(x): r=x%N; return min(r,N-r)
def dvec(v): return tuple(depth((v*a)%N) for a in U)
def margin(ri,rj):
    Di,Dj=dvec(ri),dvec(rj)
    return sum((Di[c]>Dj[c])-(Di[c]<Dj[c]) for c in range(len(U)))
def is_primitive(S): return reduce(gcd,S)==1
def score(adj,k): return tuple(sorted(sum(adj[i][j] for j in range(k) if j!=i) for i in range(k)))
def is_tour(adj,k): return all(adj[i][j]!=adj[j][i] for i in range(k) for j in range(i+1,k))

def m3_speed(S):  # actual M3 with speed tie-break
    k=len(S); adj=[[False]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i==j: continue
            m=margin(S[i]%N,S[j]%N)
            if m>0: adj[i][j]=True
            elif m<0: adj[i][j]=False
            else: adj[i][j]=S[i]<S[j]
    return adj if is_tour(adj,k) else None

print("### Q1: push vmax for primitive 5-speed sets, M3 N=14, look for regular T_5 ###")
def has_regular(vmax):
    for S in combinations(range(1,vmax+1),5):
        if not is_primitive(S): continue
        adj=m3_speed(S)
        if adj and score(adj,5)==(2,2,2,2,2):
            return S
    return None
for vmax in (25,30,40):
    r=has_regular(vmax)
    n=sum(1 for S in combinations(range(1,vmax+1),5) if is_primitive(S))
    print(f"  vmax={vmax} ({n} primitive sets): regular T_5 found = {r}")

print("\n### Q2: residue+tiebreak witness for regular T_5, and is it speed-realizable? ###")
# enumerate residue multisets (incl 0) + tie-breaks that give regular; report
# the residue multiset and required tie-break pattern.
witnesses=[]
for combo in cwr(range(0,N),5):
    k=5
    base=[[None]*k for _ in range(k)]
    ties=[]
    for i in range(k):
        for j in range(i+1,k):
            m=margin(combo[i],combo[j])
            if m>0: base[i][j]=True; base[j][i]=False
            elif m<0: base[i][j]=False; base[j][i]=True
            else: ties.append((i,j))
    for tb in product([True,False],repeat=len(ties)):
        adj=[[False]*k for _ in range(k)]
        for i in range(k):
            for j in range(k):
                if base[i][j] is not None: adj[i][j]=base[i][j]
        for idx,(i,j) in enumerate(ties):
            adj[i][j]=tb[idx]; adj[j][i]=not tb[idx]
        if is_tour(adj,k) and score(adj,k)==(2,2,2,2,2):
            witnesses.append((combo,ties,tb)); break  # one tb per combo enough
print(f"  residue multisets (with some tie-break) giving regular T_5: {len(witnesses)}")
# how many require residue 0 (parked)? how many require duplicated residues?
need0=sum(1 for w in witnesses if 0 in w[0])
distinct=sum(1 for w in witnesses if len(set(w[0]))==5)
have_tie=sum(1 for w in witnesses if len(w[1])>0)
print(f"   ... require a parked residue 0: {need0}/{len(witnesses)}")
print(f"   ... have all-distinct residues: {distinct}/{len(witnesses)}")
print(f"   ... rely on a margin-0 tie pair (tie-break dependent): {have_tie}/{len(witnesses)}")
if witnesses:
    print("   sample witnesses (residue multiset, tie pairs):")
    for w in witnesses[:8]:
        print("     ", w[0], "ties:", w[1])

# Now: for a DISTINCT-residue, NO-TIE witness (if any), a primitive speed set
# realizing those residues mod 14 would realize regular T_5. Check if such a
# witness exists and find a primitive speed set:
clean=[w for w in witnesses if len(set(w[0]))==5 and len(w[1])==0 and 0 not in w[0]]
print(f"\n  CLEAN witnesses (distinct nonzero residues, no tie dependence): {len(clean)}")
if clean:
    res=clean[0][0]
    print("   e.g. residues", res, "-> any primitive speed set with these residues mod 14 gives regular T_5")
    # build smallest speed set: speeds = residues themselves (if primitive)
    S=tuple(sorted(res))
    print("   speeds = residues:", S, "primitive?", is_primitive(S), "M3 score:", score(m3_speed(S),5) if m3_speed(S) else None)
else:
    print("   => EVERY regular-T_5 witness needs EITHER a parked residue 0, OR repeated residues, OR a tie-break.")
    print("   This explains why generic primitive speed sets (distinct residues, no parking) avoid it.")

print("\n### Q3: M5 forbidden classes over primitive sets, push vmax ###")
HALF=(1,3,5)
def m5_speed(S):
    k=len(S); adj=[[False]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1,k):
            c2=sum(1 for a in HALF if (S[i]*a)%N < (S[j]*a)%N)
            if c2%2==1: adj[i][j]=True
            else: adj[j][i]=True
    return adj if is_tour(adj,k) else None
def num3(adj,k):
    c=0
    for i,j,l in combinations(range(k),3):
        if adj[i][j] and adj[j][l] and adj[l][i]:c+=1
        if adj[i][l] and adj[l][j] and adj[j][i]:c+=1
    return c
# track score+c3 realized
import collections
real=set()
for S in combinations(range(1,40),5):
    if not is_primitive(S): continue
    adj=m5_speed(S)
    if adj: real.add((score(adj,5),num3(adj,5)))
print("  M5 (vmax=40) realized (score,c3):")
for x in sorted(real): print("    ",x)
# free classes' (score,c3):
freelist=[((0,1,2,3,4),0),((0,1,3,3,3),1),((0,2,2,2,4),1),((0,2,2,3,3),2),
          ((1,1,1,3,4),1),((1,1,2,2,4),2),((1,1,2,3,3),3),((1,2,2,2,3),4),
          ((2,2,2,2,2),5)]
miss=[x for x in freelist if x not in real]
print("  M5 (score,c3) NEVER realized at vmax=40:", miss)
