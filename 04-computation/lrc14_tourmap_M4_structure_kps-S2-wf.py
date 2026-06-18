"""
lrc14_tourmap_M4_structure_kps-S2-wf.py

STRUCTURAL analysis of the M4 'band-depth-majority' map (the most promising
non-trivial residue-orbit tournament). GOAL: explain WHY it forbids the
regular class and prove the forbidden-class claim rather than just sample it.

M4 DEFINITION (mod 14). Vertices = runners (speeds v_1..v_m). For each unit
a in (Z/14)* = {1,3,5,9,11,13}, runner i has section r_i(a)=v_i*a mod 14 and
DEPTH d_i(a) = min(r_i(a), 14 - r_i(a)) in {0,1,2,...,7} (7 = the midpoint, the
loneliest band; 0 = sitting on the observer). Arc i->j iff
  #{a : d_i(a) > d_j(a)}  >  #{a : d_j(a) > d_i(a)}     (orbit depth-majority).

KEY STRUCTURAL FACT: the DEPTH PROFILE multiset {d_i(a) : a in U} depends ONLY on
gcd(v_i, 14) and the orbit of v_i under (Z/14)*. We classify speeds by residue mod 14:
  - residue r and 14-r give the SAME depth profile (since depth is symmetric).
  - the group (Z/14)* acts on residues; runner depth profile = a class function.

So M4 is really a tournament on the QUOTIENT TYPES. We enumerate the finite set of
profile-types and show the resulting comparison is a (weak) TOTAL PREORDER -> the
tournament is a blow-up of a transitive order EXCEPT where ties are broken by speed.
This pins down exactly which iso classes are reachable.
"""
from math import gcd
from itertools import combinations, permutations
from collections import Counter

MOD = 14
U = [a for a in range(1, MOD) if gcd(a, MOD) == 1]

def depth_profile(v):
    """Sorted tuple of depths d(a)=min(va mod14, 14-(va mod14)) over units a."""
    prof = []
    for a in U:
        r = (v * a) % MOD
        prof.append(min(r, MOD - r))
    return tuple(prof)  # order matches U = [1,3,5,9,11,13]

def depth_profile_multiset(v):
    return tuple(sorted(depth_profile(v)))

print("="*70)
print("DEPTH PROFILES by residue mod 14 (units act; depth is the band 0..7)")
print(f"Units U = {U}")
profiles = {}
for r in range(MOD):
    p = depth_profile(r)
    profiles[r] = p
    print(f"  residue {r:2d}: depth-profile over U = {p}   multiset={tuple(sorted(p))}")

print()
print("Group of distinct depth-MULTISETS (these are the only 'types' M4 can see):")
types = {}
for r in range(MOD):
    ms = depth_profile_multiset(r)
    types.setdefault(ms, []).append(r)
for ms, rs in sorted(types.items()):
    print(f"  type {ms}: residues {rs}")

print()
print("="*70)
print("M4 comparison between TYPES: for residues r,s define")
print("  win(r,s) = #{a: depth_r(a) > depth_s(a)} - #{a: depth_s(a) > depth_r(a)}")
print("Build the TYPE tournament; check transitivity / 3-cycles.")
# Use representative residues, one per type
reps = sorted({rs[0]: ms for ms, rs in types.items()}.keys())
reps = [rs[0] for ms, rs in sorted(types.items())]
def win(r, s):
    w = 0
    pr = depth_profile(r); ps = depth_profile(s)
    for x, y in zip(pr, ps):
        if x > y: w += 1
        elif y > x: w -= 1
    return w
print(f"  representative residues (one per type): {reps}")
print("  pairwise win matrix (row beats col if >0):")
for r in reps:
    row = []
    for s in reps:
        row.append(f"{win(r,s):+d}" if r != s else "  .")
    print(f"   r={r:2d}: " + " ".join(f"{x:>4}" for x in row))

# Is the relation 'win(r,s)>0' transitive among types?
def beats(r, s):
    w = win(r, s)
    if w > 0: return True
    if w < 0: return False
    return r < s  # tie-break by residue (proxy; real M4 ties break by speed)
trans = True
threecyc = []
for a, b, c in permutations(reps, 3):
    if beats(a, b) and beats(b, c) and beats(c, a):
        trans = False
        threecyc.append((a, b, c))
print(f"\n  Type relation transitive (no 3-cycle among types)? {trans}")
if threecyc:
    print(f"  3-cycles among types: {threecyc[:6]}")

print()
print("="*70)
print("CONSEQUENCE for regular class. A REGULAR tournament on 5 vertices needs")
print("score (2,2,2,2,2): every vertex beats exactly 2. Check: can 5 speeds chosen")
print("from these types ever give a regular M4 tournament? Exhaustive over residue")
print("multisets (with multiplicity, since speeds can share a residue mod 14).")

def m4_tournament_from_residues(res, speeds):
    m = len(res)
    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            w = win(res[i], res[j])
            if w > 0: adj[i][j] = 1
            elif w < 0: adj[j][i] = 1
            else:
                adj[i][j] = 1 if speeds[i] < speeds[j] else 0
                adj[j][i] = 1 - adj[i][j]
    return adj

def score_seq(adj, m):
    return tuple(sorted(sum(adj[i][j] for j in range(m) if j != i) for i in range(m)))
def num_3cycles(adj, m):
    c=0
    for a,b,cc in combinations(range(m),3):
        if adj[a][b] and adj[b][cc] and adj[cc][a]: c+=1
        if adj[a][cc] and adj[cc][b] and adj[b][a]: c+=1
    return c
def canon_key(adj,m):
    best=None
    for perm in permutations(range(m)):
        bits=tuple(adj[perm[i]][perm[j]] for i in range(m) for j in range(m) if i!=j)
        if best is None or bits<best: best=bits
    return best
def h_count(adj,m):
    return sum(1 for perm in permutations(range(m)) if all(adj[perm[k]][perm[k+1]] for k in range(m-1)))

# enumerate all residue choices (0..13) for 5 runners, with speeds = residue+14*tiebreak
# to realize all possible tie-break orders, we let speeds be a strictly increasing
# sequence consistent with chosen residues. To capture tie-break freedom, we test BOTH
# possible speed orders within equal-residue groups by trying all permutations of a
# distinct speed assignment.
print("\nEnumerating ALL residue-5-multisets mod 14 with ALL tie-break orders...", flush=True)
reached = set()
score_counter = Counter()
regular_found = False
seen_resmulti = set()
import itertools
# precompute pairwise win(r,s) for all residue pairs (small)
WIN = {(r,s): win(r,s) for r in range(MOD) for s in range(MOD)}
def m4_fast(res, speeds):
    m=len(res); adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            w=WIN[(res[i],res[j])]
            if w>0: adj[i][j]=1
            elif w<0: adj[j][i]=1
            else:
                adj[i][j]=1 if speeds[i]<speeds[j] else 0; adj[j][i]=1-adj[i][j]
    return adj
allperms = list(set(permutations(range(5))))
for res in itertools.combinations_with_replacement(range(MOD), 5):
    for perm in allperms:
        speeds = [0]*5
        for rank, pos in enumerate(perm):
            speeds[pos] = rank + 1
        adj = m4_fast(list(res), speeds)
        k = canon_key(adj, 5)
        reached.add(k)
        sc = score_seq(adj, 5)
        score_counter[sc] += 1
        if sc == (2,2,2,2,2):
            regular_found = True
print(f"  distinct iso classes reachable by M4 over ALL residue-multisets+tiebreaks: {len(reached)}")
print(f"  REGULAR (2,2,2,2,2) ever reached? {regular_found}")
print("  score-sequence reachability (M4, all residues mod 14):")
for sc in sorted(score_counter):
    print(f"     {sc}{'  <-- REGULAR' if sc==(2,2,2,2,2) else ''}")

# Compare to free set
def all_iso_classes(m):
    from itertools import product as pr
    seen={}
    pairs=list(combinations(range(m),2))
    for bits in pr([0,1],repeat=len(pairs)):
        adj=[[0]*m for _ in range(m)]
        for (i,j),b in zip(pairs,bits):
            if b: adj[i][j]=1
            else: adj[j][i]=1
        kk=canon_key(adj,m)
        if kk not in seen: seen[kk]=(h_count(adj,m),num_3cycles(adj,m),score_seq(adj,m))
    return seen
free5=all_iso_classes(5)
forb = set(free5)-reached
print(f"\n  FORBIDDEN under M4 (residue-exhaustive, n=5): {len(forb)} of 12")
for k in sorted(forb, key=lambda k: free5[k]):
    print(f"     H={free5[k][0]}, #3cyc={free5[k][1]}, score={free5[k][2]}")
