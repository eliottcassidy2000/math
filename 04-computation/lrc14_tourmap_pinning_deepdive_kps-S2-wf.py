#!/usr/bin/env python3
"""
Deep dive on the two NON-TRIVIAL pinning maps M3 (section majority vote)
and M5 (half-grid parity XOR). Goals:
  1. Widen the input family (larger vmax, more speeds) to test whether the
     k=5 forbidden classes survive (or were artifacts of small vmax).
  2. Characterize the forbidden classes structurally.
  3. For M3: prove/explain why the regular tournament (c3=5,H=15) is forbidden.
  4. Test the GRID parameter N: M3/M5 hardwire N=14. Is the regular-tournament
     forbiddance specific to N=14, or general?
EXACT integers only (these maps use only mod-N residues, no Fractions needed
for the arc decisions, but we keep Fraction infra available).
"""
from math import gcd
from itertools import combinations, permutations
from functools import reduce

def is_primitive(S): return reduce(gcd, S) == 1

# ---- canonical tournament infra (k<=6) ----
def tour_key(adj, k):
    best = None
    for perm in permutations(range(k)):
        bits = 0
        for i in range(k):
            for j in range(k):
                if i == j: continue
                bits = (bits << 1) | (1 if adj[perm[i]][perm[j]] else 0)
        if best is None or bits < best: best = bits
    return best
def is_tournament(adj, k):
    for i in range(k):
        for j in range(i+1,k):
            if adj[i][j]==adj[j][i]: return False
    return True
def num_3cycles(adj,k):
    c=0
    for i,j,l in combinations(range(k),3):
        if adj[i][j] and adj[j][l] and adj[l][i]: c+=1
        if adj[i][l] and adj[l][j] and adj[j][i]: c+=1
    return c
def num_hampaths(adj,k):
    h=0
    for perm in permutations(range(k)):
        if all(adj[perm[a]][perm[a+1]] for a in range(k-1)): h+=1
    return h
def score_seq(adj,k):
    return tuple(sorted(sum(1 for j in range(k) if i!=j and adj[i][j]) for i in range(k)))

A000568={3:2,4:4,5:12,6:56}
def enumerate_free(k):
    classes={}
    edges=list(combinations(range(k),2)); m=len(edges)
    for bits in range(1<<m):
        adj=[[False]*k for _ in range(k)]
        for idx,(i,j) in enumerate(edges):
            if (bits>>idx)&1: adj[i][j]=True
            else: adj[j][i]=True
        key=tour_key(adj,k)
        if key not in classes:
            classes[key]={'score':score_seq(adj,k),'c3':num_3cycles(adj,k),'H':num_hampaths(adj,k)}
    return classes
FREE={k:enumerate_free(k) for k in (3,4,5)}
for k in (3,4,5): assert len(FREE[k])==A000568[k]
def label(k,key):
    i=FREE[k][key]; return f"[score={i['score']},c3={i['c3']},H={i['H']}]"

# ---- METHOD 3: section majority vote, parametrized by N and unit set ----
def units(N): return [a for a in range(1,N) if gcd(a,N)==1]
def method3(S, N=14):
    k=len(S); U=units(N)
    adj=[[False]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i==j: continue
            tally=0
            for a in U:
                ri=(S[i]*a)%N; rj=(S[j]*a)%N
                di=min(ri,N-ri); dj=min(rj,N-rj)
                if di>dj: tally+=1
                elif di<dj: tally-=1
            if tally>0: adj[i][j]=True
            elif tally<0: adj[i][j]=False
            else: adj[i][j]=S[i]<S[j]
    return adj if is_tournament(adj,k) else None

# ---- METHOD 5: half-grid parity XOR ----
def method5(S, N=14, half=(1,3,5)):
    k=len(S)
    adj=[[False]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1,k):
            c2=sum(1 for a in half if (S[i]*a)%N < (S[j]*a)%N)
            if c2%2==1: adj[i][j]=True
            else: adj[j][i]=True
    return adj if is_tournament(adj,k) else None

def realized_classes(fn, sets, k):
    rs=set()
    for S in sets:
        adj=fn(S)
        if adj is None: continue
        rs.add(tour_key(adj,k))
    return rs

def speed_sets(nspeed, vmax, primitive=True):
    out=[]
    for S in combinations(range(1,vmax+1),nspeed):
        if (not primitive) or is_primitive(S):
            out.append(S)
    return out

print("### WIDEN vmax for k=5, N=14 ###")
for vmax in (10, 13, 16, 20, 25):
    sets=speed_sets(5,vmax)
    r3=realized_classes(lambda S: method3(S,14), sets, 5)
    r5=realized_classes(lambda S: method5(S,14), sets, 5)
    f3=set(FREE[5])-r3; f5=set(FREE[5])-r5
    print(f" vmax={vmax} (#sets={len(sets)}): M3 realizes {len(r3)}/12 forbids {len(f3)} | M5 realizes {len(r5)}/12 forbids {len(f5)}")

print("\n### Surviving forbidden classes at vmax=25, N=14 ###")
sets=speed_sets(5,25)
r3=realized_classes(lambda S: method3(S,14), sets,5)
r5=realized_classes(lambda S: method5(S,14), sets,5)
print(" M3 forbidden:")
for key in sorted(set(FREE[5])-r3): print("   ", label(5,key))
print(" M5 forbidden:")
for key in sorted(set(FREE[5])-r5): print("   ", label(5,key))

print("\n### Does the forbiddance depend on N? Test M3 regular-tournament for various N ###")
# regular tournament on 5 vertices = the rotational T_5 (i->i+1,i+2 mod5).
# It is the UNIQUE iso class with score (2,2,2,2,2). Check if M3 EVER realizes it.
reg_key = None
for key,info in FREE[5].items():
    if info['score']==(2,2,2,2,2): reg_key=key
print(" regular T_5 canonical key:", reg_key, label(5,reg_key))
for N in (7,9,11,13,14,15,16,18,20):
    sets=speed_sets(5,min(2*N,26))
    r3=realized_classes(lambda S: method3(S,N), sets,5)
    has_reg = reg_key in r3
    f3=set(FREE[5])-r3
    print(f"  N={N}: M3 realizes {len(r3)}/12, regular-T5 realized={has_reg}, forbids {len(f3)}")

print("\n### M3: which iso classes are realizable at k=6, N=14 (sample) ###")
sets6=speed_sets(6,16)
r36=realized_classes(lambda S: method3(S,14), sets6, 6)
print(f"  M3 k=6 realizes {len(r36)}/{A000568[6]} (sampled vmax=16, #sets={len(sets6)})")

print("\n### CHARACTERIZE M3 score sequences realized vs free (k=5,N=14,vmax=25) ###")
sets=speed_sets(5,25)
realized_scores=set()
for S in sets:
    adj=method3(S,14)
    if adj: realized_scores.add(score_seq(adj,5))
free_scores=set(info['score'] for info in FREE[5].values())
print("  free score seqs:", sorted(free_scores))
print("  M3 realized score seqs:", sorted(realized_scores))
print("  score seqs NEVER realized by M3:", sorted(free_scores-realized_scores))
