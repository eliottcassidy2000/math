"""
Final adversarial stress test for the M4_condorcet_diffdist forbidden-regular
claim at n=5. Push higher vmax (sampled, since exhaustive C(vmax,5) explodes),
include very spread / sporadic / covering-flavored sets, and test ALL optimal
taus. Also tabulate the SCORE sequences seen at optimum (the deeper invariant):
if (2,2,2,2,2) never appears, regular is excluded a fortiori.

Sampled-large + exhaustive-medium. Exact rationals.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
import random

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def all_opt(S):
    c=cand(S); vals={t:g(S,t) for t in c}; Mv=max(vals.values())
    opt=[t for t in c if vals[t]==Mv]
    ext=[1-t for t in opt if 0<1-t<1 and (1-t) not in c]
    return Mv, opt+ext

def adj_from_arcfn(verts, arcfn):
    n = len(verts); A = [[0]*n for _ in range(n)]; valid = True
    for a in range(n):
        for b in range(a+1, n):
            i, j = verts[a], verts[b]
            ij = arcfn(i, j); ji = arcfn(j, i)
            if ij == ji: valid = False
            if ij: A[a][b] = 1
            else:  A[b][a] = 1
    return A, valid
def score_seq(A): n=len(A); return tuple(sorted(sum(A[i]) for i in range(n)))
def count_3cycles(A):
    n=len(A); c=0
    for i in range(n):
        for j in range(n):
            if i==j or not A[i][j]: continue
            for k in range(n):
                if k==i or k==j: continue
                if A[j][k] and A[k][i]: c+=1
    return c//3
def canon(A):
    n=len(A); best=None
    for p in permutations(range(n)):
        flat=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or flat<best: best=flat
    return best
def H_hp(A):
    n=len(A); cnt=0
    for p in permutations(range(n)):
        ok=True
        for i in range(n-1):
            if not A[p[i]][p[i+1]]: ok=False; break
        if ok: cnt+=1
    return cnt
def method4(S, tau):
    S=sorted(set(S)); verts=list(S)
    dd={(a,b):nrm((a-b)*tau) for a in S for b in S}
    def arcfn(i,j):
        wi=wj=0
        for k in S:
            if k==i or k==j: continue
            di=dd[(i,k)]; dj=dd[(j,k)]
            if di<dj: wi+=1
            elif dj<di: wj+=1
        if wi!=wj: return wi>wj
        return i<j
    return verts,arcfn
def m4(S,t):
    v,f=method4(S,t); A,ok=adj_from_arcfn(v,f)
    return (canon(A),score_seq(A),count_3cycles(A),H_hp(A)) if ok else None
def regular5():
    n=5; A=[[0]*n for _ in range(n)]
    for i in range(n):
        A[i][(i+1)%n]=1; A[i][(i+2)%n]=1
    return canon(A)
REG5=regular5()
def is_primitive(S):
    g0=0
    for v in S: g0=gcd(g0,v)
    return g0==1

random.seed(1234567)
scores_seen=set()
reg_hits=0
reg_witness=None
total_sets=0
total_opt_taus=0

# (A) sampled large vmax
for vmax in (40, 60, 100, 200, 500, 1000):
    nsamp = 4000
    seen=0
    while seen < nsamp:
        S=tuple(sorted(random.sample(range(1,vmax+1),5)))
        if not is_primitive(S): continue
        seen+=1; total_sets+=1
        Mv,opts=all_opt(S)
        for t in opts:
            total_opt_taus+=1
            r=m4(S,t)
            if r is None: continue
            c,sc,c3,H=r
            scores_seen.add(sc)
            if c==REG5:
                reg_hits+=1
                if reg_witness is None: reg_witness=(S,t,Mv)
    print(f"  vmax<={vmax:5d} sampled {nsamp}: cumulative reg_hits={reg_hits} "
          f"distinct-scores={len(scores_seen)} regular-score-seen={(2,2,2,2,2) in scores_seen}", flush=True)

# (B) structured adversarial families: AP-with-one-shifted, geometric, covering-flavored
fams=[]
for d in range(1,8):
    for a in range(1,6):
        S=tuple(a+d*i for i in range(5))
        if is_primitive(S): fams.append(S)
for base in range(2,7):
    S=tuple(base**i for i in range(5))
    if is_primitive(S): fams.append(S)
# covering-flavored: lots of small factors
for S in [(2,3,4,6,12),(2,3,5,30,7),(6,10,15,3,5),(1,2,3,6,12),(4,6,9,12,18),
          (2,6,10,14,35),(1,2,4,7,14),(3,5,15,7,35),(1,2,3,4,84)]:
    Ss=tuple(sorted(set(S)))
    if len(Ss)==5 and is_primitive(Ss): fams.append(Ss)
fams=sorted(set(fams))
for S in fams:
    Mv,opts=all_opt(S)
    for t in opts:
        total_opt_taus+=1
        r=m4(S,t)
        if r is None: continue
        c,sc,c3,H=r
        scores_seen.add(sc)
        if c==REG5:
            reg_hits+=1
            if reg_witness is None: reg_witness=(S,t,Mv)
print(f"  structured families ({len(fams)} sets): cumulative reg_hits={reg_hits} "
      f"regular-score-seen={(2,2,2,2,2) in scores_seen}", flush=True)

print()
print("="*70)
print("STRESS-TEST SUMMARY")
print("="*70)
print(f"  total speed sets probed (sampled-large + structured): {total_sets+len(fams)}")
print(f"  total optimal taus tested: {total_opt_taus}")
print(f"  REGULAR (c3=5,H=15) hits at optimum: {reg_hits}")
print(f"  distinct score-sequences at optimum: {len(scores_seen)}")
print(f"  regular SCORE (2,2,2,2,2) ever at optimum: {(2,2,2,2,2) in scores_seen}")
print(f"  all score-sequences seen: {sorted(scores_seen)}")
if reg_witness:
    print(f"  *** REFUTED: regular witness S={reg_witness[0]} tau={reg_witness[1]} M={reg_witness[2]} ***")
else:
    print(f"  CLAIM SURVIVES stress test (no regular at any optimum).")
