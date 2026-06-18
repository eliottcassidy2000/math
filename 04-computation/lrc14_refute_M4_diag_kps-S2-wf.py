"""
Diagnostic companion to lrc14_refute_pairwise-gap-load_M4condorcetdiffd_kps-S2-wf.py

Goals:
 1. List WHICH iso classes M4 realizes at the optimal lonely tau (the 5/12).
 2. CONFIRM the regular n=5 class IS reachable at NON-OPTIMAL crossing times
    (validates the pipeline can produce it at all -> forbiddenness at optimum
    is a real loneliness constraint, not a coding artifact that never makes it).
 3. Structural probe: at the optimal tau, the BINDING runners (dist == M) have a
    special role; check the score-multiset distribution of M4 tournaments at
    optimum to see whether the (2,2,2,2,2) regular score is excluded at the
    SCORE level (necessary condition) or only at the iso level.
 4. Broaden once more: ALL optimal taus, vmax up to 23, but ALSO include tau and
    1-tau AND restrict to LRC-relevant (loneliness gap >= 1/13, the "easy core"
    regime) to mimic the N=14 covering reduction.

Exact rationals throughout.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations

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
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at

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
def score_seq(A):
    n = len(A); return tuple(sorted(sum(A[i]) for i in range(n)))
def count_3cycles(A):
    n = len(A); c = 0
    for i in range(n):
        for j in range(n):
            if i == j or not A[i][j]: continue
            for k in range(n):
                if k == i or k == j: continue
                if A[j][k] and A[k][i]: c += 1
    return c // 3
def canon(A):
    n = len(A); best = None
    for p in permutations(range(n)):
        flat = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or flat < best: best = flat
    return best
def H_hampaths(A):
    n = len(A); cnt = 0
    for p in permutations(range(n)):
        ok = True
        for i in range(n-1):
            if not A[p[i]][p[i+1]]: ok = False; break
        if ok: cnt += 1
    return cnt

def method4(S, tau):
    S = sorted(set(S)); verts = list(S)
    dd = {(a, b): nrm((a - b) * tau) for a in S for b in S}
    def arcfn(i, j):
        wi = wj = 0
        for k in S:
            if k == i or k == j: continue
            di = dd[(i, k)]; dj = dd[(j, k)]
            if di < dj: wi += 1
            elif dj < di: wj += 1
        if wi != wj: return wi > wj
        return i < j
    return verts, arcfn
def m4_class(S, tau):
    v, f = method4(S, tau); A, ok = adj_from_arcfn(v, f)
    return (canon(A), A) if ok else None

def is_primitive(S):
    g0 = 0
    for v in S: g0 = gcd(g0, v)
    return g0 == 1
def gen(n, vmax):
    return [c for c in combinations(range(1, vmax+1), n) if is_primitive(c)]

# regular target
def regular5():
    n=5; A=[[0]*n for _ in range(n)]
    for i in range(n):
        A[i][(i+1)%n]=1; A[i][(i+2)%n]=1
    return canon(A)
REG5 = regular5()

# all 12 free classes signatures for reference
def free5():
    n=5; pairs=list(combinations(range(n),2)); cls={}
    for bits in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for idx,(a,b) in enumerate(pairs):
            if (bits>>idx)&1: A[a][b]=1
            else: A[b][a]=1
        c=canon(A)
        if c not in cls: cls[c]=(score_seq(A),count_3cycles(A),H_hampaths(A))
    return cls
FREE5 = free5()
print(f"Free n=5 classes: {len(FREE5)}")
print()

# ---- (1) which classes realized at OPTIMAL tau, with all optimal taus, vmax=23
def all_opt(S):
    c=cand(S); vals={t:g(S,t) for t in c}; Mv=max(vals.values())
    opt=[t for t in c if vals[t]==Mv]
    ext=[1-t for t in opt if 0<1-t<1 and (1-t) not in c]
    return Mv, opt+ext
realized={}
score_at_opt=set()
for S in gen(5,23):
    Mv,opts=all_opt(S)
    for t in opts:
        r=m4_class(S,t)
        if r is None: continue
        c,A=r
        score_at_opt.add(score_seq(A))
        if c not in realized: realized[c]=(score_seq(A),count_3cycles(A),H_hampaths(A),tuple(S),t)
print("="*70)
print("(1) Classes realized at OPTIMAL tau (all optimal taus, vmax=23):")
print("="*70)
for c,(sc,c3,H,S,t) in sorted(realized.items(), key=lambda kv: kv[1][:3]):
    star = "  <-- REGULAR (forbidden target)" if c==REG5 else ""
    print(f"   score={sc} c3={c3} H={H}{star}   e.g. S={S} tau={t}")
forbidden=[ (FREE5[c]) for c in FREE5 if c not in realized]
print(f"\n   realized {len(realized)}/12; FORBIDDEN signatures:")
for f in sorted(set(forbidden)):
    star=" <-- REGULAR" if f==(REG5 and (2,2,2,2,2),) else ""
    print(f"     {f}{' <-- REGULAR' if f[:3]==((2,2,2,2,2),5,15) else ''}")
print(f"\n   distinct SCORE sequences seen at optimum: {sorted(score_at_opt)}")
print(f"   regular score (2,2,2,2,2) present at optimum? {(2,2,2,2,2) in score_at_opt}")
print()

# ---- (2) confirm regular IS reachable at NON-optimal crossings (sanity that map can make it)
print("="*70)
print("(2) Is REGULAR reachable at NON-optimal crossing times? (pipeline sanity)")
print("="*70)
found_off=None; checked=0
for S in gen(5,17):
    Mv,opts=all_opt(S)
    optset=set(opts)
    for t in cand(S):
        if t in optset: continue   # non-optimal only
        r=m4_class(S,t)
        checked+=1
        if r is None: continue
        c,A=r
        if c==REG5:
            found_off=(tuple(S),t,g(S,t),Mv); break
    if found_off: break
if found_off:
    S,t,gv,Mv=found_off
    print(f"   YES regular reached OFF-optimum: S={S} tau={t} gap-there={gv} (optimal M={Mv})")
    print(f"   -> map CAN produce regular; its absence at OPTIMUM is a genuine loneliness constraint.")
else:
    print(f"   regular NOT found off-optimum either (checked {checked} non-opt crossings up to vmax=17)")
print()

# ---- (3) score-level exclusion check: is (2,2,2,2,2) excluded at the SCORE level at optimum?
print("="*70)
print("(3) SCORE-level exclusion of regular (2,2,2,2,2) at optimum?")
print("="*70)
print(f"   (2,2,2,2,2) appears at optimum: {(2,2,2,2,2) in score_at_opt}")
print(f"   If False, regular iso class is excluded a fortiori (no regular-SCORE tournament at optimum).")
