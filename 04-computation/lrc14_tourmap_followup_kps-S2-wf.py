"""
Follow-up analysis for the two NON-TRIVIAL methods (M2 diff-speed-sign,
M4 condorcet-diffdist). Goals:

(A) Understand WHY M2/M3 produce only "balanced" antisymmetric classes:
    these arc rules are based on snrm((v_i +/- v_j) tau), which is an
    ANTISYMMETRIC function of the pair under a relabel symmetry. Prove the
    structural constraint: M2 tournaments are exactly those orientable as
    "sign of a skew form" -> they are forced to be COMPLEMENT-INVARIANT under
    the involution v -> (something). Check empirically: are all M2-realized
    classes self-complementary (SC)? Are non-SC classes FORBIDDEN structurally?

(B) For M4, is the regular class (2,2,2,2,2) [the rotational/Paley T5, c3=5,
    H=15] robustly FORBIDDEN even over all crossings AND larger speed ranges?
    Extend vmax and sample to stress-test.

(C) For M2: prove forbidden = non-SC. The diff-speed-sign rule:
    i->j iff snrm((v_i - v_j) tau) > 0. Claim: replacing every v_i by -v_i
    (mod 1, i.e. reflecting tau) reverses every arc => the tournament is
    isomorphic to its COMPLEMENT via the identity-ish map only if symmetric.
    Test: is M2 ALWAYS self-complementary? (SC count at n=5 is 12 classes? no)
    Recall A000568(5)=12, of which SC count = 12? Let's compute #SC classes.

All exact rationals.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
import sys

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def snrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else r - 1
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

def score_seq(A):
    n = len(A); return tuple(sorted(sum(A[i]) for i in range(n)))
def count_3cycles(A):
    n = len(A); c = 0
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
def complement(A):
    n=len(A); B=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i!=j: B[i][j]=A[j][i]
    return B
def is_SC(A):
    return canon(A)==canon(complement(A))
def H_hampaths(A):
    n=len(A); cnt=0
    for p in permutations(range(n)):
        if all(A[p[i]][p[i+1]] for i in range(n-1)): cnt+=1
    return cnt

def free_set(n):
    classes={}; pairs=list(combinations(range(n),2))
    for bits in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for idx,(a,b) in enumerate(pairs):
            if (bits>>idx)&1: A[a][b]=1
            else: A[b][a]=1
        c=canon(A)
        if c not in classes:
            classes[c]=(score_seq(A),count_3cycles(A),H_hampaths(A),is_SC(A))
    return classes

def is_primitive(S):
    g0=0
    for v in S: g0=gcd(g0,v)
    return g0==1
def gen_speedsets(n,vmax):
    return [c for c in combinations(range(1,vmax+1),n) if is_primitive(c)]

def adj_from_arcfn(verts,arcfn):
    n=len(verts); A=[[0]*n for _ in range(n)]; valid=True
    for a in range(n):
        for b in range(a+1,n):
            i,j=verts[a],verts[b]
            ij=arcfn(i,j); ji=arcfn(j,i)
            if ij==ji: valid=False
            if ij: A[a][b]=1
            else: A[b][a]=1
    return A,valid

# Methods
def method2(S,tau):
    S=sorted(set(S)); verts=list(S)
    def arcfn(i,j):
        s=snrm((i-j)*tau)
        if s>0: return True
        if s<0: return False
        return i<j
    return verts,arcfn
def method4(S,tau):
    S=sorted(set(S)); verts=list(S)
    dd={}
    for a in S:
        for b in S: dd[(a,b)]=nrm((a-b)*tau)
    def arcfn(i,j):
        wi=0;wj=0
        for k in S:
            if k==i or k==j: continue
            if dd[(i,k)]<dd[(j,k)]: wi+=1
            elif dd[(j,k)]<dd[(i,k)]: wj+=1
        if wi!=wj: return wi>wj
        return i<j
    return verts,arcfn

print("="*70)
print("(0) SC-class census of the FREE set")
print("="*70)
for n in (3,4,5):
    fs=free_set(n)
    nsc=sum(1 for v in fs.values() if v[3])
    print(f"  n={n}: {len(fs)} classes, {nsc} self-complementary (SC)")
    # list SC signatures
    scsig=sorted(set((v[0],v[1],v[2]) for v in fs.values() if v[3]))
    print(f"    SC signatures (score,c3,H): {scsig}")

FREE={n:free_set(n) for n in (3,4,5)}

print()
print("="*70)
print("(A) M2 (diff-speed-sign): is EVERY realized class self-complementary?")
print("    And is every NON-SC class forbidden?")
print("="*70)
for n_base,vmax in [(4,14),(5,12)]:
    realized={}
    for S in gen_speedsets(n_base,vmax):
        _,tau=M(S)
        verts,arcfn=method2(S,tau)
        A,valid=adj_from_arcfn(verts,arcfn)
        if not valid: continue
        c=canon(A)
        if c not in realized:
            realized[c]=(score_seq(A),count_3cycles(A),H_hampaths(A),is_SC(A))
    free=FREE[n_base]
    allSC = all(v[3] for v in realized.values())
    realized_SC = set(c for c in realized)
    free_SC = set(c for c,v in free.items() if v[3])
    free_nonSC = set(c for c,v in free.items() if not v[3])
    print(f"  base n={n_base} vmax={vmax}: realized {len(realized)}/{len(free)}; "
          f"ALL realized are SC = {allSC}")
    print(f"    SC classes realized: {len(realized_SC & free_SC)}/{len(free_SC)}; "
          f"any non-SC realized? {len(realized_SC & free_nonSC)>0}")
    missingSC=free_SC-realized_SC
    if missingSC:
        ms=sorted(set((free[c][0],free[c][1],free[c][2]) for c in missingSC))
        print(f"    SC classes NOT realized: {ms}")

print()
print("="*70)
print("(A') M2 STRUCTURAL: prove arc-reversal under tau->-tau (reflection).")
print("    snrm((i-j)*(-tau)) = -snrm((i-j)*tau) (generically). So reflecting")
print("    the lonely time reverses ALL arcs => complement tournament.")
print("    Since g(S,tau)=g(S,1-tau) (||.|| symmetric), the lonely set is")
print("    reflection-symmetric: if tau optimal then 1-tau is ALSO optimal.")
print("    Check: at BOTH optimal taus we get a class and its complement; the")
print("    REALIZED SET is thus closed under complement. Verify empirically.")
print("="*70)
for n_base,vmax in [(5,11)]:
    closed=True; examples=[]
    realized=set()
    for S in gen_speedsets(n_base,vmax):
        Mv,tau=M(S)
        for t in (tau, 1-tau):
            if t<=0 or t>=1: continue
            verts,arcfn=method2(S,t)
            A,valid=adj_from_arcfn(verts,arcfn)
            if not valid: continue
            realized.add(canon(A))
    # check closure under complement
    for c in list(realized):
        # rebuild A from canon c
        n=n_base; A=[[c[i*n+j] for j in range(n)] for i in range(n)]
        if canon(complement(A)) not in realized:
            closed=False
    print(f"  base n={n_base}: realized set closed under complement (incl 1-tau)? {closed}; "
          f"#classes={len(realized)}")

print()
print("="*70)
print("(B) M4 regular-class forbiddenness STRESS TEST")
print("    Is (2,2,2,2,2) [c3=5,H=15, the rotational T5] forbidden at optimal")
print("    lonely tau over LARGER speed ranges? Also over ALL crossings?")
print("="*70)
target_sig=((2,2,2,2,2),5,15)
for vmax in [11,14,17,20]:
    found_opt=False; found_all=False; ex=None
    for S in gen_speedsets(5,vmax):
        Mv,tau=M(S)
        verts,arcfn=method4(S,tau)
        A,valid=adj_from_arcfn(verts,arcfn)
        if valid and (score_seq(A),count_3cycles(A),H_hampaths(A))==target_sig:
            found_opt=True; ex=(S,tau); break
    # all crossings
    for S in gen_speedsets(5,vmax):
        if found_all: break
        for t in cand(S):
            verts,arcfn=method4(S,t)
            A,valid=adj_from_arcfn(verts,arcfn)
            if valid and (score_seq(A),count_3cycles(A),H_hampaths(A))==target_sig:
                found_all=True; break
    ns=len(gen_speedsets(5,vmax))
    print(f"  vmax={vmax} ({ns} speedsets): regular class at OPTIMAL tau? {found_opt}; "
          f"over ALL crossings? {found_all}" + (f"  ex={ex}" if found_opt else ""))

print()
print("="*70)
print("(B') M4 full realized-class profile at optimal tau, larger range")
print("="*70)
for vmax in [11,14,17]:
    realized={}
    for S in gen_speedsets(5,vmax):
        _,tau=M(S)
        verts,arcfn=method4(S,tau)
        A,valid=adj_from_arcfn(verts,arcfn)
        if not valid: continue
        c=canon(A)
        if c not in realized:
            realized[c]=(score_seq(A),count_3cycles(A),H_hampaths(A))
    free=FREE[5]
    forb=set(free)-set(realized)
    fs=sorted(set((free[c][0],free[c][1],free[c][2]) for c in forb))
    print(f"  vmax={vmax}: realized {len(realized)}/12; forbidden sigs: {fs}")

print("\nDONE.")
