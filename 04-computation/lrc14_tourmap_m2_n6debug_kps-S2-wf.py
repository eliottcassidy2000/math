"""
Debug the n=6 M2 non-SC report. Is M2 universally SC for ALL n, or does it
break at n=6? Find the exact counterexample (S, tau) and verify by hand whether
the tournament is genuinely a valid tournament and genuinely non-SC.

Hypothesis: the snrm reflection argument proves T(tau)^op = T(1-tau) labelwise
ALWAYS. That alone does NOT prove T(tau) is SC. The SC observation at n<=5 may
be a low-n coincidence that breaks at n=6. Settle it exactly.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def snrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else r-1
def g(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def M(S):
    b=F(0); at=None
    for t in cand(S):
        v=g(S,t)
        if v>b: b=v; at=t
    return b,at
def score_seq(A):
    n=len(A); return tuple(sorted(sum(A[i]) for i in range(n)))
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
def complement(A):
    n=len(A); B=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i!=j: B[i][j]=A[j][i]
    return B
def is_SC(A): return canon(A)==canon(complement(A))
def is_primitive(S):
    g0=0
    for v in S: g0=gcd(g0,v)
    return g0==1
def gen_speedsets(n,vmax):
    return [c for c in combinations(range(1,vmax+1),n) if is_primitive(c)]
def adj_m2_strict(S,tau):
    """Strict: only return a tournament if NO boundary tie (snrm==0)."""
    S=sorted(set(S)); n=len(S); A=[[0]*n for _ in range(n)]; valid=True
    for a in range(n):
        for b in range(a+1,n):
            i,j=S[a],S[b]
            s=snrm((i-j)*tau)
            if s==0: valid=False
            elif s>0: A[a][b]=1
            else: A[b][a]=1
    return A,valid

print("Find the FIRST valid non-SC M2 tournament at n=6 base (optimal tau):")
found=False
for S in gen_speedsets(6,12):
    _,tau=M(S)
    A,valid=adj_m2_strict(S,tau)
    if not valid: continue
    if not is_SC(A):
        print(f"  COUNTEREXAMPLE: S={S}, optimal tau={tau}")
        print(f"    score={score_seq(A)}, c3={count_3cycles(A)}, SC={is_SC(A)}")
        # show the signed residues of all pairwise differences
        Sl=sorted(S)
        print(f"    pairwise snrm((vi-vj)*tau):")
        for a in range(len(Sl)):
            for b in range(a+1,len(Sl)):
                print(f"      ({Sl[a]},{Sl[b]}): snrm({(Sl[a]-Sl[b])}*{tau})={snrm((Sl[a]-Sl[b])*tau)}")
        found=True
        break
if not found:
    print("  NONE found at n=6 optimal tau -> M2 IS universally SC at n=6 too.")

print()
print("Now check ALL crossings at n=6 (broader) for ANY valid non-SC M2:")
cnt_valid=0; cnt_nonsc=0; first=None
for S in gen_speedsets(6,10):
    for tau in cand(S):
        A,valid=adj_m2_strict(S,tau)
        if not valid: continue
        cnt_valid+=1
        if not is_SC(A):
            cnt_nonsc+=1
            if first is None: first=(S,tau,score_seq(A),count_3cycles(A))
print(f"  n=6 vmax=10: valid M2 tournaments over all crossings={cnt_valid}, non-SC={cnt_nonsc}")
if first:
    print(f"  first non-SC: S={first[0]}, tau={first[1]}, score={first[2]}, c3={first[3]}")
    print("  => M2 is NOT universally SC; the n<=5 SC property is a LOW-N coincidence.")
else:
    print("  => M2 IS universally SC even over all crossings at n=6.")

print("\nDONE.")
