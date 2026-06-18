"""
M2 STRUCTURAL PROOF CHECK.

CLAIM (to prove): the M2 difference-speed-sign tournament, at ANY single time
tau (not requiring optimality), is SELF-COMPLEMENTARY whenever it is a valid
tournament (no boundary ties). Hence M2 can NEVER realize a non-SC iso class.

Mechanism: arc i->j iff snrm((v_i - v_j) tau) > 0. Consider the relabeling
sigma that REVERSES the speed order is not the point; the point is the map on
VERTICES induced by tau -> 1 - tau. We have snrm((v_i-v_j)(1-tau))
= snrm((v_i-v_j) - (v_i-v_j)tau) = snrm(-(v_i-v_j)tau) = -snrm((v_i-v_j)tau)
(since (v_i-v_j) is an integer, snrm(integer - x) = snrm(-x) = -snrm(x) when
snrm(x) != 1/2). So the tournament at time (1-tau) is the COMPLEMENT (every arc
reversed) of the tournament at time tau, with the SAME vertex labels (identity
map). Therefore T(tau)^op = T(1-tau).

That shows realized-set closure, but to show EACH T(tau) is SC we need an
ISOMORPHISM T(tau) -> T(tau)^op. Test empirically whether every individual
valid M2 tournament is SC (canon(A)==canon(complement(A))) across many tau,
many S. If YES universally, M2 is a tournament map that PROVABLY only hits SC
classes -> all non-SC classes forbidden (4 of 12 at n=5).

Also: does the constraint persist to n=6? Compute realized at n=6 base (A000568(6)=56,
SC count?) at optimal tau, and confirm all realized are SC.
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
def adj_m2(S,tau):
    S=sorted(set(S)); n=len(S); A=[[0]*n for _ in range(n)]; valid=True
    for a in range(n):
        for b in range(a+1,n):
            i,j=S[a],S[b]
            s=snrm((i-j)*tau)
            if s==0: valid=False; A[a][b]=1  # boundary tie
            elif s>0: A[a][b]=1
            else: A[b][a]=1
    return A,valid

print("="*70)
print("(1) Is EVERY individual valid M2 tournament self-complementary?")
print("    (over ALL candidate crossing times, not just optimal)")
print("="*70)
for n_base,vmax in [(3,16),(4,13),(5,11)]:
    total=0; sc=0; nonsc_examples=[]
    for S in gen_speedsets(n_base,vmax):
        for tau in cand(S):
            A,valid=adj_m2(S,tau)
            if not valid: continue
            total+=1
            if is_SC(A): sc+=1
            else:
                if len(nonsc_examples)<3: nonsc_examples.append((S,tau,score_seq(A),count_3cycles(A)))
    print(f"  base n={n_base} vmax={vmax}: valid tournaments over all crossings={total}, "
          f"self-complementary={sc}  -> ALL SC? {sc==total}")
    if nonsc_examples:
        print(f"    NON-SC counterexamples: {nonsc_examples}")

print()
print("="*70)
print("(2) Direct involution check: is the identity map an anti-automorphism")
print("    after reflecting tau->1-tau? Verify T(tau)^op == T(1-tau) labelwise.")
print("="*70)
mismatch=0; checked=0
for S in gen_speedsets(5,11):
    for tau in cand(S):
        A,vA=adj_m2(S,tau)
        if not vA: continue
        B,vB=adj_m2(S,1-tau if (1-tau)>0 else tau)
        if not vB: continue
        checked+=1
        # complement of A should equal B labelwise
        if complement(A)!=B: mismatch+=1
print(f"  checked={checked}, labelwise mismatches (complement(T(tau)) != T(1-tau)): {mismatch}")
print("  (0 mismatches confirms T(tau)^op = T(1-tau) exactly.)")

print()
print("="*70)
print("(3) But is T(tau) ITSELF SC via a NON-identity relabeling? The above")
print("    gives T(tau)^op=T(1-tau); SC needs T(tau)~T(tau)^op. The empirical")
print("    (1) result answers this. Now find the ISOMORPHISM explicitly for one")
print("    example and report the anti-automorphism permutation.")
print("="*70)
S=(1,2,4,7,11); _,tau=M(S)
A,valid=adj_m2(S,tau)
print(f"  S={S}, optimal tau={tau}, valid={valid}, SC={is_SC(A)}, score={score_seq(A)}, c3={count_3cycles(A)}")
Bc=complement(A); n=len(A)
found=None
for p in permutations(range(n)):
    if all(A[p[i]][p[j]]==Bc[i][j] for i in range(n) for j in range(n) if i!=j):
        found=p; break
print(f"  anti-automorphism permutation (T -> T^op): {found}")
# verify it's an involution-free or involution
if found:
    sq=tuple(found[found[i]] for i in range(n))
    print(f"  square of permutation: {sq} (identity => involution)")

print()
print("="*70)
print("(4) Persistence to n=6: census + realized SC-only at optimal tau.")
print("="*70)
# SC census of free set n=6 is expensive (56 classes, 6! canon). Do it.
def free_set(n):
    classes={}; pairs=list(combinations(range(n),2))
    for bits in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for idx,(a,b) in enumerate(pairs):
            if (bits>>idx)&1: A[a][b]=1
            else: A[b][a]=1
        c=canon(A)
        if c not in classes:
            classes[c]=is_SC(A)
    return classes
print("  computing n=6 free set (this enumerates 2^15=32768 tournaments)...")
fs6=free_set(6)
nsc6=sum(1 for v in fs6.values() if v)
print(f"  n=6: {len(fs6)} classes (A000568(6)=56 expected), SC count={nsc6}")

realized6={}
allSC6=True
for S in gen_speedsets(6,12):
    _,tau=M(S)
    A,valid=adj_m2(S,tau)
    if not valid: continue
    c=canon(A)
    if c not in realized6:
        scflag=is_SC(A)
        realized6[c]=scflag
        if not scflag: allSC6=False
print(f"  M2 at optimal tau, n=6 base, vmax=12: realized {len(realized6)}/{len(fs6)} classes; "
      f"ALL realized SC? {allSC6}")
print(f"    => forbidden (non-SC) classes: {len(fs6)-nsc6} non-SC are ALL unreachable by construction; "
      f"plus {nsc6-len(realized6)} SC classes also not realized in this range.")

print("\nDONE.")
