"""
Adversarial refutation test for the 'character-spectral / M5 binding-parity'
forbidden-class claim.

CLAIM (n=5): the M5 binding-parity tournament, evaluated at the LRC gap-optimal
time tau* of each primitive speed set S, NEVER realizes the iso class
  TARGET = (H=15, #3cyc=4, score(1,2,2,2,3))  [unique max-H NON-regular
  strongly-connected tournament on 5 vertices], while realizing the other 11.

M5 binding-parity (per the claim): vertices = the n speeds. At tau*:
  theta_i=frac(v_i tau*), w_i=floor(v_i tau*). Arc x->y iff
  (theta_x>theta_y) XOR ((w_x+w_y) odd), equal-phase tie-break x>y (by speed).

We adversarially try to REALIZE the forbidden class with this map over a broad,
exact (Fraction) search:
  - all primitive 5-subsets of {1..18},{1..24},{1..30}
  - many construction-detail variants (base dir, tie-break, twist parity, which
    optimal tau) so we do not accidentally test a strictly weaker map
  - ALL gap-optimal taus, plus the all-odd peak tau=1/2
  - off-grid taus and a dense exact sweep of candidate crossing times
  - random/sporadic large-speed sets, covering-flavored and AP-flavored sets
If TARGET is realized once -> REFUTED.  (H,3cyc,score) is a COMPLETE iso-class
separator at n=5 -- verified by full 2^10 enumeration.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
from collections import Counter
import random, sys

def flush(*a): print(*a); sys.stdout.flush()

# ---------- exact LRC gap tool ----------
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
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
def M_all_opt(S):
    b=F(0); ats=[]
    for t in sorted(cand(S)):
        v=g(S,t)
        if v>b: b=v; ats=[t]
        elif v==b: ats.append(t)
    return b,ats
def ffloor(x): return x.numerator//x.denominator

# ---------- M5 tournament ----------
def build(S,t,basedir='gt',tiebreak='speedhi',twistparity=1):
    n=len(S)
    theta=[S[i]*t - ffloor(S[i]*t) for i in range(n)]
    w=[ffloor(S[i]*t) for i in range(n)]
    A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if theta[i]!=theta[j]:
                base=(theta[i]>theta[j]) if basedir=='gt' else (theta[i]<theta[j])
            else:
                if tiebreak=='speedhi': base=(S[i]>S[j])
                elif tiebreak=='speedlo': base=(S[i]<S[j])
                elif tiebreak=='idxhi': base=(i>j)
                else: base=(i<j)
            tw=((w[i]+w[j])%2==twistparity)
            if base^tw: A[i][j]=1
            else: A[j][i]=1
    return A

# ---------- invariants (fast) ----------
def is_tournament(A,n):
    for i in range(n):
        for j in range(i+1,n):
            if A[i][j]+A[j][i]!=1: return False
    return True
def score_seq(A,n): return tuple(sorted(sum(A[i]) for i in range(n)))
def count_ham_paths(A,n):
    out=[[j for j in range(n) if A[i][j]==1] for i in range(n)]
    size=1<<n
    dp=[[0]*n for _ in range(size)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(size):
        row=dp[mask]
        for last in range(n):
            c=row[last]
            if not c: continue
            for nxt in out[last]:
                b=1<<nxt
                if mask&b: continue
                dp[mask|b][nxt]+=c
    return sum(dp[size-1])
def count_3cycles(A,n):
    c=0
    for i,j,k in combinations(range(n),3):
        if A[i][j] and A[j][k] and A[k][i]: c+=1
        if A[i][k] and A[k][j] and A[j][i]: c+=1
    return c
def inv5(A): return (count_ham_paths(A,5),count_3cycles(A,5),score_seq(A,5))
def primitive(S):
    g0=0
    for v in S: g0=gcd(g0,v)
    return g0==1

TARGET=(15,4,(1,2,2,2,3))
VARIANTS=[]
for bd in ['gt','lt']:
    for tb in ['speedhi','speedlo','idxhi','idxlo']:
        VARIANTS.append((bd,tb,1))   # twist parity is a global relabel; only need one

def test_set_all_variants_all_opt(S, witnesses):
    realized=set()
    gap,ats=M_all_opt(list(S))
    for topt in ats:
        for (bd,tb,tp) in VARIANTS:
            A=build(list(S),topt,bd,tb,tp)
            if not is_tournament(A,5): continue
            k=inv5(A); realized.add(k)
            if k==TARGET:
                witnesses.append((tuple(S),topt,bd,tb,tp,gap))
    return realized

def main():
    flush("TARGET (forbidden) =",TARGET)
    for R in [18,24,30]:
        configs=[s for s in combinations(range(1,R+1),5) if primitive(s)]
        allrealized=set(); witnesses=[]
        for S in configs:
            allrealized |= test_set_all_variants_all_opt(S,witnesses)
        flush(f"\n[PHASE A] R={R}: {len(configs)} primitive 5-sets, ALL {len(VARIANTS)} variants x ALL optimal taus")
        flush(f"  distinct classes realized: {len(allrealized)}")
        flush(f"  TARGET realized: {len(witnesses)} times")
        flush(f"  classes: {sorted(allrealized)}")
        if witnesses:
            flush("  *** REFUTED *** first witnesses:")
            for w in witnesses[:5]: flush("   ",w)
            return True
    flush("\n[PHASE B] random sporadic large-speed primitive 5-sets")
    random.seed(12345)
    witnesses=[]; allr=set(); ntested=0
    for trial in range(200000):
        S=tuple(sorted(random.sample(range(1,400),5)))
        if not primitive(S): continue
        ntested+=1
        allr |= test_set_all_variants_all_opt(S,witnesses)
        if witnesses: break
    flush(f"  tested {ntested} random sets; distinct classes {len(allr)}; TARGET {len(witnesses)}")
    if witnesses:
        flush("  *** REFUTED ***",witnesses[:5]); return True
    flush("\n[PHASE C] structured 5-sets (AP, near-AP, covering-flavored, lacunary)")
    structured=set()
    for a in range(1,8):
        for d in range(1,12):
            structured.add(tuple(sorted(a+d*k for k in range(5))))
    for base in [2,3,5,7]:
        structured.add(tuple(sorted(base**k for k in range(5))))
    for w in [12,24,30,42,60,84,120,210,420,840]:
        for combo in combinations([1,2,3,4,5,6,7,8,9,10,11,13],4):
            S=tuple(sorted(combo+(w,)))
            if len(set(S))==5: structured.add(S)
    structured={s for s in structured if primitive(s)}
    witnesses=[]; allr=set()
    for S in structured:
        allr |= test_set_all_variants_all_opt(S,witnesses)
        if witnesses: break
    flush(f"  tested {len(structured)} structured sets; distinct classes {len(allr)}; TARGET {len(witnesses)}")
    if witnesses:
        flush("  *** REFUTED ***",witnesses[:5]); return True
    flush("\n[PHASE D] dense exact tau sweep (all crossing times + extra off-grid), sample of sets")
    sample=[s for s in combinations(range(1,16),5) if primitive(s)]
    witnesses=[]; allr=set(); seen=0
    extra_taus=[F(p,q) for q in range(2,40) for p in range(1,q) if gcd(p,q)==1]
    for S in sample:
        seen+=1
        taus=set(cand(list(S))) | set(extra_taus)
        for t in taus:
            for (bd,tb,tp) in VARIANTS:
                A=build(list(S),t,bd,tb,tp)
                if not is_tournament(A,5): continue
                k=inv5(A); allr.add(k)
                if k==TARGET:
                    witnesses.append((tuple(S),t,bd,tb,tp))
            if witnesses: break
        if witnesses: break
    flush(f"  swept {seen} sets x ({len(extra_taus)}+crossings) taus x {len(VARIANTS)} variants")
    flush(f"  distinct classes realized over ALL taus (off-grid included): {len(allr)}")
    flush(f"  classes: {sorted(allr)}")
    flush(f"  TARGET realized: {len(witnesses)}")
    if witnesses:
        flush("  *** REFUTED (even off-grid) ***",witnesses[:5]); return True
    flush("\n==> TARGET never realized in any phase. Claim CONFIRMED within search bound.")
    return False

if __name__=='__main__':
    main()
