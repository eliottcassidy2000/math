#!/usr/bin/env python3
"""
ANGLE B (followup): pin down the BINDING-PAIR LAW precisely, and test
the M=1/2 single-runner exception + the others-clear independence.
"""
from fractions import Fraction as F
from itertools import combinations
import random

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1,2) else 1 - r
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
def all_optima(S):
    Mv,_=M(S); return Mv, sorted(t for t in cand(S) if g(S,t)==Mv)
def binding(S,t,gap): return [v for v in sorted(set(S)) if nrm(v*t)==gap]
def has_pair(S,t,gap):
    bind=binding(S,t,gap)
    for a,b in combinations(bind,2):
        for d in (a+b,abs(b-a)):
            if d>0 and (t*d).denominator==1: return (a,b,d)
    return None

print("="*72)
print("EXCEPTION CHARACTERIZATION: M=1/2 single-binding-runner case.")
print("="*72)
print("M=1/2 means every runner is at distance 1/2 at tau=1/2: v*1/2 ~ 1/2")
print("=> v ODD for all v. So M=1/2 <=> all speeds odd; tau=1/2 pins it with")
print("ALL runners binding (the whole set), trivially a pair too.")
random.seed(7)
allodd_M_half=0; allodd_total=0; mixed_half=0
for _ in range(3000):
    n=random.randint(1,8); S=random.sample(range(1,40),n)
    Mv,_=M(S)
    if all(v%2==1 for v in S):
        allodd_total+=1
        if Mv==F(1,2): allodd_M_half+=1
    if Mv==F(1,2) and not all(v%2==1 for v in S): mixed_half+=1
print(f"  all-odd sets: {allodd_total}, of which M=1/2: {allodd_M_half}")
print(f"  M=1/2 with a NON-odd speed present: {mixed_half}")
print("  => M=1/2 <=> all speeds odd (the single-runner/half-int case is")
print("     subsumed: tau=1/2 makes EVERY odd runner bind, so >=2 bind).")

print()
print("="*72)
print("BINDING-PAIR LAW for COVERING CORES with a single multiple m*14.")
print("="*72)
print("Test set: {1..13} plus one extra speed w=m*14, and {1..k,13,84} style.")
print("Question: does the small-speed/multiple-of-14 pair (v, w) bind, with")
print("tau = k/(v+w) for the SMALLEST gap-determining crossing?")
def covering_probe(S):
    Mv,opts=all_optima(S)
    w=[v for v in sorted(set(S)) if v%14==0]
    rows=[]
    for t in opts:
        b=binding(S,t,Mv); hp=has_pair(S,t,Mv)
        w_in=[v for v in b if v%14==0]
        rows.append((t,b,w_in,hp))
    return Mv,w,rows
tests = {
  "{1..13,14}":     list(range(1,14))+[14],
  "{1..13,28}":     list(range(1,14))+[28],
  "{1..13,42}":     list(range(1,14))+[42],
  "{1..13,56}":     list(range(1,14))+[56],
  "{1..13,70}":     list(range(1,14))+[70],
  "{1..13,84}":     list(range(1,14))+[84],
  "{1..11,13,84}":  list(range(1,12))+[13,84],
  "{1..11,13,98}":  list(range(1,12))+[13,98],
}
for name,S in tests.items():
    Mv,w,rows=covering_probe(S)
    print(f"\n{name}: M={Mv} ({float(Mv):.5f}), mult-of-14 = {w}")
    for t,b,w_in,hp in rows:
        d_used = hp[2] if hp else None
        print(f"   tau={t}: binding={b}, mult14-binds={w_in or 'NO'}, "
              f"pair={hp}, denom={t.denominator}")
    # is denominator == w + small? characterize
    for t,b,w_in,hp in rows:
        if w_in and hp:
            a,bb,d=hp
            other=a if bb%14==0 else bb
            print(f"   --> LAW CHECK: pair=({other},{w[0]}) sum={other+w[0]}, "
                  f"tau denom={t.denominator}, sum==denom? {other+w[0]==t.denominator}")

print()
print("="*72)
print("OTHERS-CLEAR INDEPENDENCE: can a pair-crossing have gap>=1/14 but be")
print("KILLED by another runner being closer? Show the constraint is REAL")
print("(i.e. binding-pair gap alone is NOT sufficient; need clearance).")
print("="*72)
def pair_gap_only(S):
    """Best gap among pair-crossings IGNORING other runners (just the pair),
       vs the true M (which respects all runners)."""
    Sset=sorted(set(S)); best_naive=F(0); best_true=F(0); naive_at=None
    for i in range(len(Sset)):
        for j in range(i+1,len(Sset)):
            a,b=Sset[i],Sset[j]
            for d in (a+b,b-a):
                if d<=0: continue
                k=1
                while F(k,d)<=F(1,2):
                    t=F(k,d)
                    pg=min(nrm(a*t),nrm(b*t))   # gap of JUST the pair
                    tg=g(S,t)                     # true gap (all runners)
                    if pg>best_naive: best_naive=pg; naive_at=(a,b,t)
                    if tg>best_true: best_true=tg
                    k+=1
    return best_naive,naive_at,best_true
demo = {
  "{1,2,3,...,13}": list(range(1,14)),
  "{1,5,7,11}": [1,5,7,11],
  "{2,3,11}": [2,3,11],
  "{1,6,10,15}": [1,6,10,15],
}
for name,S in demo.items():
    bn,nat,bt=pair_gap_only(S)
    Mv,_=M(S)
    print(f"\n{name}: best pair-only gap = {bn} at {nat}")
    print(f"   true M (others respected) = {bt} = {Mv}")
    print(f"   others-clear MATTERS here? {bn>bt}  "
          f"(pair-only overstates M by {bn-bt} when True)")

print()
print("="*72)
print("DECISIVE STRESS: is true M ALWAYS = max over pair-crossings of the")
print("FULL gap g(S,tau) (i.e. min over ALL runners)? This is the actual")
print("reduction we claim: M = max_{pair crossings tau} g(S,tau).")
print("(plus tau=1/2 for the all-odd M=1/2 case.) Counterexamples?")
print("="*72)
random.seed(2024)
viol=0; tested=0; ex=[]
for _ in range(8000):
    n=random.randint(2,9); S=sorted(set(random.sample(range(1,60),n)))
    if len(S)<2: continue
    Mv,_=M(S)
    # max over pair crossings of full gap, plus 1/2
    best=F(0)
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            a,b=S[i],S[j]
            for d in (a+b,b-a):
                if d<=0: continue
                k=1
                while F(k,d)<=F(1,2):
                    val=g(S,F(k,d))
                    if val>best: best=val
                    k+=1
    best=max(best,g(S,F(1,2)))  # all-odd safety
    tested+=1
    if best!=Mv:
        viol+=1
        if len(ex)<6: ex.append((S,Mv,best))
print(f"\nTested {tested} configs (n=2..9).")
print(f"Configs where max-over-pair-crossings(+1/2) != true M: {viol}")
if ex:
    for S,Mv,best in ex: print(f"  S={S}: M={Mv}, pair-max={best}")
else:
    print("==> ZERO violations. The reduction")
    print("    M(S) = max( g(S,1/2),  max_{pairs(a,b),k} g(S, k/(a+-b)) )")
    print("    holds on all sampled configs. LRC(14) <=> this max >= 1/14,")
    print("    checkable over O(n^2) pairs x O(d) crossings, each an O(n)")
    print("    clearance evaluation.")
