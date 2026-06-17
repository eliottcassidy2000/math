#!/usr/bin/env python3
"""
ANGLE A (boundary): the simple certificate "DIP <= SLACK => M(S)>=1/14" is
SUFFICIENT but NOT NECESSARY. Random search found a covering 13-set where
DIP > SLACK yet M(S)=4/37 >= 1/14. Goal here:

(1) Confirm the gap-monotonicity FACT: M(S) <= M(A) always (adding a runner can
    only lower the gap), so DIP = M(A)-M(S) >= 0 always. The certificate's real
    content is: how big can the dip be? It is bounded by M(A)-1/14, and the
    EASY core has M(A) >= 1/q. The certificate as literally stated (dip<=1/q-1/14)
    is just the rearrangement of M(S) >= 1/14 GIVEN M(A) = 1/q EXACTLY.
    => When M(A) > 1/q (core gap STRICTLY exceeds the trivial witness), the slack
    is larger than 1/q-1/14 and the certificate-as-stated understates the slack.

(2) The CORRECT certificate: M(S) >= 1/14  <=>  DIP <= M(A) - 1/14.
    Test THIS across all covering 13-sets. The "slack" should be M(A)-1/14, not
    1/q-1/14. With q the SMALLEST uncovered modulus, 1/q is only a LOWER bound on
    M(A); the true slack M(A)-1/14 is >=.

(3) Decompose: is the breaking case explained by M(A) > 1/q? Report
    M(A) vs 1/q for the random breaker.

(4) Re-run the worst DIP/(M(A)-1/14) ratio: this should NEVER reach 1 if LRC14
    holds (it equals 1 exactly when M(S)=1/14). Find the max and the minimizing
    M(S).

stdlib only.
"""
from fractions import Fraction as F
from itertools import combinations
import random, functools
print = functools.partial(print, flush=True)

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
def covers(S, q): return any(v % q == 0 for v in S)
def uncovered(S, qmax=14): return [q for q in range(2, qmax+1) if not covers(S, q)]
def smallest_uncovered(A, qmax=13):
    u = uncovered(A, qmax); return u[0] if u else None
ONE14 = F(1, 14)

# ---- (1)/(3) the random breaker, dissected ----
print("="*78)
print("DISSECTION of the random 'breaker' S where DIP>1/q-1/14 but M(S)>=1/14")
print("="*78)
S = [1, 4, 6, 8, 10, 13, 14, 15, 18, 22, 24, 25, 27]
A = [x for x in S if x != 14]
MS, tS = M(S); MA, tA = M(A)
q = smallest_uncovered(A,13)
print(f"S = {S}")
print(f"A = S\\{{14}} = {A}")
print(f"M(S) = {MS} ({float(MS):.6f}) at tau*={tS};  >=1/14 ? {MS>=ONE14}")
print(f"M(A) = {MA} ({float(MA):.6f}) at tau*={tA}")
print(f"smallest uncovered q of A = {q};  1/q = {F(1,q)} ({float(F(1,q)):.6f})")
print(f"  Is M(A) > 1/q (core gap beats trivial witness)?  {MA > F(1,q)}  "
      f"[M(A)={MA} vs 1/q={F(1,q)}]")
print(f"  stated slack 1/q-1/14 = {F(1,q)-ONE14} ({float(F(1,q)-ONE14):.6f})")
print(f"  TRUE slack M(A)-1/14  = {MA-ONE14} ({float(MA-ONE14):.6f})")
print(f"  DIP = M(A)-M(S) = {MA-MS} ({float(MA-MS):.6f})")
print(f"  => DIP <= 1/q-1/14 (stated)? {MA-MS <= F(1,q)-ONE14}")
print(f"  => DIP <= M(A)-1/14 (true)?  {MA-MS <= MA-ONE14}  (this is just M(S)>=1/14)")
print("CONCLUSION: the stated certificate fails ONLY because 1/q UNDERSTATES M(A).")
print("The core is EASIER than the trivial witness shows; true slack is bigger.")

# ---- (4) exhaustive + random: worst DIP/(M(A)-1/14), and whether it hits 1 ----
def scan(label, sets_iter):
    worst = (F(0), None); minMS = (F(1), None); broke = []; n=0
    stated_fail = 0
    for A, S in sets_iter:
        MS,tS = M(S); MA,tA = M(A)
        n += 1
        if MS < minMS[0]: minMS = (MS,(A,S,MS,tS,MA))
        if MS < ONE14: broke.append((A,S,MS,tS))
        true_slack = MA - ONE14
        dip = MA - MS
        # certificate-as-stated check
        q = smallest_uncovered(A,13)
        if q is not None and dip > F(1,q)-ONE14:
            stated_fail += 1
        if true_slack > 0:
            r = dip/true_slack
            if r > worst[0]: worst = (r,(A,S,MS,MA,dip,true_slack))
    print(f"\n[{label}] tested {n} covering 13-sets")
    print(f"  GLOBAL MIN M(S) = {minMS[0]} ({float(minMS[0]):.6f}) >=1/14 ? {minMS[0]>=ONE14}")
    if minMS[1]:
        A,S,ms,ts,ma = minMS[1]
        print(f"     at S={S}, tau*={ts}, M(A)={ma}")
    print(f"  WORST DIP/(M(A)-1/14) = {worst[0]} ({float(worst[0]):.6f})  "
          f"[=1 exactly iff M(S)=1/14]")
    if worst[1]:
        A,S,ms,ma,dip,sl = worst[1]
        print(f"     at S={S}")
        print(f"     M(A)={ma}, M(S)={ms}, dip={dip}, true_slack={sl}")
    print(f"  certificate-AS-STATED (dip<=1/q-1/14) FAILED on {stated_fail}/{n} "
          f"(but all still M(S)>=1/14 if breaks=0)")
    print(f"  sets actually breaking LRC14 (M(S)<1/14): {len(broke)}")
    for c in broke[:5]:
        print(f"     S={c[1]} M(S)={c[2]} tau={c[3]}")
    return worst, minMS, broke, stated_fail

# exhaustive {1..13}+14m
def gen_exhaustive():
    base = list(range(1,14))
    for A in combinations(base,12):
        A=list(A)
        for m in range(1,25):
            w=14*m; S=sorted(set(A)|{w})
            if len(S)!=13: continue
            if uncovered(S,14): continue
            if smallest_uncovered(A,13) is None: continue
            yield A,S

# random covering 13-sets, broad
def gen_random(seed, count, hi):
    random.seed(seed); got=0; tries=0
    while got<count and tries<400000:
        tries+=1
        parked=14*random.randint(1,hi//14)
        pool=[x for x in range(1,hi+1) if x%14!=0]
        rest=random.sample(pool,12)
        S=sorted(set(rest)|{parked})
        if len(S)!=13: continue
        if uncovered(S,14): continue
        A=sorted(set(S)-{parked})
        if len(A)!=12: continue
        if smallest_uncovered(A,13) is None: continue
        got+=1
        yield A,S

scan("EXHAUSTIVE {1..13}+14m", gen_exhaustive())
scan("RANDOM hi=28", gen_random(1, 3000, 28))
scan("RANDOM hi=42", gen_random(2, 1500, 42))
scan("RANDOM hi=56", gen_random(3, 800, 56))
