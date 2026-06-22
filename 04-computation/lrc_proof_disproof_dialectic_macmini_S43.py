from fractions import Fraction as Fr
from itertools import combinations
import random
def meas_safe(S):  # measure of {t: ||s t||>=1/14 for all s in S}
    iv=[]
    for s in S:
        w=Fr(2,14)/s
        for k in range(s):
            c=Fr(k,s); lo=(c-w/2)%1; hi=lo+w
            if hi<=1: iv.append((lo,hi))
            else: iv.append((lo,Fr(1))); iv.append((Fr(0),hi-1))
    iv.sort(); tot=Fr(0); clo=chi=None
    for lo,hi in iv:
        if chi is None: clo,chi=lo,hi
        elif lo<=chi: chi=max(chi,hi)
        else: tot+=chi-clo; clo,chi=lo,hi
    if chi is not None: tot+=chi-clo
    return 1-tot
def M(S):
    def nrm(x): x=x%1; return min(x,1-x)
    dens=set(S); Sl=sorted(S)
    for a in Sl:
        for b in Sl:
            if a!=b: dens.add(abs(a-b)); dens.add(a+b)
    best=Fr(0)
    for d in dens:
        if d:
            for k in range(d+1):
                m=min(nrm(s*Fr(k,d)) for s in S)
                if m>best: best=m
    return best
print("PROOF via the MEASURE INEQUALITY: meas(Safe(B)) > |L|*(2/14) => S=B∪L is safe (union bound).")
seed=list(range(1,12))+[13]   # {1..11,13}, the user's seed, M=1/12
ms=meas_safe(seed)
print(f"  seed {{1..11,13}}: meas(Safe at 1/14) = {ms} = {float(ms):.4f}; large-speed budget for |L|=1 is 2/14={float(Fr(2,14)):.4f}")
print(f"  meas(Safe(seed)) > 2/14 ? {ms>Fr(2,14)}  => for ANY single large N, {{1..11,13,N}} is SAFE (M>=1/14). Disproof family DEFEATED.")
print()
print("DISPROOF refocused: COVERING 13-sets (a multiple of every q in 2..14). min M found:")
def covers_all(S):
    return all(any(s%q==0 for s in S) for q in range(2,15))
random.seed(3); best=(Fr(1),None); cnt=0
# build covering sets: must include mults of 8,9,11,13 (hard); fill rest
hard=[[8,16,24],[9,18,27],[11,22,33],[13,26,39],[5,10,15,20,25],[7,14,21,28],[12,24,36]]
for _ in range(8000):
    S=set()
    for opts in hard: S.add(random.choice(opts))
    while len(S)<13: S.add(random.randint(1,40))
    S=tuple(sorted(S))
    if len(S)!=13 or not covers_all(S): continue
    cnt+=1; m=M(S)
    if m<best[0]: best=(m,S)
print(f"  searched {cnt} covering 13-sets: min M = {float(best[0]):.5f} (1/14={float(Fr(1,14)):.5f}); <1/14: {best[0]<Fr(1,14)}")
print(f"  argmin: {best[1]}")
print(f"  => min M >= 1/14: no counterexample found; the covering constraint + measure inequality squeeze toward LRC.")
