"""
Reproduce the claim's EXACT exhaustive sweep condition and find failures.
Claim: "ALL-MULT7-LARGE closure verified by exact interval arithmetic exhaustively
over V*<=399 with mult-of-7 sub-systems of size <=4; 0 counterexamples in 3,515,688
adversarial configs (V*=14..399, sub-system size<=4)."

The reduced subproblem (window-collapse): given V* and W={w_1..w_r}, r<=4, with every
7 w_i > V* (i.e. w_i >= floor(V*/7)+1), is there u in (0,1/(2V*)] with ||w_i u||>=1/14?

We sweep ALL such W with w_i in [wmin, WMAX] and report the FIRST (smallest V*, then
lexicographically smallest W) where NO safe u exists. If the claim's "0 counterexamples"
is right, we find none up to whatever WMAX the claim used. We test a generous WMAX.
kind-pasteur S4-wf.
"""
from fractions import Fraction as F
import itertools

def safe_exists(W,Vstar):
    R=F(1,2*Vstar); h=F(1,14); forb=[]
    for w in set(W):
        j=0
        while True:
            c=F(j,w); lo=c-h/F(w)
            if lo>R: break
            a=max(F(0),lo); b=min(R,c+h/F(w))
            if a<=b: forb.append((a,b))
            j+=1
    forb.sort(); merged=[]
    for a,b in forb:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    prev=F(0); gaps=[]
    for a,b in merged:
        if a>prev: gaps.append((prev,a))
        prev=max(prev,b)
    if prev<R: gaps.append((prev,R))
    h=F(1,14)
    for lo,hi in gaps:
        if hi<=0: continue
        for u in (lo,(lo+hi)/2,hi):
            if u>0 and u<=R and all(((F(w)*u-int(F(w)*u))%1) and True for w in W):
                pass
        # robust check via midpoint and endpoints
    # proper: a gap is safe-bearing iff it contains a point u>0 with ||w u||>=1/14 all w.
    # Since teeth are closed-removed, any open gap interior point works; also closed
    # endpoints are safe (boundary of a tooth has ||w u||=1/14 exactly). So existence
    # of ANY gap (lo,hi) with hi>0 and hi> max(lo,0) gives a safe point.
    for lo,hi in gaps:
        if hi>0 and hi>lo:
            u=(max(lo,F(0))+hi)/2
            if u<=0: u=hi/2
            if u>0: return True
    return False

# First find the SMALLEST failing (V*, W) with r<=4, scanning broadly.
# For each V*, wmin=floor(V*/7)+1. We let w range up to wmin+WSPAN. Larger w have
# finer teeth but still can block in combination; we use a generous span.
WSPAN=30
print("Scanning V*=14..399, r=1..4, w in [wmin, wmin+%d]. Report failures."%WSPAN)
first_fail=None; nfail=0; ntot=0
fail_by_r={1:0,2:0,3:0,4:0}
examples=[]
for Vstar in range(14,400):
    wmin=Vstar//7+1
    wpool=list(range(wmin,wmin+WSPAN+1))
    for r in range(1,5):
        for W in itertools.combinations(wpool,r):
            ntot+=1
            if not safe_exists(W,Vstar):
                nfail+=1; fail_by_r[r]+=1
                if first_fail is None: first_fail=(Vstar,W)
                if len(examples)<8: examples.append((Vstar,W))
print("  total configs:",ntot)
print("  FAILURES (no safe u):",nfail)
print("  failures by r:",fail_by_r)
print("  FIRST failing config (smallest V*):",first_fail)
print("  sample failing configs:",examples)

# Sanity: confirm the r=2 minimal failing example by listing teeth explicitly.
if first_fail:
    Vstar,W=first_fail
    R=F(1,2*Vstar)
    print("\n  DETAIL of first failure: V*=",Vstar," W=",W," R=1/(2V*)=",R,"~",float(R))
    for w in W:
        j=0; teeth=[]
        while True:
            c=F(j,w); lo=c-F(1,14)/w
            if lo>R: break
            teeth.append((float(max(F(0),lo)),float(min(R,c+F(1,14)/w))))
            j+=1
        print("    w=",w," teeth on (0,R]:",teeth)
