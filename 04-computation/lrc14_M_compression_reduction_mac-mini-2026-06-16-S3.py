#!/usr/bin/env python3
"""
LRC(14) PROVE route, part 4: the COMPRESSION / bounded-max reduction.

The single blocker to closing LRC(14) via quantization is: M is scale-invariant
(so WLOG primitive) but max(S) is unbounded among primitive reps.  To get a
finite certificate we need: "a putative M-minimizer can be taken with max <= B"
for an explicit B = B(14).

This script tests the structural claim that drives such a reduction:

  CLAIM (compression):  if S has a 'large gap' (a speed v far above the rest),
  then we can REPLACE v by a smaller value without DECREASING M -- i.e. the
  minimizer is 'compressed'.  Concretely test the operation:
      replace the largest speed v_max by v_max - lcm-of-others-step, or by the
      smallest value keeping distinctness/primitivity, and see if M does not drop.

  If M is monotone under such compression, every minimizer is equivalent to a
  COMPRESSED config of bounded max, reducing LRC(14) to a finite check.

We ALSO directly measure: over random primitive sets, how does M correlate with
max?  Does large max force large M (away from 1/14)?  If min-M configs always
have SMALL max, that is strong evidence the extremal problem is compact.

DISPROVE-seed: the counterexample search is exactly the compressed-config search.
If compression holds, enumerating compressed configs (small max) is BOTH the
proof certificate AND the complete disproof search.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, random

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
def M(S):
    b=F(0); at=None
    for t in cand(S):
        v=g(S,t)
        if v>b: b=v; at=t
    return b,at
def prim(S): return reduce(gcd,S,0)==1

print("="*70)
print("PART A: do MIN-M configs have small max? (compactness evidence)")
print("="*70)
print("""
Sample primitive 13-sets with max ranging over several scales.  For each scale
band, record the MINIMUM M found.  If min-M does NOT keep dropping as max grows
(stays pinned at 1/14), the extremal problem is effectively compact: large-max
configs cannot beat the small-max tight AP.
""")
random.seed(7)
bands=[(13,20),(20,40),(40,80),(80,200),(200,1000),(1000,10000)]
print(f"{'max-band':>14} {'#prim':>8} {'min M':>14} {'float':>11} {'argmin (sorted)':>30}")
for lo,hi in bands:
    minM=None; argm=None; cnt=0
    trials = 3000 if hi<=80 else 1500
    for _ in range(trials):
        S=sorted(random.sample(range(1,hi+1),13))
        if max(S)<lo: continue
        if not prim(S):
            # divide out gcd to make primitive (scale-invariant: same M)
            gg=reduce(gcd,S,0); S=[x//gg for x in S]
        cnt+=1
        m,_=M(S)
        if minM is None or m<minM: minM=m; argm=S[:]
    print(f"{f'[{lo},{hi}]':>14} {cnt:>8} {str(minM):>14} {float(minM):>11.7f} {str(argm[:6])+'...':>30}")
print()
print("If every band's min M == 1/14 (or >), large max never helps => compact.")
print()

print("="*70)
print("PART B: compression operator -- replace v_max by smaller, track M")
print("="*70)
print("""
Take a config, repeatedly try to LOWER its maximum speed (to the smallest unused
value keeping 13 distinct positive primitive speeds) and check whether M stays
>= the original.  If M never drops below 1/14 under compression toward a small
window, the minimizer lives in a bounded window.
""")
def compress_max(S):
    """try to reduce max(S) to smallest unused value > 0, distinct."""
    S=sorted(set(S))
    used=set(S)
    # candidate replacements for the top element: any positive int not used,
    # below current max, keeping set distinct
    top=S[-1]
    rest=S[:-1]
    best=None
    for nv in range(1, top):
        if nv in used: continue
        cand_set=sorted(rest+[nv])
        if len(set(cand_set))!=13: continue
        best=cand_set; break  # smallest such
    return best

tests=[
 list(range(1,14)),
 [1,2,3,4,5,6,7,8,9,10,11,13,24],
 [1,2,3,4,5,6,7,8,9,10,11,13,36],
 sorted(random.sample(range(1,500),13)),
 sorted(random.sample(range(1,2000),13)),
]
for S0 in tests:
    S=S0[:]
    m0,_=M(S)
    chain=[max(S)]
    mvals=[m0]
    for _ in range(40):
        nx=compress_max(S)
        if nx is None or max(nx)>=max(S): break
        S=nx
        m,_=M(S)
        chain.append(max(S)); mvals.append(m)
        if max(S)<=14: break
    print(f"start max={S0[-1] if False else max(S0)} M0={m0}={float(m0):.5f} -> "
          f"compressed to max={max(S)} M={mvals[-1]}={float(mvals[-1]):.5f}  "
          f"{'M PRESERVED/UP' if mvals[-1]>=m0 else 'M DROPPED'}")
print()
print("Note: naive 'lower the top' compression is NOT guaranteed monotone; this")
print("is a heuristic probe.  The honest statement: a PROVABLE compression lemma")
print("(M non-decreasing under some explicit move toward bounded max) would close")
print("LRC(14).  Below we report whether ANY single-step lowering raised M<1/14.")
print()

print("="*70)
print("PART C: the forced floor M >= 1/(2 max) -- attempt at a proof certificate")
print("="*70)
print("""
PROOF ATTEMPT of M(S) >= 1/(2 max(S)):
Let V=max(S).  Consider tau ranging in [0, 1/(2V)].  Pick tau* = 1/(2V) ... no.
Better: the runners' positions v*tau for v in S at small tau are all in (0,1/2)
as long as tau < 1/(2V).  At tau = 1/(2V), the fastest runner V has ||V tau||=1/2,
and every slower v<V has v*tau in (0,1/2) so ||v tau|| = v/(2V) >= 1/(2V) (since
v>=1).  Hence min_v ||v tau|| >= 1/(2V) at tau=1/(2V).  Therefore
        M(S) = max_tau min_v ||v tau|| >= 1/(2V) = 1/(2 max(S)).   QED (clean!)
This is a GENUINE forced floor.  Verify it numerically as a sanity check, and
note: it gives M>=1/14 immediately whenever max(S) <= 7.  But 13 distinct
positive speeds need max>=13>7, so this exact floor never reaches 1/14 alone.
""")
# verify the certificate tau = 1/(2 max) gives min >= 1/(2 max)
ok=True
random.seed(55)
for _ in range(500):
    S=sorted(random.sample(range(1,300),13))
    V=max(S)
    t=F(1,2*V)
    val=g(S,t)
    if val < F(1,2*V):
        ok=False
        print("  CERTIFICATE FAILS at", S, "val", val); break
print(f"certificate tau=1/(2max) gives min_v||v t|| >= 1/(2max) on all 500 samples: {ok}")
# and M >= that
ok2=True
for _ in range(200):
    S=sorted(random.sample(range(1,150),13)); V=max(S)
    m,_=M(S)
    if m < F(1,2*V): ok2=False; print("M<1/(2max) at",S); break
print(f"M(S) >= 1/(2max) on all 200 samples: {ok2}")
print()
print("CONCLUSION: M >= 1/(2max) is PROVED (tau=1/(2max) certificate).  It reduces")
print("LRC(14) to configs with max>=14, but to REACH 1/14 we need max<=7 which is")
print("impossible for 13 distinct speeds.  So a SHARPER per-runner floor is needed:")
print("the floor should use the SECOND-largest etc., i.e. a Newman-style bound.")
print()
print("DONE part 4.")
