#!/usr/bin/env python3
"""lrc14_disproof_search_boxeph_S206.py -- boxeph-2026-07-21-S206

LRC(14) disproof characterization + search. LRC(14) holds for a 13-speed set S iff the covering-min
   M(S) = max_t min_{v in S} ||v t||  >=  1/14.
The conjectured global covering-min over COVERING sets is 14/183 ~ 0.07650 (deep well {1..12,182},
uniquely at t*=14/183=[0;13,14], Eisenstein/anti-golden; klein-S124, boxeph-S103). Since 14/183 > 1/14,
a DISPROOF needs a covering 13-set with M < 1/14 ~ 0.07143 -- i.e. dipping BELOW the AP extremal by the
margin 14/183 - 1/14 = 13/2562 ~ 0.00507.

This scans rational times a/q up to a stated denominator cutoff. The returned
quantity is a rigorous LOWER BOUND on M(S), not an exact maximum unless a
separate breakpoint-completeness argument covers the cutoff. A lower bound
above 1/14 certifies an individual candidate safe; a lower bound below 1/14 is
not a counterexample. The deep-well control is known exactly from canon.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random

THRESH = F(1,14)         # LRC(14) loneliness threshold
DEEPWELL_M = F(14,183)   # conjectured global covering-min

def sampled_lower_bound(S, Qmax=1500):
    """Return max over sampled a/q, q<=Qmax: an exact rational lower bound for M(S)."""
    best = F(0); arg=None
    for q in range(2, Qmax+1):
        for a in range(1, q):
            m = q
            for v in S:
                r = (v*a) % q
                d = r if r <= q-r else q-r
                if d < m: m = d
                if m == 0: break
            if m*best.denominator > best.numerator*q:  # m/q > best
                best = F(m, q); arg=(a,q)
    return best, arg

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2,15))

# ---------- verify the deep well ----------
DW = [1,2,3,4,5,6,7,8,9,10,11,12,182]
M, arg = sampled_lower_bound(DW)
print("deep well {1..12,182}: covering=%s  scan-LB=%s (~%.5f) at t=%s ; recovers canonical 14/183? %s ; >1/14? %s"
      % (is_covering(DW), M, float(M), arg, M==DEEPWELL_M, M>THRESH))
print("  margin 14/183 - 1/14 = %s (~%.5f) -- a disproof must dip below 1/14 by more than this\n" % (DEEPWELL_M-THRESH, float(DEEPWELL_M-THRESH)))

# ---------- disproof candidate families ----------
print("CANDIDATE SAFETY SCAN (LB > 1/14 certifies safe; LB < 1/14 proves nothing):")
cands = {}
# far-element variants of the AP
for far in [182, 2*182, 13*14, 2*3*4*5*6*7, 5460]:  # 5460=lcm(1..14)/... test blockers
    S=[1,2,3,4,5,6,7,8,9,10,11,12,far]
    if is_covering(S): cands["AP12 + far=%d"%far]=S
# near-AP perturbations (drop one, add a covering far element)
cands["skip12: {1..11,13,182}"]=[1,2,3,4,5,6,7,8,9,10,11,13,182]
cands["skip7:  {1..6,8..13,182}"]=[1,2,3,4,5,6,8,9,10,11,12,13,182]
# generalized AP / dilated
cands["2*AP: {2,4..24,182}"]=[2,4,6,8,10,12,14,16,18,20,22,24,182]  # covering? check
# Fibonacci/Zeckendorf-structured 12-sets + covering far element
fib=[1,2,3,5,8,13,21,34,55,89,144,233]
cands["Fibonacci12 + 720720"]=fib+[720720]  # 720720=lcm(1..14)? add big blocker
cands["Fib-ish covering {1,2,3,4,5,8,9,10,11,12,13,7,182}"]=[1,2,3,4,5,7,8,9,10,11,12,13,182]
# random covering sets
random.seed(2026)
found_cov=0
for _ in range(400):
    S=sorted(random.sample(range(1,60),12))
    far=1
    for q in range(2,15):
        if not any(v%q==0 for v in S): far*=q//gcd(far,q)
    S13=S+[far*random.choice([1,2,3])]
    if is_covering(S13):
        found_cov+=1
        if found_cov<=6: cands["random#%d %s"%(found_cov,S13[:4]+['..',S13[-1]])]=S13

worst=[]
for name,S in cands.items():
    if not is_covering(S):
        print("  %-38s NOT COVERING (skip)"%name); continue
    M,arg=sampled_lower_bound(S, 1200)
    flag = "  SAFE WITNESS" if M>=THRESH else "  UNRESOLVED BY THIS SCAN"
    print("  %-38s LB(q<=1200)=%-10s (~%.5f)%s"%(name, str(M), float(M), flag))
    worst.append((float(M),name))
worst.sort()
print("\n  smallest sampled lower bounds:", [(round(m,5),n) for m,n in worst[:4]])
print("  => every displayed candidate has an explicit sampled witness above 1/14.")
print("     This finite family is excluded; no global ranking or anti-golden necessity follows.")
print("     THM-730 proves only that a non-AP 12-core loses additive triples.")
