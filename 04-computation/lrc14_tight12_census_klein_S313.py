#!/usr/bin/env python3
"""
lrc14_tight12_census_klein_S313.py
==================================
klein-2026-07-14-S313 (owner: prove the n=12 sporadic branch empty).

DIRECT census of tight (M=1/13) primitive 12-speed families = the extremal LRC instances at n=12.
The n=12 sporadic branch is EMPTY iff the AP {1..12} (+ dilations) is the unique primitive tight 12-set
(= Tao's optimistic uniqueness conjecture; Goddyn-Wong "Tight instances of the lonely runner", Integers 6
(2006) A38 constructs non-AP tight instances at OTHER n via v->2v acceleration under an arithmetic condition).

RESULTS (evidence, NOT proof; the exhaustive proof is codex's H6 owner-feasibility program, scale by scale):
 - EXHAUSTIVE window {1..16}: the AP {1,..,12} is the UNIQUE primitive tight (M=1/13) 12-set.
 - Goddyn-Wong single accelerations {1..12} with v->2v: ALL loose (M in {1/10,2/21,1/11,2/23,1/12,2/25}) --
   the GW arithmetic condition ("v shares a common factor with every integer in [13-v, 25-2v]") holds for NO
   v at n=12, so GW produces no new tight instance here.
 - ADVERSARIAL ~256k primitive 12-sets (double accel, H1 faces<=60, multi-replace, random<=40): 0 non-AP tight.
 - The nearest M above 1/13 is 2/25 (v=12->24) => the Farey gap (1/13, 2/25) = codex CRUX (C), confirmed.
Cross-checks codex's H6 ledger (scales 1..28 empty) from a direct-enumeration angle. Uniform emptiness OPEN.
"""
import itertools, random
from math import gcd
from functools import reduce
from fractions import Fraction as Fr
def g(xs): return reduce(gcd,xs)
def Mex(S,Qs):
    best=Fr(0)
    for q in range(2,Qs):
        for a in range(1,(q+2)//2):
            if gcd(a,q)!=1: continue
            m=min(min((a*s)%q,q-(a*s)%q) for s in S); v=Fr(m,q)
            if v>best:
                best=v
                if best>Fr(1,13): return best   # early-out: loose
    return best
one=Fr(1,13)
print("== exhaustive window {1..K}: primitive tight (M=1/13) 12-sets ==")
for K in [13,14,15,16]:
    tt=[list(S) for S in itertools.combinations(range(1,K+1),12) if g(S)==1 and Mex(S,2*13*K)==one]
    print("  K=%2d: %d tight -> %s"%(K,len(tt),tt))
print("== Goddyn-Wong v->2v of {1..12} ==")
base=list(range(1,13))
for v in range(1,13):
    S=sorted(set([2*v if x==v else x for x in base]))
    if len(S)==12: print("  v=%2d %s M=%s"%(v,S,Mex(S,200)))
print("== adversarial non-AP primitive tight hunt ==")
random.seed(5); found=[]; n=0
def test(S):
    global n
    S=sorted(set(S))
    if len(S)!=12 or g(S)!=1 or min(S)<1: return
    n+=1
    if Mex(S,60)==one and S!=base: found.append(S)
for v in range(1,13):
    for w in range(1,13): test([(2*x if x in (v,w) else x) for x in base])
for j in range(12):
    for nv in range(13,61): test(base[:j]+[nv]+base[j+1:])
for _ in range(150000):
    S=base[:]
    for _ in range(random.randint(1,3)): S[random.randint(0,11)]=random.randint(1,60)
    test(S)
for _ in range(150000): test(random.sample(range(1,41),12))
print("  tested ~%d primitive 12-sets; non-AP tight found: %d %s"%(n,len(found),found[:5]))
