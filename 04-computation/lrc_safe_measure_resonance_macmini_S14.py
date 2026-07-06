#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S14 -- reasoning about (G) gap-emptiness: the safe-measure /
Newman-Fourier decomposition and WHY the 2/25 gap is harder than the 1/14 threshold.

SETUP. For speeds S = {v_1,...,v_m} and level beta, the beta-SAFE SET is
  Safe(S,beta) = { t in [0,1) : ||v_i t|| >= beta for all i }.
Then M(S) = max_t min_i ||v_i t|| satisfies:
  M(S) >= beta  <=>  Safe(S,beta) nonempty ;   M(S) < beta  <=>  beta-arcs COVER.
So safe_measure(S,beta) = |Safe(S,beta)| = INT_0^1 prod_i (1 - g(v_i t)) dt,
where g(u)=1[||u||<beta] is the danger arc (measure 2beta).

NEWMAN-FOURIER IDENTITY (exact): expanding the product,
  safe_measure = (1 - 2beta)^m  +  SUM over nontrivial resonances (k: sum k_i v_i=0)
                 of  prod_i ghat(k_i),      ghat(k)=sin(2 pi k beta)/(pi k).
The FIRST term (1-2beta)^m is the INDEPENDENT baseline (what you'd get with no
additive relations).  The AP has MAXIMAL resonances -> drives safe_measure to 0
(its arcs cover) -> M(AP)=1/13 < 2/25.  (G)+rigidity = "only the AP has enough
resonance to zero the safe measure at beta=2/25."

THIS SCRIPT computes safe_measure(S,beta) EXACTLY (arc union over breakpoints),
compares to the baseline (1-2beta)^m, and reports the Bonferroni partial sums
S_1,S_2,S_3,... to show WHICH ORDER the covering needs -- pairs for the 1/14
threshold (the 'seven wall'), higher order for the 2/25 gap.
"""
from fractions import Fraction as F
from math import gcd, sin, pi
from functools import reduce
from itertools import combinations
import sys
sys.path.insert(0,'04-computation')
from lonely_profile import profile

def log(m=""): print(m, flush=True)

def M_exact(S):
    S=sorted(set(abs(x) for x in S))
    for cap in (14,11,8,6,4,3,2):
        p=profile(S,F(1,cap)); m=p.M()
        if m is not None: return m

def safe_measure(S, beta):
    """|{t in [0,1): ||v_i t||>=beta for all i}| EXACTLY, beta a Fraction.
       Danger of v_i = union over j of ( (j-beta)/v_i, (j+beta)/v_i ), j=0..v_i.
       Safe = complement.  Compute via merging danger intervals (exact Fractions)."""
    intervals=[]
    for v in S:
        # ||v t||<beta <=> exists integer j with |v t - j|<beta <=> t in ((j-beta)/v,(j+beta)/v)
        j=0
        while F(j-beta, v) < 1:
            lo=F(j-beta, v); hi=F(j+beta, v)
            lo=max(lo,F(0)); hi=min(hi,F(1))
            if lo<hi: intervals.append((lo,hi))
            j+=1
    # merge
    intervals.sort()
    danger=F(0); cur_lo=None; cur_hi=None
    for lo,hi in intervals:
        if cur_hi is None:
            cur_lo,cur_hi=lo,hi
        elif lo<=cur_hi:
            if hi>cur_hi: cur_hi=hi
        else:
            danger+=cur_hi-cur_lo; cur_lo,cur_hi=lo,hi
    if cur_hi is not None: danger+=cur_hi-cur_lo
    return F(1)-danger

def bonferroni(S, beta, upto=4):
    """S_k = sum over k-subsets T of measure(intersection of dangers)."""
    def inter_measure(T):
        # {t: ||v_i t||<beta for all i in T}: intersect danger interval-unions
        segs=[(F(0),F(1))]
        for v in T:
            arcs=[]
            j=0
            while F(j-beta,v)<1:
                lo=max(F(j-beta,v),F(0)); hi=min(F(j+beta,v),F(1))
                if lo<hi: arcs.append((lo,hi))
                j+=1
            newsegs=[]
            for (a,b) in segs:
                for (c,d) in arcs:
                    lo=max(a,c); hi=min(b,d)
                    if lo<hi: newsegs.append((lo,hi))
            segs=newsegs
            if not segs: break
        return sum((b-a for a,b in segs), F(0))
    Ss=[]
    m=len(S)
    for k in range(1,upto+1):
        tot=F(0)
        for T in combinations(S,k):
            tot+=inter_measure(T)
        Ss.append(tot)
    return Ss

# ---- families ----
AP12=list(range(1,13))                       # the AP, M=1/13
apex12=list(range(1,12))+[24]                # doubled apex {1..11,24}, M=2/25
block=[1,2,3,5,7,8,9,10,11,12,17,19]         # block lift, M=2/25
generic=[1,2,3,4,5,6,7,8,9,10,11,23]         # covering-ish, likely loose

log("="*70)
log("SAFE MEASURE vs INDEPENDENT BASELINE (1-2beta)^12,  beta near the gap")
log("="*70)
for name,S in [("AP {1..12}",AP12),("doubled-apex {1..11,24}",apex12),
               ("block lift",block),("generic {1..11,23}",generic)]:
    M=M_exact(S)
    log(f"\n{name}:  M(S) = {M} = {float(M):.5f}")
    for beta in [F(1,13),F(3,38),F(2,25),F(1,12),F(2,23)]:
        sm=safe_measure(S,beta); base=(1-2*beta)**12
        cov = "COVER (M<beta)" if sm==0 else f"safe={float(sm):.5f}"
        log(f"   beta={str(beta):>6}={float(beta):.5f}: {cov:22s} baseline(1-2b)^12={float(base):.5f}")

log("\n"+"="*70)
log("BONFERRONI ORDER: threshold 1/14 (7-wall) vs gap 2/25")
log("="*70)
for beta,label in [(F(1,14),"threshold 1/14"),(F(2,25),"gap 2/25")]:
    log(f"\nbeta = {beta} ({label}); danger per runner = {float(2*beta):.4f}; "
        f"S_1 = 12*{float(2*beta):.4f} = {float(12*2*beta):.4f}")
    for name,S in [("AP {1..12}",AP12),("doubled-apex",apex12)]:
        Ss=bonferroni(S,beta,upto=4)
        # partial alternating sums 1 - S1 + S2 - S3 + S4 (inclusion-exclusion partials)
        partials=[]
        acc=F(1); sign=-1
        for k,Sk in enumerate(Ss,1):
            acc+=sign*Sk; partials.append(float(acc)); sign*=-1
        log(f"   {name}: S1={float(Ss[0]):.3f} S2={float(Ss[1]):.3f} S3={float(Ss[2]):.3f} S4={float(Ss[3]):.3f}")
        log(f"      IE partials 1-S1(+S2)(-S3)(+S4): " + " ".join(f"{p:.3f}" for p in partials)
            + f"   [true safe={float(safe_measure(S,beta)):.4f}]")

log("\nREADING: at beta=1/14 the pair term S_2 already dominates (7-wall closes ~c<=7);")
log("at beta=2/25 the alternating IE partials swing WIDE -- covering needs high-order")
log("resonance (the AP's additive energy), so no low-order Bonferroni certificate exists.")
log("This is WHY the gap (G) is strictly harder than the LRC(14) threshold.")
