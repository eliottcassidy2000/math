#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S32 (HYP-4622) -- the TWO-MODULUS completion of kps's mod-25 reduction.

kps-S41 (HYP-4567): (G) reduces to "every >=3-defect 12-family w/o a multiple of 25 admits a
clearing rotation c in (Z/25)* with all c*v_i mod 25 in [2,23] (=> M>=2/25)".
CHARACTERIZATION (this script): NO clearing rotation mod 25 <=> the unit-speeds hit ALL 10
+/- classes of (Z/25)* (a FULL TRANSVERSAL mod 25).  So kps's fact FAILS exactly for full
transversals -- and my S7 showed gap-member => full transversal.  So the full transversals
are the GENUINE RESIDUAL; mod-25 rotation cannot clear them.

TEST: among near-tight >=3-defect families that ARE full transversals mod 25 (mod-25 fails),
  (a) is M ever strictly in the gap (1/13, 2/25)?  (b) what witness clears each to M>=2/25?
If the two moduli {25 (miss-a-class), and small-q/13 for the transversals} TOGETHER clear all,
(G) closes as a two-modulus rigidity: 25 necessary (full transversal) + 13 pinning => AP only.
Uses the S20-FIXED exact fast-M.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations, product

GAP_LO, GAP_HI = F(1,13), F(2,25)

def Mfast(S):
    S=sorted(set(S)); Q=set()
    for v in S: Q.add(2*v)
    for a,b in combinations(S,2): Q.add(a+b); Q.add(abs(a-b))
    Q.discard(0); best=F(0)
    for q in Q:
        for a in range(1,q):
            mn=min(min((v*a)%q,q-((v*a)%q)) for v in S)
            val=F(mn,q)
            if val>best: best=val
    return best

def clears_mod(S,p,mu):
    """exists c in (Z/p)* with all c*v mod p in [mu,p-mu] => M>=mu/p. returns c or None."""
    for c in range(1,p):
        if gcd(c,p)!=1: continue
        if all(mu <= (c*v)%p <= p-mu for v in S): return c
    return None

def pmclass(x,p): return frozenset({x%p,(-x)%p})

def transversal_mod25(S):
    hit=set(pmclass(v,25) for v in S if v%5!=0)
    return len(hit)==10

def longest_ap_subset(S):
    """size of the largest arithmetic-progression subset (any spacing)."""
    S=sorted(set(S)); best=1; n=len(S)
    idx={v:i for i,v in enumerate(S)}
    for i in range(n):
        for j in range(i+1,n):
            g=S[j]-S[i]; L=2; nxt=S[j]+g
            while nxt in idx: L+=1; nxt+=g
            best=max(best,L)
    return best

def analyze(S):
    S=[v//reduce(gcd,S) for v in S]           # primitive
    if len(set(S))<len(S): return None
    M=Mfast(S)
    return dict(S=S, M=M, ingap=(GAP_LO<M<GAP_HI), defects=len(S)-longest_ap_subset(S),
               transv=transversal_mod25(S), mult25=any(v%25==0 for v in S),
               c25=clears_mod(S,25,2))

def log(m=""): print(m,flush=True)

if __name__=="__main__":
    log("TWO-MODULUS analysis of >=3-defect near-tight 12-families (N=12).")
    log("Goal: are full-transversal-mod-25 families (mod-25 FAILS) ever in the gap? what clears them?\n")
    fams=set()
    # (1) dilated sub-APs (len 9, spacing g, dilation c) + 3 defects near the span
    for g in range(1,5):
        for c in range(1,6):
            base=[c*(1+g*j) for j in range(9)]         # 9-term sub-AP (>=3 defects)
            span=max(base)
            pool=[x for x in range(1,span+8) if x not in base]
            for defs in combinations(pool[:16],3):
                S=tuple(sorted(base+list(defs)))
                if len(S)==12: fams.add(S)
                if len(fams)>60000: break
            if len(fams)>60000: break
        if len(fams)>60000: break
    # (2) explicit full-transversal seeds: one rep per +/- class, scaled, + fill
    reps=[1,2,3,4,6,7,8,9,11,12]     # 10 reps, all 10 classes, no mult of 5
    for c in range(1,4):
        for extra in combinations([13,14,16,17,18,19,21,22,23,24,26,27],2):
            S=tuple(sorted([c*r for r in reps]+list(extra)))
            if len(set(S))==12: fams.add(S)
    log(f"generated {len(fams)} candidate families")
    ingap=[]; transv_all=[]; transv_nomult=[]; hard=[]
    for S in fams:
        r=analyze(list(S))
        if r is None: continue
        if r['ingap']: ingap.append(r)
        if r['transv']:
            transv_all.append(r)
            if not r['mult25']:
                transv_nomult.append(r)
                if r['M']<GAP_HI:      # candidate that mod-25 CANNOT clear and isn't obviously loose
                    hard.append(r)
    log(f"\n  IN-GAP families found: {len(ingap)}")
    for r in ingap[:8]:
        log(f"    {r['S']}  M={r['M']}  defects={r['defects']} transversal={r['transv']} mult25={r['mult25']}")
    log(f"\n  full-transversal mod 25 (mod-25 rotation FAILS, c25=None): {len(transv_all)}")
    log(f"    of these, NO multiple of 25 (kps's supposed-clearable set): {len(transv_nomult)}")
    log(f"    all have M >= gap-top 2/25?  {all(r['M']>=GAP_HI for r in transv_nomult)}")
    log(f"    min M among transversal-no-mult-25: {min((r['M'] for r in transv_nomult), default='n/a')}")
    log(f"\n  HARD residual (full transversal, no mult 25, M < 2/25) count: {len(hard)}")
    for r in sorted(hard,key=lambda r:-r['M'])[:12]:
        # what clears it above 1/13? scan small denominators for the actual M-witness
        log(f"    {r['S']}  M={r['M']} ({float(r['M']):.4f})  defects={r['defects']}  ingap={r['ingap']}")
    log("\n  => if HARD=0 (or all HARD have M>=2/25 via small-q), the transversal residual is loose the EASY way")
    log("     and (G) = [miss-a-class => mod-25 clears] + [full transversal => loose via small q / AP boundary].")
