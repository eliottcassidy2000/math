#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S32 (HYP-4622) -- ADVERSARIAL hunt for a full-transversal-mod-25 gap member,
and the clearing-witness denominator for non-AP full transversals.

The S32 dichotomy: v not full transversal mod 25 (no mult 25) => M>=2/25 (mod-25 rotation, a clean
witness = kps LRCMod25Floor).  RESIDUAL = full transversals.  Computation (structured) said the AP
is the UNIQUE full transversal with M<2/25 (boundary 1/13).  Here: (1) ADVERSARIALLY build full
transversals across heights, hunt for one with M in (1/13,2/25); (2) for non-AP transversals report
the M-witness q (is it bounded/small => "loose the easy way"?).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

GAP_LO, GAP_HI = F(1,13), F(2,25)
CLASSES=[[1,24],[2,23],[3,22],[4,21],[6,19],[7,18],[8,17],[9,16],[11,14],[12,13]]

def Mfull(S):
    """exact M and its witness denominator q (S20-fixed; witness q | pairwise sum/diff or 2v)."""
    S=sorted(set(S)); Q=set()
    for v in S: Q.add(2*v)
    for a,b in combinations(S,2): Q.add(a+b); Q.add(abs(a-b))
    Q.discard(0); best=F(0); bq=None
    for q in Q:
        for a in range(1,q):
            mn=min(min((v*a)%q,q-((v*a)%q)) for v in S)
            if F(mn,q)>best: best=F(mn,q); bq=q
    return best,bq

def transversal_mod25(S):
    hit=set(frozenset({v%25,(-v)%25}) for v in S if v%5!=0)
    return len(hit)==10

def prim(S):
    g=reduce(gcd,S); return tuple(sorted(x//g for x in S))

def rng(seed):
    x=seed
    while True:
        x=(1103515245*x+12345)&0x7fffffff; yield x

def log(m=""): print(m,flush=True)

if __name__=="__main__":
    log("ADVERSARIAL full-transversal-mod-25 search (no mult of 25), N=12, hunting for M in gap.\n")
    g=rng(20260706)
    seen=set(); ingap=[]; nonap=[]; total=0
    # build: for each of 10 classes pick a signed rep lifted by residue+25*k (k in 0..H), + 2 fillers
    for trial in range(120000):
        S=[]
        for cls in CLASSES:
            rep=cls[next(g)%2]                 # +a or -a residue
            k=next(g)%6                         # height layer 0..5
            val=rep+25*k
            if val==0: val=rep+25
            S.append(val)
        # 2 fillers: multiples of 5 (stay safe under mod-25) or extra units
        for _ in range(2):
            f=5*((next(g)%10)+1)                # 5..50, multiple of 5 (not nec 25)
            S.append(f)
        S=[v for v in S if v>0]
        S=prim(S)
        if len(set(S))!=12: continue
        if any(v%25==0 for v in S): continue
        if not transversal_mod25(S): continue
        if S in seen: continue
        seen.add(S); total+=1
        M,bq=Mfull(list(S))
        if GAP_LO<M<GAP_HI: ingap.append((S,M,bq))
        elif M>=GAP_HI: nonap.append((S,M,bq))
        if total>=4000: break
    log(f"distinct full-transversal families tested: {total}")
    log(f"  IN-GAP (1/13 < M < 2/25): {len(ingap)}   <-- if >0, kps's mod-25 reduction is INCOMPLETE")
    for S,M,bq in ingap[:10]: log(f"    {S} M={M} q={bq}")
    log(f"  M >= 2/25 (loose): {len(nonap)}")
    if nonap:
        qs=[bq for _,_,bq in nonap]
        import statistics
        log(f"    M-witness denominator q: min={min(qs)} max={max(qs)} median={int(statistics.median(qs))}")
        log(f"    fraction with q<=25 (small-q 'easy' witness): {sum(1 for q in qs if q<=25)/len(qs):.3f}")
        log(f"    fraction with q<=50: {sum(1 for q in qs if q<=50)/len(qs):.3f}")
        log("    sample (highest M):")
        for S,M,bq in sorted(nonap,key=lambda t:-t[1])[:6]:
            log(f"      {S} M={M} ({float(M):.4f}) witness-q={bq}")
    log("\n  => IN-GAP=0 supports: full transversal & M<2/25 => AP (boundary). Non-AP transversals loose;")
    log("     if witness-q is bounded, they are 'loose the easy way' at a small denominator (two-modulus complete).")
