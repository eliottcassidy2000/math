# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont28: FLESH OUT THE ARCHITECTURE. An exact identity collapses the depth-5 Bonferroni
# B5 into a transparent "live minus deep-coverage penalty", reducing hB5 to a clean CHECKABLE condition.
#
# KEY IDENTITY (partial alternating binomial sum): T5(n) := sum_{d=0}^5 (-1)^d C(n,d) = -C(n-1,5).
#   T5(0)=1, T5(n)=0 for 1<=n<=5, T5(n)=-C(n-1,5)<0 for n>=6.
# THEREFORE:  B5(v,q) = sum_p T5(bandCount(p)) = liveCount(q) - PENALTY(q),
#   PENALTY(q) := sum_{p: bandCount(p)>=6} C(bandCount(p)-1, 5).
# So B5>0  <=>  liveCount(q) > PENALTY(q).  And a CLEAN RULER (max bandCount<=5) has PENALTY=0, giving
# B5 = liveCount exactly. hB5 then reduces to: every residual family has a ruler with a live multiplier
# and NO multiplier covering >=6 runners.  THIS is the fleshed-out architecture's base lemma.
from math import comb

def bandCount(v,q,p): return sum(1 for vi in v if not (q <= 14*((vi*p)%q) <= 13*q))
def T5(n): return sum((-1)**d * comb(n,d) for d in range(6))
def B5_direct(v,q):
    return sum(T5(bandCount(v,q,p)) for p in range(1,q))
def liveCount(v,q): return sum(1 for p in range(1,q) if bandCount(v,q,p)==0)
def penalty(v,q): return sum(comb(bandCount(v,q,p)-1,5) for p in range(1,q) if bandCount(v,q,p)>=6)
def B5_reform(v,q): return liveCount(v,q) - penalty(v,q)
def maxclean(v,qmax=320):
    # best q by B5, and whether a CLEAN ruler (penalty 0, live>=1) exists
    bestB5=(-10**9,0,None); cleanq=None
    for q in range(14,qmax+1):
        lc=liveCount(v,q); pen=penalty(v,q); b=lc-pen
        if b>bestB5[0]: bestB5=(b,q,(lc,pen))
        if pen==0 and lc>=1 and cleanq is None: cleanq=(q,lc,max((bandCount(v,q,p) for p in range(1,q)),default=0))
    return bestB5, cleanq

def main():
    print("(1) IDENTITY T5(n) = -C(n-1,5):  n:  T5(n)  -C(n-1,5)")
    for n in range(0,10):
        c = -comb(n-1,5) if n>=1 else 1   # T5(0)=1
        print(f"        n={n}:  {T5(n):4d}   {c:4d}   {'OK' if T5(n)==c else 'MISMATCH'}")
    print("   => T5 = +1 at n=0, ZERO for 1<=n<=5, negative (=-C(n-1,5)) for n>=6.\n")

    fams = {
        'AP {1..13}':            list(range(1,14)),
        'near-AP {1..12,26} (BINDING)': list(range(1,13))+[26],
        'covering {1..11,13,26}':list(range(1,12))+[13,26],
        'covering {1..11,13,36}':list(range(1,12))+[13,36],
        'resid A [.. ,23,25]':   [1,2,3,4,5,6,7,8,9,10,11,23,25],
        'far {1..12,120}':       list(range(1,13))+[120],
    }
    print("(2) B5 = liveCount - PENALTY  (verify reformulation == direct), and the CLEAN RULER:")
    print("     family                          bestB5 @q  (live,pen)   clean-ruler q (live, maxBand)")
    for name,v in fams.items():
        (b,q,lp),cleanq = maxclean(v)
        # verify reform == direct at the best q
        assert B5_reform(v,q)==B5_direct(v,q), (name,q)
        cq = f"q={cleanq[0]} (live={cleanq[1]}, maxBand={cleanq[2]}<=5)" if cleanq else "NONE FOUND (penalty unavoidable)"
        print(f"     {name:31s} {b:4d} @{q:4d}  {str(lp):9s}   {cq}")
    print()
    print("(3) THE BASE LEMMA (clean ruler): if max_p bandCount(v,q,p) <= 5 and liveCount(q) >= 1,")
    print("    then PENALTY=0 so B5(v,q) = liveCount(q) > 0.  Reduces the depth-5 obligation to:")
    print("    'a live multiplier with no multiplier covering >= 6 runners' -- transparent + Lean-friendly.")
    # verify the binding near-AP is discharged by a clean ruler
    v = list(range(1,13))+[26]
    (b,q,lp),cleanq = maxclean(v)
    print(f"    BINDING near-AP {{1..12,26}}: clean ruler {cleanq[0]}, maxBand {cleanq[2]}<=5, B5=liveCount={cleanq[1]}>0. DISCHARGED.")

if __name__ == '__main__':
    main()
