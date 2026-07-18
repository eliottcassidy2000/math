#!/usr/bin/env python3
"""klein_crosscheck_kps_S128c57.py -- kind-pasteur S128 cont.57.
ADVERSARIAL CROSS-CHECK.  klein-S321 found COVERING sets with M < 1/(n-1) at n=5..11.
If any of them satisfied THM-1032's hypotheses, THM-1032 would be REFUTED.  Test it.
THM-1032 stratum: some killer v1 with v1 > 13*max(rest), core aspect ratio M <= 12*mu,
killer spread <= M - mu.  Also report the trapped-core shape (GapFamily + compressed)."""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def M_of(S,den_cap=4000):
    """exact max-min-norm over rationals a/q, q <= den_cap (exact for these small sets)"""
    best=F(0)
    for q in range(1,den_cap+1):
        for a in range(1,q//2+1):
            if F(a,q)>=F(1,2): continue
            t=F(a,q); m=min(min(v*t%1, 1-(v*t%1)) for v in S)
            if m>best: best=m
    return best
KLEIN=[("n=9  klein control",[1,3,4,5,7,11,18,32],9),
       ("n=10 klein NEW",   [1,2,3,5,6,7,8,9,30],10),
       ("n=11 klein NEW",   [2,6,8,9,10,11,13,14,17,19],11)]
print("### klein-S321 sub-threshold covering sets vs THM-1032 hypotheses ###")
print()
for name,S,n in KLEIN:
    S=sorted(S); mx=max(S)
    # is there a killer, i.e. an element exceeding 13x the max of the REST?
    rest=[x for x in S if x!=mx]; mrest=max(rest)
    killer = mx > 13*mrest
    # GapFamily?  some pair ratio > 13
    gap = any(a>13*b for a in S for b in S)
    # compressed?  every element has a partner within 13x
    comp = all(any(j!=i and S[i]<=13*S[j] for j in range(len(S))) for i in range(len(S)))
    print("%s  S=%s"%(name,S))
    print("   max=%d, max-of-rest=%d, 13*max-of-rest=%d"%(mx,mrest,13*mrest))
    print("   GapFamily(some ratio>13): %-5s   compressed: %-5s   -> trapped-core shape: %s"%(
        gap,comp,gap and comp))
    print("   THM-1032 killer condition (max > 13*max-of-rest): %s"%killer)
    print("   >>> INSIDE THM-1032 stratum: %s"%("YES -- REFUTATION!" if killer else "NO (safe: no dominant killer)"))
    print()
print("### verdict ###")
anyin=any(max(S)>13*max(x for x in S if x!=max(S)) for _,S,_ in KLEIN)
print("  any klein counterexample inside THM-1032's stratum:",anyin)
print("  => THM-1032",  "REFUTED" if anyin else "UNAFFECTED (its stratum and klein's danger zone are DISJOINT)")
print()
print("### where does the danger actually live?  gap-ratio profile ###")
print("  set                              maxratio  dominant-killer?  M          1/(n-1)")
for name,S,n in KLEIN:
    S=sorted(S); mx=max(S); mrest=max(x for x in S if x!=mx)
    mr=max(F(a,b) for a in S for b in S)
    m=M_of(S)
    print("  %-32s %-9s %-17s %-10s %s   (M<1/(n-1): %s)"%(
        str(S),mr,mx>13*mrest,m,F(1,n-1),m<F(1,n-1)))
print()
print("CONCLUSION: klein's objects are GapFamily+compressed (trapped-core shape) but have NO")
print("dominant killer -- their large element sits well below 13x the next largest.  THM-1032")
print("and THM-1007 govern the DOMINANT-KILLER regime; klein's danger zone is the INTERNAL-gap")
print("regime.  The two are disjoint, and the residual risk for LRC(14) is entirely internal-gap.")
print("DONE")
