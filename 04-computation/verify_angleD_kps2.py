#!/usr/bin/env python3
# Part 2: maxgap-sufficiency pointwise, AP tie, no-coupling, no-FOSD claims.
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.path.insert(0,'04-computation')
from verify_angleD_kps import N_at, maxgap, breakpoints, dist_p, consec

print("=== (4b) POINTWISE: maxgap<1/7 ==> N=0 at every breakpoint midpoint ===")
for k in [8,9,10]:
    E=consec(k); bps=breakpoints(E); one7=Fraction(1,7)
    bad=0; checked=0
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; checked+=1
        if maxgap(E,mid)<one7 and N_at(E,mid)!=0: bad+=1
    print(f"k={k}: intervals={checked}, violations(maxgap<1/7 but N>0)={bad}")

print("\n=== (3b) AP class TIES consec exactly (k=8, d=1 is consec; test general E -> only AP ties) ===")
# Confirm: a dilated AP has identical distribution (so ties p0). Already shown identical in part1.
# Now confirm UNIQUENESS more strongly: among ALL valid k=8 sets with spread<=11 none tie.
# (AP {0,d,..7d} needs 7d<=11 => d=1 only.) So uniqueness holds in that window. Confirmed in part1.
print("   AP {0,2,...,14} (k=8) p0 ==", dist_p([2*i for i in range(8)])[0], " consec p0 ==", dist_p(consec(8))[0])

print("\n=== (a) NO pointwise coupling: N_consec(x) <= N_E(x) FAILS ===")
for k in [8,9]:
    C=consec(k)
    if k==8: comps=[[0,1,2,3,4,5,6,8],[0,1,2,3,4,5,7,8],[0,2,3,4,5,6,7,9]]
    else: comps=[[0,1,2,3,4,5,6,7,9],[0,1,2,3,4,5,6,8,9]]
    for E in comps:
        # common refinement of breakpoints
        bps=sorted(set(breakpoints(C))|set(breakpoints(E)))
        viol=0; tot=0
        for i in range(len(bps)-1):
            lo,hi=bps[i],bps[i+1]
            if hi==lo: continue
            mid=(lo+hi)/2; tot+=1
            if N_at(C,mid) > N_at(E,mid): viol+=1
        print(f"k={k} E={E}: intervals where N_consec>N_E (violates coupling)={viol}/{tot}")

print("\n=== (b) NO first-order stochastic dominance: consec not FOSD-smallest ===")
# FOSD-smallest in N means P(N>=t) smallest for all t. Check cdf tail.
def tailcdf(p):  # P(N>=t) for t=0..7
    return [sum(p[t:]) for t in range(8)]
for k in [8,9]:
    C=consec(k); pc=dist_p(C); tc=tailcdf(pc)
    nondom=0; tested=0
    elems=list(range(1, (12 if k==8 else 11)))
    for combo in itertools.combinations(elems,k-1):
        E=[0]+list(combo)
        if reduce(gcd,E)!=1: continue
        tested+=1
        pe=dist_p(E); te=tailcdf(pe)
        # consec FOSD-smaller than E iff tc[t]<=te[t] for all t
        if not all(tc[t]<=te[t] for t in range(8)):
            nondom+=1
    print(f"k={k}: among {tested} competitors, consec NOT tail-dominated by it in {nondom} (i.e. consec fails to be FOSD-smallest vs {nondom})")
