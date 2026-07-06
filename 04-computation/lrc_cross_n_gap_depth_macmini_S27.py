#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S27 -- CROSS-N understanding of the LRC gap, toward opus-S117's
(O-depth-monotone): map the achievable gap ORDER/depth across N=4..13 and confirm it
decreases -> 0 at N=12/13.

The second gap at N speeds is (1/(N+1), 2/(2N+1)); a member has value s/((N+1)s+... )
-- opus: below-1/(N+1) is a rung s/(ns+k), n=N+1; ORDER k>=2 for interior; k<s<2k.
Gap members are BORDERED / defected dilated APs (mac-mini construction): a spine
a+d*{0..m-1} plus 'border' elements (spine_i +- small).  This search builds those
(including INTERIOR borders opus's search missed) and records, per N, the achievable
gap members and their MIN ORDER k (= q - (N+1)*s).  Uses the S20-FIXED fast-M.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations, product
import time

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

def order_k(M, N):   # N speeds; M=s/q; k = q - (N+1)*s
    return M.numerator, M.denominator - (N+1)*M.numerator

def log(m=""): print(m, flush=True)

log("cross-N gap map (bordered/defected dilated APs), FIXED fast-M:\n")
log(f"  {'N':>2} {'gap window':>18} {'#gap-members':>13} {'min order k':>12} {'example (M, order)':>26}")
for N in range(4,14):
    lo,hi=F(1,N+1),F(2,2*N+1)
    found={}   # M -> example set
    t0=time.time(); tries=0
    # bordered dilated APs: spine a + d*{0..m-1}, plus borders = spine_i +- e (e in 1..3), total N speeds
    for d in range(1,9):
        for a in range(1,4):
            for m in range(max(3,N-4), N+1):     # spine length
                spine=[a+d*i for i in range(m)]
                nb=N-m                            # number of border/defect elements
                if nb<0: continue
                if nb==0:
                    cand=[tuple(spine)]
                else:
                    # choose nb border elements: each = some spine element +- e
                    borders=set()
                    for i in range(m):
                        for e in (1,2,3):
                            for s in (1,-1):
                                b=spine[i]+s*e
                                if b>0 and b not in spine: borders.add(b)
                    borders=sorted(borders)
                    cand=[]
                    for combo in combinations(borders, nb):
                        cand.append(tuple(sorted(spine+list(combo))))
                        if len(cand)>4000: break
                for S in cand:
                    tries+=1
                    if len(set(S))<N: continue
                    Sp=[x//reduce(gcd,S) for x in S]
                    if len(set(Sp))<N: continue
                    M=Mfast(Sp)
                    if lo<M<hi and M not in found:
                        found[M]=Sp
                    if time.time()-t0>20: break
                if time.time()-t0>20: break
            if time.time()-t0>20: break
        if time.time()-t0>20: break
    if found:
        orders={order_k(M,N)[1]:M for M in found}
        mink=min(orders)
        Mex=orders[mink]
        exset=found[Mex]
        log(f"  {N:>2} ({float(lo):.4f},{float(hi):.4f}) {len(found):>13} {mink:>12} "
            f"  {str(Mex)+' k='+str(mink):>24}  {exset}")
    else:
        log(f"  {N:>2} ({float(lo):.4f},{float(hi):.4f}) {0:>13} {'-':>12}   (none found)")
log("\n=> if min achievable order k INCREASES / #members ->0 as N->13, that IS opus's (O-depth-monotone):")
log("   the window 1/(2N^2) outruns the achievable denominator; the mediant 3/(3N+2) is the last to survive.")
