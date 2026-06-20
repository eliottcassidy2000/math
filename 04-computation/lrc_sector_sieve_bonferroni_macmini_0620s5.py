#!/usr/bin/env python3
"""
lrc_sector_sieve_bonferroni_macmini_0620s5.py  (mac-mini-2026-06-20-S5)

Inclusion-exclusion over the 7 SECTORS, Bonferroni truncation.
p0(E)=P(N=0)=sum_r (-1)^r S_r,  S_r = sum_{|A|=r, A subset {1..6}} meas{all sectors in A missed by E}
     = E[C(N,r)] (factorial moments; N=#missed nonzero sectors).
Bonferroni: p0 <= 1-S1+S2 (truncate even), p0 >= 1-S1+S2-S3 (odd), etc.
TEST: does a low-order Bonferroni UPPER bound <= cap_k for the dangerous extremal clusters?
Also CUBE-ROOT organization: group the 6 nonzero sectors by C_3 orbits {1,2,4},{3,5,6}; compute
the orbit-resolved moments and see if the cube-root structure tightens the bound.
"""
import itertools, sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def sector_of(p): return int((p%1)*7)
def breakpoints(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    return sorted(bps)
def missed_measure(E, A):
    """meas{x: orbit of E misses ALL sectors in A} (A subset {1..6})."""
    Aset=set(A); bps=breakpoints(E); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        secs=set(sector_of(e*((x0+x1)/2)) for e in E)
        if not (Aset & secs): tot+=x1-x0   # none of A's sectors hit => all missed
    return tot
def S_r(E, r):
    return sum(missed_measure(E,A) for A in itertools.combinations(range(1,7), r))
def p0(E):
    bps=breakpoints(E); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        if len(set(sector_of(e*((x0+x1)/2)) for e in E))==7: tot+=x1-x0
    return tot

caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7)}
print("Sector-sieve Bonferroni partial sums vs cap (extremal/dangerous clusters):")
rows={
 'consec_8 (k=8)':( (0,1,2,3,4,5,6,7), 8),
 'consec_9 (k=9)':( (0,1,2,3,4,5,6,7,8), 9),
 'true-wide leader (k=9)':((0,4,6,8,10,12,14,15,16),9),
 'wide AP (0,2,..,16) k=9':((0,2,4,6,8,10,12,14,16),9),
}
for name,(E,k) in rows.items():
    Sr=[S_r(E,r) for r in range(0,7)]   # S0=1
    Sr[0]=F(1)
    partial=[]; acc=F(0)
    for r in range(0,7):
        acc+= (-1)**r * Sr[r]; partial.append(acc)
    actual=p0(E); cap=caps[k]
    print(f"\n{name}: p0={float(actual):.5f}  cap_{k}={float(cap):.5f}")
    print(f"   S_r = {[f'{float(x):.4f}' for x in Sr]}")
    for r in range(1,7):
        b=partial[r]; typ='UPPER' if r%2==0 else 'lower'
        mark = ('<=cap OK' if (r%2==0 and b<=cap) else ('' if r%2 else '')) 
        print(f"   Bonferroni level {r} ({typ}): {float(b):+.5f}  {('<= cap? '+('YES' if b<=cap else 'no')) if r%2==0 else ''}")
print("\nIf an even-level (UPPER) Bonferroni partial <= cap for ALL clusters, p0<=cap is proved at that level.")
