#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) C_3 reframing — DIAGNOSE the obstruction.  Why does C_3 bundling NOT crack
the signed d>=2 bound, even though orbit-level |T3|/sum|D7| ~ 0.64?

The trace collapse sum_a D7(ac)=Tr(D7) is exact ONLY if the lattice weights W_{ac}
are EQUAL across the dilation orbit.  The correction is sum_c W_c D7(c).  Bundling by
C_3 orbit gives sum_O [ sum_{c in O} W_c D7(c) ].  This equals sum_O W_bar(O) Tr3(D7)
ONLY in the equal-weight limit.  The residual is sum_O sum_{c in O}(W_c - W_bar) D7(c).

So the obstruction is the WITHIN-ORBIT WEIGHT SPREAD.  We measure:
  (a) within-C_3-orbit relative weight spread (vs within-full-orbit spread).
  (b) decompose correction = TRACE part (W_bar * Tr) + RESIDUAL part, for C_3 and full,
      and report |trace part| and |residual| separately.  If residual dominates and is
      itself only controllable by absolute bound, C_3 cannot help.
  (c) the TRUE cancellation is the lattice/Abel one (conditional convergence in n_j),
      orthogonal to the residue symmetry.  Show |corr| keeps shrinking with L while the
      trace part does NOT -> confirms residue symmetry is the wrong axis.
"""
import sys, itertools, cmath, math
from collections import defaultdict
from fractions import Fraction
sys.path.insert(0,"04-computation")
from lrc14_c3_partial_trace_macmini_0620s2 import (
    D7, T3, Tfull, znum, zadd, zscale, ZERO, C3, FULL,
    relations_support6, w_real)

def F(x): return Fraction(x).limit_denominator(10**12)

if __name__=="__main__":
    for E,lab in [(list(range(8)),"AP8 consec"),([0,1,3,5,7,9,11,13],"WIDE8 odd-AP")]:
        print(f"\n=== {lab}  E={E} ===")
        prev_exact=None
        for L in [3,4,5]:
            rels=relations_support6(E,L)
            byres=defaultdict(float)
            corr=ZERO
            for n in rels:
                c=tuple(v%7 for v in n); wv=w_real(n)
                byres[c]+=wv
                corr=zadd(corr, zscale(D7(c),F(wv)))
            exact=abs(znum(corr))

            # C_3 orbit decomposition into trace + residual
            for grp,name in [(C3,"C_3 "),(FULL,"full")]:
                orbits=defaultdict(list)
                for c,W in byres.items():
                    key=min(tuple((a*cj)%7 for cj in c) for a in grp)
                    orbits[key].append((c,W))
                trace_part=ZERO; resid_part=ZERO; spreads=[]
                for key,members in orbits.items():
                    Ws=[W for _,W in members]
                    Wbar=sum(Ws)/len(members)
                    # trace approx: Wbar * sum_{c in orbit-of-rep} D7(c) but only over members present.
                    # Use exact: trace contribution = Wbar * sum_{c in members} D7(c)
                    tr=ZERO
                    for c,_ in members: tr=zadd(tr, zscale(D7(c),F(Wbar)))
                    trace_part=zadd(trace_part,tr)
                    for c,W in members:
                        resid_part=zadd(resid_part, zscale(D7(c),F(W-Wbar)))
                    if len(members)>1:
                        mx=max(abs(x) for x in Ws)
                        if mx>0: spreads.append((max(Ws)-min(Ws))/mx)
                tp=abs(znum(trace_part)); rp=abs(znum(resid_part))
                avgspr=sum(spreads)/len(spreads) if spreads else float('nan')
                print(f"  L={L} {name}: exact={exact:.3e}  |trace_part|={tp:.3e}  "
                      f"|resid_part|={rp:.3e}  within-orbit wt-spread={avgspr:.3f}  "
                      f"#orbits={len(orbits)}")
            if prev_exact is not None:
                print(f"        exact|corr| L:{L-1}->{L}  ratio={exact/prev_exact:.3f} "
                      f"(<1 => lattice/Abel cancellation still active)")
            prev_exact=exact
