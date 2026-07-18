#!/usr/bin/env python3
"""
lrc_additive_cert_klein_S326.py -- klein-2026-07-18-S326
Owner: build the additive certificate route (the direction surviving S322-S325).

BUILT: the measure-recursion certificate. Track (mu, N) = (measure of the good set, component count)
instead of the largest interval:
      mu_{k+1} >= mu_k (1-2d) - 2 d N_k / w_k        (each component loses <= 2d*len + 2d/w)
      N_{k+1}  <= N_k + w_k                          (w disjoint arcs split at most w times)
This REMOVES the r < 1/(2d) = 7 cap of the largest-interval tail (THM-1004): each speed now removes a
FRACTION, so the decay is (1-2d)^r > 0 for every r.

RESULT: IT STILL FAILS -- 0/8 on named + random covering families -- but for a DIFFERENT and more
fundamental reason, which is the useful output.

EXACT DIAGNOSIS (deep well, base {1..7}, mu_0 = 0.334694, N_0 = 18):
   add w |  decay mu(1-2d) | boundary 2dN/w |    mu bound  |  EXACT |G|
      8  |     0.286880    |    0.321429    |  -0.034548   |  0.265816
      9  |    -0.029613    |    0.412698    |  -0.442311   |  0.181066
     10  |    -0.379124    |    0.500000    |  -0.879124   |  0.137982
    182  |    -1.554412    |    0.053375    |  -1.607788   |  0.023897
The boundary term kills it at STEP ONE: at w=8, 2dN/w = 0.321 exceeds the whole surviving measure 0.287.
Note the last row: at w=182 the boundary is only 0.053. SMALL speeds are the problem, not large ones.
And the EXACT |G| column stays healthy throughout -- the geometry is fine; the ACCOUNTING collapses.

CONSEQUENCE. Both formulations of the additive route -- largest-interval (THM-1004/1015) and
measure-recursion (here) -- require the ADDED SPEEDS TO BE LARGE. The measure version escapes the r-cap
but not the large-speed requirement. So THM-1015's large-killer regime is not a convenient special case;
it is the NATURAL LIMIT of the additive method. Small speeds have arcs of width 2d/w that are wide
relative to the component structure, and no amount of reformulation hides that.

THE MISSING INPUT, precisely. The accounting charges 2d/w per component per step, but a component of
length l >> 1/w loses only ~2d*l, proportionally. The crude bound is tight only for components of length
~1/w. A working certificate needs a bound on the MEASURE HELD IN SHORT COMPONENTS (length < c/w); with
that, mu' >= (mu - mu_short)(1-2d)(1-1/c). The component-length distribution of G_B(delta) is the
missing ingredient -- that is the concrete next object, not another reformulation of the recursion.
"""
import numpy as np
d=1.0/14
def good_meas_comps(B):
    cuts=[0.0,1.0]
    for v in B:
        iv=1.0/v; dv=d/v
        for k in range(0,v+1):
            for e in (k*iv-dv,k*iv+dv):
                if 0.0<=e<=1.0: cuts.append(e)
    cs=sorted(set(cuts)); segs=[]
    for i in range(len(cs)-1):
        a,b=cs[i],cs[i+1]
        if b<=a: continue
        mid=0.5*(a+b)
        if all(min((v*mid)%1.0, 1-((v*mid)%1.0))>=d-1e-12 for v in B): segs.append((a,b))
    m=[]
    for a,b in segs:
        if m and abs(m[-1][1]-a)<1e-15: m[-1]=(m[-1][0],b)
        else: m.append((a,b))
    return sum(b-a for a,b in m), max(1,len(m))
def cert(V,bs=7,order='asc'):
    V=sorted(V); B=V[:bs]; K=sorted(V[bs:],reverse=(order=='desc'))
    mu,N=good_meas_comps(B)
    for w in K:
        mu = mu*(1-2*d) - 2*d*N/w
        N  = N + w
        if mu<=0: return False,mu
    return mu>0,mu
if __name__=="__main__":
    print("deep well:", cert([1,2,3,4,5,6,7,8,9,10,11,12,182]))
    print("(fails at step 1: boundary 2dN/w = 0.321 > decay 0.287)")
