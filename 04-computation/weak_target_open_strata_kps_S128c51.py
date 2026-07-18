#!/usr/bin/env python3
"""weak_target_open_strata_kps_S128c51.py -- kind-pasteur S128 cont.51.
THE WEAK-TARGET RELAXATION. THM-724/726 prove covering => M >= 14/183 (sharp), with open
strata at |P| <= 8 (>=5 far outliers). But the equality horn needs only the WEAK target
M > 1/14. Since 14/183 = 0.07650 vs 1/14 = 0.07143 (7% gap), the weak target has more room.
TEST: sample covering families IN the open strata (>=5 far outliers) and measure how far
M sits above 1/14 vs above 14/183 vs above 1/13 (THM-726's multi-killer target)."""
import sys, random
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
random.seed(523726)
LAM=F(1,14); SHARP=F(14,183); MK=F(1,13)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def Mgrid(V,Q=20000):
    best=0.0; iq=1.0/Q
    for p in range(1,Q):
        t=p*iq; m=2.0
        for v in V:
            d=abs(v*t-round(v*t))
            if d<m:
                m=d
                if m<=best: break
        if m>best: best=m
    return best
def is_covering(V):
    return all(any(v%q==0 for v in V) for q in range(2,15))
def gen_open_stratum():
    """covering family with >=5 FAR outliers (|P| <= 8): small core + >=5 big speeds."""
    for _ in range(4000):
        core=sorted(random.sample(range(1,20), random.randint(6,8)))
        nfar=13-len(core)
        if nfar<5: continue
        far=sorted(random.sample(range(200,3000), nfar))
        V=sorted(set(core+far))
        if len(V)!=13: continue
        if is_covering(V): return V, len(core), nfar
    return None,None,None
res=[]
for s in range(400):
    V,nc,nf=gen_open_stratum()
    if V is None: continue
    m=Mgrid(V)
    res.append((m,tuple(V),nc,nf))
res.sort()
if not res:
    print("no covering open-stratum families generated"); sys.exit()
print("== WEAK-TARGET TEST in THM-726's open strata (|P|<=8, >=5 far outliers) ==")
print("  covering families sampled: %d"%len(res))
mn=res[0][0]
print("  min M = %.6f  (core size %d, far %d)"%(mn,res[0][2],res[0][3]))
print("  V_min = %s"%(list(res[0][1]),))
print("  thresholds: 1/14=%.6f (WEAK) | 14/183=%.6f (SHARP) | 1/13=%.6f (THM-726)"%(1/14,14/183,1/13))
below_weak=[m for m,_,_,_ in res if m<=1/14]
below_sharp=[m for m,_,_,_ in res if m<=14/183]
below_mk=[m for m,_,_,_ in res if m<=1/13]
print("  families at or below WEAK 1/14 : %d"%len(below_weak))
print("  families at or below SHARP 14/183: %d"%len(below_sharp))
print("  families at or below THM-726 1/13: %d"%len(below_mk))
print("  min margin over 1/14 = %+.6f  (%.2fx threshold)"%(mn-1/14, mn*14))
print()
if not below_weak:
    print("  >>> WEAK TARGET COMFORTABLE: every open-stratum covering family has M > 1/14.")
    if below_sharp or below_mk:
        print("  >>> AND some fall below the SHARP/THM-726 targets -- the weak target is STRICTLY EASIER here:")
        print("      the open strata can be closed for M > 1/14 even where M >= 14/183 or 1/13 is delicate.")
    else:
        print("  >>> (none below sharp either in this sample; the strata are loose at both targets)")
else:
    print("  *** some open-stratum covering family has M <= 1/14 -- investigate ***")
print("DONE")
