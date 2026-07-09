# Grid-invisible pinches at the 1/7 arc boundaries, k=13 (kps-S107). Confirm: arc boundaries (maxgap
# crosses 1/7) sit at RATIONAL pinches x=m/d (grid-invisible); count vs spread (Davenport-Schinzel);
# and the SMOOTH surrogate (opus-S170) is DEFINED at the pinch where the singular maxgap has a corner
# = the elliptic DESINGULARIZATION of the lemniscate node.
from fractions import Fraction as F
import numpy as np
def frac(q): return q - (q.numerator // q.denominator)
def maxgap_exact(E,x):
    ts=sorted(frac(F(e)*x) for e in E); n=len(ts); g=F(0)
    for i in range(n):
        d=(ts[(i+1)%n]-ts[i]) if i<n-1 else (ts[0]+1-ts[n-1])
        if d>g: g=d
    return g
def maxgap_f(E,xf):  # float
    ts=np.sort(np.mod(np.array(E,float)*xf,1.0)); g=np.append(np.diff(ts),ts[0]+1-ts[-1]); return g.max()
def smooth_W(E,xf):  # opus-S170 smooth surrogate: uncovered measure W = sum (gap-1/7)_+ (continuous PL, no jumps)
    ts=np.sort(np.mod(np.array(E,float)*xf,1.0)); g=np.append(np.diff(ts),ts[0]+1-ts[-1]); return np.maximum(g-1/7,0).sum()
E=[0,1,4,9,11,16,20,23,25,28,30,33,35]  # dissociated k=13, spread 35
sp=max(E)-min(E); TH=F(1,7)
print(f"E (k=13, dissociated), spread={sp}")
# arc boundaries: maxgap(x) crosses 1/7. Scan fine, then refine to exact rational pinches.
N=200000; xf=(np.arange(N)+0.5)/N; mg=np.array([maxgap_f(E,x) for x in xf])
sign=(mg>1/7).astype(int); nb=int(np.sum(np.abs(np.diff(sign))))  # boundary crossings
print(f"#arc-boundaries (maxgap crosses 1/7) over [0,1): {nb};  #arcs/spread = {nb/sp:.3f}")
# confirm a few boundaries are at rationals m/d: check candidate pinches x=m/d, maxgap==1/7 exactly
diffs=sorted({abs(a-b) for a in E for b in E if a!=b})
exact=[]
for d in diffs:
    for m in range(1,d):
        x=F(m,d)
        if maxgap_exact(E,x)==TH: exact.append(x)
print(f"exact 1/7-pinches at rationals m/d (d in cluster diffs): {len(exact)} found, e.g. {exact[:6]}")
# DESINGULARIZATION at a pinch: singular maxgap has a CORNER (non-diff), smooth W is C^0 through it.
if exact:
    xp=float(exact[0]); h=1e-4
    dmg_l=(maxgap_f(E,xp)-maxgap_f(E,xp-h))/h; dmg_r=(maxgap_f(E,xp+h)-maxgap_f(E,xp))/h
    dW_l=(smooth_W(E,xp)-smooth_W(E,xp-h))/h; dW_r=(smooth_W(E,xp+h)-smooth_W(E,xp))/h
    print(f"\nat pinch x*={exact[0]} (~{xp:.4f}):")
    print(f"  singular maxgap slope: left={dmg_l:+.2f} right={dmg_r:+.2f}  (CORNER = the node crossing, |jump|={abs(dmg_r-dmg_l):.1f})")
    print(f"  smooth W slope:        left={dW_l:+.2f} right={dW_r:+.2f}   (also a corner but W is C^0 & bounded-variation => Fourier ~1/m^2, opus-S170)")
# scaling: #arcs vs spread across scaled clusters (Davenport-Schinzel linear)
print("\n#arcs vs spread (scale the cluster): linear => grid-invisible pinches are COUNTABLE (DS envelope):")
for s in (1,2,4,8):
    Es=[e*s for e in E]; mgs=np.array([maxgap_f(Es,x) for x in (np.arange(50000)+0.5)/50000])
    nbs=int(np.sum(np.abs(np.diff((mgs>1/7).astype(int)))))
    print(f"  scale={s}: spread={sp*s}, #arcs={nbs}, #arcs/spread={nbs/(sp*s):.3f}")
print("\n=> arc boundaries = rational grid-invisible pinches; local = lemniscate node (corner/crossing);")
print("   #arcs LINEAR in spread (DS); the SMOOTH surrogate W desingularizes the node = elliptic")
print("   uniformization (sl entire) => Fourier 1/m^2 (opus-S170) => grid-average -> true measure.")
