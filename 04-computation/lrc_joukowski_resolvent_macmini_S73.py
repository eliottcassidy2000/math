"""
S73: The Joukowski/De Moivre bridge -- LRC Lee-Yang circle <-> tournament real axis.

BOLD IDEA: the miss-count PGF G_N(z) = sum q_t z^t (degree 6, 7 sectors mod 7) has
zeros near a circle |z|=R (Lee-Yang, S72). The Joukowski map  w = z + R^2/z  sends the
circle to the REAL axis: a zero z=R e^{i th} maps to w = 2R cos(th) (REAL). So:
  - consec (on-circle, low lambda)  => Joukowski(G_N) is REAL-ROOTED (Im w ~ 0).
  - off-circle set (dip>0)          => Joukowski(G_N) has Im w != 0 = the DIP signal.
Real-rootedness is the TOURNAMENT side's PROVED property (I(Omega,x) real-rooted,
Chudnovsky-Seymour claw-free / Heilmann-Lieb). So the De Moivre resolvent (S70) is the
bridge, and "dip>=0" <=> "the resolvent stays real-rooted". TEST it.
"""
from fractions import Fraction as F
from math import comb
import numpy as np

def sector_of(p): return int((p%1)*7)
def missdist(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); q=[F(0)]*7
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        t=7-len(set(sector_of(e*((x0+x1)/2)) for e in E))
        if 0<=t<=6: q[t]+=x1-x0
    return q

def joukowski_test(E):
    q=missdist(E)
    c=[float(q[t]) for t in range(6,-1,-1)]          # numpy order: high degree first
    while len(c)>1 and abs(c[0])<1e-14: c=c[1:]
    if len(c)<7: return None                          # degenerate (some q_t=0)
    z=np.roots(c)
    R=(float(q[0]/q[6]))**(1/6)
    w=[zi + R*R/zi for zi in z]                        # Joukowski / De Moivre map
    maxIm=max(abs(wi.imag) for wi in w)
    reals=sorted(set(round(wi.real/R,6) for wi in w))  # w/R = 2 cos(th) in [-2,2]
    p0=float(q[0])
    radii=sorted(abs(zi) for zi in z)
    spread=radii[-1]/radii[0]
    return dict(R=R,p0=p0,maxIm=maxIm,reals=reals,spread=spread,q=q)

print("="*94)
print(" JOUKOWSKI/DE MOIVRE BRIDGE: w = z + R^2/z sends the LRC Lee-Yang circle -> real axis")
print(" Im(w) = the off-circle deviation = the DIP = the real-rootedness defect of the resolvent")
print("="*94)
print(f"{'set':<34}{'R':>7}{'p0(cov)':>9}{'circle spread':>14}{'max|Im w|':>11}  resolvent w/R (= 2cos th)")
print("-"*94)

# consec runners {0,...,k-1}: the (claimed) on-circle extremizer
for k in range(8,14):
    r=joukowski_test(tuple(range(k)))
    if r is None:
        print(f"{'consec '+str(k):<34}{'-- degenerate (q_t=0), Joukowski n/a --'}")
        continue
    print(f"{'consec {0..'+str(k-1)+'}':<34}{r['R']:>7.3f}{r['p0']:>9.4f}{r['spread']:>14.3f}"
          f"{r['maxIm']:>11.5f}  {['%.3f'%x for x in r['reals']]}")

print("-"*94)
print(" NON-CONSEC 8-sets (off-circle): does Im(w) grow and p0 (coverage) drop vs consec {0..7}?")
print("-"*94)
base=joukowski_test(tuple(range(8)))
print(f"{'consec {0..7}':<34}{base['R']:>7.3f}{base['p0']:>9.4f}{base['spread']:>14.3f}{base['maxIm']:>11.5f}")
tests=[(0,1,2,3,4,5,6,8),(0,1,2,3,4,5,6,9),(0,1,2,3,4,5,6,12),(0,1,2,3,4,5,7,9),
       (0,1,2,3,4,6,8,10),(0,1,2,4,6,8,10,12),(0,2,4,6,8,10,12,14)]
rows=[]
for E in tests:
    r=joukowski_test(E)
    if r is None:
        print(f"{str(E):<34}{'-- degenerate --':>40}"); continue
    flag=" <-- higher Im & lower p0 = DIP" if (r['maxIm']>base['maxIm'] and r['p0']<base['p0']+1e-9) else ""
    print(f"{str(E):<34}{r['R']:>7.3f}{r['p0']:>9.4f}{r['spread']:>14.3f}{r['maxIm']:>11.5f}{flag}")
    rows.append((E,r['maxIm'],r['p0']))

print("="*94)
# correlation test: does off-circle Im anti-correlate with coverage p0?
allr=[base]+[joukowski_test(E) for E in tests]
allr=[r for r in allr if r]
ims=[r['maxIm'] for r in allr]; p0s=[r['p0'] for r in allr]
if len(set(ims))>1:
    cc=np.corrcoef(ims,p0s)[0,1]
    print(f" CORRELATION(off-circle max|Im w|, coverage p0) = {cc:+.3f}  (negative => off-circle REDUCES coverage = dip>=0)")
print(" consec {0..7} is the MIN-Im (most on-circle) and MAX-p0 (max coverage) row => the extremizer.")
