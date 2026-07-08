"""
mac-mini-2026-07-07-S57 (HYP-5237) -- INDEPENDENT CONFIRMATION of the k=10 conditional
composition (de-risking kps-S75 HYP-5247).

The conditional Markov on G_P, with the window mass added (windows inside G_P where q|d
pairs pay the tent exactly), is EXACTLY additive:

  meas(S∩G_P)*toll <= int_{G_P} F  and  int_{G_P} F >= toll*meas(S∩G_P) + W_F^{G_P}
  (windows in G_P) AND int_{G_P} F = avgc*k(k-1)*meas*int f = avgc*(1-floor)*toll*meas.
  => meas(S∩G_P) <= avgc*(1-floor)*meas - W_F^{G_P}/toll
  => rho* >= meas*(1 - avgc*(1-floor)) + W_F^{G_P}/toll.

So the composed bound = THM-655 conditional term + W_F^{G_P}/toll.  W_F^{G_P} = the window
mass from the q<=6 windows whose CENTER p'/q lies in G_P (the window half-width c_q/diam is
tiny, so center-in-G_P is the leading criterion; being conservative we require the whole
window in G_P: keep a window only if [p'/q - c_q/diam, p'/q + c_q/diam] ⊂ G_P).

TEST: for the k=10 residual families (block, block+outlier, 2-AP) at every shape, does the
COMPOSED bound rho*_comp >= m_P?  Report min over (shape, residual family).
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import numpy as np

MP = F(14249, 252252); TH = F(1, 7); DMAX = 30000; k = 10
BETA = F(14-k, 7*k); W = TH - BETA; INTF = W*W/2
TOLL = 1 - k*BETA
FLOOR = 1 - F(2*(k-1)*(k-7), 7*k)
def phi(q): return sum(1 for a in range(1,q+1) if gcd(a,q)==1)
CQ = {q: F(7-q, 7*q) for q in range(1,7)}
mp = float(MP)

def GP(P):
    bad=[]
    for p in P:
        w=F(1,14*p)
        for j in range(p+1): bad.append((F(j,p)-w,F(j,p)+w))
    bad=[(max(l,F(0)),min(h,F(1))) for l,h in bad if h>0 and l<1]; bad.sort(); m=[]
    for l,h in bad:
        if m and l<=m[-1][1]: m[-1]=(m[-1][0],max(m[-1][1],h))
        else: m.append((l,h))
    g=[];pv=F(0)
    for l,h in m:
        if l>pv: g.append((pv,l))
        pv=max(pv,h)
    if pv<1: g.append((pv,F(1)))
    return g

def in_GP(x, iv):
    return any(l <= x <= h for l,h in iv)
def window_in_GP(center, half, iv):
    lo, hi = center-half, center+half
    return any(l <= lo and hi <= h for l,h in iv)

def carr(P):
    bf,tf,iff_=float(BETA),float(TH),float(INTF)
    iv=GP(P);meas=sum(h-l for l,h in iv)
    def D(t): return np.where(t<=bf,-t*iff_,np.where(t<=tf,(t-bf)**2/2-t*iff_,(1.0-t)*iff_))
    da=np.arange(1,DMAX+1,dtype=np.int64);acc=np.zeros(DMAX)
    for (l,h) in iv:
        for (pt,sg) in ((h,1.0),(l,-1.0)):
            nu,de=pt.numerator,pt.denominator;acc+=sg*D(((da*nu)%de).astype(float)/de)
    c=1.0+acc/(da*float(meas)*iff_);ct=1.0+len(iv)*(6/7+bf+float(W)**2/4)/(DMAX*float(meas))
    return c,float(meas),ct,iv

def cst(meas): return float((1-MP/F(meas).limit_denominator(10**9))*F(7*k,2*(k-1)*(k-7)))
def avgc(E,c,ct):
    E=sorted(E);t=0.0;n=0
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i];t+=c[d-1] if d<=DMAX else ct;n+=1
    return t/n

def tentint(u):
    x=min(u,float(TH))-float(BETA); return (x*x/2) if x>0 else 0.0
def Int(L):
    fl=int(L); return fl*float(INTF)+tentint(L-fl)

def WF_GP(E, iv):
    """window mass from q<=6 windows whose full extent lies inside G_P."""
    E=sorted(E)
    g=0
    for i in range(len(E)):
        for j in range(i+1,len(E)): g=gcd(g,E[j]-E[i])
    g=max(g,1)
    diam=(E[-1]-E[0])//g
    diffs=[(E[j]-E[i])//g for i in range(len(E)) for j in range(i+1,len(E))]
    tot=0.0
    for q in range(1,7):
        cq=CQ[q]; half=float(cq)/diam
        for pp in range(0,q):
            if q>1 and gcd(pp,q)!=1: continue
            if q==1 and pp!=0: continue
            center=F(pp,q)
            if not window_in_GP(center, half, iv): continue
            for d in diffs:
                if d % q == 0:
                    tot += 2*(1.0/d)*Int(d*float(cq)/diam)
    return tot

def rho_comp(E, c, ct, meas, iv):
    base = meas*(1 - avgc(E,c,ct)*float(1-FLOOR))
    return base + WF_GP(E, iv)/float(TOLL)

SHAPES=list(combinations(range(1,14),3))
def resid_fams():
    fams=[]
    for cc in range(1,7): fams.append([cc*t for t in range(k)])          # dilated blocks
    for X in [9,10,11,12,15,20,30,50,100]: fams.append(list(range(9))+[X])  # block+1 outlier
    for X in [10,12,20,40]: fams.append(list(range(8))+[X,X+1])          # block+2 outliers
    fams.append([2*t for t in range(k)])                                 # 2-AP
    fams.append(list(range(5))+[X+5 for X in range(5)])                  # two-block
    for holes in [(3,),(4,),(5,),(3,6)]:                                 # perforated blocks diam ~11-13
        base=[t for t in range(k+len(holes)+1) if t not in holes][:k]
        if len(set(base))==k: fams.append(base)
    return fams
from math import gcd as _gcd
from functools import reduce
def primitive(E):
    E=sorted(E); dif=[E[i+1]-E[i] for i in range(len(E)-1)]
    return reduce(_gcd, dif) == 1
def diam(E): return max(E)-min(E)
ALL=resid_fams()
# klein-S174 THM-653 exhausts PRIMITIVE diam<=10 for k=10; our job = primitive diam>=11
FAMS=[E for E in ALL if primitive(E) and diam(E) >= 11]
print(f"=== k=10 CONDITIONAL COMPOSITION on PRIMITIVE diam>=11 (klein covers diam<=10) vs m_P={mp:.5f} ===")
print(f"    {len(FAMS)} primitive-diam>=11 residual families x {len(SHAPES)} shapes")
print(f"    (excluded: {sum(1 for E in ALL if not primitive(E))} non-primitive = dilations, reduce to primitive core;")
print(f"     {sum(1 for E in ALL if primitive(E) and diam(E)<=10)} primitive diam<=10 = klein-covered)\n")

worst=(1e9,None,None); per_shape=[]
for P in SHAPES:
    c,meas,ct,iv=carr(P)
    sw=(1e9,None)
    for E in FAMS:
        r=rho_comp(E,c,ct,meas,iv)
        if r<sw[0]: sw=(r,tuple(E))
    per_shape.append((sw[0]-mp,P,sw[1]))
    if sw[0]<worst[0]: worst=(sw[0],P,sw[1])
per_shape.sort()
nclosed=sum(1 for d,_,_ in per_shape if d>=0)
print(f"shapes CLOSED by the composed bound (all residual families >= m_P): {nclosed}/{len(SHAPES)}")
print(f"min over all (shape,family) of rho_comp - m_P = {worst[0]-mp:+.5f} at P={worst[1]}, E={worst[2]}")
print(f"\nworst 12 (deficit, shape, family):")
for d,P,E in per_shape[:12]:
    Es=f"{E[:4]}..{E[-1]}" if E and len(E)>6 else E
    print(f"  {d:+.5f}  P={P}  E={Es}")
