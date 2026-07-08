"""
mac-mini-2026-07-07-S57 (HYP-5237) -- k=10: does max(conditional avgc, global window)
close the residual?

Two independent lower bounds on rho*(P,E) (both proved):
 (a) CONDITIONAL average-form tent (THM-655 machinery, mac-mini):
       rho* >= meas(G_P) * (1 - avgc(E,P) * (1 - floor_k)),   floor_k = 1 - 2(k-1)(k-7)/(7k).
 (b) GLOBAL window composition (THM-653, klein) via the union bound rho* >= meas + mu - 1:
       mu(E) >= 1 - (E[F] - W_F(E)) / toll,  toll = 1 - k beta,
       W_F(E) = sum_{q<=6} phi(q) sum_{unordered pairs, q|d} 2*(1/d)*Int(d c_q/diam),
       c_q = (7-q)/(7q),  Int(L) = floor(L)(1/7-beta)^2/2 + tentint(L-floor(L)),
       tentint(u) = ((min(u,1/7)-beta)_+)^2/2.
     => rho*_b = meas(G_P) + mu - 1.

For each (shape, family), rho* >= max(a, b).  k=10 CLOSES if max(a,b) >= m_P for every
residual family.  We test the residual class (block, block+outlier at several X, 2-AP,
two-block, dilated blocks) at every shape; report the min over residual families of
[max(a,b) - m_P] -- if >= 0, k=10 is closed by the existing two bricks; else the exact
deficit is the spread-floor (klein-S175) target.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import numpy as np

MP = F(14249, 252252); TH = F(1, 7); DMAX = 30000; k = 10
BETA = F(14-k, 7*k); W = TH - BETA; INTF = W*W/2
TOLL = 1 - k*BETA                       # (k-7)/7 = 3/7
EF = F(k*(k-1),1) * INTF                # E[F] = k(k-1) int f
FLOOR = 1 - F(2*(k-1)*(k-7), 7*k)
def phi(q): return sum(1 for a in range(1,q+1) if gcd(a,q)==1)
CQ = {q: F(7-q, 7*q) for q in range(1,7)}

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

def carr(P):
    bf,tf,iff_=float(BETA),float(TH),float(INTF)
    iv=GP(P);meas=sum(h-l for l,h in iv)
    def D(t): return np.where(t<=bf,-t*iff_,np.where(t<=tf,(t-bf)**2/2-t*iff_,(1.0-t)*iff_))
    da=np.arange(1,DMAX+1,dtype=np.int64);acc=np.zeros(DMAX)
    for (l,h) in iv:
        for (pt,sg) in ((h,1.0),(l,-1.0)):
            nu,de=pt.numerator,pt.denominator;acc+=sg*D(((da*nu)%de).astype(float)/de)
    c=1.0+acc/(da*float(meas)*iff_);ct=1.0+len(iv)*(6/7+bf+float(W)**2/4)/(DMAX*float(meas))
    return c,float(meas),ct

def cst(meas): return float((1-MP/F(meas).limit_denominator(10**9))*F(7*k,2*(k-1)*(k-7)))
def avgc(E,c,ct):
    E=sorted(E);t=0.0;n=0
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i];t+=c[d-1] if d<=DMAX else ct;n+=1
    return t/n

def tentint(u):
    x = min(u, float(TH)) - float(BETA)
    return (x*x/2) if x > 0 else 0.0
def Int(L):
    fl = int(L)
    return fl*float(INTF) + tentint(L-fl)
def WF(E):
    E=sorted(E); diam=E[-1]-E[0]
    g=0
    for i in range(len(E)):
        for j in range(i+1,len(E)): g=gcd(g,E[j]-E[i])
    diam=diam//max(g,1)                 # primitive diam
    diffs=[(E[j]-E[i])//max(g,1) for i in range(len(E)) for j in range(i+1,len(E))]
    tot=0.0
    for q in range(1,7):
        cq=float(CQ[q]); ph=phi(q)
        for d in diffs:
            if d % q == 0:
                tot += ph * 2*(1.0/d)*Int(d*cq/diam)
    return tot

def rho_a(E,c,ct,meas):     # conditional average-form
    return meas*(1 - avgc(E,c,ct)*float(1-FLOOR))
def rho_b(E,meas):          # global window composition via union
    mu = 1 - (float(EF) - WF(E))/float(TOLL)
    return meas + mu - 1

SHAPES=list(combinations(range(1,14),3))
mp=float(MP)
print(f"=== k=10: max(conditional-a, window-b) vs m_P={mp:.5f}; E[F]={float(EF):.4f}, toll={float(TOLL):.4f} ===\n")

def residual_families(diam_cap=200):
    fams=[]
    # block, dilated blocks
    for cc in range(1,20): fams.append([cc*t for t in range(k)])
    # block+outlier X
    for X in [9,10,11,12,15,20,30,50,80,120,180]:
        fams.append(list(range(9))+[X])
    for X in [9,10,12,20,40,100]:
        fams.append(list(range(8))+[X,X+1])
    # compact diam 9..20 blocks-with-holes (sample: all diam<=14 exhaustive)
    for Wd in range(9,15):
        for e in combinations(range(1,Wd),k-2):
            fams.append((0,)+e+(Wd,))
    return fams
FAMS=residual_families()
print(f"testing {len(FAMS)} candidate residual families per shape\n")

worst_overall=(1e9,None,None)
per_shape=[]
for P in SHAPES:
    c,meas,ct=carr(P); cs=cst(meas)
    worst=(1e9,None)
    a_fail_b_ok=0; both_fail=0; nfail_a=0
    for E in FAMS:
        ra=rho_a(E,c,ct,meas); rb=rho_b(E,meas); r=max(ra,rb)
        if ra < mp: nfail_a+=1
        if r < worst[0]: worst=(r,tuple(E))
        if ra<mp and rb>=mp: a_fail_b_ok+=1
        if ra<mp and rb<mp: both_fail+=1
    per_shape.append((worst[0]-mp, P, worst[1], nfail_a, both_fail))
    if worst[0]<worst_overall[0]: worst_overall=(worst[0],P,worst[1])
per_shape.sort()
nclosed=sum(1 for d,_,_,_,_ in per_shape if d>=0)
print(f"shapes CLOSED by max(a,b) (all tested families >= m_P): {nclosed}/{len(SHAPES)}")
print(f"min over all (shape,family) of max(a,b) - m_P = {worst_overall[0]-mp:+.5f}  "
      f"at P={worst_overall[1]}, E={worst_overall[2]}")
print(f"\nworst 12 shapes (deficit, shape, witness family, #a-fails, #both-fail):")
for d,P,E,na,bf in per_shape[:12]:
    Es = f"{E[:4]}..{E[-1]}" if E and len(E)>6 else E
    print(f"  {d:+.5f}  P={P}  E={Es}  a-fails={na}  both-fail={bf}")
