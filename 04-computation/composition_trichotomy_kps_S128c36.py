#!/usr/bin/env python3
"""composition_trichotomy_kps_S128c36.py -- kind-pasteur S128 cont.36.
(A) referee the relation-mass identity B_m = B_m^eq + sum_s (-1)^s E_s^(m) M_s (n=3 exact-ish).
(B) THE COMPOSITION ALGORITHM: certify (exact B5>0) OR surrender small relations (stratum witness)
    OR flag MIDDLE. Battery: structured + small-scale + random packets. Exhaustion = no MIDDLE."""
import sys, random
from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd, sin, pi
sys.stdout.reconfigure(line_buffering=True)
random.seed(360360)
LAM=1/14
def c(h): return 2*LAM if h==0 else sin(2*pi*h*LAM)/(pi*h)
def depth_mu(speeds):
    ev=[]
    for v in speeds:
        for j in range(v):
            lo=F(14*j-1,14*v); hi=F(14*j+1,14*v)
            if lo<0: ev.append((F(0),1)); ev.append((hi,-1)); ev.append((lo+1,1)); ev.append((F(1),-1))
            else: ev.append((lo,1)); ev.append((hi,-1))
    ev.sort()
    mu=[F(0)]*(len(speeds)+1); depth=0; last=F(0)
    for x,d in ev:
        if x>last: mu[depth]+=x-last; last=x
        depth+=d
    if F(1)>last: mu[depth]+=F(1)-last
    return mu
def B5exact(V):
    mu=depth_mu(V); n=len(V)
    return mu[0]-sum(comb(d-1,5)*mu[d] for d in range(6,n+1))

print("== (A) identity referee, V={1,2,3}, H=4000 ==")
V=[1,2,3]; n=3; H=4000
mu=depth_mu(V)
S=[sum(comb(d,k)*mu[d] for d in range(4)) for k in range(4)]
M2=0.0
for (i,j) in combinations(range(3),2):
    a,b=V[i],V[j]; g=gcd(a,b); pa,pb=b//g,a//g
    m=1
    while max(pa*m,pb*m)<=H:
        M2+=2*c(pa*m)*c(pb*m); m+=1
M3=0.0
for h1 in range(-H,H+1):
    if h1==0: continue
    for h2 in range(-H,H+1):
        if h2==0: continue
        r=h1+2*h2
        if r%3==0:
            h3=-r//3
            if h3!=0 and abs(h3)<=H: M3+=c(h1)*c(h2)*c(h3)
for m_ in (2,3):
    Bm=float(sum((-1)**k*S[k] for k in range(m_+1)))
    Bmeq=float(sum((-1)**k*comb(3,k)*(1/7)**k for k in range(m_+1)))
    pred=Bmeq
    for s,Ms in ((2,M2),(3,M3)):
        if s<=m_:
            Es=sum((-1)**j*comb(3-s,j)*(1/7)**j for j in range(0,m_-s+1))
            pred+=((-1)**s)*Es*Ms
    print("  m=%d: B_m exact %.8f | eq %.8f | eq + sum (-1)^s E_s M_s = %.8f  err %.1e"%(m_,Bm,Bmeq,pred,abs(Bm-pred)))

print("== (B) THE COMPOSITION ALGORITHM ==")
def small_relations(V,H2=30,H34=8,H5=3):
    """return witness list of small relations: (support-size, subset, coeffs)"""
    n=len(V); wit=[]
    for i in range(n):
        for j in range(i+1,n):
            a,b=V[i],V[j]; g=gcd(a,b)
            if b//g<=H2 and a//g<=H2: wit.append((2,(V[i],V[j]),(b//g,-(a//g))))
    for A in combinations(range(n),3):
        a,b,cc=[V[x] for x in A]
        for h1 in range(-H34,H34+1):
            for h2 in range(-H34,H34+1):
                if h1==0 or h2==0: continue
                r=h1*a+h2*b
                if r!=0 and r%cc==0 and 0<abs(r//cc)<=H34:
                    wit.append((3,(a,b,cc),(h1,h2,-r//cc)))
    for A in combinations(range(n),4):
        vs=[V[x] for x in A]
        for h1 in range(-H34,H34+1):
            for h2 in range(-H34,H34+1):
                for h3 in range(-H34,H34+1):
                    if 0 in (h1,h2,h3): continue
                    r=h1*vs[0]+h2*vs[1]+h3*vs[2]
                    if r!=0 and r%vs[3]==0 and 0<abs(r//vs[3])<=H34:
                        wit.append((4,tuple(vs),(h1,h2,h3,-r//vs[3])))
        if len(wit)>200: break
    return wit
def compose(name,V):
    b5=B5exact(V)
    if b5>0:
        print("  %-22s B5=%+.5f  CERTIFIED (level-5)"%(name,float(b5))); return 'C'
    wit=small_relations(V)
    if wit:
        s2=[w for w in wit if w[0]==2]; s34=[w for w in wit if w[0]>2]
        route="ratio/lacunary" if s2 else "linear-forms/near-AP"
        print("  %-22s B5=%+.5f  STRUCTURED -> %s (witnesses: %d ratio, %d linform; e.g. %s)"%(
            name,float(b5),route,len(s2),len(s34),(s2+s34)[0][1:]))
        return 'S'
    print("  %-22s B5=%+.5f  *** MIDDLE (no small witness!) ***"%(name,float(b5))); return 'M'
battery=[]
battery.append(("tight",list(range(1,14))))
battery.append(("deepwell",list(range(1,13))+[182]))
battery.append(("geometric x2",[5*2**i for i in range(13)]))
battery.append(("geometric x3",[7*3**i for i in range(13)]))
battery.append(("opus30Z",[420,450,510,570,690,870,1230,1770,2370,3210,4170,5190,7230]))
battery.append(("AP d=17",[17*i for i in range(1,14)]))
battery.append(("AP+stranger",[17*i for i in range(1,13)]+[493]))
battery.append(("near-AP pert",[17,34,51,68,85,102,119,136,153,170,187,204,225]))
for trial in range(12):
    lo=random.choice([13,20,40,80])
    xs=sorted(random.sample(range(lo,lo*14),13))
    battery.append(("small-scale [%d,%d]"%(lo,lo*14),xs))
for trial in range(6):
    xs=sorted(random.sample(range(300,4001),13))
    battery.append(("random [300,4000]",xs))
counts={'C':0,'S':0,'M':0}
for name,V in battery:
    counts[compose(name,V)]+=1
print("VERDICT: certified %d | structured-with-witness %d | MIDDLE %d"%(counts['C'],counts['S'],counts['M']))
print("EXHAUSTION %s on this battery"%("HOLDS (zero middle)" if counts['M']==0 else "FAILS"))
print("== (C) the s=2 tail lemma constant + pigeonhole floor ==")
print("  |E2| * M2-tail(>H) <= (24/343)*13/H: at H=30: %.5f  (vs B5_eq = 2052/16807 = %.5f)"%((24/343)*13/30,2052/16807))
print("  pigeonhole: max(V) < 40 => 78 pairwise sums in [3,2*40-1=79] => Sidon collision forced => +-1-coeff support<=4 relation")
print("DONE")
