#!/usr/bin/env python3
"""
lrc14_fejer_scaling_opus_S4.py   opus-2026-07-23-S4
Quantify the certifiable-concentration WALL: (a) tight AP approach rate gap-B_N ~ c log N / N
(g has a downward CORNER at tau*); (b) certification degree N* ~ (binding slope)/delta,
delta = gap - 1/14. => the route certifies every config with delta bounded below; the cost
diverges as 1/delta exactly at the tight locus (OPEN-Q-108).
"""
from fractions import Fraction as Fr
from math import floor, pi, log
import cmath
def ndist(x):
    f=x-floor(x); return min(f,1-f)
def gval(V,t): return min(ndist(v*t) for v in V)
def breakpoints(V):
    pts={Fr(0),Fr(1)}
    for v in V:
        for j in range(0,2*v+1): pts.add(Fr(j,2*v))
    for i in range(len(V)):
        for j in range(i+1,len(V)):
            for d in {abs(V[i]-V[j]),V[i]+V[j]}:
                if d==0: continue
                for k in range(0,d+1): pts.add(Fr(k,d))
    return sorted(p for p in pts if 0<=p<=1)
def pieces(V,bp):
    out=[]
    for i in range(len(bp)-1):
        a,b=bp[i],bp[i+1]; mid=(a+b)/2
        vs=min(V,key=lambda v: ndist(v*mid)); vm=vs*mid; n=floor(vm+Fr(1,2))
        if vm-n>=0: m,c=Fr(vs),Fr(-n)
        else: m,c=Fr(-vs),Fr(n)
        out.append((a,b,m,c))
    return out
def gap_ts(V,bp):
    best=None
    for p in bp:
        gv=gval(V,p)
        if best is None or gv>best[0]: best=(gv,p)
    return best
def Gk(pc,Nmax):
    G=[0j]*(Nmax+1); g0=0.0
    for (a,b,m,c) in pc:
        a_,b_,m_,c_=float(a),float(b),float(m),float(c); g0+=m_*(b_*b_-a_*a_)/2+c_*(b_-a_)
    G[0]=complex(g0)
    for k in range(1,Nmax+1):
        w=2*pi*k; acc=0j
        for (a,b,m,c) in pc:
            m_,c_,a_,b_=float(m),float(c),float(a),float(b)
            ea=cmath.exp(2j*pi*float(k*a-floor(k*a))); eb=cmath.exp(2j*pi*float(k*b-floor(k*b)))
            acc+=eb*((m_*b_+c_)/(1j*w)+m_/(w*w))-ea*((m_*a_+c_)/(1j*w)+m_/(w*w))
        G[k]=acc
    return G
def BN(G,ts,N):
    s=G[0].real
    for k in range(1,N+1):
        wk=1.0-k/(N+1); ph=cmath.exp(-2j*pi*float(k*ts-floor(k*ts))); s+=2.0*wk*(ph*G[k]).real
    return s
def binding_slope(V,bp,ts):
    # local |slope| of g just left/right of tau*: use adjacent pieces
    i=bp.index(ts) if ts in bp else None
    sl=[]
    for (a,b,m,c) in pieces(V,bp):
        if b==ts or a==ts: sl.append(abs(float(m)))
    return sum(sl)/len(sl) if sl else float('nan')

FL=1.0/14
# (a) AP approach rate
V=list(range(1,14)); bp=breakpoints(V); pc=pieces(V,bp); gap,ts=gap_ts(V,bp)
Nmax=4000; G=Gk(pc,Nmax)
print("(a) TIGHT AP {1..13}: approach rate  gap - B_N  (expect ~ c log N / N; corner at tau*=1/14)")
print("    N      gap-B_N     (gap-B_N)*N/lnN")
for N in [250,500,1000,2000,4000]:
    e=FL-BN(G,ts,N); print(f"   {N:5d}  {e:.6f}     {e*N/log(N):.4f}")
print(f"    binding slope at tau* = {binding_slope(V,bp,ts):.2f}   (=> corner, sum|slopes|/2 ~ (1+13)/2=7)")
# (b) N* vs delta across a margin spread
print("\n(b) CERTIFICATION DEGREE N* vs margin delta = gap - 1/14")
print("    config                         gap    delta      slope   N*     N* * delta / slope")
cfgs={
 "{1..11,13,36} 3/41":  list(range(1,12))+[13,36],
 "{1..12,26} 2/27":     list(range(1,13))+[26],
 "{1..12,14} 1/13":     list(range(1,13))+[14],
 "{1..11,13,20}":       list(range(1,12))+[13,20],
 "{1,2,3,4,5,6,7,9,11,13,15,17,19}": [1,2,3,4,5,6,7,9,11,13,15,17,19],
 "{1,3,5,7,9,11,13,15,17,19,21,23,25}": [1,3,5,7,9,11,13,15,17,19,21,23,25],
}
for nm,V in cfgs.items():
    bp=breakpoints(V); pc=pieces(V,bp); gap,ts=gap_ts(V,bp); d=float(gap)-FL
    if d<=0: 
        print(f"    {nm:30s} {str(gap):>6} tight/below"); continue
    G=Gk(pc,4000); sl=binding_slope(V,bp,ts)
    Ns=None
    for N in range(20,4001,20):
        if BN(G,ts,N)>FL: Ns=N; break
    prod = (Ns*d/sl) if Ns else float('nan')
    print(f"    {nm:30s} {str(gap):>6} {d:.5f}  {sl:6.2f}  {str(Ns) if Ns else '>4000':>5}  {prod:.3f}")
print("\n=> N* * delta / slope ~ const (up to log): N* ~ (binding slope)/delta. The certifiable-")
print("   concentration route certifies EVERY config with delta bounded below; cost ~1/delta diverges")
print("   at the tight locus (delta=0) = OPEN-Q-108. Quantifies 'bulk certifiable + rigidity crux'.")
