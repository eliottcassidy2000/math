#!/usr/bin/env python3
"""t4_strip_referee_kps_S128c41.py -- kind-pasteur S128 cont.41.
THE T4 RESONANCE STRIP: decompose the support-4 tail by k = v3 u + v4 t:
STRIP (|k| <= K0 = v1'v2'): inner two-pole separation Delta_k = |k|/(v1'v2') <= 1 (no decay);
the strip outer points lie on <= 2K0+1 coset LINES (one per k).  COMPLEMENT: Delta >= |k|/(v1'v2') decays.
Referee: (a) exact strip vs complement masses on dissociated quadruples at varying H;
(b) the envelope (strip+complement)*H/L^3 bounded; (c) the per-line structure (points on k-lines)."""
import sys
from math import gcd, sin, pi, log
sys.stdout.reconfigure(line_buffering=True)
LAM=1/14
def c(h): return 2*LAM if h==0 else abs(sin(2*pi*h*LAM)/(pi*h))
def t4_split(v1,v2,v3,v4,H,BOX=250):
    """full-support tail split into strip (|k| <= v1v2/g12-scale) and complement; exact enumeration in box."""
    g12=gcd(v1,v2); K0=(v1//g12)*(v2//g12)
    K0=min(K0, 4000)
    strip=0.0; comp=0.0; npts_strip=0; npts_comp=0
    for u in range(-BOX,BOX+1):
        if u==0: continue
        for t in range(-BOX,BOX+1):
            if t==0: continue
            k=v3*u+v4*t
            # inner solutions: v1 h1 + v2 h2 = -k ; need g12 | k
            if k% g12!=0: continue
            # enumerate inner h1 line: h1 = h1_0 + m*(v2/g12)
            v1r,v2r=v1//g12,v2//g12
            # ext gcd inverse
            def inv(x,m):
                x%=m; r0,r1,s0,s1=m,x,0,1
                while r1: q=r0//r1; r0,r1=r1,r0-q*r1; s0,s1=s1,s0-q*s1
                return s0%m if r0==1 else None
            ia=inv(v1r,v2r) if v2r>1 else 0
            kk=-k//g12
            h10=(ia*kk)%v2r if v2r>1 else 0
            # walk m so |h1|,|h2|<=BOX
            for m in range(-(BOX//max(1,v2r))-2,(BOX//max(1,v2r))+3):
                h1=h10+m*v2r
                if h1==0 or abs(h1)>BOX: continue
                num=kk-v1r*h1
                if v2r==0: continue
                if num % v2r!=0: continue
                h2=num//v2r
                if h2==0 or abs(h2)>BOX: continue
                if max(abs(h1),abs(h2),abs(u),abs(t))<=H: continue
                val=c(h1)*c(h2)*c(u)*c(t)
                if abs(k)<=K0: strip+=val; npts_strip+=1
                else: comp+=val; npts_comp+=1
    return strip,comp,npts_strip,npts_comp,K0
quads=[(307,425,541,671),(800,944,1413,2147),(541,1087,1943,3310)]
print("T4 strip/complement split (box 250):")
for (a,b,cc,d) in quads:
    for H in (10,40,160):
        s,co,ns,nc,K0=t4_split(a,b,cc,d,H)
        L=1+log(2+d)
        env=(s+co)*H/L**3
        print("  (%d,%d,%d,%d) H=%d K0=%d: strip=%.2e (%d pts) comp=%.2e (%d pts)  env=%.3e"%(
            a,b,cc,d,H,K0,s,ns,co,nc,env))
print("DONE")
