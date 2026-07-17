#!/usr/bin/env python3
"""t5_slab_referee_kps_S128c42.py -- kind-pasteur S128 cont.42.
THE T5 SLAB: k = w3 u + w4 t + w5 r; slab {|k| <= K0} = coset PLANES.  Referee:
(a) slab vs complement masses on dissociated quintuples; (b) envelope T5*H/L^4 bounded;
(c) within-plane structure: the three outer coordinates' vanishing lines (sub-resonances)
carry the mass concentration -- count near-line points vs plane bulk."""
import sys
from math import gcd, sin, pi, log
sys.stdout.reconfigure(line_buffering=True)
LAM=1/14
def c(h): return 2*LAM if h==0 else abs(sin(2*pi*h*LAM)/(pi*h))
def inv(x,m):
    x%=m; r0,r1,s0,s1=m,x,0,1
    while r1: q=r0//r1; r0,r1=r1,r0-q*r1; s0,s1=s1,s0-q*s1
    return s0%m if r0==1 else None
def t5_split(w,H,BOX=60):
    w1,w2,w3,w4,w5=w
    g12=gcd(w1,w2); K0=min((w1//g12)*(w2//g12),3000)
    w1r,w2r=w1//g12,w2//g12
    ia=inv(w1r,w2r) if w2r>1 else 0
    slab=0.0; comp=0.0; ns=0; nc=0; nearline=0.0
    for u in range(-BOX,BOX+1):
        if u==0: continue
        for t in range(-BOX,BOX+1):
            if t==0: continue
            for r in range(-BOX,BOX+1):
                if r==0: continue
                k=w3*u+w4*t+w5*r
                if k%g12!=0: continue
                kk=-k//g12
                h10=(ia*kk)%w2r if w2r>1 else 0
                for m in range(-(BOX//max(1,w2r))-2,(BOX//max(1,w2r))+3):
                    h1=h10+m*w2r
                    if h1==0 or abs(h1)>BOX: continue
                    num=kk-w1r*h1
                    if num%w2r!=0: continue
                    h2=num//w2r
                    if h2==0 or abs(h2)>BOX: continue
                    if max(abs(h1),abs(h2),abs(u),abs(t),abs(r))<=H: continue
                    val=c(h1)*c(h2)*c(u)*c(t)*c(r)
                    if abs(k)<=K0:
                        slab+=val; ns+=1
                        if min(abs(u),abs(t),abs(r))<=2: nearline+=val
                    else: comp+=val; nc+=1
    return slab,comp,ns,nc,nearline,K0
quints=[(307,425,541,671,800),(541,800,1087,1943,2570),(671,944,1413,2147,3310)]
print("T5 slab/complement split (box 60):")
for w in quints:
    for H in (10,40):
        s,co,ns,nc,nl,K0=t5_split(w,H)
        L=1+log(2+w[4])
        env=(s+co)*H/L**4
        print("  %s H=%d K0=%d: slab=%.2e (%d pts, near-line frac %.0f%%) comp=%.2e (%d)  env=%.3e"%(
            w,H,K0,s,ns,100*nl/s if s>0 else 0,co,nc,env))
print("DONE")
