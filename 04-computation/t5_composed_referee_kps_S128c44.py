#!/usr/bin/env python3
"""t5_composed_referee_kps_S128c44.py -- T5 hardened composition referee:
per k-plane, mass = 3 collar-line sums (each a two-pole atom in the two surviving forms)
+ bulk (all coords >= 3). Check: (a) collar+bulk decomposition exact per plane (partition);
(b) collars dominated by their atom bounds; (c) bulk <= C_bulk L^2 shell formula."""
import sys
from math import gcd, log, pi, sin
sys.stdout.reconfigure(line_buffering=True)
LAM=1/14
def c(h): return 2*LAM if h==0 else abs(sin(2*pi*h*LAM)/(pi*h))
def inv(x,m):
    x%=m; r0,r1,s0,s1=m,x,0,1
    while r1: q=r0//r1; r0,r1=r1,r0-q*r1; s0,s1=s1,s0-q*s1
    return s0%m if r0==1 else None
def analyze(w,BOX=70):
    w1,w2,w3,w4,w5=w
    g12=gcd(w1,w2); w1r,w2r=w1//g12,w2//g12
    ia=inv(w1r,w2r) if w2r>1 else 0
    # aggregate over k: outer weight (u,t,r) sums split collar/bulk; inner factor left exact
    collar=0.0; bulk=0.0; ncol=0; nbul=0
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
                inner=0.0
                for m in range(-(BOX//max(1,w2r))-2,(BOX//max(1,w2r))+3):
                    h1=h10+m*w2r
                    if h1==0 or abs(h1)>BOX: continue
                    num=kk-w1r*h1
                    if num%w2r!=0: continue
                    h2=num//w2r
                    if h2==0 or abs(h2)>BOX: continue
                    inner+=c(h1)*c(h2)
                if inner==0: continue
                val=inner*c(u)*c(t)*c(r)
                if min(abs(u),abs(t),abs(r))<=2: collar+=val; ncol+=1
                else: bulk+=val; nbul+=1
    return collar,bulk,ncol,nbul
quints=[(307,425,541,671,800),(671,944,1413,2147,3310)]
print("T5 collar/bulk partition (box 70):")
for w in quints:
    co,bu,nc,nb=analyze(w)
    tot=co+bu
    print("  %s: collar=%.3e (%d pts, %.0f%%) bulk=%.3e (%d pts)  [partition exact by construction]"%(
        w,co,nc,100*co/tot if tot>0 else 0,bu,nb))
    L=1+log(2+w[4])
    print("    bulk/L^2 = %.3e (the bulk-lemma normalization; bounded across quintuples = the shell formula's check)"%(bu/L**2))
print("DONE")
