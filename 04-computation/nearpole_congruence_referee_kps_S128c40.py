#!/usr/bin/env python3
"""nearpole_congruence_referee_kps_S128c40.py -- referee the punctured near-pole congruence lemma:
for dissociated triples, the near-pole portion of the full-support tail obeys <= C(1+ln Vmax)/H.
Near-pole point on line t: the (h1,h2) with minimal |h1| (and the one with minimal |h2|).
Verify: (a) r(t) = least|kappa t mod b'| congruence structure; (b) dissociation floor kicks in
(surviving near-pole points have a coordinate > H); (c) the summed near-pole mass * H stays bounded."""
import sys
from math import gcd, log
sys.stdout.reconfigure(line_buffering=True)
def nearpole_sum(a,b,cc,H,TMAX=6000):
    g=gcd(a,b); g0=gcd(g,cc); d=g//g0
    a2,b2=a//g,b//g; c2=cc//g0
    # kappa: h1 = -c2 t * inv(a2) mod b2
    def inv(x,m):
        x%=m
        r0,r1,s0,s1=m,x,0,1
        while r1: q=r0//r1; r0,r1=r1,r0-q*r1; s0,s1=s1,s0-q*s1
        return s0%m if r0==1 else None
    ia=inv(a2,b2)
    tot=0.0; cnt=0
    for t in range(1,TMAX+1):
        # near-pole in h1: minimal |h1| with a2 h1 = -c2 t mod b2
        if ia is not None and b2>1:
            r=(-c2*t*ia)%b2
            r=min(r,b2-r)
            if r>0:
                h1=r  # (sign choice irrelevant for magnitudes)
                # h2 from a2 h1 + b2 h2 = -c2 t (choose the branch minimizing |h2| consistent)
                num=-c2*t-a2*h1
                if num%b2==0:
                    h2=num//b2
                else:
                    h1=-r; num=-c2*t-a2*h1
                    h2=num//b2 if num%b2==0 else None
                if h2 is not None and h2!=0:
                    h3=d*t
                    if max(abs(h1),abs(h2),abs(h3))>H:
                        tot+=1.0/(abs(h1)*abs(h2)*abs(h3))
                        cnt+=1
        # near-pole in h2 (symmetric): minimal |h2| with b2 h2 = -c2 t mod a2
        if a2>1:
            ib=inv(b2,a2)
            if ib is not None:
                r=(-c2*t*ib)%a2
                r=min(r,a2-r)
                if r>0:
                    h2=r; num=-c2*t-b2*h2
                    if num%a2==0: h1=num//a2
                    else:
                        h2=-r; num=-c2*t-b2*h2
                        h1=num//a2 if num%a2==0 else None
                    if h1 is not None and h1!=0:
                        h3=d*t
                        if max(abs(h1),abs(h2),abs(h3))>H:
                            tot+=1.0/(abs(h1)*abs(h2)*abs(h3))
                            cnt+=1
    return tot,cnt
triples=[(307,425,541),(671,944,1413),(541,800,2147),(1087,1943,3310),(313,741,1531),(420,451,873)]
print("near-pole mass NP(H) and the lemma envelope NP*H/(1+ln Vmax):")
for (a,b,cc) in triples:
    line="(%d,%d,%d):"%(a,b,cc)
    for H in (5,20,80,320):
        np_,cnt=nearpole_sum(a,b,cc,H)
        env=np_*H/(1+log(cc))
        line+="  H=%d: NP=%.2e env=%.2e"%(H,np_,env)
    print(line)
print("DONE")
