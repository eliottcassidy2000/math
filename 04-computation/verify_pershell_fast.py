#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fast per-shell signed/abs sums for k=8 AP, plus L=9 box-cum, to confirm the
NEGATIVE decay claim independently. Precompute D7 over all 6^6 cosets once,
then a single enumeration. Uses the SAME factorization K(n)=D7(nmod7)/prod n.
"""
import sys, itertools, math, cmath
from collections import defaultdict
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def A_coef(r): return (1 - cmath.exp(-2j*math.pi*r/7.0))/(2j*math.pi)
def hT(T,r): return A_coef(r)*sum(cmath.exp(-2j*math.pi*r*j/7.0) for j in T)
Ts=[T for cnt in range(7) for T in itertools.combinations(range(1,7),cnt)]
sgn=[(-1)**len(T) for T in Ts]
hc={(ti,r):hT(T,r) for ti,T in enumerate(Ts) for r in range(1,7)}

# precompute D7 for all residue tuples (real part only needed)
D7re={}
for c in itertools.product(range(1,7),repeat=6):
    tot=0j
    for ti in range(len(Ts)):
        p=1+0j
        for r in c: p*=hc[(ti,r)]
        tot+=sgn[ti]*p
    D7re[c]=tot.real   # correction is real -> use Re

E=[0,1,2,3,4,5,6,7]; k=8
L=9
shell_s=defaultdict(float); shell_a=defaultdict(float)
total=0.0
for idxs in itertools.combinations(range(k),6):
    es=[E[i] for i in idxs]
    dep=max(range(6),key=lambda i:abs(es[i])); e_dep=es[dep]
    free=[i for i in range(6) if i!=dep]; efree=[es[i] for i in free]
    for vf in itertools.product(range(-L,L+1),repeat=5):
        bad=False
        for c in vf:
            if c==0 or c%7==0: bad=True;break
        if bad: continue
        s=0
        for c,e in zip(vf,efree): s+=c*e
        if e_dep==0:
            if s!=0: continue
            for vd in range(-L,L+1):
                if vd==0 or vd%7==0: continue
                combo=[0]*6
                for i,c in zip(free,vf): combo[i]=c
                combo[dep]=vd
                res=tuple(v%7 for v in combo); ip=1.0
                for v in combo: ip/=v
                kv=ip*D7re[res]
                ninf=max(abs(x) for x in combo)
                total+=kv; shell_s[ninf]+=kv; shell_a[ninf]+=abs(kv)
        else:
            if s%e_dep!=0: continue
            vd=-s//e_dep
            if vd==0 or vd%7==0 or abs(vd)>L: continue
            combo=[0]*6
            for i,c in zip(free,vf): combo[i]=c
            combo[dep]=vd
            res=tuple(v%7 for v in combo); ip=1.0
            for v in combo: ip/=v
            kv=ip*D7re[res]
            ninf=max(abs(x) for x in combo)
            total+=kv; shell_s[ninf]+=kv; shell_a[ninf]+=abs(kv)

target=0.302731
print(f"L=9 box-cum corr = {total:+.6f}  frac of target = {total/target:+.4f}")
print("per-shell signed / abs:")
csig=0.0
for s in sorted(shell_s):
    csig+=shell_s[s]
    print(f"  |n|inf={s}: signed={shell_s[s]:+.5f}  abs={shell_a[s]:.5f}  cumsigned={csig:+.5f}")
