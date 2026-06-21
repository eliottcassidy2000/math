#!/usr/bin/env python3
"""
lrc_gk8_correction_uniform_macmini_0622s24.py  (mac-mini-2026-06-22-S24)
The gK8 concentration L_yK8(E) <= 10cap, via L_yK8 = DECORRELATED main + correction.
Decorrelated miss-dist: k-1 nonzero runners land in iid uniform sectors {0..6}; q_t^dec = P(exactly t
of the 6 inner sectors empty). L_yK8^dec(k) = 10 q0^dec + q3^dec + 10 q6^dec.
If correction(E)=L_yK8(E)-L_yK8^dec(k) is UNIFORMLY bounded < (10cap - L_yK8^dec), the concentration
follows WITHOUT 'consec is max' -- decorrelated main + bounded correction.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(1)
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
        if t<=6: q[t]+=x1-x0
    return q
def Lyk8(q): return 10*q[0]+q[3]+10*q[6]
def decorr_missdist(k):
    """k-1 nonzero runners iid uniform on {0..6}; q_t = P(exactly t of inner sectors {1..6} empty)."""
    m=k-1; q=[F(0)]*7
    # enumerate occupancy: each runner picks a sector; count empty inner sectors. Use exact prob.
    # P(set A of inner sectors all empty) = ((7-|A|)/7)^m ; incl-excl for exactly t.
    from math import comb
    for t in range(7):
        # P(exactly t inner empty) = C(6,t) sum_{i} (-1)^i C(6-t,i) ((7-t-i)/7)^m
        s=F(0)
        for i in range(0,6-t+1):
            s+=(-1)**i*comb(6-t,i)*F(7-t-i,7)**m
        q[t]=comb(6,t)*s
    return q
caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7),11:F(5,7),12:F(6,7)}
print(f"{'k':>3}{'Lyk8_dec':>10}{'10cap':>9}{'room=10cap-dec':>15}{'consec corr':>13}{'max|corr| sample':>18}")
for k in (8,9,10,11,12):
    qd=decorr_missdist(k); Ld=Lyk8(qd); cap10=10*caps[k]; room=cap10-Ld
    cons=tuple(range(k)); corr_c=Lyk8(missdist(cons))-Ld
    # sample configs, max correction
    mxc=corr_c; argmx='consec'
    for _ in range(50):
        kk=k; E=tuple(sorted([0]+random.sample(range(1,40),kk-1)))
        c=Lyk8(missdist(E))-Ld
        if c>mxc: mxc=c; argmx=str(E)[:30]
    print(f"{k:>3}{float(Ld):>10.4f}{float(cap10):>9.4f}{float(room):>15.4f}{float(corr_c):>13.4f}{float(mxc):>18.4f}")
    print(f"      => correction <= room? consec corr={float(corr_c):.4f} vs room {float(room):.4f}: {'OK' if corr_c<room else 'OVER'}; max sample corr {float(mxc):.4f} @ {argmx}")
print("\nIf max correction over ALL configs < room (=10cap - Lyk8_dec) UNIFORMLY => gK8 concentration PROVED")
print("(decorrelated main is base-independent; correction = the relation-lattice/R-tail object).")
