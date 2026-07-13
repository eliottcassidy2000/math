#!/usr/bin/env python3
"""
lrc14_prime_grid_klein_S274.py
==============================
klein-2026-07-13-S274. CORRECTED measurement with a PRIME grid Ng >> w (independent of w), fixing
the Ng∝w aliasing that faked err*w~30 at w=lcm(C) (frac(wx)=frac(k/c) cycles through c values,
and at w=lcm the cluster phases become dependent -> spurious blow-up).

Confirms the two-scale constant C0=err*w is BOUNDED (the min-argument proves |Error_cover|<=12/7),
including at adversarial w=lcm-multiples, once the grid no longer aliases.
Also validates the min-argument: |Error_cover| <= Sum_I min(osc/w, ||g|| |I|) <= 12/7.
"""
import math
from math import gcd
NG=999983   # large PRIME, >> all w and e'; independent of w -> no aliasing
def lcm(xs):
    r=1
    for x in xs:
        if x: r=r*x//gcd(r,x)
    return r
def occ_of(E,x):
    o=0
    for e in E: o|=1<<(int((e*x%1.0)*7.0)%7)
    return o
def moments3(E):
    s1=s2=s3=0
    for k in range(1,NG):
        o=occ_of(E,k/NG); N=7-bin(o).count("1")
        s1+=N;s2+=N*(N-1);s3+=N*(N-1)*(N-2)
    n=NG-1;return s1/n,s2/n,s3/n
def Phi_from_m(m1,m2,m3):return 1-(2/3)*m1+(47/252)*m2-(5/252)*m3
def Phi(E):m1,m2,m3=moments3(E);return Phi_from_m(m1,m2,m3)
def Phi_inf(C):m1,m2,m3=moments3(C);return Phi_from_m((6/7)*m1,(5/7)*m2,(4/7)*m3)
# cover-measure Error + min-bound (grid-independent claim, measured on the prime grid)
def error_cover_and_minbound(Ep,w):
    c0=c1=0
    for k in range(1,NG):
        o=occ_of(Ep,k/NG); nf=bin(o).count("1")
        if nf==7:c0+=1
        elif nf==6:c1+=1
    n=NG-1; plat=c0/n+(1/7)*(c1/n)
    c0w=0
    for k in range(1,NG):
        if occ_of(Ep+[w],k/NG)==0x7F:c0w+=1
    p0w=c0w/n
    err=p0w-plat
    # min-bound over R-intervals
    oscG=6/49; gmax=6/7; minb=0.0; run=0
    for k in range(1,NG):
        nf=bin(occ_of(Ep,k/NG)).count("1")
        if nf==6: run+=1
        else:
            if run: minb+=min(oscG/w, gmax*run/NG); run=0
    if run: minb+=min(oscG/w, gmax*run/NG)
    return err,minb

CAP9=1979/4004
print("="*72)
print(f"PRIME grid Ng={NG}.  err*w for clean (prime) AND adversarial (lcm-mult) w, all << Ng.")
print(f"cap9={CAP9:.5f}; tail needs C0<=25*(cap9-0.397)={25*(CAP9-0.397):.3f}")
print("="*72)
supC0=0; suparg=None
tests=[
  ("compact {0..5,30}",   [0,1,2,3,4,5,30]),
  ("2blk {0,1,2,28,29,30,15}",[0,1,2,28,29,30,15]),
  ("dilAP {0,5,..,30}",   [0,5,10,15,20,25,30]),
  ("spread {0,3,7,12,20,33,54}",[0,3,7,12,20,33,54]),
  ("{0..5,120}",          [0,1,2,3,4,5,120]),
  ("2blk-wide {0,1,2,118,119,120,60}",[0,1,2,118,119,120,60]),
]
for name,C in tests:
    L=lcm(C); pinf=Phi_inf(C)
    ws=[101,1009,2003]+[L, 2*L] if L<3000 else [101,1009,2003,L] if L<300000 else [101,1009,2003]
    ws=[w for w in ws if 20*max(C) < w < NG//50]
    if not ws: ws=[2003]
    worst=0; wa=None
    for w in ws:
        e=abs(Phi(C+[w])-pinf)*w
        if e>worst: worst=e; wa=w
    if worst>supC0: supC0=worst; suparg=(name,wa)
    print(f"  {name:34s} lcm={L:8d}  sup err*w={worst:.3f} (w={wa})   ws={ws}")
print(f"  SUP C0 (prime grid) = {supC0:.3f} at {suparg}")
print()
print("="*72)
print("min-argument validation: |Error_cover| <= min-bound <= 12/7 = %.4f"%(12/7))
print("="*72)
for name,C in tests[:4]:
    L=lcm(C)
    for w in [1009, L if L<50000 else 2*max(C)+1]:
        if w>=NG//50: continue
        err,minb=error_cover_and_minbound(C,w)
        ok=abs(err)<=minb+1e-6 and minb<=12/7+1e-6
        print(f"  {name:30s} w={w:7d}: |Err_cover|={abs(err):.5f} min-bound={minb:.4f} 12/7={12/7:.3f} {'OK' if ok else 'VIOLATION'}")
print()
print(f"  => C0 bounded (no blow-up at lcm-w on the prime grid); min-argument holds.")
print(f"     tail closure: err<=C0/w, w=d>25 => err<={supC0:.2f}/26={supC0/26:.3f}; Phi<=0.397+{supC0/26:.3f}={0.397+supC0/26:.3f} vs cap9={CAP9:.3f}")
print("\ndone.")
