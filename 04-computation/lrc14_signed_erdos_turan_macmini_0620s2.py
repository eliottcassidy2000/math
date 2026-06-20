#!/usr/bin/env python3
"""
lrc14_signed_erdos_turan_macmini_0620s2.py  (mac-mini-2026-06-20-S2)

SIGNED Erdos-Turan packet estimate for the far-element deviation Delta_w (sharpens THM-546).
THM-546 absolute bound |Delta_w| <= kappa V(E')/(pi^2 w) is ~5-30x loose. Sharpen using:
 (1) the QR REALITY (HYP-2657): Delta_w = Sum_j Sum_{n>=1} 2 Re[ shat_j(n) * 1hat_{Bj}(-n w) ]
     (imag parts cancel via n<->-n since 6=-1 is NQR mod 7) -> work with real parts only.
 (2) an Erdos-Turan CUTOFF M: |Delta_w| <= |head_{n<=M}| + tail, tail <= (V/pi^2 w)(1/M).
 (3) the 7-vanishing shat_j(7m)=0 (apex prime) thins the sum.
Compare: exact Delta_w (rational engine), the absolute THM-546 bound, and the signed
head+tail bound at cutoff M. Report the SHARPENING factor and the implied smaller cutoff.
shat_j(n) = exp(-2pi i n j/7) * (1-exp(-2pi i n/7))/(2pi i n)   [sector j Fourier coeff].
1hat_{Bj}(m) = sum over arcs (a,b) of (exp(-2pi i m a)-exp(-2pi i m b))/(2pi i m).
"""
import sys, itertools, cmath, math
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
SEV=F(1,7); TWO_PI=2*math.pi

def sector_of(p): return int((p%1)*7)
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        if len(set(sector_of(e*xm) for e in E))==7: tot+=x1-x0
    return tot
def Bj_arcs(E):
    """exact arcs (as float (lo,hi)) of B_j = {x: E misses EXACTLY sector j}, j=0..6; and p1."""
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); arcs=[[] for _ in range(7)]; cur=[None]*7; p1=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set(sector_of(e*xm) for e in E)
        missed=[j for j in range(7) if j not in secs]
        mj=missed[0] if len(missed)==1 else None
        for j in range(7):
            if j==mj:
                if cur[j] is None: cur[j]=x0
            else:
                if cur[j] is not None: arcs[j].append((cur[j],x0)); cur[j]=None
        if mj is not None: p1+=x1-x0
    for j in range(7):
        if cur[j] is not None: arcs[j].append((cur[j],F(1)))
    return arcs, p1
def delta_w_exact(Ep,w):
    E=sorted(set(list(Ep)+[w])); p0E=measS7(Ep); p0Ew=measS7(E)
    _,p1=Bj_arcs(Ep); return float(p0Ew-(p0E+SEV*p1))

def shat(n,j):
    if n==0: return 1/7
    return cmath.exp(-1j*TWO_PI*n*j/7)*(1-cmath.exp(-1j*TWO_PI*n/7))/(1j*TWO_PI*n)
def Bjhat(m, arcs_j):
    if m==0: return sum(float(b-a) for a,b in arcs_j)
    s=0j
    for a,b in arcs_j:
        af,bf=float(a),float(b)
        s+=(cmath.exp(-1j*TWO_PI*m*af)-cmath.exp(-1j*TWO_PI*m*bf))/(1j*TWO_PI*m)
    return s

def signed_head(Ep,w,M):
    arcs,_=Bj_arcs(Ep); tot=0.0; V=sum(len(arcs[j]) for j in range(1,7))
    for j in range(1,7):
        for n in range(1,M+1):
            if n%7==0: continue
            tot+=2*(shat(n,j)*Bjhat(-n*w,arcs[j])).real
    return tot, V

kappa=2*sum(abs(math.sin(math.pi*n/7))/n**2 for n in range(1,100000) if n%7!=0)
print(f"kappa={kappa:.5f}")
print("="*100)
print(f"{'E_prime':<26}{'w':>4}{'Delta_w':>11}{'absBound':>10}{'signedHead(M)':>14}{'tail<=':>9}{'sgnBnd':>9}{'sharpen':>8}")
print("="*100)
cases=[([0,1,2,3,4,5,6,7],[20,50]),([0,1,2,4,6,7,8,10],[24,60]),
       ([0,1,2,4,8,12,16,20],[24,48]),([0,3,5,16,28,30,33],[70])]
M=400
for Ep,ws in cases:
    for w in ws:
        dw=delta_w_exact(Ep,w); head,V=signed_head(Ep,w,M)
        absB=kappa*V/math.pi**2/w; tail=V/math.pi**2/w/M
        sgnB=abs(head)+tail; sharp=absB/sgnB if sgnB>0 else float('inf')
        print(f"{str(Ep):<26}{w:>4}{dw:>11.6f}{absB:>10.4f}{head:>14.6f}{tail:>9.5f}{sgnB:>9.5f}{sharp:>8.1f}x")
print("="*100)
print("signedHead = actual signed sum truncated at M (the head of Erdos-Turan).")
print("sgnBnd = |head| + tail. sharpen = absBound/sgnBnd (how much tighter the signed bound is).")
print("If sgnBnd tracks |Delta_w| (much < absBound), the signed estimate shrinks the finite base.")
print("\nDONE.")
