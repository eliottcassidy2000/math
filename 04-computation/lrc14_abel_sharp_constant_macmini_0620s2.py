#!/usr/bin/env python3
"""
lrc14_abel_sharp_constant_macmini_0620s2.py  (mac-mini-2026-06-20-S2)

Closed-form SIGNED (Abel-summation) sharpening of THM-546.
term_j = Sum_{n!=0} shat_j(n) 1hat_{Bj}(-n w)
       = (1/(-2pi i w)) Sum_{arcs (a,b)} [ H_j(w a) - H_j(w b) ],
  H_j(y) = Sum_{n!=0} (shat_j(n)/n) e(n y) = 2pi i * F_j(y),
  F_j(y) = centered antiderivative of (1_{sector_j} - 1/7)  (a sawtooth, slope +6/7 on sector j,
           -1/7 elsewhere; periodic, mean-zero-adjusted).
=> |Delta_w| <= C_H * V(E') / (pi w),  C_H = max_j sup_y |H_j(y)| = 2pi * max_j (sup-to-... |F_j|).
Compute C_H exactly (F_j is piecewise linear), compare C_H/pi to the THM-546 absolute kappa/pi^2,
and verify |Delta_w| <= C_H V/(pi w) on the test cores.  Also report the EXACT closed form of C_H.
"""
import sys, itertools, math
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
SEV=F(1,7)

# F_j: centered antiderivative of (1_{[j/7,(j+1)/7)} - 1/7). Sawtooth.
# On [0,1): g(t)=1_{sector j}(t)-1/7. Antideriv A(y)=int_0^y g. A is piecewise linear:
#   slope -1/7 on [0,j/7], +6/7 on [j/7,(j+1)/7], -1/7 on [(j+1)/7,1]. A(0)=0,A(1)=0 (mean-zero g).
# F_j = A - mean(A). sup|F_j| = max|A-mean(A)|. H_j=2pi i F_j so |H_j|=2pi|F_j|, sup|H_j|=2pi sup|F_j|.
def supF_exact(j):
    # vertices of A at y=0, j/7, (j+1)/7, 1
    ys=[F(0),F(j,7),F(j+1,7),F(1)]
    A={}
    A[F(0)]=F(0)
    A[F(j,7)]=A[F(0)]+(-F(1,7))*(F(j,7)-F(0))
    A[F(j+1,7)]=A[F(j,7)]+(F(6,7))*(F(1,7))
    A[F(1)]=A[F(j+1,7)]+(-F(1,7))*(F(1)-F(j+1,7))
    # mean of A (piecewise linear): integral of A over [0,1] = sum of trapezoids
    mean=F(0)
    for i in range(len(ys)-1):
        y0,y1=ys[i],ys[i+1]; mean+=(A[y0]+A[y1])/2*(y1-y0)
    # sup|A-mean| at the vertices (piecewise linear => extrema at vertices)
    return max(abs(A[y]-mean) for y in ys)

print("EXACT sup|F_j| and C_H = 2*pi*max_j sup|F_j|:")
supFs=[supF_exact(j) for j in range(1,7)]
for j in range(1,7): print(f"   j={j}: sup|F_j| = {supFs[j-1]} = {float(supFs[j-1]):.5f}")
maxsupF=max(supFs); C_H=2*math.pi*float(maxsupF)
print(f"   max_j sup|F_j| = {maxsupF} = {float(maxsupF):.5f};  C_H = 2*pi*that = {C_H:.5f}")
kappa=2*sum(abs(math.sin(math.pi*n/7))/n**2 for n in range(1,100000) if n%7!=0)
print(f"\n   THM-546 absolute constant: kappa/pi^2 = {kappa/math.pi**2:.5f}")
print(f"   Abel signed constant:      C_H/pi      = {C_H/math.pi:.5f}   "
      f"(sharper by {(kappa/math.pi**2)/(C_H/math.pi):.2f}x)")

# verify |Delta_w| <= C_H V/(pi w) on cores
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
def V_and_p1(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); arccount=[0]*7; inBj=[False]*7; p1=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set(sector_of(e*xm) for e in E)
        missed=[j for j in range(7) if j not in secs]; mj=missed[0] if len(missed)==1 else None
        for j in range(7):
            if j==mj:
                if not inBj[j]: arccount[j]+=1; inBj[j]=True
            else: inBj[j]=False
        if mj is not None: p1+=x1-x0
    return sum(arccount[1:7]), p1
print(f"\n{'E_prime':<26}{'w':>4}{'|Delta_w|':>11}{'V':>5}{'C_H V/(pi w)':>14}{'OK?':>5}")
for Ep,ws in [([0,1,2,3,4,5,6,7],[20,50]),([0,1,2,4,8,12,16,20],[24,48]),([0,3,5,16,28,30,33],[70])]:
    V,p1=V_and_p1(Ep)
    for w in ws:
        dw=abs(float(measS7(sorted(set(Ep+[w])))-(measS7(Ep)+SEV*p1)))
        bnd=C_H*V/math.pi/w
        print(f"{str(Ep):<26}{w:>4}{dw:>11.6f}{V:>5}{bnd:>14.5f}{'yes' if dw<=bnd+1e-9 else 'NO':>5}")
print("\nDONE.")
