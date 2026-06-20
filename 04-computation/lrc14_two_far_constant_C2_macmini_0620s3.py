#!/usr/bin/env python3
"""
lrc14_two_far_constant_C2_macmini_0620s3.py  (mac-mini-2026-06-20-S3)  -> THM-548 Part B

The TWO-FAR curvature constant C_2 (parabolic analogue of the one-far sup|F_j|=3/49).

One-far:  Delta_w = -(1/w) sum_arcs [F_j(wa)-F_j(wb)],  F_j = 1st centered antiderivative of
          psi_j = 1_{sector j}-1/7,  sup|F_j| = 3/49  =>  |Delta_w| <= (6/49) V/w.
Two-far interaction term:  T_{jj'}(u,v) = INT 1_A(x) psi_j(ux) psi_{j'}(vx) dx
          = sum_{m,n!=0} shat_j(m) shat_{j'}(n) 1hat_A(-(mu+nv)).
A SECOND Abel summation (antiderivative in each frequency) introduces G_j = 2nd centered
antiderivative of psi_j (piecewise QUADRATIC).  C_2 := sup_y |G_j(y)| is the parabolic constant.

This script:
 (1) computes sup|G_j| EXACTLY (piecewise quadratic; extrema at vertices AND segment interiors
     where G_j' = F_j = 0) for all j -> the rational constant C_2;
 (2) verifies the QR-REALITY of the product kernel: shat_j(m)shat_{j'}(n) is real under the joint
     (m,n)->(-m,-n) pairing (the 6=-1 NQR mod 7 argument, extended to the product);
 (3) measures the actual two-far interaction sup_{u,v} |T_{jj'}| for a bounded base and its decay
     with the resonance distance dist((u,v)) = min_{small (m,n)} |mu+nv|.
"""
import sys, math, cmath, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
TWO_PI=2*math.pi

# ---------- (1) sup|G_j| exact: G_j = 2nd centered antiderivative of psi_j=1_{sector j}-1/7 ----------
# psi_j on [0,1): value 6/7 on [j/7,(j+1)/7), -1/7 elsewhere.  F_j(y)=int_0^y psi_j, centered.
# G_j(y)=int_0^y F_j, centered.  Both periodic mean-zero.  Compute F_j at nodes, then G_j piecewise-quadratic.
def supG_exact(j):
    # nodes where psi_j changes: 0, j/7, (j+1)/7, 1
    nodes=[F(0),F(j,7),F(j+1,7),F(1)]
    # slope of F_j on each piece = psi_j value there
    def psi_val(seg):  # seg index 0:[0,j/7],1:[j/7,(j+1)/7],2:[(j+1)/7,1]
        return F(6,7) if seg==1 else F(-1,7)
    # F_j(node) values (raw antiderivative from 0), then center
    Fraw={F(0):F(0)}
    segpts=[(F(0),F(j,7),0),(F(j,7),F(j+1,7),1),(F(j+1,7),F(1),2)]
    for (a,b,seg) in segpts:
        Fraw[b]=Fraw[a]+psi_val(seg)*(b-a)
    meanF=F(0)  # integral of F_j over [0,1] (trapezoid since piecewise linear)
    for (a,b,seg) in segpts:
        meanF+=(Fraw[a]+Fraw[b])/2*(b-a)
    Fc=lambda y_node: Fraw[y_node]-meanF  # centered F at a node
    # G_j(y)=int_0^y Fc.  Fc is piecewise linear; G_j piecewise quadratic.
    # Build G at nodes (raw), find interior extrema where Fc=0.
    Graw={F(0):F(0)}; cand=[(F(0),F(0))]  # (y, G value) candidates for sup
    for (a,b,seg) in segpts:
        slope=psi_val(seg)  # Fc'(y)=psi on this seg
        Fa=Fc(a)            # Fc at a
        # G(b)=G(a)+int_a^b Fc = G(a)+Fa*(b-a)+slope*(b-a)^2/2
        Graw[b]=Graw[a]+Fa*(b-a)+slope*(b-a)**2/2
        # interior extremum where Fc(y)=0: Fa+slope*(y-a)=0 -> y=a-Fa/slope
        if slope!=0:
            ystar=a-Fa/slope
            if a<ystar<b:
                Gstar=Graw[a]+Fa*(ystar-a)+slope*(ystar-a)**2/2
                cand.append((ystar,Gstar))
        cand.append((b,Graw[b]))
    meanG=F(0)  # center G: integral of Graw over [0,1] (piecewise quadratic -> Simpson exact)
    for (a,b,seg) in segpts:
        slope=psi_val(seg); Fa=Fc(a); Ga=Graw[a]
        # int_a^b [Ga+Fa(t-a)+slope(t-a)^2/2] dt = Ga*h + Fa*h^2/2 + slope*h^3/6, h=b-a
        h=b-a; meanG+=Ga*h+Fa*h**2/2+slope*h**3/6
    sup=max(abs(g-meanG) for (_,g) in cand)
    return sup

print("(1) sup|G_j| (second centered antiderivative; the parabolic constant C_2):")
sups=[supG_exact(j) for j in range(1,7)]
for j in range(1,7): print(f"    j={j}: sup|G_j| = {sups[j-1]} = {float(sups[j-1]):.6f}")
C2=max(sups)
print(f"    C_2 = max_j sup|G_j| = {C2} = {float(C2):.6f}")
print(f"    (compare one-far sup|F_j| = 3/49 = {float(F(3,49)):.6f})")

# ---------- (2) QR-reality of the product kernel ----------
def shat(n,j):
    if n==0: return 1/7
    return cmath.exp(-1j*TWO_PI*n*j/7)*(1-cmath.exp(-1j*TWO_PI*n/7))/(1j*TWO_PI*n)
print("\n(2) QR-REALITY of product kernel: shat_j(m)shat_{j'}(n) + shat_j(-m)shat_{j'}(-n) real?")
maxim=0.0
for j in range(1,7):
    for jp in range(1,7):
        for m in range(1,8):
            for n in range(1,8):
                pair=shat(m,j)*shat(n,jp)+shat(-m,j)*shat(-n,jp)
                maxim=max(maxim,abs(pair.imag))
print(f"    max |Im( shat_j(m)shat_j'(n) + shat_j(-m)shat_j'(-n) )| over j,j',1<=m,n<=7 = {maxim:.2e}")
print(f"    => {'REAL (joint (m,n)<->(-m,-n) pairing kills imaginary part; product QR-reality HOLDS)' if maxim<1e-12 else 'NOT real!'}")

# ---------- (3) interaction sup vs resonance distance ----------
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
def IB(B,u,v): return measS7(list(B)+[u,v])-measS7(list(B)+[u])-measS7(list(B)+[v])+measS7(list(B))
def p1p2(B):
    B=sorted(set(B)); bps=set([F(0),F(1)])
    for e in B:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); p1=F(0);p2=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        nmiss=7-len(set(sector_of(e*((x0+x1)/2)) for e in B))
        if nmiss==1:p1+=x1-x0
        elif nmiss==2:p2+=x1-x0
    return p1,p2
def resdist(u,v,H=7):
    return min(abs(m*u+n*v) for m in range(-H,H+1) for n in range(-H,H+1) if (m,n)!=(0,0))
print("\n(3) interaction |I_B - Phi_2| vs resonance distance (B=(0,4,6,8,10,12,14)):")
B=(0,4,6,8,10,12,14); p1,p2=p1p2(B); Phi2=(2*p2-p1)/49
print(f"    Phi_2={float(Phi2):.6f}")
print(f"    {'(u,v)':>12}{'resdist':>9}{'|IB-Phi2|':>12}{'*resdist':>10}")
for (u,v) in [(15,16),(15,23),(20,21),(31,47),(40,41),(61,97),(80,81),(101,211),(200,401)]:
    dev=abs(float(IB(B,u,v)-Phi2)); rd=resdist(u,v)
    print(f"    {f'({u},{v})':>12}{rd:>9}{dev:>12.6f}{dev*rd:>10.4f}")
print("\nIf |IB-Phi2|*resdist is bounded, the bound |I_B-Phi_2| <= C2'*V/resdist holds (resonance-gated).")
print("\nDONE.")
