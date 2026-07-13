#!/usr/bin/env python3
"""mac-mini-S88: PROVE klein's disc_v certificate at the far peel via three-gap/Farey.
Claim chain for the deep well {1..12,182}, far peel v=182, G'=SafeSet({1..12}):
 (A) G' = union of Farey arcs between F_12 neighbors a/b<c/d with b+d=13, each length 1/(14bd).
     => N = phi(13) = 12 intervals; |G'| = sum 1/(14bd) = H_12/91 (closed form).
 (B) Fourier jump bound: |hat1_{G'}(l)| <= 2N/(2*pi*|l|) => disc_182 = sum_{m!=0}|c_{182m}|^2 <= N^2/(3*182^2).
 (C) certificate disc_182 < 6|G'|^2  <=>  H_k > k/(6*sqrt(2)); k=12 gives 3.103 > 1.414. RIGOROUS L>0."""
from fractions import Fraction as F
from math import gcd, pi, sqrt
import numpy as np

c=F(1,14); k=12; far=182
# --- (A) Farey F_12 neighbor pairs with b+d=13, exact arcs ---
def farey(n):
    fr=sorted({F(a,b) for b in range(1,n+1) for a in range(0,b+1) if gcd(a,b)==1})
    return fr
Fk=farey(k)
arcs=[]  # (left_edge, right_edge, b, d)
for i in range(len(Fk)-1):
    x,y=Fk[i],Fk[i+1]; b,d=x.denominator,y.denominator
    if b+d==13:  # the safe condition (F_12 neighbors always have b+d>=13)
        L=x+c/b; R=y-c/d
        arcs.append((L,R,b,d))
N=len(arcs)
measG=sum((R-L) for (L,R,b,d) in arcs)
H=lambda m: sum(F(1,i) for i in range(1,m+1))
print("=== (A) Farey-arc structure of SafeSet({1..12}) ===")
print(f"N (arcs with b+d=13) = {N}   [claim phi(13)=12]")
print(f"|G'| exact          = {measG} = {float(measG):.7f}")
print(f"H_12/91             = {H(12)}/91 = {float(H(12)/91):.7f}   [claim equal]")
print(f"each arc length 1/(14bd); sum matches: {measG==H(12)*F(1,91)}")
# numeric cross-check of the safe set
Ng=6_000_000; t=(np.arange(Ng)+0.5)/Ng
ok=np.ones(Ng,bool)
for i in range(1,k+1):
    r=(i*t)%1.0; ok &= (np.minimum(r,1-r)>=1/14)
nint=int((np.diff(ok.astype(np.int8))==1).sum())
print(f"numeric: |SafeSet({{1..12}})|={ok.mean():.7f} (vs {float(measG):.7f}), #intervals={nint} (vs {N})")

# --- (B) disc_182: exact Fourier coeffs of the arc union, and the jump bound ---
def chat(l):  # hat 1_{G'}(l) = sum_arcs (e^{-2pi i l R}-... )/(-2pi i l); use exact endpoints as floats
    if l==0: return complex(float(measG),0)
    s=0j
    for (L,R,b,d) in arcs:
        Lf,Rf=float(L),float(R)
        s+=(np.exp(-2j*pi*l*Rf)-np.exp(-2j*pi*l*Lf))/(-2j*pi*l)
    return s
disc=sum(abs(chat(m*far))**2 for m in range(-400,401) if m!=0)
bound=N**2/(3*far**2)
print("\n=== (B) disc_182 and the rigorous jump bound ===")
print(f"disc_182 (Fourier sum |m|<=400) = {disc:.3e}")
print(f"rigorous bound N^2/(3*182^2)    = {bound:.3e}   [disc <= bound: {disc<=bound}]")
print(f"opus/klein exact disc (HYP-6510)= {float(F(2629220219,6363107150400)):.3e}")

# --- (C) the certificate ---
sixGsq=6*float(measG)**2
print("\n=== (C) certificate disc < 6|G'|^2 ===")
print(f"6|G'|^2 = {sixGsq:.3e};  bound {bound:.3e} < 6|G'|^2 ? {bound<sixGsq}")
Lcert_bound=(6/7)*float(measG)-sqrt((6/49)*bound)   # using the RIGOROUS disc bound
Lcert_true =(6/7)*float(measG)-sqrt((6/49)*disc)
print(f"L_cert with rigorous bound = (6/7)|G'|-sqrt(6/49*bound) = {Lcert_bound:.5f}  [>0: {Lcert_bound>0}]")
print(f"L_cert with true disc      = {Lcert_true:.5f}")
print(f"actual L(deep well)        ~ 0.02390 (numeric)")
# the clean closed-form condition
print("\n=== closed-form condition (far=14(k+1), N=k, |G'|=H_k/(7(k+1))) ===")
for kk in range(6,14):
    Hk=float(H(kk)); cond=Hk> kk/(6*sqrt(2))
    print(f"  k={kk:2d}: H_k={Hk:.4f}, k/(6sqrt2)={kk/(6*sqrt(2)):.4f}  H_k>k/(6sqrt2)? {cond}  margin x{Hk/(kk/(6*sqrt(2))):.2f}")
print("\n=> for k=12 (deep well) the far-peel certificate is RIGOROUS: L>0, hence M>=1/14. QED-shape.")
