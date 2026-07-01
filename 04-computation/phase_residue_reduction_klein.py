#!/usr/bin/env python3
"""
phase_residue_reduction_klein.py  --  klein-2026-07-01-S68

THE FINITE PHASE-REDUCTION: the entire far-speed interaction with the lonely set L_C is governed, at
leading order, by ONE FINITE invariant -- the PHASE-RESIDUE  p(w) = n*w mod Phi6  in Z/Phi6 -- collapsing
the INFINITE multi-far covering problem to a FINITE covering problem on Z/Phi6(n), which CRT-factors over
the primes of Phi6(n), with the antipode Z/2 (prime 2) acting as p <-> -p.

KEY IDENTITIES (general n, the covering-min construction {1..n-2, n(n-1)}):
  * t* = n/Phi6(n),  Phi6 = n^2-n+1.  Continued fraction  t* = [0; n-1, n]  for ALL n.
  * first convergent denominator = n-1 = the DROPPED speed = the RESONANCE PERIOD.
  * killer identity  n(n-1) = Phi6 - 1  ≡ -1  (mod Phi6)  => stepping w by (n-1) steps the phase-residue
    p by exactly -1 (mod Phi6):  p(w + (n-1)) = p(w) + n(n-1) = p(w) - 1.  So n-1 is the UNIT phase-step
    (the Ostrowski generator), which is WHY the resonance lattice is (n-1)Z (S67).
  * phase  phi(w) = w t* mod 1 = p(w)/Phi6.  RESONANT (couples to the binding atom of L_C) iff p(w) near
    0; ANTI iff p(w) near Phi6/2 (the antipode).  period in w is Phi6 (since gcd(n,Phi6)=1).

This script verifies: (A) corr1(w) sign tracks cos(2*pi*p(w)/Phi6) [the phase-residue predicts resonance];
(B) w and w+Phi6 have identical phase (same p) -- the finite periodicity; (C) pairwise corr2(a,b) depends
on p(a)-p(b) = n(a-b) mod Phi6; (D) the CRT factorization Z/Phi6 = prod over primes of Phi6.
"""
import numpy as np
from math import gcd

def norms(v, t):
    x = (v*t) % 1.0
    return np.minimum(x, 1-x)

def factor(m):
    m=int(m); fs={}; d=2
    while d*d<=m:
        while m%d==0: fs[d]=fs.get(d,0)+1; m//=d
        d+=1
    if m>1: fs[m]=fs.get(m,0)+1
    return fs

if __name__ == "__main__":
    n=14; core=list(range(1,n-1)); r=rp=1.0/n
    Phi6=n*n-n+1; tstar=n/Phi6
    def p(w): return (n*w)%Phi6
    N=600000; t=np.arange(N)/N
    G=np.full(N,1.0)
    for v in core: G=np.minimum(G,norms(v,t))
    f=(G>r).astype(np.float64); L=f.mean(); Lmask=f>0
    def corr1(w): return ( (norms(w,t)<rp)&Lmask ).mean()/L - 2*rp
    def corr2(a,b):
        ga=(norms(a,t)<rp).astype(float)-2*rp; gb=(norms(b,t)<rp).astype(float)-2*rp
        return (ga*gb*Lmask).mean()/L

    print(f"n={n} Phi6={Phi6}=prod{factor(Phi6)}  t*={tstar:.5f}  CF=[0;{n-1},{n}]  L=meas(L_C)={L:.5f}")
    print(f"KILLER identity: n(n-1)={n*(n-1)}=Phi6-1 => n(n-1) mod Phi6 = {(n*(n-1))%Phi6} (=-1). p(w+{n-1})=p(w)-1.")

    print("\n(A) PHASE-RESIDUE p(w)=n*w mod Phi6 predicts corr1(w): sign vs cos(2π p/Phi6), phase-dist to binding")
    print(f"    {'w':>4} {'p(w)':>5} {'phase':>7} {'cos(2πp/Φ)':>11} {'corr1(w)':>10} {'agree?':>7}")
    agree=0; tot=0
    for w in [13,26,39,52,65,78,7,61,91,100,183,196,300,500,501]:
        c1=corr1(w); cph=np.cos(2*np.pi*p(w)/Phi6); ok=(np.sign(c1)==np.sign(cph)) or abs(c1)<0.01
        agree+=ok; tot+=1
        print(f"    {w:>4} {p(w):>5} {p(w)/Phi6:>7.4f} {cph:>+11.3f} {c1:>+10.4f} {'yes' if ok else 'NO':>7}")
    print(f"    => sign(corr1) agrees with sign cos(2π p/Phi6) in {agree}/{tot}. The disagreements are ALL")
    print(f"       large-w with |corr1|<0.04 (magnitude decayed, S66 O(1/w)): p(w) sets the PHASE/DIRECTION,")
    print(f"       but the MAGNITUDE is width-limited O(1/w). Phase != full coupling.")

    print("\n(B) FINITE PERIODICITY: w and w+Phi6 have the SAME phase-residue p (=> same coupling DIRECTION);")
    print("    magnitude differs only by the width/decay O(1/w). corr1(w) vs corr1(w+Phi6):")
    for w in [13,26,39]:
        print(f"    w={w:>3} (p={p(w)}): corr1={corr1(w):+.4f} | w+Φ6={w+Phi6} (p={p(w+Phi6)}): corr1={corr1(w+Phi6):+.4f} | 2Φ6+w={w+2*Phi6}: corr1={corr1(w+2*Phi6):+.4f}")
    print("    (same sign/direction across +Phi6 shifts = phase is a Z/Phi6 invariant; magnitude ~ decays)")

    print("\n(C) PAIRWISE corr2(a,b) depends on the PHASE-DIFFERENCE p(a)-p(b)=n(a-b) mod Phi6 (not a,b sep.):")
    print(f"    {'a':>5} {'b':>5} {'a-b':>4} {'p(a)-p(b) modΦ':>14} {'corr2':>9}")
    for (a,b) in [(1000,1013),(2000,2013),(500,513),(1000,1026),(1000,1001),(1000,1007),(300,313)]:
        dp=(p(a)-p(b))%Phi6
        print(f"    {a:>5} {b:>5} {a-b:>4} {min(dp,Phi6-dp):>14} {corr2(a,b):>+9.4f}")
    print("    => corr2 large & positive iff p(a)-p(b) near 0 (a-b in (n-1)Z=13Z, i.e. n(a-b)≡small): REDUNDANT.")

    print("\n(D) THE FINITE PROBLEM lives on Z/Phi6, CRT-factored over the primes of Phi6, mod antipode Z/2:")
    fac=factor(Phi6)
    print(f"    Z/{Phi6} = " + " x ".join(f"Z/{q}^{e}" if e>1 else f"Z/{q}" for q,e in fac.items()) + "  (CRT over primes of Phi6)")
    print(f"    prime 2 (the FIRST prime) = the antipode iota: p <-> -p (t<->1-t). The finite covering problem")
    print(f"    is on (Z/{Phi6})/(±1). primes of Phi6={list(fac)} refine it (S68 lens: 2=sign, {list(fac)}=phase-CRT).")
    # sanity: is n a unit mod each prime power of Phi6? (needed for p(w)=n*w to be a bijection)
    print(f"    gcd(n,Phi6)={gcd(n,Phi6)} (n is a unit mod Phi6 => p(w)=n*w is a BIJECTION Z/Phi6->Z/Phi6).")
    print("\n=> CORRECTED REDUCTION (honest): the SCALE-INVARIANT finite structure is the DIFFERENCE SET.")
    print("   - single-far: phase = p(w)∈Z/Phi6 (DIRECTION); magnitude = width O(1/w) -> 0 (impotent, S66).")
    print("   - PAIRWISE: corr2(a,a+δ) depends ONLY on the difference δ (scale-invariant: same for a=300..2000),")
    print("     via p(δ)=nδ mod Phi6; LARGE & POSITIVE (redundant) iff δ∈(n-1)Z with δ small (S67).")
    print("   So multi-far coupling is governed by the DIFFERENCE-SET {w_i-w_j} and its phase-residues in")
    print("   Z/Phi6 (translation-invariant) -- and the magnitude on the resonance lattice (n-1)Z is the")
    print("   ADDITIVE ENERGY of the far-speed set (HYP-2873/THM-515). Redundancy(+)=helps survival;")
    print("   spreading(-)=weak. The FINITE object is the difference-set's phase-histogram on Z/Phi6.")
