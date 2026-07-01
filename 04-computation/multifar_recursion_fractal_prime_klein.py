#!/usr/bin/env python3
"""
multifar_recursion_fractal_prime_klein.py  --  klein-2026-07-01-S67

THE NEXT LEVER (r-far multi-patch correction), viewed RECURSIVELY / FRACTALLY, with 2 as the first prime.

Core = {1..n-2} (the small speeds that DEFINE the lonely set); L_C = {t: ||v t||>r for all v in core}.
Far speeds H={w_1..w_r} (large). Survival_r = meas(L_C ∩ ∩_i safe(w_i))/L. mac-mini-S74: survival ~
(1-2r')^r = (6/7)^r (independence). kind-pasteur-S3: survival_r = (6/7)^r - [signed resonance correction].
S66 (klein): single-far (r=1) correction = O(1/w), an EXACT Fourier sum on hat1 = FT of 1_{L_C}.

THIS SESSION: the r-far correction as a RECURSION (cumulant/Mobius expansion of pairwise 2-far), its
SELF-SIMILAR (fractal) structure from the continued fraction of t*, and the PRIME lens:
  * PRIME 2 (the first prime) = the SIGN. ||.|| has a ± (nearest integer both ways) => the danger Fourier
    coeff s(t)=sin(2πt r')... is ODD (ι: t->1-t, a->-a). The r-far correction is a product of r mean-zero
    ι-odd factors => its SIGN is (-1)^r (parity). This IS why the correction is "signed" (the (-1)^|T| in
    THM-501's L=(6/7)^13 + Σ(-1)^|T|(6/7)^{13-|T|}∏s). Prime 2 = parity = the sign grading.
  * PRIME 7 = n/2 (the OTHER prime of n=14): the 7-VANISHING (THM-503). s(t)=sin(πt/7)/(πt)=0 when 7|t.
    Sign period = 14 = 2·7. Prime 2 gives the sign, prime 7 gives the ZEROS. => hat1 suppressed at 7|k.
  * PRIMES 3,61 = Φ6(14)=183: the RESONANCE denominators. t*=n/Φ6=14/183, CF=[0;13,14], convergents
    1/13 (denom 13=n-1=the missing speed) and 14/183. Resonances at multiples of 13; fractal/Ostrowski.
So the same prime that is "first" (2 = the sign/±) has echoes: 7 = the zeros, 3·61 = the resonances --
each prime plays the analogous structural role one level finer. FRACTAL = self-similar combs at CF scales.
"""
import numpy as np
from fractions import Fraction as F

def norms(v, t):
    x = (v*t) % 1.0
    return np.minimum(x, 1-x)

def factor(m):
    m=int(abs(m)); fs=[]; d=2
    while d*d<=m:
        while m%d==0: fs.append(d); m//=d
        d+=1
    if m>1: fs.append(m)
    return fs

if __name__ == "__main__":
    n=14; core=list(range(1,n-1)); r=1.0/n; rp=1.0/n     # canonical: 2r'=1/7 (the prime 7)
    Phi6=n*n-n+1; tstar=n/Phi6                            # 183; t*=14/183=[0;13,14]
    N=600000; t=np.arange(N)/N
    G=np.full(N,1.0)
    for v in core: G=np.minimum(G,norms(v,t))
    f=(G>r).astype(np.float64); L=f.mean()
    d=np.diff((f>0).astype(int)); I=int((d==1).sum())+(1 if f[0]>0 else 0)
    Lmask=f>0
    print(f"n={n} core={core}  r=r'=1/{n}  2r'=1/7={2*rp:.4f}  Phi6={Phi6}  t*={tstar:.5f}=[0;13,14]")
    print(f"L=meas(L_C(core))={L:.5f}  #intervals I={I}  (Cantor set at the binding)")
    hat1=(np.fft.rfft(f)/N).real
    Kmax=N//2

    print("\n===== PART A: FRACTAL Fourier structure (self-similar combs at CF convergent denominators) =====")
    print("CF of t*=14/183 = [0;13,14]; convergent denominators q1=13 (=n-1, the MISSING speed), q2=183=Phi6.")
    print(" two-atom law hat1(k) ~ L cos(2π k t*); peaks where k t* ~ integer i.e. k ~ multiples of 1/t*=13.07")
    peaks=[13,26,39,52,65,183,196]  # comb-1 (denom 13) and near the second convergent 183
    for k in peaks:
        if k<len(hat1):
            print(f"   k={k:>4}: hat1={hat1[k]:+.5f}  hat1/L={hat1[k]/L:+.3f}  cos(2πk t*)={np.cos(2*np.pi*k*tstar):+.3f}")
    print(" SELF-SIMILARITY: the comb around k=0 (DC) repeats around each resonance k~13m. Ratio of combs:")
    for m in [1,2,3]:
        base=abs(hat1[13*m]); nbr=abs(hat1[13*m+13]) if 13*m+13<len(hat1) else 0
        print(f"   comb at 13*{m}={13*m}: |hat1|={base:.5f}, next rung 13*{m+1}: |hat1|={nbr:.5f}, ratio={nbr/base if base>0 else 0:.3f}")
    print(" PRIME-7 GATING (7-vanishing, THM-503: s(t)=0 at 7|t): |hat1(k)| at k = multiples of 7 vs neighbors:")
    for k in [7,14,21,28,49,61]:
        nb=(abs(hat1[k-1])+abs(hat1[k+1]))/2
        print(f"   k={k:>3} ({'7|k' if k%7==0 else 'prime '+str(k) if len(factor(k))==1 else 'x'.join(map(str,factor(k)))}): |hat1|={abs(hat1[k]):.5f}  nbr-avg={nb:.5f}  {'SUPPRESSED' if abs(hat1[k])<0.6*nb else ''}")

    print("\n===== PART B: MULTI-FAR survival vs independence; PRIME-2 lives in the ODD kernel (not survival sign) =====")
    def survival(H):
        m=Lmask.copy()
        for w in H: m=m&(norms(w,t)>rp)
        return m.mean()/L
    def block(W,rr): return list(range(W, W+rr))   # near-equal combs (kind-pasteur worst case)
    print(" H = consecutive far block {W..W+r-1} (near-equal combs).  W=1000:")
    print(f"   {'r':>2} {'survival':>9} {'(6/7)^r':>9} {'correction':>11} {'sign':>5}")
    W=1000
    for rr in range(1,6):
        H=block(W,rr); s=survival(H); indep=(1-2*rp)**rr; corr=s-indep
        print(f"   {rr:>2} {s:>9.5f} {indep:>9.5f} {corr:>+11.5f} {'+' if corr>0 else '-':>5}")
    print(" NOTE: survival-correction sign is NOT (-1)^r (data: +,+,+,+,-). The prime-2 (ι-odd) structure")
    print(" is in the DANGER KERNEL s(j)=sin(2πj r')/(πj) (ODD in j = ι-odd = the ± band = why it's SIGNED),")
    print(" NOT in the survival parity. hat1 (FT of 1_{L_C}) is real & EVEN (cos) -- the ι-EVEN partner.")

    print("\n===== PART C: PAIRWISE 2-far corr2(a,b): resonance = |a-b| or a+b hits the 13-lattice (fractal) =====")
    def corr2(a,b):
        ga=(norms(a,t)<rp).astype(float)-2*rp; gb=(norms(b,t)<rp).astype(float)-2*rp
        return (ga*gb*Lmask).mean()/L
    print(" corr2 for a=1000, b=1000+delta (near-equal combs): resonance vs the difference delta:")
    a=1000
    for delta in [1,2,3,6,7,12,13,14,26]:
        c=corr2(a,a+delta)
        fs=factor(delta)
        print(f"   delta={delta:>3} (={'x'.join(map(str,fs)) if fs else '1'}): corr2={c:+.5f}  {'RESONANT (delta ~ 13-lattice? '+str(delta%13==0)+')' if abs(c)>0.01 else ''}")
    print(" => corr2 peaks when the comb-DIFFERENCE delta resonates with the 13-lattice (=1/first-convergent).")
    print("    The prime structure of delta (2,3,7,13...) controls which combs co-resonate: PRIME-graded.")
    print(" SIGN CENSUS (scan delta=1..40): are the STRONG resonances POSITIVE (redundant, bad for beater)")
    print(" while anti-correlations (negative, coverage-SPREADING, needed to beat) stay WEAK?")
    scan=[(corr2(a,a+dd),dd) for dd in range(1,41)]
    pos=sorted(scan,reverse=True)[:4]; neg=sorted(scan)[:4]
    print(f"   most POSITIVE (redundant): {[(round(c,4),d,'13|d='+str(d%13==0)) for c,d in pos]}")
    print(f"   most NEGATIVE (spreading): {[(round(c,4),d) for c,d in neg]}")
    print("   => if strong = positive/13-lattice (redundant) and negatives are weak, far combs cannot")
    print("      SPREAD coverage over L_C (only pile up redundantly) => no multi-far beater (heuristic).")

    print("\n===== PART D: does the 13-lattice resonance SURVIVE to all scales? (the OPEN-Q-108 crux) =====")
    print(" THE decisive test: a 13-SPACED pair (W, W+13) resonates (PART C). Does corr2(W,W+13) DECAY in W,")
    print(" or persist at all scales? If it decays like single-far O(1/W), multi-far is ALSO impotent.")
    def corr1(w): return ( (norms(w,t)<rp)&Lmask ).mean()/L - 2*rp
    print(f"   {'W':>7} {'|corr1(W)|':>11} {'|corr2(W,W+1)|':>15} {'|corr2(W,W+13)|':>16} {'x |corr1|':>10}")
    for W in [200,1000,5000,20000,50000]:
        c1=abs(corr1(W)); c2nr=abs(corr2(W,W+1)); c2r=abs(corr2(W,W+13))
        print(f"   {W:>7} {c1:>11.5f} {c2nr:>15.5f} {c2r:>16.5f} {c2r/c1 if c1>0 else 0:>10.1f}")
    print(" KEY: if |corr2(W,W+13)| ALSO -> 0 as W->inf, the 13-resonance is a FINITE-SCALE effect: far-far")
    print(" combs cannot resonate at large scale => the multi-far correction vanishes => far element impotent")
    print(" at ALL orders r (recursive extension of S66). If it PERSISTS, that is the OPEN-Q-108 danger locus.")
    print(" (grid N=6e5 resolves W up to ~N/20; W=50000 is near the resolution edge -- read the trend.)")
