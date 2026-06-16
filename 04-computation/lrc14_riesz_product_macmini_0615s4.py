#!/usr/bin/env python3
"""
lrc14_riesz_product_macmini_0615s4.py  (mac-mini-2026-06-15-S4, THM-515 / HYP-2540)

APPLY THE RIESZ-PRODUCT METHOD (arXiv:2511.16636, "A wider gap of loneliness", 2025)
to the LRC(14) singular series. The lonely set {τ: M(τ)=0}, M(τ)=Σ_v 1_{||vτ||≤1/14}(τ),
has Lebesgue measure L(S) (THM-515). S is LOOSE iff L>0 iff the covering Φ=M fails.

RIESZ-PRODUCT CERTIFICATE OF LOOSENESS: build a probability density R(τ)=∏_j(1+a_j cos
2π m_j τ) ≥ 0 (dissociated frequencies m_j). Then ∫M·R = Σ_v Σ_k s(k) R̂(v k) (s=sinc
band coeff). The main term (k=0) is Σ_v s(0) = 13/7 ≈ 1.857 > 1; the SIGNED Riesz
coefficients R̂(vk)<0 on the relation frequencies must subtract ≥ 0.857 to push ∫MR < 1
= ∫R, which (M≥1 on the cover) CERTIFIES the lonely set has positive R-mass ⟹ LOOSE.

This script:
 (1) multiplicity moments ∫M, ∫M², ∫M³ (Tao's second/third-moment method) for the
     extremizer cores — confirm ∫M ≈ 13/7 and the AP-forcing structure;
 (2) the ADDITIVE-ENERGY predictor: L vs E(S)=#short relations (corrects λ_1);
 (3) a RIESZ-PRODUCT pairing ∫M·R for a hand-built R adapted to a core's relations —
     does the signed cancellation pull ∫MR toward / below 1 (proof-of-concept)?
"""
import sys, itertools, math
from math import gcd, sin, pi

sys.stdout.reconfigure(line_buffering=True)

def s(k):
    return sin(pi*k/7)/(pi*k) if k != 0 else 1/7

def lonely_measure(S, Q):
    rad = Q//14; cnt = 0
    for a in range(Q):
        if all(not ((v*a) % Q <= rad or (v*a) % Q >= Q-rad) for v in S): cnt += 1
    return cnt/Q

def mult_moments(S, Q):
    rad = Q//14
    m1=m2=m3=z=0
    for a in range(Q):
        M = sum(1 for v in S if ((v*a)%Q)<=rad or ((v*a)%Q)>=Q-rad)
        m1+=M; m2+=M*M; m3+=M*M*M
        if M==0: z+=1
    return m1/Q, m2/Q, m3/Q, z/Q

def short_rel_count(S, B):
    n=len(S); c=0
    for sz in range(2,5):
        for T in itertools.combinations(range(n),sz):
            for co in itertools.product([x for x in range(-B,B+1) if x],repeat=sz):
                if sum(x*S[i] for x,i in zip(co,T))==0:
                    g=0
                    for x in co: g=gcd(g,x)
                    if g==1: c+=1
    return c

def riesz_pairing(S, Q, freqs_amps):
    """∫ M(τ) R(τ) dτ with R(τ)=∏(1+a cos 2π m τ), evaluated on Z/Q grid (R normalized to ∫R=1).
    Returns (∫MR, ∫R) so we see if ∫MR < ∫R (certificate)."""
    rad=Q//14
    sMR=0.0; sR=0.0
    for a in range(Q):
        tau=a/Q
        R=1.0
        for (m,amp) in freqs_amps:
            R*= (1+amp*math.cos(2*pi*m*tau))
        M=sum(1 for v in S if ((v*a)%Q)<=rad or ((v*a)%Q)>=Q-rad)
        sMR+=M*R; sR+=R
    return sMR/Q, sR/Q

cores = {
  "extremizer {1..13}\\{6}∪{56}": sorted(set([x for x in range(1,14) if x!=6]+[56])),
  "extremizer {1..13}\\{12}∪{84}": sorted(set([x for x in range(1,14) if x!=12]+[84])),
  "evader 7*{1..12}∪{13}": sorted(set([7*i for i in range(1,13)]+[13])),
  "generic primitive {3,5,9,14,17,33,65,129,257,513,1025,2049,4097}": [3,5,9,14,17,33,65,129,257,513,1025,2049,4097],
}

print("="*74)
print("(1) MULTIPLICITY MOMENTS ∫M, ∫M², ∫M³  (Tao second/third-moment) + lonely measure")
print("="*74)
print(f"   main term ∫M should ≈ 13·(2⌊Q/14⌋+1)/Q ≈ 13/7 = {13/7:.4f}")
for name,S in cores.items():
    if len(S)!=13: print(f"   {name}: |S|={len(S)} skip"); continue
    m1,m2,m3,L = mult_moments(S, 14000)
    print(f"   {name[:40]:40s}: ∫M={m1:.4f} ∫M²={m2:.4f} ∫M³={m3:.4f}  L=μ{{M=0}}={L:.5f}")

print("\n" + "="*74)
print("(2) ADDITIVE-ENERGY PREDICTOR: L vs E(S)=#short 7-primitive relations (B=4)")
print("="*74)
rows=[]
for name,S in cores.items():
    if len(S)!=13: continue
    L=lonely_measure(S,11000); E=short_rel_count(S,4)
    rows.append((L,E,name)); print(f"   {name[:40]:40s}: L={L:.5f}  E(short rels,B=4)={E}")
rows.sort()
print(f"   --> ordered by L: {[(round(L,4),'E='+str(E)) for L,E,_ in rows]}  (low L ⟺ high E)")

print("\n" + "="*74)
print("(3) RIESZ-PRODUCT PAIRING ∫M·R  (certificate test: can signed R̂ pull ∫MR < 1?)")
print("="*74)
Sx = sorted(set([x for x in range(1,14) if x!=6]+[56]))   # the hardest extremizer
print(f"   target core: {Sx}, L≈0.0056 (the inf extremizer)")
print("   R=∏(1+a cos 2π m τ): tune frequencies m to the speeds (danger bands) + amplitudes:")
# baseline R=1 gives ∫M = 13/7
mr0, r0 = riesz_pairing(Sx, 7000, [])
print(f"     R=1 (no product): ∫MR={mr0:.4f}, ∫R={r0:.4f}  (= ∫M = 13/7 baseline)")
# try R with a few frequencies at the speeds themselves (each cos at freq v subtracts from danger_v overlap)
for amps in [[(v, -0.3) for v in (1,2,3,4,5)],
             [(v, -0.5) for v in Sx[:6]],
             [(v, -0.6) for v in (7,14,21)]]:
    mr, r = riesz_pairing(Sx, 7000, amps)
    print(f"     R with {len(amps)} factors (freqs {[m for m,_ in amps][:4]}...): "
          f"∫MR={mr:.4f}, ∫R={r:.4f}, ∫MR/∫R={mr/r:.4f}  {'< 1 CERT!' if mr/r<1 else '(≥1)'}")
print("   NOTE: a full looseness certificate needs an OPTIMIZED dissociated Riesz product")
print("   (the arXiv:2511.16636 construction); this is a feasibility probe of the pairing.")
print("\nDONE.")
