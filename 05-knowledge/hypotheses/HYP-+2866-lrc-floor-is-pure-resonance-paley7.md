---
id: HYP-+2866
title: The LRC(14) witness floor R' is EXACTLY 1 off-resonance (perfect decorrelation) -- the deviation is PURELY a finite low-height resonance phenomenon at the apex prime 7 (Paley P_7 / QR mod 7), real by HYP-2657's 6=-1 NQR; closes the floor to a finite QR/Gauss-sum resonance bound
status: STRUCTURE VERIFIED. R'=1 exact off-resonance; R' in [0.81,1.0] for LRC14 (resonant). The remaining floor = bounding the finite low-height resonances via Paley-P_7/Gauss-sum (kps Node-3 + my sqrt-cancellation + HYP-2657 reality).
source: mac-mini-2026-06-22-S31
related:
  - HYP-2657   # QR/Gauss-sum coset kernel D7: reality via 6=-1 NQR (Paley P_7 structure) PROVED
  - HYP-2863   # the floor rho* = R'*meas(GOOD)*meas(G_P); R' the quasi-independence
  - HYP-2840   # sqrt-cancellation (the resonance magnitude) = kps Node-3 spectrum
---

# HYP-+2866 -- the LRC floor is pure resonance at Paley P_7

## The finding (think tournaments: the apex prime 7 = Paley P_7 = QR mod 7)
The witness-floor quasi-independence R' = meas(GOOD cap G_P)/(meas(GOOD)*meas(G_P)) is
**EXACTLY 1 off-resonance**: for scale-separated/dilated clusters E = D*base (GOOD = cluster
maxgap event) and small part P, R' = 1.00000 (verified, grid 180180, D=7,14,100,101; common
factors 7,14,15 ALL give R'=1). So GOOD and G_P are PERFECTLY independent (disjoint Fourier
supports) unless their frequency lattices RESONATE.
- R' != 1 ONLY for OVERLAPPING/interleaved scales (E={0..8}, P={1,2,3}: R'=0.915). LRC(14) is
  in this regime (cluster co-offsets {Vmax-u: u>13} and small-part speeds {<=13} both small =>
  overlap), so R' in [0.81,1.0] -- the RESONANT case.

## Consequence: the floor reduces to a FINITE resonance bound
R' - 1 = SPEC/baseline, SPEC = Sum_{shared n} c^(n) g^(n). Off-resonance SPEC = 0 (exact). The
deviation comes ONLY from the finitely many LOW-HEIGHT shared frequencies (the overlap + the
apex prime 7). So the floor R' >= c reduces to BOUNDING these finite resonances, NOT a diffuse
decorrelation estimate. This is exactly the user's "route finitely many low-height resonances
to the resonant-nbhd patch" + kps's Node-3 spectrum-intersection sum.

## The tournament structure at 7 (Paley P_7 / QR / Gauss sum)
The dominant resonance is the apex prime 7 (seven sectors, gap 1/7). Its kernel is the D7 coset
kernel (HYP-2657), governed by QR mod 7 = the **Paley tournament P_7**:
- **Reality (PROVED, HYP-2657):** 6 = -1 is a NON-residue mod 7, so D7(-c) = conj(D7(c)) and
  the n<->-n pairing cancels ALL imaginary parts -- SPEC is REAL. (This is the QR(-1) dichotomy.)
- **Magnitude:** each resonance term is bounded by the Gauss sum |g_7| = sqrt(7) (the Paley
  eigenvalue scale); with the 1/n decay of g^ (G_P) the tail converges. The low-height
  resonances are a finite check.
- **q-uniform:** ports to LRC(2q) via the Paley P_q / Gauss-sum sqrt(q) at each apex prime q.

## Net (the floor-closure structure)
R' >= c > 0 via: [R' = 1 EXACT off-resonance, decorrelation] + [finite low-height resonances at
apex 7, REAL (HYP-2657 Paley/QR) + bounded by Gauss-sum sqrt(7) (my sqrt-cancellation magnitude
= kps Node-3 spectrum)] + [V<=b regime -> proven LRC(<=13)]. The remaining analytic content is
the finite apex-7 (+ common-factor) resonance magnitude -- a Paley-P_7/Gauss-sum bound, q-uniform.
-> HYP-2657, HYP-2863, HYP-2840, OPEN-Q-108.

Script: lrc14_floor_resonance_paley7_macmini_S31.py (to save).


## ⚠️ CORRECTION (mac-mini-S31, same session): the floor resonance is LOW-FREQUENCY, NOT apex-7
FFT decomposition of the floor covariance SPEC = Sum_n c^(n)conj(g^(n)) for the resonant config
(E={0..8}, P={1,2,3}, R'=0.915) shows the DOMINANT terms are at LOW frequencies n=+-3,+-6,+-4,+-12
(from the small part P={1,2,3} + the cluster), with only **9.8% of |covariance| at multiples of 7**.
So the apex-7/Paley-P_7 framing above is WRONG -- the floor deviation is NOT a 7-resonance; it is a
LOW-FREQUENCY overlap between GOOD (cluster maxgap) and G_P (small part). HYP-2657's apex-7 QR/Gauss
reality is about the K(n) sector-correction KERNEL, a DIFFERENT object, not the floor covariance.

## What's VERIFIED (the honest, clean findings)
1. **R' = 1 EXACTLY off-resonance** (scale-separated/dilated clusters: D=7,14,100,101 + common
   factors all give R'=1.00000). GOOD and G_P perfectly decorrelate when their frequency lattices
   don't overlap.
2. **R' in [0.81,1.0] for LRC(14)** (resonant: cluster co-offsets + small-part speeds overlap in
   low frequencies). The floor rho* = R'*meas(GOOD)*meas(G_P) > 0 since R'>=c~0.8.
3. **The floor deviation is FINITE LOW-FREQUENCY** (small n, the small-part + cluster overlap) =
   exactly the user's "finitely many low-height resonances." Bound via the spectrum-intersection
   sum Sum_n |c^(n)||g^(n)| with g^'s 1/n decay (kps Node-3). The reality of SPEC (it is real,
   no imaginary part) still holds via the n<->-n pairing.
4. **The TOURNAMENT connection that DOES hold:** the small part P (<=13 speeds) routes to PROVEN
   LRC(<=13) (giving meas(G_P)>0) -- the Rédei-adjacent proven-case routing (user's guidance).

Lesson: verify before claiming a QR/Paley specialization (3rd such over-claim this arc).
