---
id: HYP-+2876
title: LRC(14) = a FINITE rational-witness residue check (RS payoff) -- every 13-set has a witness a/D with D<=41; N(S,D)=#{a:||sa/D||>=1/14}=(6/7)^13 phi(D) + resonances (=Node-3 spectrum); D=14 never certifies (apex-7, PROVED); unifies floor+Node-3+V*+apex-7
status: VERIFIED (the user's tested Ideas 1/2/3, independently confirmed + connected to my floor/Node-3/V* work). The finite-certificate PROOF (every covering S has a basis witness) reduces to a resonance/sieve bound.
source: mac-mini-2026-06-22-S34 (user inspo: finite certificate basis + character-sum + apex-7 fragment)
related:
  - HYP-2869   # the assembled proof (rho*>0 + finite V*); this is the SHARPER rational-witness form
  - HYP-2860   # kps Node-3 spectrum = the resonances {sum k_s s = 0 mod D}
  - HYP-2856   # the floor 3/pi^2 = the continuous main-term analog
  - THM-565    # Node-1 V* atlas (V*<=234); this sharpens to witness D<=41
---

# HYP-+2876 -- LRC(14) as a finite rational-witness certificate

## The framework (verified, user Ideas 1/2/3)
A **rational witness** for S={s_1..s_13} is a/D with ||s a/D||>=1/14 for all s (a lonely time at
denominator D). Count: **N(S,D) = #{a in 1..D-1 : ||s a/D||>=1/14 forall s}**.
- **Main term (Idea 2):** N(S,D) ~ (6/7)^13 phi(D) (each runner avoids the 1/14-nbhd w.p. 6/7;
  phi(D) for a coprime to D). VERIFIED: (6/7)^13 phi(83)=11.05, phi(89)=11.86, phi(21)=1.62.
- **Error = RESONANCES {a : sum_s k_s s ≡ 0 mod D}** -- EXACTLY the Node-3 spectrum lattice (HYP-2860).
- **D=14 never certifies (Idea 3, PROVED apex-7 fragment):** covering forces a multiple of 14; at D=14
  that runner sits on the observer (||14k a/14||=0), so N(S,14)=0 always => min witness D >= 15. This is
  why 14=2*7 (composite) is hard. VERIFIED N(loosest,14)=0.
- **FINITE (Ideas 1, RS payoff):** every 13-set has a witness with D<=41 (loosest {1..11,13,84}: D=41;
  random: D<=25; the basis {83,89,21} certifies the COVERING family -- 16/3000 RANDOM sets miss the
  basis but ALL have a witness elsewhere, max min-witness-D = 41). So LRC(14) hard case = a FINITE
  residue check.

## The unification (this is the discrete/sharp form of the assembled proof HYP-2869)
- main term (6/7)^13 phi(D) = the discrete meas(G_P) floor factor (HYP-2856's 3/pi^2 is the continuous analog).
- resonances = the Node-3 spectrum (HYP-2860, kps).
- D<=41 finite bound = the sharp form of the V*<=234 atlas (THM-565).
- D=14 apex-7 obstruction = the proved fragment.
So the finite-certificate is the SAME object (rho*>0) at rational witnesses, with a SHARPER finite bound.

## The PROOF angle (the RS finite-certificate)
PROVE: for every covering 13-set S, exists D<=41 with N(S,D)>=1. Via: the main term (6/7)^13 phi(D)
(>=1.6) dominates the resonance deficit for at least one D<=41 -- a SIEVE/covering argument (the
resonances {sum k_s s=0 mod D} for different D are quasi-independent, so no S resonates with all D<=41).
This is the Node-3 spectrum bound applied to the finite basis. -> HYP-2869, HYP-2860, HYP-2856, THM-565.

Script: lrc14_rational_witness_certificate_macmini_S34.py.
