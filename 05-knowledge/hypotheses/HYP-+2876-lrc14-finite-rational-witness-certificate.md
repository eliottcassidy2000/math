---
id: HYP-+2876
title: LRC(14) finite rational-witness residue check lead -- character-count and apex fragments survive, fixed finite denominator closure is refuted
status: REFUTED as a global finite-certificate claim by THM-566/HYP-2876; retained as a sampled residue-atlas lead
source: mac-mini-2026-06-22-S34 (user inspo: finite certificate basis + character-sum + apex-7 fragment)
related:
  - HYP-2869   # the assembled proof (rho*>0 + finite V*); this is the SHARPER rational-witness form
  - HYP-2860   # kps Node-3 spectrum = the resonances {sum k_s s = 0 mod D}
  - HYP-2856   # the floor 3/pi^2 = the continuous main-term analog
  - THM-565    # Node-1 V* atlas (V*<=234); this sharpens to witness D<=41
  - HYP-2876   # correction: finite denominator bases are atlases, not closures
---

# HYP-+2876 -- finite rational-witness certificate lead

## Correction after codex S98

The global finite-certificate reading below is false.  HYP-2876 proves the
stronger finite-basis obstruction:

```text
for any finite denominator list B,
S_B = {1,2,...,11,13,84*lcm(B)}
```

is primitive and covering, but has `N(S_B,D)=0` for every `D in B`, because
the tail speed is divisible by each such `D` and sits at the observer for
every numerator.  In particular, no fixed `D<=41` certificate and no fixed
finite basis such as `{83,89,21}` can prove LRC(14) globally.

What survives from this note is valuable but local:

- the exact count `N(S,D)` and its character/resonance expansion;
- the apex obstruction, strengthened in HYP-2876 from `D=14` to every reduced
  denominator `2..14` for covering rows;
- the idea of finite residue bases as adaptive or scaled atlases for sampled
  and non-divisor-loaded families.

The text below is retained as a historical sampled lead, not as a current
theorem-level claim.

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


## SIEVE REFINEMENT (mac-mini-S34): tight in [15,41], COMFORTABLE at large D
The proof "exists D with N(S,D)>=1" is RAZOR-THIN in the small range [15,41]: consec/loosest/dense-AP
have EXACTLY 1 good D of 27 (resonances kill 26/27). 0/2000 random sets are bad for ALL [15,41] (max
20 bad). So the small-D certificate is tight. BUT at LARGE D the main term (6/7)^13 phi(D) ~ 11
(D=83,89) gives a COMFORTABLE margin: N(S,83)>=1 unless the resonance deficit exceeds 11. So the
proof route is the LARGE-D basis {83,89}: bound the resonance deficit (Node-3 spectrum {sum k_s s=0
mod 83}) < 11 = the main term. For a prime D=83, resonances are sparse => deficit small => N>=1
(comfortable). The basis includes 89 (+21) so if S resonates with 83, use 89 -- no covering S resonates
with all. This is the Node-3 spectrum bound at the finite basis: deficit(S,D) < (6/7)^13 phi(D) for some basis D.
