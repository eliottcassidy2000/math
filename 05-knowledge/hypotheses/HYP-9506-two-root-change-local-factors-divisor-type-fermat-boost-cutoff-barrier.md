---
id: HYP-9506
title: "Local factors of the two-root-change multiplicity: (A) the residue automaton of THM-4916 predicts divisor-type local factors L_p(0) = 1 + c_p/p with c_p -> about 1.5 and uniform nonzero classes, the mod-3 limit (34/27, 47/54, 47/54); (B) at census scale the Fermat-type multiplier primes 3, 5, 17 (and 11) carry the largest boosts, which decay like 1/log N; (C) the census factors are height-weighted and exceed the letter-count automaton at those primes; (D) near-misses of two root changes are rough (no factor 3, few small factors); (E) two root changes cover every odd N > 3,780,967, with no proof route yet (product/multiplier barrier)"
status: "(A) EXPLORATORY: the automaton is exact for its letter-count model; its limits are numerical (Richardson to 4 digits); the census agreement is in (C). (B) CONFIRMED by pre-registered R3 (rho_5 > rho_17 > rho_11 > the rest). (C) Pre-registered R2 FAILED: the census rho_p exceeds the automaton band at p = 5, 11, 13, 17; inside at 7, 23, 31. R4 (uniform nonzero classes) and R5 (27/27 classes mod 27) CONFIRMED. (D) FINITE-EXACT observation. (E) OPEN, FINITE-EXACT below 10^9; R1 (to 2^32) pending. AUDIT: pending."
session: kind-pasteur-2026-11-07-S15
date: 2026-11-07
related: [THM-4916, THM-4904, HYP-9494, HYP-9480, THM-4888, 04-computation/experiments/chanpak_mass_exploratory_20261106/local_factor_mod3_prediction.py, 05-knowledge/results/session223_root_changes_and_local_factors_20261107/README.md]
script: 05-knowledge/results/root_change_local_factors_20261107/resdp.c; 05-knowledge/results/root_change_local_factors_20261107/rc_lf.c; 05-knowledge/results/root_change_local_factors_20261107/eval_s15.py; 05-knowledge/results/root_change_local_factors_20261107/decomp.c
---

# HYP-9506 - local factors of the two-root-change multiplicity, the Fermat boost, and the cutoff barrier

Packet: [README](../results/root_change_local_factors_20261107/README.md),
[PROOF](../results/root_change_local_factors_20261107/PROOF.md),
[ATLAS](../results/root_change_local_factors_20261107/ATLAS.md). The exact identities are THM-4916.

## (A) The residue automaton

Count the prefixes with n letters O/Z and at most two J, with the M_2 weights, by class of the value mod q. This gives
local factors L_n(r).
- The nonzero classes become uniform.
- At n = 1000, L(0) is 1.2125, 1.1687, 1.1219, 1.1023, 1.0877, 1.0754, 1.0644, 1.0511, 1.0479 for
  p = 5, 7, 11, 13, 17, 19, 23, 29, 31. So (L(0) - 1) p tends to about 1.5, a divisor-type factor consistent with
  the product structure THM-4916 (1e).
- For p = 3, Richardson extrapolation of n = 4000-16000 gives L(0) = 1.2592 and L(1) = L(2) = 0.8703. The conjectured
  exact values are 34/27 and 47/54.

## (B) The Fermat boost

The cheap multiplier tails of THM-4916 (1c) (multipliers 3, 5, 9, 17, 33, ...) boost the primes dividing 2^i + 1 for
small i at finite scale. The census ratios rho_p = mean(p | N)/mean(p ∤ N) on [2^23, 2^24) are:

| p | 5 | 17 | 11 | 7 | 13 | 23 | 31 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| rho_p | 1.614 | 1.534 | 1.374 | 1.240 | 1.171 | 1.124 | 1.102 |

The ordering 5 > 17 > 11 > the rest is pre-registered R3, CONFIRMED.

## (C) Height weighting

The letter-count automaton underestimates rho_p at census scale for 5, 11, 13 and 17 (pre-registered R2 FAILED). It
matches 7, 23 and 31, the uniformity of the nonzero classes (R4), and all 27 classes mod 27 to within 5% (R5).
Conjecture: weighting words by value instead of letter count closes the gap. Session223 independently names a
"height-weighted 3-adic factor".

## (D) Rough near-misses

The 15 misses above 10^5 and the 12 values with M_2 <= 6 in (2^22, 2^28) are never divisible by 3, and are mostly 1
(mod 3). Many are prime (5 of 15, and 3 of 12). The low points of the local factor are the rough numbers.

## (E) The cutoff (OPEN)

Two root changes cover every odd N > 3,780,967. This is FINITE-EXACT to 10^9; R1 tests it to 2^32. No proof route is
known.
- Product and multiplier families cover only composites, so primes need the non-product part of the K-form
  (THM-4916 (1d)).
- Single seeds and finite banks fail.
- A proof needs one of:
  - a state-cone blow-up in the style of Gerbicz's square-sum induction, plus the census;
  - a second-moment (collision) bound for K(v,w) over two independent states.
