---
id: HYP-9506
title: "Local factors of the two-root-change multiplicity beyond THM-4916: (A) the census (height-weighted) factors tend to THM-4916's letter-count closed forms L_p(0) = 1 + (3p-2)(p-1)/(2p^3); (B) at census scale the Fermat-type primes 3, 5, 17 (then 11) carry the largest boosts, which decay with N; (C) the census is more contrasted than the letter-count automaton at every class; (D) near-misses avoid 3 and are depleted in 5, 7, 11, 17, 19; (E) two root changes cover every odd N > 3,780,967, with no proof route yet"
status: "(A) OPEN. The letter-count limits are PROVED in THM-4916; the census factor at N = 0 (mod 3) is 1.444 at [2^23, 2^24), against 34/27 = 1.259 in the limit. (B) The ordering is CONFIRMED at one scale (pre-registered R3). The decay rate is untested: rho_5 is 1.667 on [2^19, 2^20) and 1.614 on [2^23, 2^24). The attribution to the multiplier family THM-4916 (1c) is partial: about a quarter of the excess, with root-rectangle words another fifth and seed confinement about half (audit M3). (C) Pre-registered R2 FAILED: the census exceeds the automaton band at p = 5, 11, 13, 17 and is inside at 7, 23, 31. R5 passes only with its 5% widening (9/27 classes inside the raw band). R4 CONFIRMED. (D) FINITE-EXACT observation. (E) OPEN; FINITE-EXACT below 10^9, with R1 (to 2^32) recorded in the packet. AUDIT: independent, PASS WITH CORRECTIONS, MISTAKE-734."
session: kind-pasteur-2026-11-07-S15
date: 2026-11-07
related: [THM-4916, THM-4904, HYP-9494, HYP-9480, THM-4888, 05-knowledge/results/session223_root_changes_and_local_factors_20261107/README.md, 04-computation/experiments/chanpak_mass_exploratory_20261106/local_factor_mod3_prediction.py]
script: 05-knowledge/results/root_change_local_factors_20261107/rc_lf.c; 05-knowledge/results/root_change_local_factors_20261107/eval_s15.py; 05-knowledge/results/root_change_local_factors_20261107/decomp.c; 05-knowledge/results/root_change_local_factors_20261107/audit/census_mechanism.c
---

# HYP-9506 - local factors at census scale, the Fermat boost, and the cutoff

Packet: [README](../results/root_change_local_factors_20261107/README.md),
[PROOF](../results/root_change_local_factors_20261107/PROOF.md),
[ATLAS](../results/root_change_local_factors_20261107/ATLAS.md),
[AUDIT](../results/root_change_local_factors_20261107/audit/AUDIT.md). The proved limits and identities are THM-4916.

## (A) The census factors tend to the letter-count limits (OPEN)

THM-4916 proves, for the letter-count model, L_p(0) = 1 + (3p-2)(p-1)/(2p^3) with uniform nonzero classes. For
p = 3 this is 34/27, which Session223 also proved at every 3-adic precision.

The census counts words by value, not by letter count. At [2^23, 2^24) the census factor at N = 0 (mod 3) is 1.444.
Whether the height-weighted factors converge to the same limits is open. Session223 section 2 shows that residues
correlate with real-state history, which is a reason for doubt.

## (B) The Fermat boost

The cheap multiplier tails of THM-4916 (1c) (multipliers 3, 5, 9, 17, 33, ...) boost the primes dividing 2^i + 1 for
small i at finite scale. They carry about a quarter of the zero-class excess at 3, 5 and 17, and an eighth at 11 and
13. The rest is root-rectangle words and seed confinement, and the R3 ordering survives the removal of both
Fermat-type families. The census ratios rho_p = mean(p | N)/mean(p ∤ N) on [2^23, 2^24) are:

| p | 5 | 17 | 11 | 7 | 13 | 23 | 31 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| rho_p | 1.614 | 1.534 | 1.374 | 1.240 | 1.171 | 1.124 | 1.102 |

The Fermat primes are exactly the primes that are themselves multipliers 2^i + 1. Every prime with even ord_p(2)
divides some multiplier.

## (C) Height weighting

The letter-count automaton underestimates rho_p at census scale (R2 FAILED at 5, 11, 13, 17). At every class mod 27 the
census is more contrasted than the automaton (R5 passes only with its widening). Conjecture: weighting words by value
instead of by letter count closes the gap.

## (D) Near-misses

The 130 values with M_2 <= 12 in (2^22, 2^28) avoid 3 and are depleted in 5, 7, 11, 17 and 19, but not in 13. The 15
misses above 10^5 are 14 times 1 (mod 3), never 0, five prime, and six have a prime factor of at most 19.

## (E) The cutoff (OPEN)

Two root changes cover every odd N > 3,780,967 (FINITE-EXACT to 10^9; R1 is recorded in the packet).

A single O/Z seed misses a positive proportion (Session206). For finite O/Z banks this is OPEN (HYP-9480 (C)).

Product and multiplier families cover only composites, so every representation of a prime has a nonzero non-product
part (THM-4916 (1d)). A proof needs one of:
- a state-cone blow-up in the style of Gerbicz's square-sum induction, plus the census;
- a second-moment (collision) bound for the bilinear form.
