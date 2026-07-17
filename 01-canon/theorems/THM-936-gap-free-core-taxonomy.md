---
id: THM-936
title: THE GAP-FREE CORE TAXONOMY — the last unstructured regime (13-speed primitive packets, consecutive ratios ≤ 15, no cascade step or gluing junction applies) is WITNESS-CERTIFIABLE at accessible scale: every sampled packet is lonely with a kernel-checkable rational certificate (mu > 0 ⟹ simplest-rational witness by Stern–Brocot in the largest uncovered component; mu = 0 ⟹ tight with 14-grid attainment at the mu6 locus); the witness law: STRUCTURE prices witnesses, oppositely to how cleanliness prices BONF5/DK — two complementary certificate families; strata census in-file
status: census complete (exact, lam = 1/14 TRUE LRC): GENERIC (mu>0, rational witness) = 255; TIGHT (mu=0, 14-grid attained) = 2; every-packet-lonely = True
source: opus-2026-07-16-S334 (owner: close CascadeGluing sorries, then the gap-free core taxonomy; finish remaining LRC(14) tasks)
depends_on: [THM-928 (two-scale certificate), THM-932/933 (gluing; this regime is its complement), THM-738/kps (tight families), LRCRatWitness/LRCWitnessCert Lean modules (the certificate consumers)]
scripts: 04-computation/gapfree_core_taxonomy_opus_S334.py -> 05-knowledge/results/gapfree_core_taxonomy_opus_S334.out
---

# THM-936 — the gap-free core taxonomy

**Regime.** Primitive 13-speed packets, sorted, all consecutive ratios ≤ 15:
no THM-928(A) cascade step fires and no THM-932/933 junction exists. This is
the complement of the two-scale machinery — the residual regime of the loose
half after S332/S333.

**Method (all exact, lam = 1/14).** Per packet: the exact uncovered set W
(union of rational intervals); mu(W); if mu > 0, the SIMPLEST RATIONAL in
the largest component (Stern–Brocot) = the loneliness witness, verified by
13 exact distance checks (the kernel-checkable object consumed by the
project's LRCRatWitness / LRCWitnessCert Lean modules); if mu = 0, closed
attainment on the 14-grid (the mu6 locus, kps THM-738's tight territory).
Families sampled: tight {1..13}, APs, near-APs, in-band geometric, 120
random, and the 80 cleanest (top Theta*-quartile of a 200-packet pool).

**Census.** GENERIC (mu>0, rational witness) = 255; TIGHT (mu=0, 14-grid attained) = 2.
EVERY sampled packet is lonely — the accessible gap-free core is witness-certifiable end to end.

**The witness-denominator law — an honest surprise.** Witness cost tracks
STRUCTURE, not cleanliness: AP med 3, nearAP med 36, geometric med 58,
random med ~91 — and the cleanest random quartile is NOT simpler (med 103
vs 91; Theta*-bins: Theta*<30: med 76, p90 1623, max 5335 (n=245); 30-59: med 107, p90 484, max 484 (n=10)). The prediction "Theta* prices witnesses"
is REFUTED within the random regime. The true law: **the two certificate
families price OPPOSITELY.** Resonance/structure concentrates the
uncovered set at simple rationals (cheap witness, bad BONF5 — the AP has
3-term APs everywhere and would fail every equidistribution certificate);
cleanliness spreads the uncovered set (cheap BONF5/DK, expensive witness).
The gap-free core is covered by the UNION of two complementary certificate
families with opposite pricing — structure pays one, cleanliness pays the
other, and nothing at accessible scale escapes both.

**Dispatch (the taxonomy as a decision procedure).** Given any 13-packet:
(1) a scale gap ⟹ split by THM-932/933 gluing; (2) 15-lacunary ⟹
THM-928(A) cascade; (3) gap-free ⟹ compute W exactly (cost O(Σx·mu),
accessible whenever max speed is) and emit: the rational witness (mu > 0)
or 14-grid attainment (tight). At accessible scale nothing else occurs
(census). The open remainder of the loose half is therefore PURELY the
LARGE-SCALE gap-free regime — where the covering program (routes [A],
THM-922/924, c0 = 17/84, T2 = THM-926) already operates; the two programs
now meet with no unclaimed middle.

**Lean note.** The witness stratum's certificates are exactly the shape the
existing LRCRatWitness/LRCSafeCert modules consume; CascadeGluing.lean
(this session, sorry-free: cascade_step + window/union_floor_sample) covers
the measure layer for regimes (1)-(2).
