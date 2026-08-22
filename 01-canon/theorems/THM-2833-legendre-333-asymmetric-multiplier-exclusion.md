---
id: THM-2833
title: "Legendre pairs of length 333: asymmetric multiplier exclusion (orders >= 7, <= 25 orbits)"
status: >
  FINITE-EXACT with positive control.  No Legendre pair of length 333 exists
  in which the two sequences are separately invariant under (possibly
  DIFFERENT) multiplier subgroups H_A, H_B of (Z/333)^x of order >= 7 having
  at most 25 orbits each.  Extends the common-subgroup obstruction of
  arXiv:2607.20765 (which leaves the asymmetric case explicitly open) to
  two-sided structured pairs.  The unrestricted Legendre-pair existence
  problem remains OPEN.  Hadamard order 668 is now PROVED independently by
  THM-3393 through a length-166 four-sequence construction; no converse from
  that matrix to a length-333 Legendre pair is available.
source: mac-mini-2026-07-28-S171 (external open-problem raid; Epoch
  FrontierMath "Hadamard Matrices", order 668 = 2*333 + 2)
depends_on: []
related:
  - THM-3393-hadamard-order-668-explicit-certificate
script: 04-computation/legendre333_asym_multiplier_mim_macmini_S171.py
output: 05-knowledge/results/legendre333_asym_multiplier_macmini_S171.out
script_sha256: 0389a1cb412faee47d5c2bee574f92e6095649d7893378b0b3f36c8dfc37abce
output_sha256: bd22115d574593e67edabce3b99dd639c664eab4f26878d923d4337344fe45af
hash_basis: LF-normalized bytes
---

# THM-2833 — no asymmetric multiplier-invariant Legendre pair at length 333

## Statement

Let `L = 333`.  There is **no** pair of `±1` sequences `A, B` on `Z_L` with

* `PAF_A(s) + PAF_B(s) = -2` for all `s != 0` (Legendre pair condition), and
* `A` invariant under a multiplier subgroup `H_A <= (Z/333)^x`, `|H_A| >= 7`,
  with at most 25 orbits on `Z_333`, and likewise `B` under some `H_B`
  (independent of `H_A`).

The July 2026 obstruction paper (arXiv:2607.20765) covers only `H_A = H_B`
(common subgroup) and proves `|H| <= 6` there; pairs with different
invariances and small intersection were open.  All such pairs in the stated
enumerable range are now excluded.

## Method

Row sums satisfy `a^2 + b^2 = 2` (from summing the PAF condition), so each
side has row sum `±1`; many subgroups admit no invariant sequence at all
(orbit sizes divisible by 3 force row sum ≡ 0 mod 3).  For each of the ~30
subgroups of `(Z/333)^x ≅ C6 x C36` with `|H| >= 7` and `<= 25` orbits, all
invariant row-sum-`±1` sequences were enumerated (up to `2^25`), exact PAF
vectors computed, and all subgroup pairs hash-joined on
`PAF_A = -2 - PAF_B` (full 332-coordinate match).  Hits are re-verified by
direct `O(L^2)` integer autocorrelation.  Result: **0 pairs**.

Positive control: at `L = 31` the same pipeline (with common subgroups
allowed) finds and exactly verifies the classical Legendre-symbol pair,
invariant under the order-15 QR subgroup.

## Boundary / loss ledger

* NOT covered: subgroups of order >= 7 with more than 25 orbits acting on one
  side only (e.g. |H| = 9 or 12 one-sided with an unstructured partner); the
  fully unstructured problem; affine (translation-composed) symmetries;
  two-block or Yang-multiplication structures.
* Historically this left Hadamard order 668 untouched.  THM-3393 now proves
  that order by a bordered four-sequence construction of length 166.  The
  present exclusion remains valid and the unrestricted length-333 Legendre
  pair problem remains open; the two statements are no longer to be conflated.

## One-sided boundary (same session)

For A invariant and B fully UNSTRUCTURED, the necessary condition is
`PSD_A(k) <= 2L+2 = 668` for all `k != 0` (since `PSD_B = 668 - PSD_A >= 0`).
A vectorized scan over every enumerable invariant side
(`04-computation/legendre333_onesided_psd_scan_macmini_S171.py`, output
`05-knowledge/results/legendre333_onesided_psd_scan_macmini_S171.out`,
sha256 6385e0a888b3574c0f95ed816504c5d596e6624332afaac54549c4a55f7f1956) gives:

* every side of order >= 24 with any row-sum-valid classes is FULLY
  PSD-excluded (e.g. |H| = 54: min max-PSD 986.8 > 668);
* several order-18 sides are also fully excluded (min max-PSD 693.3);
* **12,048** invariant sequences remain PSD-admissible, forming **1,536
  distinct PAF classes** (a second, slower per-class scan deduplicating by
  PAF vector gives the class count; the two scans agree on every min
  max-PSD, e.g. 512.8 on the 25-orbit order-18 side — class table in
  `05-knowledge/results/legendre333_onesided_psd_classes_macmini_S171.out`).
  Each PAF class poses ONE prescribed-PAF problem for the partner sequence,
  so 1,536 is the count of genuinely distinct remaining attack targets
  (best class: max-PSD 512.8, order-18, 25 orbits).

So the one-sided Legendre-pair route is NOT closed; it is reduced to 12,048 explicit
candidate classes, each posing a prescribed-PAF problem for the partner sequence
(`PAF_B = -2 - PAF_A`), best attacked per-candidate by SAT/pseudo-Boolean
with XOR-cardinality encodings, flattest PSD first.  It is a structured route
to a length-333 Legendre pair, no longer a necessary route to Hadamard order
668 after THM-3393.
