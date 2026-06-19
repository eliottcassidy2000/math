---
id: HYP-2640
title: LRC(14) the gcd-3 stratum mod 27 is NOT the carrier of the near-AP danger - the archimedean low-height integer relations are
status: REFUTED (the specific stratum-danger split) / CLARIFIED (the danger is archimedean, matching codex HYP-2634)
source: kind-pasteur-2026-06-19-S13
related:
  - HYP-2637   # Freiman dimension pockets / dimension penalty
  - HYP-2607   # consec maximizes L_y (the crux)
  - HYP-2083   # antipodal summand-unit bridge: V* in gcd-3 blind spot
  - HYP-2632   # codex finite chi_7/affine character kernel
  - HYP-2633   # codex residue-lift equidistribution crux
  - HYP-2634   # codex: opposition is low-height integer relation DEFECTS, not finite QR packets
  - THM-401    # pair-sum sieve modulus C=2n-1=27
---

# HYP-2640: the summand-graph STRATUM mod C=27=3^3 is NOT where the near-AP danger lives

## The tested ANGLE
Per the task: classify the relation lattice Lambda(E) by stratum mod C=2n-1=27=3^3
(unit / gcd-3 / gcd-9), pair antipodal shells {a, 27-a}={+a,-a}, and test whether
the gcd-3 stratum (where the tight sporadic V* lives, HYP-2083) carries the
DANGEROUS near-AP corrections while the unit stratum carries the harmless ones.

## RESULT: REFUTED at the binding (small-spread) scale; danger is ARCHIMEDEAN

Scripts:
- `04-computation/lrc14_summand_strata_relation_lattice_kps_S13.py` (+.out)
- `04-computation/lrc14_stratum_correction_sensitivity_kps_S13.py` (+.out)

Three converging findings:

1. **No modular wraparound at the binding scale.** For the AP `consec` and all
   near-AP / d=2 GAP competitors at k=8,9,10, every element sits in residues
   `0..k-1 < 27`. The relation lattice is dominated by genuine *integer*
   relations (`mixed`/`pure-unit`), and the count of *antipodal* relations
   (using a true `{a,27-a}` pair) is **0**. The modulus-27 structure is INERT
   for the sets that actually approach the cap — it only activates once
   elements exceed ~C/2 (e.g. the full V*, where `24 = -3` makes `{3,24}` a
   real gcd-3 antipodal shell). So the stratification cannot be the mechanism
   separating AP from near-AP.

2. **Residue-preserving shifts (+27, +81) destroy L_y regardless of stratum.**
   Shifting one consec element by a multiple of C keeps its residue and stratum
   identical but breaks the small-integer relations. L_y collapses by ~0.2 in
   ALL cases (k=8: 0.358 -> ~0.10-0.21; k=9: 0.493 -> ~0.29-0.41). The danger
   travels with the **archimedean integer position**, not the residue class.

3. **CONTROLLED test kills the stratum signal.** Fixing position (always
   replacing the central element 4 of k=9 consec) and varying ONLY the
   stranger's residue/stratum, the L_y drop is FLAT across strata:
   `dL in [-0.224, -0.201]`, with gcd-3 values (-0.213, -0.216, -0.204)
   interleaving unit values, gcd-9 (-0.221) in the middle. The apparent
   "gcd-3 is more load-bearing" effect in the uncontrolled test
   (gcd3 mean drop 0.247 vs unit 0.229 at k=8) was a **position confound**:
   the gcd-3 elements 3,6 sit mid-progression, and middle elements are more
   load-bearing than the endpoint (element 7/8 has by far the smallest drop,
   0.043-0.149). Stratum explains nothing once position is held fixed.

## This CONFIRMS codex's HYP-2634 from the empty-sector side

Codex (HYP-2634) found the opposite-sign / dangerous corrections are exact
**low-height INTEGER relation defects** (e.g. `2a-8=0`, `7a-28=0`), localized to
the start-aligned archimedean shape and removed by shifting one period up — i.e.
an integer-lift phenomenon, NOT a finite chi_7 / residue-class phenomenon. My
empty-sector experiment independently reaches the same conclusion: the residue
stratum mod 27 is the wrong coordinate for the binding tail; the archimedean
low-height integer relation lattice is the carrier.

## CONSEQUENCE for the proof program (HYP-2637 step 3a + HYP-2633)
The "summand-shell stratum" decomposition should NOT be expected to separate the
dangerous near-AP corrections analytically. The correct separating coordinate is
the **archimedean height** of the integer relations (the `+27/+81` experiment
shows L_y is governed by which small integer relations survive). This aligns the
two routes:
- HYP-2633's "residue-lift equidistribution after finite low-height wall
  deletion" must do its work in the **integer height** filtration, and the
  finite wall it deletes is precisely the archimedean low-height relation
  defects (HYP-2634), not a gcd-stratum selection.
- The gcd-3 stratum DOES matter for V*'s *existence as a lattice point* mod 27
  (HYP-2083 antipodal {3,24} shell) — but that is the witness/realizability
  question, distinct from the L_y-tail danger. Two different roles for the
  modulus: realizability (gcd-stratum, V*) vs cap-approach danger (archimedean
  height). Conflating them is the error this hypothesis corrects.

## Honest status
LRC(14) NOT proved. This is a clean NEGATIVE result that prunes a plausible-looking
route (stratum-based danger separation) and reinforces the archimedean-height
route already opened by codex. The remaining HYP-2637 step 3a (GAP dimension
penalty, margin >=0.21) and HYP-2633 (height-filtered residue-lift equidistribution)
are unaffected and remain the live crux.
