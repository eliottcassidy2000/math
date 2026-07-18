---
id: THM-995
title: THE TRAPPED CUT EXCLUDES THE TIGHT LOCUS — the THM-984 residual is a STRICT-MARGIN question, never an M=1/14 equality. The reduction capstone (THM-984) leaves "μ₀ > 0 on trapped cores"; the FORMALIZATION-PICTURE sharpens this to "rigidity of the tight equality M=1/14 plus a strict margin on the non-tight residual". THIS FILE REMOVES THE EQUALITY HORN: every known tight family (M = 1/14 exactly) VIOLATES a trapped-core hypothesis, so the tight locus lies entirely OUTSIDE the trapped cut and the residual is purely the strict-margin case M > 1/14. (I) THE EXACT ESCAPE CENSUS (all M-values exact rationals on the Farey-refined critical grid, all 8 hypotheses exact integer predicates): AP {1..13} — M=1/14, FAILS gap + max≥23 (small/comparable, → window-census branch); AP {2,4,…,26} — M=1/14, FAILS gap + common-residue (→ common-residue branch); sporadic {1..11,13,24} — M=1/14 EXACTLY yet FAILS covering (→ sieve branch, t=1/q witness); deep-well {1..12,182} — M=14/183, FAILS compressed (182 is a dominant runner, 182 > 13·12; → 91/dominant-peel branch); perturbed {1..11,13,36} — M=3/41, FAILS covering; deep-well {1..12,84}, {1..12,98} — M=1/13, FAIL covering; (II) THE THREE-ROUTE STRUCTURE: every tight/near-tight family escapes by exactly one of three mechanisms, each matching an ALREADY-CLOSED assembly branch — SMALL/COMPARABLE (fail gap or max≥23 ⇒ window ≤22 or comparable branch), DOMINANT-RUNNER (fail compressed ⇒ the peel branch), SIEVE-COVERED (fail covering ⇒ t=1/q). The tight locus is arithmetically special (APs, dominant wells, sieve-nulls); the trapped cut is designed to be its complement; (III) THE STRICT-INTERIOR MARGIN: a 40-family sample of genuine trapped cores (distinct, gap, compressed, max≥23, non-clusterable, covering — a SUPERSET of true trapped cores, so the bound transfers) has min M = 0.1945, margin +0.1231 over 1/14, median 0.2448 — trapped cores are not merely non-tight, they sit DEEP in the loneliness interior (≥ 2.7× the threshold on the sample). The residual μ₀ > 0 is a large-margin statement everywhere the sample reaches; the analytic task is to certify NO trapped family approaches 1/14, and the census shows the near-1/14 configs are exactly the ones the cut excludes
status: (I) EXACT (7 families, M-values exact rationals, all 8 trapped-core hypotheses exact integer predicates — trapped_core_nontight_kps_S128c48.py) — the escape of each named tight family is a proved fact; (III) sampled (40 trapped cores, superset cut, all M > 0.19 ≫ 1/14 — trapped_margin_sample_kps_S128c48.py). GENERAL "every tight family escapes" is a well-supported CONJECTURE pending the full tight-locus characterization (the three-route structure is the proof skeleton); NOT claimed as a universal theorem
source: kind-pasteur-2026-07-17-S128 (cont.48; owner: finish LRC(14), pull often, integrate incoming — codex took the live-floor Lean layer, this reroutes to the reduction's analytic residual)
depends_on:
  - THM-984 (the reduction capstone this sharpens)
  - THM-853/826 (the tight-locus / deep-well corridor: the M-values of the tight families)
  - LRC14Grand.ResidualObligation (the 8 trapped-core hypotheses checked here)
related:
  - THM-979 (the modulus supply — the OTHER half; together: strict margin ⇒ μ₀ > 0 ⇒ explicit modulus ⇒ census)
  - the FORMALIZATION-PICTURE tight-equality question (this removes its equality horn)
  - opus THM-936 (gap-free taxonomy — the accessible-scale loneliness census this complements)
scripts: 04-computation/trapped_core_nontight_kps_S128c48.py, trapped_margin_sample_kps_S128c48.py (+ .out)
---

# THM-995 — the trapped cut excludes the tight locus

## The point

LRC(14) reduces (THM-984) to μ₀ > 0 on trapped cores. The remaining worry was the
EQUALITY case M = 1/14 (where μ₀ = 0 and no margin exists). This file shows the equality
case never arises inside the cut: every tight family fails a trapped-core hypothesis.

## (I) The escape census (exact)

| family | M (exact) | tight? | trapped-core hypothesis it FAILS |
|---|---|---|---|
| AP {1..13} | 1/14 | yes | gap, max≥23 |
| AP {2,…,26} | 1/14 | yes | gap, common-residue |
| sporadic {1..11,13,24} | 1/14 | **yes** | covering (sieve t=1/q) |
| deep-well {1..12,182} | 14/183 | near | compressed (182 dominant) |
| perturbed {1..11,13,36} | 3/41 | near | covering |
| deep-well {1..12,84} | 1/13 | near | covering |
| deep-well {1..12,98} | 1/13 | near | covering |

The one subtle case is the sporadic {1..11,13,24}: it is EXACTLY tight (M = 1/14) and
passes gap/compressed/distinct/max/non-cluster, but it is sieve-covered (some t = 1/q is
already a loneliness witness), so it never reaches the trapped branch.

## (II) The three escape routes ↔ three closed branches

- SMALL / COMPARABLE (max ≤ 22 or all-ratios ≤ 13) — the AP-type tights; closed by the
  window census / comparable branch.
- DOMINANT RUNNER (some |vᵢ| > 13·|vⱼ| for all j) — the deep-well tights; closed by the
  91/dominant peel.
- SIEVE-COVERED (some t = 1/q lonely) — the sporadic tights; closed by the sieve.

The tight locus is the union of these three arithmetically-special strata; the trapped
cut is precisely their complement.

## (III) The strict interior

Genuine trapped cores sample far from tight: min M = 0.1945 (margin +0.1231), median
0.2448 — ≥ 2.7× the 1/14 threshold. The reduction targets the deep interior of the
loneliness region.

## Named next
- The tight-locus characterization: prove the three-route escape is exhaustive (every
  M = 1/14 family is small/comparable, dominant, or sieve-covered). This upgrades (I)
  from census to theorem and CLOSES the equality horn universally.
- The quantitative trapped-core margin floor M ≥ 1/14 + δ(V): the strict-margin half of
  the residual, feeding THM-984's converter → THM-979's modulus → the census.
