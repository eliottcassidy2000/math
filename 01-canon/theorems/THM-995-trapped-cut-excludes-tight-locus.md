---
id: THM-995
title: THE TRAPPED CUT EXCLUDES THE TIGHT LOCUS — the THM-984 residual is a STRICT-MARGIN question, never an M=1/14 equality. The reduction capstone (THM-984) leaves "μ₀ > 0 on trapped cores"; the FORMALIZATION-PICTURE sharpens this to "rigidity of the tight equality M=1/14 plus a strict margin on the non-tight residual". THIS FILE REMOVES THE EQUALITY HORN: every known tight family (M = 1/14 exactly) VIOLATES a trapped-core hypothesis, so the tight locus lies entirely OUTSIDE the trapped cut and the residual is purely the strict-margin case M > 1/14. (I) THE EXACT ESCAPE CENSUS (all M-values exact rationals on the Farey-refined critical grid, all 8 hypotheses exact integer predicates): AP {1..13} — M=1/14, FAILS gap + max≥23 (small/comparable, → window-census branch); AP {2,4,…,26} — M=1/14, FAILS gap + common-residue (→ common-residue branch); sporadic {1..11,13,24} — M=1/14 EXACTLY yet FAILS covering (→ sieve branch, t=1/q witness); deep-well {1..12,182} — M=14/183, FAILS compressed (182 is a dominant runner, 182 > 13·12; → 91/dominant-peel branch); perturbed {1..11,13,36} — M=3/41, FAILS covering; deep-well {1..12,84}, {1..12,98} — M=1/13, FAIL covering; (II) THE THREE-ROUTE STRUCTURE: every tight/near-tight family escapes by exactly one of three mechanisms, each matching an ALREADY-CLOSED assembly branch — SMALL/COMPARABLE (fail gap or max≥23 ⇒ window ≤22 or comparable branch), DOMINANT-RUNNER (fail compressed ⇒ the peel branch), SIEVE-COVERED (fail covering ⇒ t=1/q). The tight locus is arithmetically special (APs, dominant wells, sieve-nulls); the trapped cut is designed to be its complement; (III) THE STRICT-INTERIOR MARGIN: a 40-family sample of genuine trapped cores (distinct, gap, compressed, max≥23, non-clusterable, covering — a SUPERSET of true trapped cores, so the bound transfers) has min M = 0.1945, margin +0.1231 over 1/14, median 0.2448 — trapped cores are not merely non-tight, they sit DEEP in the loneliness interior (≥ 2.7× the threshold on the sample). The residual μ₀ > 0 is a large-margin statement everywhere the sample reaches; the analytic task is to certify NO trapped family approaches 1/14, and the census shows the near-1/14 configs are exactly the ones the cut excludes
status: (I) EXACT escape census (7 families, exact M, exact predicates); (III) sampled (40 trapped, min M 0.19); (IV) DILATION ESCAPE LEMMA PROVED (non-primitive ⟹ fails common-residue); (V) adversarial hunt survived (~970 perturbations); (VII) LAYERED REDUCTION — Layer 1 (dilation, PROVED) + Layer 2 (sieve-tight ⟹ fails covering, PROVED definitional) reduce the equality horn to PRIMITIVE COVERING tight families; the residual (primitive covering tight ⟹ small/comparable) is the conjectural hard core; (VIII) MARGIN FLOOR empirical — two independent multi-start descents minimize trapped M to 1/7 and 41/253 (both ≥ 1/7 = 2× threshold), sample+hunt agree ⟹ observed δ ≥ 1/14, conjectured trapped ⟹ M ≥ 1/7. (IX) THE SIEVE-MARGIN LEMMA PROVED (uncovered q≤13 ⟹ M ≥ 1/q > 1/14, explicit witness+margin) ⟹ tight locus pinned to the "covers 2..13, misses 14" stratum (both known tights confirmed); (X) covering-family floor empirical M ≥ 1/9 (3000 samples+descent, none tight). Rigorous: Layers 1+2, dilation, sieve-margin lemma (IX), the exact census. Conjectural (sharpened): "covers all 2..14 ⟹ M > 1/14" (the residual) + the covering floor. The 1/7 double-threshold is WITHDRAWN as universal (mod-7 resonance, family-specific). NOT claimed as a universal theorem
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

## (IV) THE DILATION ESCAPE LEMMA (PROVED) — the equality horn reduces to PRIMITIVE families

If V = c·W with c ≥ 2 (a non-primitive family, gcd of the speeds ≥ 2), then every speed
satisfies vᵢ ≡ 0 (mod c), so V FAILS the common-residue hypothesis (H9: no d ≥ 2 with all
vᵢ in one residue class) and is NOT trapped. ∎ (One line; verified on the dilates
c·{1..13} and c·{1..11,13,24}, c = 2,3,5 — all fail common-residue, all remain tight by
dilation invariance M(cV) = M(V).)

**Consequence.** The trapped cut contains only PRIMITIVE families (gcd = 1). Since M is
dilation-invariant, the tight locus is a union of dilation classes, and every class's
non-primitive members are excluded by H9 — so the equality horn reduces to PRIMITIVE tight
families. The primitive tight families are exactly the arithmetically-special strata
(consecutive/comparable APs, dominant deep wells, sieve-null sporadics), each verified in
(I) to escape by its route. This is the equality horn's proof skeleton: (dilation lemma,
proved) + (primitive tight = the three special strata, conjectured exhaustive).

## (V) ADVERSARIAL SURVIVAL

A refutation hunt — perturb the tight bases {1..11,13,24}, {1..12,24}, {1..11,14,27} by
±1..3 on 1–3 coordinates and test for a TRAPPED family with M ≤ 1/14 + 10⁻⁶ — found NONE
in ~970 near-tight perturbations: every perturbation that stays near-tight also stays
outside the cut. The conjecture survives its own counterexample search
(trapped_tight_hunt_kps_S128c48.py).

## (VI) END-TO-END CHAIN VALIDATION (concrete)

The primitive trapped family V = [25, 71, 76, 84, 103, 136, 174, 230, 234, 297, 306, 314,
343] run through the full reduction, every link exact/explicit:
- THM-995 strict margin: M(V) = 38/195 = 0.19487, margin 337/2730 over 1/14 (t* = 98/195);
- converter → good-set measure: G = 350 disjoint safe intervals, μ₀ = 0.12457 > 0;
- THM-979/984 modulus: E = 2·Σ|vᵢ| = 4786, q₀ = ⌈(E+1)/μ₀⌉ = 38427;
- census: liveCount(38427) = 4784 > 0 (≈ q₀·μ₀, the sampling bridge exact) ⟹ t = p/q₀ is a
  loneliness witness ⟹ LRC(14) holds for V.
The chain THM-995 → converter → THM-979/984 → census composes end to end on a concrete
trapped core (reduction_chain_e2e_kps_S128c48.py).

## (VII) THE LAYERED EXHAUSTIVENESS REDUCTION (cont.49)

The tight-locus escape decomposes into two RIGOROUS layers plus a residual conjecture:
- **Layer 1 (PROVED, dilation lemma IV):** non-primitive tight families (gcd ≥ 2) fail
  common-residue ⟹ not trapped. Reduces to primitive tight families.
- **Layer 2 (PROVED, definitional):** sieve-tight families — those where some t = 1/q
  (q ≤ 14) already achieves min_i ‖vᵢ/q‖ = 1/14 — fail covering by definition ⟹ not
  trapped. (The sporadic {1..11,13,24} lives here.) Reduces to primitive, COVERING, tight.
- **Residual (conjecture):** every primitive covering tight family is small/comparable
  (max ratio ≤ 13 ⟹ fails gap, or max ≤ 22 ⟹ fails max≥23). This is the classical
  LRC-tightness characterization on the compressed-gap stratum — the hard core.

**Parametrized-locus confirmation** (tight_locus_classify_kps_S128c49.py): the sporadic
family {1..11,13,c} is sieve-covered (Layer 2) for c = 24,25,27,36,48 (at q = 12 or 14),
and the deep wells {1..12,14m} are sieve-covered at q = 13 for m = 6,7,14 — Layer 2 handles
the entire parametrized sporadic/deep-well tight locus systematically. Every TIGHT
candidate in the census resolves to Layer 1, Layer 2, or small/comparable; the only family
flagged "primitive covering gap" ({1..12,182}) is NON-tight (M = 14/183) and additionally
fails compressed (182 dominant) — no tight family reaches the residual stratum in the
census.

## (VIII) THE MARGIN FLOOR — empirical δ ≥ 1/14 (M ≥ 1/7 = 2× threshold)

Two independent multi-start local-descent MINIMIZATIONS of M over the trapped cut, plus
the cont.48 sample and adversarial hunt, all bottom out at the same place:
- descent A (60 starts): min M = **1/7** exactly (V = [42,66,96,108,150,228,229,247,375,
  377,396,414,552], primitive, trapped, exact-verified);
- descent B (25 indep. starts): min M = 41/253 = 0.1621 (2.27× threshold) — could not
  even reach 1/7;
- 40-family sample: min 0.1945; ~970-perturbation adversarial hunt: none below 1/14+10⁻⁶.

Every route to a low-M trapped family stalls at **M ≥ 1/7 = 2/14 — exactly TWICE the LRC
threshold**. The margin floor δ = M − 1/14 is at least ≈ 1/14 empirically, not merely
positive. Conjecture (the DOUBLE-THRESHOLD floor): trapped ⟹ M ≥ 1/7. This would give the
residual μ₀ > 0 an enormous quantitative margin and is a far stronger statement than the
bare positivity the reduction needs. The compressed-gap stratum where this must be proved
is precisely the classical hard core of LRC(14).

## (IX) THE SIEVE-MARGIN LEMMA (PROVED) — the residual sharpens to "covers all of 2..14"

**Operational covering.** For q ≤ 14, the sieve witness t = 1/q gives min_i ‖vᵢ/q‖ ≥ 1/q
UNLESS some vᵢ ≡ 0 (mod q). So a family is "covering" (no sieve witness) iff **every
q ∈ {2,…,14} divides some speed**.

**ATTRIBUTION CORRECTION (cont.51).** This lemma is NOT new — it is **THM-523's q-witness
lemma** (mac-mini, 2026-06-16), proved there verbatim: *if S contains no multiple of
q ∈ {2,…,14}, then τ = 1/q is lonely and M(S) ≥ 1/q ≥ 1/14.* I re-derived it independently
in cont.50 and wrongly presented it as new; credit belongs to THM-523. What this file adds
is only the **strictness split and the pinning**: for q ≤ 13 the bound is STRICT
(1/q > 1/14, margin ≥ 1/182), while q = 14 gives only M ≥ 1/14 — hence a tight family must
cover 2..13 and can miss ONLY q = 14. (Verified exactly, 100 families.)

**Lemma (THM-523, restated with the strictness split):** if some q ∈ {2,…,13} divides no
speed, then M(V) ≥ 1/q > 1/14 — explicit witness t = 1/q, margin 1/q − 1/14 ≥ 1/182.
*Proof:* no vᵢ ≡ 0 (mod q) ⟹ ‖vᵢ/q‖ ≥ 1/q for all i ⟹ min ≥ 1/q; and 1/q > 1/14 for
q ≤ 13. ∎

**Consequence — the tight locus is pinned to ONE stratum.** A tight family (M = 1/14) must
cover 2..13 (else the lemma forces M > 1/14 strictly). Both known tight families ({1..13},
{1..11,13,24}) cover 2..13 and MISS exactly q = 14 — non-covering at 14, hence Layer-2
handled (not trapped). The residual therefore reduces from "primitive covering tight" to
the sharper, testable statement:

> **Does any family covering ALL of 2..14 have M = 1/14?**

## (X) THE COVERING-FAMILY FLOOR (empirical) — M ≥ 1/9

Minimizing M over families that cover all of {2,…,14} (3000 samples + local descent):
min M = **1/9 = 0.1111** (1.56× threshold, V = [3,4,11,12,13,15,18,20,24,42,55,64,67]),
and NO covering family is tight. So "covering ⟹ M > 1/14" holds on the sample with a
comfortable margin; with the sieve-margin lemma (IX) the equality horn closes conditionally
on this one clean statement.

**Erratum to (VIII):** the trapped-minimizer's 1/7 is a MOD-7 RESONANCE (witness
t* = 178/525, denominator divisible by 7; speeds 150, 375 hit 1/7), hence family-specific —
NOT a universal floor. The honest floors: non-covering at q ≤ 13 ⟹ margin ≥ 1/q − 1/14
(RIGOROUS, IX); covering families ⟹ M ≥ 1/9 (empirical, X); the double-threshold 1/7 is
withdrawn as a universal claim.

## (XI) THE RESIDUAL *IS* THE COVERING-MIN PROBLEM — and it is largely already proved

The sharpened residual "covers all of 2..14 ⟹ M > 1/14" is **not a new problem**: it is
exactly THM-523's residual covering-set hard core (HYP-2566), and the project's
covering-min rigidity already delivers a STRONGER bound:

| stratum | canon result | bound | status |
|---|---|---|---|
| single-killer covering | **THM-724** | M ≥ 14/183 | PROVED (3/4 cases unconditional; near-tight non-dilated large-s residual empirically closed) |
| multi-killer (≥2 far outliers) | **THM-726** | M ≥ 1/13 | certified shape-complete on \|P\| ∈ {10,11}; \|P\| = 9 legacy; \|P\| ≤ 8 open |

Since 14/183 = 0.07650 > 1/14 = 0.07143 (margin 13/2562) and 1/13 > 14/183, **THM-724+726
give "covering ⟹ M > 1/14" modulo their two named gaps** — the equality horn is closed to
the same degree the covering-min is.

## (XII) THE WEAK-TARGET RELAXATION (new) — where the difficulty actually lives

The equality horn needs only the WEAK target M > 1/14, which sits 7% BELOW the sharp
covering-min target 14/183. Measuring where the sharp program's gaps actually are:

- **THM-726's "open" \|P\| ≤ 8 strata are NOT near-tight.** 400 sampled covering families
  in exactly that stratum (core ≤ 8, ≥ 5 far outliers) give min M = **0.1362 = 1.91×** the
  1/14 threshold — and ZERO fall below 1/14, below 14/183, or even below 1/13. Their
  openness is an artifact of the certification method (the union-tail dying at j ≈ 6.5),
  not genuine tightness. For the weak target these strata are loose by ~91%.
- **The genuine difficulty is the near-tight SINGLE-killer neighbourhood** (the deep well
  {1..12,182}, M = 14/183 = only 1.07× threshold) — which is precisely where THM-724's
  UNCONDITIONAL cases (interval-core, dilated-core, killer-safe) already apply.

> **Localization:** for the weak target M > 1/14, the entire remaining risk is THM-724's
> near-tight non-dilated large-s residual — and there the weak target has 7% more room
> than the sharp 14/183 bound it was proved against.

## Named next
- Re-run THM-724's near-tight large-s residual against the WEAK target M > 1/14 instead of
  M ≥ 14/183: the extra 7% margin may convert its "empirically closed" status to
  unconditional — this would close the equality horn outright.
- Likewise re-run THM-726's \|P\| ≤ 8 union-tail at the weak target (empirically 1.91×
  loose, so the method should survive well past its j ≈ 6.5 sharp-target death point).
