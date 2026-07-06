# The Two-Band Theorem: a height-uniform witness mechanism for the 1d crux

**opus-2026-07-06-S96 (HYP-4246).** The finish-list route "(a) … or generalize
opus-S81's descent off the lift stratum," executed. The mechanism that makes witnesses
height-independent on CRT-frozen rays is not descent and not fees: it is the
**scale-invariance of the free fraction** plus elementary period counting.

## Statement (candidate theorem, verified numerically end to end)

Let the family split as CORE ∪ S·P where:
* CORE = bounded small speeds (the census's bottom band), with a good component
  J ⊆ [0,1] on which every core runner has margin ≥ 1/14 + ε (existence: cite
  LRC(|CORE|+1) at a point, Lipschitz-widen; |J| ≥ 2(1/13 − 1/14)/B_core ≥ 1/1092
  in the worst case);
* P = a fixed pattern of distinct integers (the top band's shape), S = the ray scale;
* φ(P) = 1 − |union of the 1/14-teeth of P| = the pattern's FREE FRACTION.

Then for every S > 2/|J|:
    |free(S·P) ∩ J| ≥ φ(P)·|J| − 2φ(P)/S > 0,
and any point of the intersection is a strict 1/14-witness for the whole family.

**Proof of the displayed bound (period counting, fully elementary):** free(S·P) is
(1/S)·free(P) tiled with period 1/S. J contains at least S·|J| − 2 complete periods;
each contributes φ(P)/S of free measure. Points of free(S·P) have dist > 1/14 to every
tooth of every S·P-runner (open complement of closed teeth); points of J have
margin ≥ 1/14 + ε for the core. ∎

## Verification (two_band_H1_opus_S96.out)

The pole-necessity family (mac-mini S55): CORE = {1,2,3,4,5}, P = {20,21,24,25,45,46,66}
(the floating 7-cluster that defeats every profile filter — no C_7 exists, profile
periodic on frozen rays):
* φ(P) = 1632811/5578650 ≈ 0.2927, K = 110 free intervals per period;
* the largest core component J = (0.0734, 0.1853), |J| = 0.1119;
* exact clear measure inside J: 0.0192 (S=1) → 0.03275 = φ·|J| (S large) — measured at
  S = 1, 3, 7, 20, 143, 1001, 30030: **witness verified at every scale**;
* uniformity threshold 2/|J| = 18; S ≤ 18 checked directly. THE FLOATING CLUSTER IS
  DEAD AT ALL SCALES, elementarily.

## Why φ > 0 is robust (and the one lemma to nail)

* |P| = 7 at band 1/14 sits EXACTLY on the fee-mean ceiling (Σ 2ρ = 1): union bounds
  give nothing; the free measure is pure overlap credit.
* Pairwise overlap densities are exactly (2ρ)² regardless of gcds (dilation reduces the
  pair to the coprime case), so second-moment credit is pattern-independent: 21/49 = 3/7.
* Pairwise-coprime patterns: full CRT independence gives φ = (1 − 1/7)^7 = (6/7)^7 ≈ 0.3399
  EXACTLY.
* COMP scan (300 random + adversarial survivor-shaped patterns in [19,70]):
  φ ∈ [0.2996, 0.402]; consecutive blocks minimize (0.2926 for {21..27}) — the same
  extremality as the desert catalog's.
* **The φ > 0 lemma (open, Newman-shaped):** distinct-frequency interval-combs never
  tile [0,1). Analogy: Newman's disjoint-covering-systems rigidity (exact covers of ℤ
  by APs force repeated largest modulus). A tiling here would be an exact cover by
  1/14-teeth at distinct frequencies; boundary-matching at the fastest frequency should
  force a duplicate. Until proven, φ > 0 is verified per pattern (a one-line exact
  computation — table-able over the census anatomy).

## What this does to the crux (the reduction, refined)

Gap emptiness on a CRT-frozen ray = (i) the ray's family is two-band (census anatomy:
w_min ≤ 4, bottom ≤ 25, top S-scaled); (ii) φ(top pattern) > 0 (computable per ray
shape; conjecturally always); (iii) S > 2/|J|: witness by THIS theorem, uniform in S;
(iv) S ≤ 2/|J|: heights ≤ ~2184·max(P) — inside the verified census/CRT-lift ranges.
The template lane (Q0, kps/mac-mini) handles single-scale checks; the two-band theorem
handles the rays BETWEEN scales — the pair is exactly the "periodicity reduction"
(mac-mini S55) with the witness-side mechanism now NAMED and PROVED-modulo-φ.

Multi-band families recurse (CORE ∪ band₂ becomes the new core at the next scale gap).

## Formalization notes

Period counting is Finset arithmetic. Free-set membership at scale S is the toothMiss
integer schema at denominator 14·Q·S (kernel decide; the ClearCert pattern transplants
with 13 → 14 — a sibling file, mechanical). The per-pattern φ > 0 checks are exact
rational interval sweeps (decide-able as list computations if wanted).


## S97 UPGRADE: the theorem goes EXACT and FORMAL (HYP-4256)

The measure/period-counting version above is superseded by the EXACT TRANSPORT — and
the open lemma disappears:

**Theorem (two_band_transport, GREEN in LRCTwoBand.lean).** Core clear interval
(a, a+L) (pointwise > 1/14, certified by toothMiss tables), pattern witness t_P with
all |p·t_P − m| ≥ 1/13 (supplied by CITING LRC(13) on the pattern: |P| ≤ 12 speeds),
and any scale S with S·L > 1. Set k = ⌊S·a − t_P⌋ + 1 and t = (t_P + k)/S. Then
t ∈ (a, a+L) and for every p ∈ P: (S·p)·t = p·t_P + p·k ≡ p·t_P (mod 1) — the top
band's margins are EXACTLY the pattern's, at every scale. Corollary two_band_lonely14:
all runners strictly above 1/14.

* No φ, no measure, no Lipschitz on the top band. The S96 free-fraction machinery
  remains true but is no longer needed; the Newman-shaped φ > 0 lemma is BYPASSED
  (and is separately a corollary of the citation: the pattern's LRC(13) witness has a
  free neighborhood of width ≥ 2(1/13 − 1/14)/max(P)).
* The scale threshold improves to S > 1/L (from 2/|J|); L comes from the citation on
  the core + Lipschitz, so the gap condition between bands is explicit and small.
* Verified end-to-end at S up to 10^12 + 7 (two_band_exact_opus_S97.out): top margins
  exactly 2/13 (the pattern witness's margin, transported without loss).
* Height-uniformity on CRT-frozen rays is now a FOUR-LINE formal theorem consuming two
  citation instances (core and pattern) plus integer arithmetic.

REMAINING for the crux: the GAPLESS residual — families where no band split satisfies
S·L > 1 (all consecutive scale ratios below ~84·B_core-ish). These cannot form frozen
rays (scaling any sub-block of a ray creates a gap), so they are height-bounded per the
census's reach or fall to the pinning-sparsification of near-unit constraints — the
single-scale template lane (Q0, kps/mac-mini) is exactly their home. The division of
labor is now clean and total: gapped = two-band transport (formal); gapless = single
scale (template/census).
