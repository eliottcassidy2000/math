---
id: THM-732
title: The certificate discrepancy disc_v (THM-731) is an EXACT Bernoulli edge-pair sum — a generalized Dedekind sum over the good-set edges — hence (i) exactly computable in ℚ (per-family certificates become finite rational-arithmetic PROOFS of L>0 through the proved THM-731 chain), and (ii) uniformly small in the peel v (|B₂|≤1/6 ⟹ disc_v ≤ r²/(3v²)), closing entire FAR-ELEMENT RAYS (all large dilates of the 13th speed) with one inequality + finitely many exact checks
status: PROVED (identity (i): 3-line Fourier proof below; tail (iii): |B₂|≤1/6, trivial) + VERIFIED EXACTLY (MISTAKE-136 rule: definition == pair form as exact rationals on 8 random toy unions × v∈{5,7,12} AND all four covering families' leave-one-out good sets; script lrc14_disc_bernoulli_dedekind_kps_S128.py). CONSEQUENCES DELIVERED: exact-ℚ certificates prove L>0 for deep well {1..12,182} (L=4637/194040, margin 6|G'|²−disc = 5966922061/909015307200 > 0), near-AP residue {1..11,13,84} (L=563/105105, margin 1815833/3156303150 > 0), {2..14} (disc≡0 EXACTLY — perfect collapse, certificate exact), variant {1,3..13,182}; peeling identity and (*) hold EXACTLY on all four (independent validation of THM-731's chain); exact values match klein-S287's grid values to grid precision. RAY COROLLARIES (all v): {1..12,v} (v₀≈82.9<182 ⟹ zero exact checks), {1..11,13,v} (v₀≈77.5<84 ⟹ zero exact checks — the WORST-BODY ray closes on the crude tail alone), {2..13,v}, {1,3..13,v} (one exact check v=14, margin>0)
source: kind-pasteur-2026-07-13-S128
depends_on:
  - THM-731   # the rigorous chain |ε_v|² ≤ (6/49)disc_v + peeling; this makes its disc_v exact/rational and bounds it in v
related:
  - HYP-6495  # (this session) block/distribution-relation collapse structure toward the GENERAL analytic disc_v bound
  - HYP-6485  # klein-S287 (the x-integral construction)
  - HYP-6248  # kps cont.70 ({1..11,13,84} = worst |core|=1 body — the residue ray this closes)
  - MISTAKE-078 / HYP-2645  # the convergent finite x-cell form the relation-lattice sum kept pointing to — this IS its closed form on the covering side
  - dedekind_e2_lrc_residual_klein.py / dedekind_e2_b2_margin_macmini_20260630.py  # the dormant June-30 Dedekind/E2 thread: covering-min MARGIN = Dedekind sum s(n,Φ₆(n)); now the CERTIFICATE DISCREPANCY is a generalized Dedekind sum too — same technology (distribution relation, reciprocity) now aimed at the live analytic residue
---

# THM-732 — disc_v is a Bernoulli edge-pair (generalized Dedekind) sum

## Statement

Let `G ⊂ ℝ/ℤ` be a finite union of `r` closed intervals with indicator `g`, and let
`E(G)` be its `2r` endpoints ("edges"), signed `σ_e = +1` at left endpoints, `σ_e = −1`
at right endpoints. For every positive integer `v`, with `B₂(t) = t² − t + 1/6` the second
Bernoulli polynomial:

**(i) — exact identity.**
```
disc_v(G) := Σ_{m≠0} |ĝ(mv)|²  =  (1/(2v²)) · Σ_{e,e'∈E(G)} σ_e σ_e' · B₂({v(e−e')}).
```
This is the `disc_v` of THM-731 (Poisson: `disc_v = (1/v)Σ_j A(j/v) − |G|²`, `A` = autocorrelation).

**(ii) — rationality / decidability.** For LRC leave-one-out good sets `G'_{~v}` every edge is a
rational `(14k±1)/(14w)` (w a speed ≠ v, k < w), so every argument `{v(e−e')}` is an explicit
rational with denominator dividing `14·w·w'`, and `disc_v ∈ ℚ` **exactly**. THM-731's certificate
condition `disc_v < 6·|G'_{~v}|²` is then decidable in **exact integer arithmetic** — through the
PROVED THM-731 chain, `L > 0` for any explicit covering family becomes a finite rational check
(no grids, no analysis).

**(iii) — far-element tail.** `|B₂({x})| ≤ 1/6` gives, for every `v`,
```
disc_v(G) ≤ (2r)²/(12 v²) = r²/(3v²),
```
so the certificate holds whenever `v > r / (3√2 · |G|)` with `r, |G|` those of the FIXED
leave-one-out set. **Corollary (rays):** for a fixed 12-speed base `B`, ALL far elements
`v ≥ v₀(B) = r/(3√2·|G'|)` certify at once; the finitely many covering `v < v₀` are exact checks
per (ii). Target rays: the deep-well ray `{1..12} ∪ {182j}` (all j≥1) and the **near-AP residue
ray** `{1..11,13} ∪ {84j}` (all j≥1) — the latter lives in the OPEN separate-13/14 class and
contains the worst known |core|=1 body (kps cont.70).

## Proof of (i) (3 lines)

`ĝ(k) = ∫ g e(−kx) dx = (1/(2πik)) Σ_e σ_e e(−ke)` (endpoint evaluation of the integral over each
interval). Hence `|ĝ(mv)|² = (1/(4π²m²v²)) Σ_{e,e'} σ_e σ_e' e(−mv(e−e'))`. Sum over `m≠0` and use
the Fourier expansion `Σ_{m≠0} e(mθ)/m² = 2π² B₂({θ})` (absolutely convergent; interchange valid):
`disc_v = (1/(2v²)) Σ_{e,e'} σ_e σ_e' B₂({v(e−e')})`. ∎

(Diagonal pairs give `+2r·(1/6)`; positivity of the whole sum is automatic from the LHS.)

## Why this matters

- THM-731 left ONE item: an analytic upper bound on `disc_v`. (i) converts `disc_v` from a grid
  quantity into a **closed-form finite arithmetic sum on the good-set edges** — a generalized
  Dedekind sum with small explicit denominators. (iii) is the first uniform-in-`v` consequence:
  it closes infinite one-parameter slices (rays) of the open class that NO per-family
  computation (grid or exact) could reach.
- The June-30 dormant thread (mac-mini-S64/klein-S56) proved the covering-min MARGIN is a
  Dedekind sum `s(n,Φ₆(n))`; the pentagonal/η²⁴ memory thread supplies the same toolbox
  (distribution relation, reciprocity). This theorem puts the live analytic residue in that
  toolbox's domain: full arc-families collapse via `Σ_{k mod q} B₂({x+k/q}) = B₂({qx})/q`
  (→ HYP-6495).

## Evidence log (this session — all done)

- [x] identity (i) verified EXACTLY (def == pair as rationals): 8 random toy unions × v∈{5,7,12}, plus
      all 4 covering families. Zero discrepancies.
- [x] exact values match klein-S287's grid numbers (L_cert 0.022117/0.004243/0.061224/0.024941 vs klein's
      0.0221/0.0042/exact/0.0249; L_full 0.023897/0.005357 vs 0.0239/0.0054).
- [x] exact-ℚ certificates all POSITIVE: L>0 PROVED for all four families (through the proved THM-731
      chain — first rigorous L>0 for the two extremal bodies, replacing VERIFIED-numerically).
- [x] rays closed: v₀({1..12})≈82.9, v₀({1..11,13})≈77.5, v₀({2..13})≈6.6, v₀({1,3..13})≈15.6.
      Both extremal rays need ZERO exact checks (all covering v exceed v₀).
- Peeling identity `L_full=(6/7)|G'_{~v}|−ε_v` and (*) `|ε_v|²≤(6/49)disc_v` verified as EXACT rational
  identities/inequalities on all four families.

## Notable structure found (feeds HYP-6495)

- `{2..13}` at v=14: **disc ≡ 0 exactly** — ĝ vanishes on 14ℤ\{0} entirely (two good intervals of
  length exactly 1/28 at distance 1/2: parity kills odd m, sinc zeros kill even m). The certificate is
  EXACT there (L_cert = L = 3/49). Perfect Bernoulli collapse exists in nature.
- Exposure is extreme: deep-well base {1..12} has only r=12 good intervals, exactly 2 exposed edges per
  speed (of 2w endpoints). Residue base {1..11,13}: r=4, edges generated only by w∈{5,7,11,13}.
- Off-diagonal blocks carry ~6× the diagonal at the deep well (near-pair structure of narrow intervals),
  NOT equidistributed-B₂ noise: |ĝ(mv)|² ≈ width-weighted Weyl sum of midpoints — the DENSITY route's
  object appears verbatim, making the [A]/[B] unification concrete at formula level.
