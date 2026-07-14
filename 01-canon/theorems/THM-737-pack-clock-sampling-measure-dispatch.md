---
id: THM-737
renumber_note: "RENUMBERED from THM-735 by kind-pasteur-S128 cont.4 (routine ID housekeeping; precedent kps-cont.51 / mac-mini-S89 / klein-S289): three concurrent THM-735 claims on 2026-07-13 — kps claim-checkpoint 17:51:08 (first-pusher, keeps THM-735 = Bonferroni simultaneous multi-peel), mac-mini 17:54 (self-renumbered to THM-736 in S89), opus-S272 17:59 (this file → THM-737). Content unchanged; the S272 script filename lrc14_packclock_sampling_thm735_opus_S272.py retains the pre-renumber id."
title: The pack-clock sampling lemma (measure form of the detuned dispatch) — every 13-family c·R ∪ D (c ≥ 2, d = |D| detuned) has 1/14-safe measure L ≥ |G_R|·(1 − d/7 − Σᵢ gcd(uᵢ,c)/c); uniform-in-c tower closure with the EXACT constant (tower bound-exact for c < 7), and the d ≥ 2 generalization of THM-668 (its open item 4) for all c > 7d/(7−d), with NO congruence condition on the detuned elements
status: PROVED (5-line counting proof below) + VERIFIED-EXACT (all bounds checked as Fraction inequalities; tower c = 2,4,6 EQUALITY; λ-form equality at c=2 recovers M = 1/13 exactly)
source: opus-2026-07-13-S272 (owner prompt: prove the uniform-in-c closure of the compressed tower; prior-art check found THM-668 already gives the tower M ≥ 1/13 — this theorem is the measure-level form and the d ≥ 2 extension)
depends_on: []   # at λ = 1/14 the pack input |G_R| > 0 is a finite exact computation (no LRC citation); the λ→1/13 corollary for |R| = 12 packs uses LRC(13) (named citation per policy)
related:
  - THM-668 (monad-S3)         # the point-witness detuned dispatch g·H ∪ {δ} ⟹ M ≥ 1/13; mechanism precedent; its "Scope" item 4 (d ≥ 2) is partially closed here
  - LRCClusterGcd (kps-S18)    # 1/g-periodic margins + pigeonhole, the shared engine (THM-668 credits it too)
  - death-star-S14 V_L          # the compressed-stratum floor 1/13 at every diameter (context)
  - MISTAKE-140/141             # the compressed near-dilate class; M(2·{1..12}∪{13}) = 1/13 exact — this lemma reproduces both the exact L AND the exact M at c=2
  - THM-731/732, HYP-6520/6525 (opus S270-271)  # the peel/perspective certificates; this lemma closes the coherent-pack sector of the non-isolated residual they map
---

# THM-735 — the pack-clock sampling lemma (measure form of the detuned dispatch)

## Statement

Let `c ≥ 2`, let `R` be a finite set of nonzero integers (the **pack**), `D = {u₁, …, u_d}`
nonzero integers (the **detuned elements**), and consider the speed family `v = c·R ∪ D`.
Let `G = G_R^λ = {s ∈ [0,1) : ‖w s‖ ≥ λ ∀ w ∈ R}` be the pack good set at threshold `λ`,
and `L_λ(v)` the measure of the `λ`-safe times of `v`. Then

> **`L_λ(v) ≥ |G_R^λ| · (1 − Σᵢ (2λc + gᵢ)/c)`, where `gᵢ = gcd(uᵢ, c)`** —

and in the coprime case (`gᵢ = 1` for all `i`), with the sharp integer count,

> **`L_λ(v) ≥ |G_R^λ| · (1 − d·(⌊2λc⌋ + 1)/c)`,  (⌊2λc⌋+1 → 2λc when 2λc ∈ ℤ).**

At `λ = 1/14` (LRC(14)): `L(v) ≥ |G_R|·(1 − d·(⌊c/7⌋+1)/c) > 0` whenever `c > 7d/(7−d)`
(d = 1: **all c ≥ 2**; d = 2: all c ≥ 3; d = 3: c ≥ 6; d = 4: c ≥ 10; d = 5: c ≥ 18; d = 6: c ≥ 43).
No covering, primitivity, ratio, or congruence hypothesis on `D` is needed; `|G_R^{1/14}| > 0`
is a finite exact computation for any concrete `R` (and holds for every `|R| ≤ 12` by LRC(≤13)
+ fattening, if a citation is preferred).

## Proof

Partition `[0,1)` into the `c` **branches** `t = (j+s)/c`, `j = 0..c−1`, `s ∈ [0,1)`.

1. **The pack shares one clock.** `‖(cw)t‖ = ‖w s‖`: the pack condition is `s ∈ G`, *the same
   for every branch j*.
2. **Branch danger count.** Fix `s`. Branch `j` is `uᵢ`-dangerous iff
   `uᵢ(j+s) ∈ (−λc, λc) mod c`. The points `{uᵢ j mod c : j = 0..c−1}` form a `gᵢ`-spaced
   lattice on `ℝ/cℤ` (each point hit `gᵢ` times, `c/gᵢ` distinct points); a `g`-spaced lattice
   meets an open arc of length `2λc` in at most `⌊2λc/g⌋ + 1` points. Hence at most
   `gᵢ·(⌊2λc/gᵢ⌋ + 1) ≤ 2λc + gᵢ` branches are `uᵢ`-dangerous (coprime: `⌊2λc⌋+1`; and `≤ 2λc`
   when `gᵢ = 1`, `2λc ∈ ℤ`, since an open arc of integer length `m` meets a unit lattice in
   `≤ m` points).
3. **Union bound and integrate.** `L_λ(v) = ∫_G (1/c)·#{j : all uᵢ safe} ds ≥
   ∫_G (1 − Σᵢ Nᵢ(s)/c) ds ≥ |G|·(1 − Σᵢ(2λc + gᵢ)/c)`. ∎

## Sharpness (verified exact, `lrc14_packclock_sampling_thm735_opus_S272.py`)

- **The compressed tower `c·{1..12} ∪ {13}` is bound-EXACT for c < 7:** `L = |G_{1..12}|·(1 − 1/c)`
  with **equality as fractions** at c = 2, 4, 6 (6617/388080, 6617/258720, 6617/232848). So for
  a.e. pack-safe `s`, **exactly one** branch is 13-dangerous — the c < 7 rigidity (window
  shorter than the lattice spacing catches ≤ 1 point, and on `G` it always catches one).
  First strict inequality at c = 8 (bound 0.025576 < true 0.029839).
- **The λ-form is tight at c = 2 up the whole threshold:** `L_λ = |G^λ|/2` exactly (checked at
  λ = 5/66); since `|G^λ_{1..12}| > 0 ⟺ λ < 1/13`, this gives **M(2·{1..12}∪{13}) = 1/13
  exactly** (lower bound here; upper = MISTAKE-141) — the lemma reproduces the known exact
  floor of the bounded-DC class, not just positivity.
- **d = 2 closures** (all verified exact, bounds positive and valid): coprime `3{1..11}∪{13,14}`,
  `4{1..11}∪{13,21}`; the **degenerate-diagonal** `5{1..11}∪{13,18}` (13 ≡ 18 mod 5 — exactly
  the branch-subgroup obstruction THM-668 item 4 names as open: the measure argument does not
  see it); a `gcd > 1` case `6{1..11}∪{13,21}`.

## Scope and context

1. **Relation to THM-668** (the point-witness dispatch, PROVED + Lean): THM-668 gives
   `M(g·H ∪ {δ}) ≥ 1/13` via LRC(13) + branch pigeonhole — a stronger *threshold* on a point
   witness, whose Lipschitz fattening yields only `L ≳ 1/(1092c)`, **decaying in c**. THM-735
   gives the **uniform measure floor** `L ≥ |G|(1 − 1/c) ≥ |G|/2` — the covering route's residue
   object (`L`) itself, non-decaying (true `L` is increasing in `c`, opus-S271 scan). The two
   are complementary: point-threshold vs measure-uniformity.
2. **The d ≥ 2 zone of THM-668 item 4** is closed for `c > 7d/(7−d)` with no congruence
   condition; the remaining d ≥ 2 zone is the finitely many small `c` per `d` (e.g. d = 2:
   only c = 2), which stay with THM-668's finite/decidable branch-walk program.
3. **What it does NOT close:** bodies with no coherent sub-pack — `gcd-incoherent` families
   (e.g. klein-S289's `{1,90..101}`: max common-factor pack too small, d would exceed 6).
   Those remain with the perspective certificates (THM-731/732 + HYP-6520: empirically
   certified at 8–13/13 peels) and the exposure-collapse program (kps HYP-6495).
4. **Perspective reading (owner's frame, opus-S270/271):** the pack is ONE perspective — all
   pack pairs are c-commensurate, so the double-counted (cycle) sector of the `(n−1)²`
   decomposition degenerates and the pack shares a single clock; the detuned runner, sampled
   from that clock, is an **equally-spaced c-point sweep** of the circle, and equal spacing
   converts equidistribution into exact counting. This is the frozen-fan/sweeper dichotomy
   with the roles inverted: from the pack's clock, the detuned element is the sweeper, and
   sweeps are countable.
5. **Lean shape:** interval measures (finite rational interval arithmetic) + the lattice-in-arc
   count + a `Finset` union bound — same infrastructure class as `LRCDetunedDispatch.lean`
   (THM-668); the `λ = 1/14` pack input for concrete `R` is a kernel-pure interval computation.
   No analysis at runtime.

## Verification & files

`04-computation/lrc14_packclock_sampling_thm735_opus_S272.py` (+ `.out`): exact-Fraction
verification of every displayed claim (tower c = 2..8 with the three equalities; λ-form; the
four d = 2 / gcd instances). Earlier exact data: opus-S270 (`lrc14_exact_rational_certificates`,
the c = 2 equality first appeared there as `L = 6617/388080`), opus-S271 tower scan.
