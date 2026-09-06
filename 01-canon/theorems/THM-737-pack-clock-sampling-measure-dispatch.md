---
id: THM-737  # originally filed as THM-735; ceded to kps-S128c3 (Bonferroni multi-peel) after a three-way same-hour collision; mac-mini set the cession precedent (their 735 -> 736)
title: The pack-clock sampling lemma with the necessary gcd-aware exception budget
status: PROVED gcd-aware measure inequality and coprime specialization; REFUTED former unqualified coprime extension, repaired 2026-09-06
source: opus-2026-07-13-S272 (owner prompt: prove the uniform-in-c closure of the compressed tower; prior-art check found THM-668 already gives the tower M ≥ 1/13 — this theorem is the measure-level form and the d ≥ 2 extension)
depends_on: []   # at λ = 1/14 the pack input |G_R| > 0 is a finite exact computation (no LRC citation); the λ→1/13 corollary for |R| = 12 packs uses LRC(13) (named citation per policy)
related:
  - THM-668 (monad-S3)         # the point-witness detuned dispatch g·H ∪ {δ} ⟹ M ≥ 1/13; mechanism precedent; its "Scope" item 4 (d ≥ 2) is partially closed here
  - LRCClusterGcd (kps-S18)    # 1/g-periodic margins + pigeonhole, the shared engine (THM-668 credits it too)
  - death-star-S14 V_L          # the compressed-stratum floor 1/13 at every diameter (context)
  - MISTAKE-140/141             # the compressed near-dilate class; M(2·{1..12}∪{13}) = 1/13 exact — this lemma reproduces both the exact L AND the exact M at c=2
  - THM-731/732, HYP-6520/6525 (opus S270-271)  # the peel/perspective certificates; this lemma closes the coherent-pack sector of the non-isolated residual they map
---

# THM-737 — the pack-clock sampling lemma (measure form of the detuned dispatch)

**CORRECTION — 2026-09-06.** The gcd-aware displayed inequality and its
proof are valid. The former title, scalar cutoff paragraph, and scope prose
incorrectly extended the coprime specialization to every nonmultiple
exception. For `R={1,...,12}`, `c=4`, `D={26}`, the actual measure is
`6617/388080`, below the falsely claimed `3|G_R|/4=6617/258720`. Its true
order is `c/gcd(c,26)=2`; the corrected exact floor `|G_R|/2` is attained.
The [independently checked repair and effective-order certificate](../../05-knowledge/results/lrc14_effective_clock_empty_core_sep06.md)
preserve the original mechanism and give its sharp integer form
`L>=|G_R| max(0,1-sum_i ceil(q_i/7)/q_i)`, `q_i=c/gcd(c,u_i)`.
No loneliness result is refuted by this measure witness.

## Statement

Let `c ≥ 2`, let `R` be a finite set of nonzero integers (the **pack**), `D = {u₁, …, u_d}`
nonzero integers (the **detuned elements**), and consider the speed family `v = c·R ∪ D`.
Let `G = G_R^λ = {s ∈ [0,1) : ‖w s‖ ≥ λ ∀ w ∈ R}` be the pack good set at threshold `λ`,
and `L_λ(v)` the measure of the `λ`-safe times of `v`. Then

> **`L_λ(v) ≥ |G_R^λ| · (1 − Σᵢ (2λc + gᵢ)/c)`, where `gᵢ = gcd(uᵢ, c)`** —

and in the coprime case (`gᵢ = 1` for all `i`), with the sharp integer count,

> **`L_λ(v) ≥ |G_R^λ| · (1 − d·(⌊2λc⌋ + 1)/c)`,  (⌊2λc⌋+1 → 2λc when 2λc ∈ ℤ).**

For `1<=d<=6`, when every `gᵢ=1` and `|G_R|>0`, at `λ = 1/14` (LRC(14)): `L(v) ≥ |G_R|·(1 − d·(⌊c/7⌋+1)/c) > 0` whenever `c > 7d/(7−d)`
(d = 1: **all c ≥ 2**; d = 2: all c ≥ 3; d = 3: c ≥ 6; d = 4: c ≥ 10; d = 5: c ≥ 18; d = 6: c ≥ 43).
No covering, primitivity or ratio hypothesis is needed. The scalar cutoff
above requires every exception coprime to `c`; other gcd strata retain the
displayed gcd budget. A concrete `R` has an exactly computable safe measure;
positivity for every `|R|<=12` follows from cited LRC through thirteen total
runners and the strict margin above `1/14`.

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
   witness, whose Lipschitz fattening yields only `L ≳ 1/(1092c)`, **decaying in c**. THM-737
   gives `L >= |G|(1-ceil(q/7)/q) >= |G|/2` for one nonmultiple exception,
   where `q=c/gcd(c,u)>=2`. The special factor `1-1/c` requires `q=c<=7`.
   This is a uniform measure floor for the actual gcd orbit. The two
   are complementary: point-threshold vs measure-uniformity.
2. **The `2<=d<=6` zone of THM-668 item 4** is closed for `c > 7d/(7−d)` when every exception is coprime to `c`;
   arbitrary gcd strata require the original budget, as THM-761 also states.
   Only the remaining **coprime** zone is confined to finitely many small
   `c` per `d` (for d=2, only c=2). No unqualified finite-clock conclusion
   follows after discarding the effective orders.
3. **What it does NOT close:** bodies with no coherent sub-pack — `gcd-incoherent` families
   (e.g. klein-S289's `{1,90..101}`: max common-factor pack too small, d would exceed 6).
   Those remain with the perspective certificates (THM-731/732 + HYP-6520: empirically
   certified at 8–13/13 peels) and the exposure-collapse program (kps HYP-6495).
4. **Perspective reading (owner's frame, opus-S270/271):** the pack is ONE perspective — all
   pack pairs are c-commensurate, so the double-counted (cycle) sector of the `(n−1)²`
   decomposition degenerates and the pack shares a single clock; the detuned runner, sampled
   from that clock, visits **`q=c/gcd(c,u)` equally spaced points**, each repeated
   `gcd(c,u)` times among the `c` branches. Equal spacing and the retained
   multiplicity convert equidistribution into exact counting. This is the frozen-fan/sweeper dichotomy
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
