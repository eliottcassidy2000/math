---
id: THM-515
title: The LRC singular series is the lonely-set measure (theta-form L(S)=Σ_{t∈Λ}∏h(t_i), h positive-definite); the relation-lattice ADDITIVE ENERGY governs L; and the Riesz-product method is the right tool for inf L>0
status: PROVED (the theta/lonely-measure form and positive-definiteness, verified n=13; the additive-energy predictor verified) + PROGRAM (the Riesz-product route to inf L>0, set up concretely, feasibility-probed). The wall (|T|≥3 conditional convergence, THM-504) is NOT broken — every clean reframe confirms it.
source: mac-mini-2026-06-15-S4
depends_on:
  - THM-501   # the singular series L(S) = lim D(q,S)/q
  - THM-503   # 7-vanishing, almost-Sidon loose class, L = archimedean integral
  - THM-504   # the |T|≥3 conditional-convergence wall (cross-level alternation)
related:
  - THM-002   # OCF as a hard-core lattice gas (the free-gas analogue)
  - HYP-2540  # the Riesz-product certificate program
  - OPEN-Q-104
---

# THM-515 — L(S) is the lonely measure; additive energy governs it; Riesz products are the key

## A. The theta-function / lonely-measure form (PROVED)

Let `S = {v_1,…,v_13}`, `Λ = {t ∈ ℤ^13 : Σ_i t_i v_i = 0}` the **relation lattice**
(rank 12), and `h: ℤ → ℝ` with `h(0)=6/7`, `h(t) = −sin(πt/7)/(πt) = −s(t)` for `t≠0`.
Then the THM-501 singular series is a **theta-like lattice sum**:

> **L(S) = Σ_{t ∈ Λ} ∏_{i=1}^{13} h(t_i)** (grouping by support `T=supp(t)` recovers the
> THM-501 form `(6/7)^13 + Σ_{exact relations}(6/7)^{13−|T|}(−1)^{|T|}∏s(t_i)`).

**`h` is positive-definite.** Its Fourier series is
`ĥ(θ) = Σ_t h(t)e(tθ) = 1 − 1_{danger}(θ) = 1_{safe}(θ) ≥ 0`, because the danger-band
indicator `1_{||θ||≤1/14}` has Fourier coefficients exactly `s(t)` (`s(0)=1/7`). Hence
by Poisson summation over the orbit `{θ_i = v_i τ}`,

> **L(S) = ∫_0^1 ∏_{i=1}^{13} 1_{safe}(v_i τ) dτ** = the **Lebesgue measure of the lonely
> set** `{τ : ||v_i τ|| > 1/14 ∀i}`.

This recovers THM-501's `L = lim D(q,S)/q` and makes **`L(S) ≥ 0` structural — from the
MEASURE side** (it is `∫` of a product of `{0,1}`-indicators). NOTE: the *lattice sum*
`Σ_{t∈Λ}∏h` genuinely ALTERNATES (≈ equal positive/negative relation terms — this IS the
THM-504 cross-level cancellation); positive-definiteness gives the safe-set FOURIER side
(`ĥ=1_safe≥0`) and hence the measure identity, NOT termwise nonnegativity of the sum. `L>0 ⟺
loose`, `L=0 ⟺ tight` (THM-501). Verified: `ĥ ≈ 1_{safe}` numerically; `L = `lonely
measure matches `D(Q)/Q` at the cores (extremizer `{1..13}\{6}∪{56}`: `L≈0.0056`).

## B. The relation-lattice ADDITIVE ENERGY governs L (verified; corrects the λ_1 guess)

The naive geometry-of-numbers predictor `λ_1(Λ)` (shortest relation, in `Σ|t_i|`) is
**uninformative**: `λ_1 = 3` for essentially every 13-set (a 3-term relation
`v_i + v_j = v_k` or `2v_i = v_j+v_k` almost always exists). The correct predictor is the
**additive energy** `E(S) = #{short 7-primitive relations}` (the relation *density*, not
the shortest vector):

> **High additive energy ⟺ small L** (strong correlation, not an exact functional).
> Generic/Sidon-like sets carry `≈ 700–10k` short relations and `L ≈ 0.12–0.15 ≈
> (6/7)^13`; near-tight / AP cores carry `≈ 17k–24k` and `L ≈ 0.005–0.03`. The exact
> value is set by the `∏s`-WEIGHTED relation sum, not the raw count: the evader
> `7·{1..12}∪{13}` has the *highest* short-relation count (`18172`) yet `L≈0.029`, while
> the lower-count interior-drop extremizers `{1..13}\{j}∪{14m}` reach `L≈0.0053` — so
> additive energy separates the *regimes* (generic vs core) cleanly but the extremizer is
> pinned by the detailed sinc weighting.

**Tao's third moment confirms the mechanism.** The multiplicity `M(τ)=Σ_v 1_{danger}`
has `∫M = 13/7` (always), and the *third* moment `∫M³` (triple danger-overlaps =
arithmetic-progression structure, Tao's 2017 method) tracks the regime: cores have
`∫M³ ≈ 35–39` (heavy triple overlaps ⟺ AP ⟺ low loneliness), generic has `∫M³ ≈ 17.8`.
This matches the mesh finding (THM-503 corr / HYP-2526): `L` is set by the balanced
*core*'s resonances, untouched by dominant strangers (dominance ⊥ `L`). Geometrically,
the t=0 term gives the ideal-gas main `(6/7)^13`; the corrections are the `∏s`-weighted
relation-lattice sum, whose suppression is driven by the relation density.

## C. The Riesz-product method is the right tool for inf L>0 (PROGRAM, HYP-2540)

The conditional-convergence wall (THM-504: `|T|≥3` absolute sum diverges, cross-level
alternation) is **not broken by any positivity reframe** — theta, free-gas, FKN all
recover `L≥0` and confirm the wall. The published state-of-the-art for the *cousin*
problem (the lonely-runner gap) is the **Riesz-product method** (Bedert, arXiv:2511.16636, *Riesz products and the Lonely Runner Conjecture*:
"A wider gap of loneliness", 2025), which handles *exactly* a signed sinc-weighted
lattice/character sum while keeping a **positive** test measure:

- S loose ⟺ the covering `M(τ) = Σ_v 1_{||vτ||≤1/14}` has a zero ⟺ `L>0`.
- Test against a **Riesz product** `R(τ)=∏_{m∈D}(1 + a_m cos 2πmτ) ≥ 0` over a
  **dissociated** set `D` (a probability density, `R̂` supported on the dissociated
  sum-set with **signed** values `(a/2)^{level}`).
- `∫ M·R = Σ_v Σ_k s(k) R̂(vk)`. The main term (`k=0`) is `Σ_v s(0) = 13/7 ≈ 1.857 > 1`;
  the **signed** Riesz coefficients on the relation frequencies must subtract `≥ 0.857`
  to force `∫MR < 1 = ∫R`, which (`M≥1` on any cover) **certifies the lonely set has
  positive measure ⟹ loose**. The dissociated set `D` is tuned to `S`'s relation lattice
  (the additive-energy structure of part B), and Bonami hypercontractivity controls the
  higher levels — exactly the 2025 construction. (Gaussian subordination, Carneiro–
  Littmann, gives the alternate finitization: a band-limited minorant summed over `Λ` via
  Poisson = a positive theta lower bound minus a `1/δ`-per-axis tail.)

This places `inf L>0` in a framework with a *working* precedent on the cousin problem,
and identifies the relation-lattice additive energy (part B) as exactly what the Riesz
product must be adapted to. **Feasibility probe** (this session): the pairing `∫M·R` is
computable, and signed cosine factors do pull it down (`1.857 → 1.41` for hand-built `R`
on the `j=6` extremizer); reaching `< 1` (the certificate) requires the *optimized
dissociated* construction (the 2025 paper's nontrivial core), not naive frequencies —
the concrete next step, HYP-2540.

## D. Closed doors (honest negatives)

- **Polymer-gas / Mayer cluster expansion CANNOT break the wall** (confirmed): standard
  Mayer/Kotecký–Preiss needs *absolutely summable* activities, but `A_3 = ∞` (THM-504).
  The almost-Sidon (`|T|=2`) class is exactly the absolute-convergence boundary.
- **No fermionic/Pfaffian rescue**: the within-level signs are *positively correlated*
  (half-period positivity, THM-504-A), not free/quasi-random, so they do not assemble
  into an antisymmetric kernel; and the project's Pfaffian/determinant objects (THM-174/
  442) are A-affine on the tournament adjacency matrix and **provably spectral** — the
  *wrong axis* (the LRC content is on the speed-relation lattice, off the spectral side
  of the Valiant det/permanent wall). A determinantal handle would need each level's
  signed sinc sum realized as a Pfaffian of a finite alternating kernel on the relations;
  no such structure exists.
- **Second-moment / Paley–Zygmund** applies only *after* positivization (it needs `G≥0`);
  it is the right *finishing* tool once a Riesz product or Gaussian-subordinated minorant
  has produced a positive measure, but not on the raw signed sum.

## Honesty

`inf L > 0` (≡ C'(14) ≡ LRC(14)) is **not proved**; it is genuinely as hard as the cousin
lonely-runner-gap problem. This theorem's content is (A) the exact theta/lonely-measure
identity and its structural positivity, (B) the verified additive-energy predictor
(correcting the λ_1 guess), and (C) the identification + concrete setup of the Riesz-
product method as the right route, with (D) the precise map of which reframes are closed
doors. The proved direction `L>0 ⟹ loose` is THM-501; the converse and the uniform lower
bound remain open.
