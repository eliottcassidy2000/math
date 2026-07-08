---
id: LEM-007 (RENUMBERED from LEM-006 by mac-mini-S57: klein-S180 factorial-moment LEM-006 pushed 10:37 < mine 11:14; collision protocol)
title: The energy–variance bridge for the covering floor — far ≤ E[W]² ⟺ Var(W) ≤ near, and Var(W) = Σ_{m≠0}|Ŵ(m)|² is governed by the reduced additive energy R2 = Σ_{d≠0} r_d² of the speed set (empirically Var(W) ≈ 0.000056·R2, 30% band), the SAME R2 that controls the tent-variance floor (THM-656). This reduces BOTH the tail decorrelation (far ≤ E[W]² for spread) AND brick (B) (R2 ≤ 614 ⟹ D3 ≥ bar) to one covering-energy inequality Var(W) ≤ c·R2; and D3 is provably decreasing in the second moment m2 (the clean half of brick B). 614 = the additive energy of the diam-16 extremizer — a value of the universal spread parameter, not a new constant
status: PARTIAL — the two REDUCTIONS and the equivalences are rigorous; the core inequality Var(W) ≤ c·R2 (with the true small c ≈ 5.6e-5, 27× below the pair-part coefficient V1_φ = 1.53e-3 because of inclusion–exclusion cancellation) is exactly klein-S179b's LEM-005 discrepancy estimate — the barely-covers wall — NOT independently proved here. Machine-verified: the far⟺Var⟺near equivalence (0 violations); Var(W)/R2 ∈ [4.9e-5, 7.0e-5]; min D3 over k=11 primitives with R2 ≤ 614 = 0.4933 ≥ bar 0.3312 (+0.162); ∂D3/∂m2 < 0 exactly.
source: mac-mini-2026-07-08-S57 (HYP-5357 program)
depends_on:
  - THM-657   # W = uncovered measure, near/far split
  - THM-660   # PZ covering floor; D3 = its degree-3 sharpening
  - THM-656   # tent-variance / additive-energy floor — SAME R2 = Σ r_d²
related:
  - LEM-005   # klein's near/far decorrelation = the core Var(W) ≤ c·R2 via discrepancy
external: Parseval; the additive-energy method (Gowers); circle covering (Stevens).
---

# LEM-007 — the energy–variance bridge

## Three rigorous reductions (the useful content)

Let `W(x)` be the uncovered measure (THM-657), `E[W²] = near + far` the near/far split (LEM-005),
`m_i = E[W^i]`, `M = 6/7`. Then:

**(1) far ≤ E[W]² ⟺ Var(W) ≤ near.** *Proof.* `far = E[W²] − near = m2 − near`; `E[W]² = m1²`;
so `far ≤ E[W]² ⟺ m2 − near ≤ m1² ⟺ Var(W) = m2 − m1² ≤ near`. ∎ (Verified: 0 violations.)
This replaces the analytic "disjoint-arc decorrelation" with a single scalar inequality on Var.

**(2) Var(W) = Σ_{m≠0} |Ŵ(m)|²** (Parseval), and `Ŵ(m) = Σ_{i<j: d_{ij}|m} φ̂(m/d_{ij}) −
Σ_{triples} … ` (coverage inclusion–exclusion; `φ = ` arc-autocorrelation triangle, width 1/7,
`φ̂(n) = sin²(πn/7)/(π²n²)`). The PAIR term alone gives `Var_pair = R2·V1_φ`,
`V1_φ = ∫φ² − (∫φ)² = 2/(3·343) − 1/49² = 1.527e-3`; but the triple+ terms cancel ~96% of it
(`Var(W)_full ≈ 0.037·Var_pair`), leaving the **covering-energy relation
`Var(W) ≈ c·R2`, `c ≈ 5.6e-5`** (band `[4.9e-5, 7.0e-5]`). So the additive energy
`R2 = Σ_{d≠0} r_d²` is the spread parameter for the covering floor — the SAME `R2` as the
tent-variance floor THM-656 (`Var(F) = R2·V1`), now running through `k ≥ 11` via the covering `W`.

**(3) D3 is decreasing in `m2` (the clean half of brick B).** With `m1, m3` fixed,
`∂D3/∂m2 = (m1 − m2/M)·[−(2/M)(m2 − m3/M) − (m1 − m2/M)]/(m2 − m3/M)² < 0`
(both bracket terms `< 0` since `m1 > m2/M` and `m2 > m3/M`). So higher variance ⇒ lower D3;
the D3-minimizer is the max-variance = max-energy family.

## Consequences (both target lemmas, reduced to ONE inequality)

- **Tail decorrelation (far ≤ E[W]² for spread):** by (1), ⟺ `Var(W) ≤ near`; by (2),
  `Var(W) ≈ c·R2`, and `near ≈ 0.015–0.023` (exact), so it holds whenever `R2 ≤ near/c ≈ 270–410`
  — i.e. for SPREAD (low-energy) families. Spread ⟺ low `R2` ⟺ `Var(W)` small ⟺ `far ≤ E[W]²`.
- **Brick (B) (R2 ≤ 614 ⟹ D3 ≥ bar, k=11):** by (3) `min D3 = D3(max Var) = D3(max R2)`; with
  `Var ≤ c·R2 ≤ c·614`, `D3 ≥ D3(Var = c·614) ≥ bar`. Verified margin: `min D3` over `R2 ≤ 614`
  primitives `= 0.4933 ≥ 0.3312` (`+0.162`).

Both reduce to the SINGLE covering-energy inequality **`Var(W) ≤ c·R2`** (`c ≈ 5.6e-5`).

## The 614 connection (what "look back" surfaces)

`614` is not a new constant: it is the value of the **universal spread parameter `R2`** at the
k=11 diameter-16 extremizer (the 1+10 split `{0,..,9,16}`). The additive energy `R2 = Σ r_d²`
is the ONE quantity governing the spread-side of the ENTIRE density floor:

| k | max R2 (block) `= k(k−1)(2k−1)/3` | role of R2 |
|---|---|---|
| 8 | 280 | THM-656 tent-variance; `R2* = 385 > 280` ⇒ every 8-set |
| 9 | 408 | THM-656; discharge for `R2 ≤ 217` (spread) |
| 10 | 570 | THM-656 tent |
| 11 | 770 | D3 covering; `maxR2@diam≥16 = 614` |
| 12 | 1012 | D3; `852` |
| 13 | 1300 | D3; `1148` |

So the tent floor (k≤10) and the D3 covering floor (k≥11) are ONE spread-side story in the single
variable `R2`; `614` is where it sits for the binding k=11 tail. This is the same additive energy
that appears in the metagraph 2nd moment (THM-589 W(n)) and the tournament-side variance — the
recurring "energy is the spread axis" motif.

## Honest scope

The core `Var(W) ≤ c·R2` with the true (cancellation-reduced) `c ≈ 5.6e-5` is exactly the
inclusion–exclusion cancellation at the barely-covers regime (`k/7 = 1.86 > 1`) — klein-S179b's
LEM-005 discrepancy estimate. The pair term overcounts 27×; deriving the true `c` needs the full
cancellation (LEM-005). What is NEW and rigorous here: the three reductions above, which turn "two
analytic lemmas" (far, brick B) into ONE scalar covering-energy inequality and expose `R2` as the
unifying spread parameter (the 614 connection). File:
`04-computation/lrc14_energy_variance_bridge_macmini_S57.py`.

## The overlap Fourier mass law + the resonance mechanism (mac-mini-S57, HYP-5357)

Deriving the resonance lemma's ingredients via THM-641's method (kps-S81's flagged target):
expand the uncovered indicator `1{y uncov} = Π_i (1 − h(y − e_i x))`, `h = 1_{[0,1/7)}`,
`ĥ_n = (1 − e(−n/7))/(2πin)`, `ĥ_0 = 1/7`. The y-integration forces `Σ n_i = 0` (balanced),
the x-character forces `Σ n_i e_i = m`:

> **`Ŵ(m) = Σ_{(n_i): Σn_i=0, Σn_i e_i=m} Π_i ĝ_{n_i}`**, `ĝ_0 = 1−1/7 = 6/7`, `ĝ_n = −ĥ_n` (n≠0).

Restricting to support `S` gives the **S-arc overlap Fourier law**
`L̂_S(m) = Σ_{balanced on S} Π_{s∈S} ĥ_{n_s}`; the **triple law** is the balanced sum over
`n_i d_{ik} + n_j d_{jk} = m` (VERIFIED exactly vs numeric Fourier, `E[L_{ijk}] = L̂(0)`).

**The resonance = the `(6/7)^{k−2}` inactive-arc damping.** In `Ŵ(m)`, the `k − |support(n)|`
zero entries each carry `ĝ_0 = 6/7`. So the leading (support-2, pair) term is
`Ŵ(m) ⊇ (6/7)^{k−2} Σ_{i<j: d_{ij}|m} φ̂(m/d_{ij})` (`φ̂ = |ĥ|² = ` arc autocorrelation), and its
variance contribution is `(6/7)^{2(k−2)} ·` (pair-resonance kernel) `≈ (6/7)^{2(k−2)}·R2·V1_φ`.
This damping is the BULK of the 96% cancellation: naive `R2·V1_φ = 1.18` (block) → damped
`(6/7)^{18}·1.18 = 0.073` (a 15× cut) → true `Var(W) = 0.047` (higher-support terms trim the
last 1.5×). So the explicit leading coefficient is

> **`c_k ≈ (6/7)^{2(k−2)}·V1_φ`**  (`= 9.5×10⁻⁵` at k=11; empirical `5.6×10⁻⁵` after the
> support-3,4 correction), `V1_φ = 2/(3·343) − 1/49² = 1.527×10⁻³`.

This turns "Var(W) ≤ c·R2 with mysterious c" into an explicit-c statement whose ONLY remaining
gap is (i) bounding the off-diagonal pair resonances (difference-multiplicity, THM-641/638) and
(ii) the sign of the support-3,4 correction (it REDUCES Var, verified `support-2/Var(W) ≈ 1.5–1.7 > 1`,
so `Var(W) ≤ (6/7)^{2(k−2)}·`(pair kernel) holds on the sampled zoo). Files:
`04-computation/lrc14_overlap_fourier_law_macmini_S57.py`.
