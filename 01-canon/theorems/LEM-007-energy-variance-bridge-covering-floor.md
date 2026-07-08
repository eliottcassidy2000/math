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

This SUGGESTED an explicit-c leading order, but the follow-up support-decomposition check (CORRECTION below) shows it is NOT a valid bound (the higher-support terms grow, not decay). The genuine (i) bounding the off-diagonal pair resonances (difference-multiplicity, THM-641/638) and
(ii) the sign of the support-3,4 correction (it REDUCES Var, verified `support-2/Var(W) ≈ 1.5–1.7 > 1`,
so `Var(W) ≤ (6/7)^{2(k−2)}·`(pair kernel) holds on the sampled zoo). Files:
`04-computation/lrc14_overlap_fourier_law_macmini_S57.py`.

## CORRECTION (mac-mini-S57, same session): the support decomposition does NOT truncate

Working the "two remaining pieces" exposed an error in the leading-order claim above. Computing
the support-EXACTLY-r Fourier parts `Ŵ_r(m)` (all r entries nonzero) shows `Σ_m|Ŵ_r(m)|²` GROWS
with r, not decays: `0.077, 0.226, 0.932, …` (block, r=2,3,4), and the truncated `Σ_m|Ŵ_2+Ŵ_3+
Ŵ_4|²` OVERSHOOTS the true `Var(W)=0.047` by 12×. So the higher-support terms are large and the
`Σ_r Ŵ_r(m)` cancellation is EXTREME and across ALL orders — NOT a "support-2 leading + small
correction". Consequently:
- The overlap Fourier mass laws (`L̂_S(m) = Σ_{balanced} Π ĥ`) and the master law
  (`Ŵ(m) = Σ_{balanced n} Π ĝ`) are CORRECT and verified.
- But **`c_k ≈ (6/7)^{2(k−2)}·V1_φ` is NOT a rigorous leading order** — it is the value of the
  support-2 part `Σ_m|Ŵ_2|²` (≈ 1.6× the true Var), a numerical correlate, not a bound: the
  claim that support-3,4 merely "trim" it is FALSE (they are large and require extreme
  cancellation). The `(6/7)^{k−2}` damping is real per-term but does not dominate.
- **`Var(W) ≤ c·R2` is therefore NOT reachable by support truncation** of the master law — the
  balanced-vector sum has the same non-truncating cancellation as the coverage inclusion-exclusion
  (the barely-covers wall, `k/7 > 1`). The honest rigorous route is klein-S179b's LEM-005
  DISCREPANCY estimate (phase-vector equidistribution), whose explicit uniform rate is the
  remaining analytic mile. The Fourier master law is exact but does not tame the cancellation.

What SURVIVES rigorously: the equivalences (`far ≤ E[W]² ⟺ Var ≤ near`, `Var(W) = Σ|Ŵ(m)|²`),
the exact overlap Fourier laws (a correct tool for the discrepancy computation), and D3 decreasing
in m2. The resonance lemma `Var(W) ≤ c·R2` remains open, its difficulty correctly located.

## The two remaining pieces, worked (mac-mini-S57): piece 1 proved, piece 2 collapses

**Piece 1 (off-diagonal pair-resonance kernel) — PROVED (in isolation).** The support-2
contribution is `Σ_m|Ŵ_2(m)|² = (6/7)^{2(k-2)} Σ_{P,P'} K(d_P,d_{P'})`,
`K(d,d') = Σ_{s≠0} φ̂((d'/g)s)φ̂((d/g)s)`, `g = gcd(d,d')`, `a=d/g, b=d'/g` coprime. Since
`φ̂(t) = sin²(πt/7)/(π²t²)` and `Σ_s sin⁴(πcs/7)/(π⁴c⁴s⁴) ≤ (2ζ(4))/(π⁴c⁴) = 1/(45c⁴)`,
Cauchy–Schwarz gives **`K(d,d') ≤ 1/(45 a²b²)`**. Then, writing `R_d = ` difference
multiplicity, AM–GM + `Σ_b 1/b² = π²/6` + `Σ_g R_{ga}² ≤ R2/2`:
`Σ_{P,P'} R_dR_{d'} K ≤ (π⁴/3240)·R2`. So **`Σ_m|Ŵ_2(m)|² ≤ (π⁴/3240)(6/7)^{2(k-2)}·R2`** —
a rigorous `R2`-linear bound on the pair part (verified: `0.079 ≤ 1.44` block; `0.014 ≤ 0.29`
spread). This is a clean rigorous result, but it bounds ONLY the support-2 part.

**Piece 2 (support-3,4 cancellation sign) — the hoped-for reduction is FALSE.** The support-
exactly-r Fourier parts have `Σ_m|Ŵ_r|² = 0.077, 0.226, 0.932, …` (block, r=2,3,4) — GROWING,
and the truncated `Σ_m|Σ_{r≤4}Ŵ_r|²` OVERSHOOTS the true `Var(W)=0.047` by **12×**. So the
higher-support terms are large, and `Var(W) ≤ Σ_m|Ŵ_2|²` (which holds numerically, `0.047 <
0.077`) is NOT a "small correction" — it requires extreme cancellation across all orders. The
master-law/support route therefore CANNOT yield `Var(W) ≤ c·R2` by truncation: the pair bound
(piece 1) does not control the exploding tail. This is the barely-covers cancellation in full.

**Net.** Piece 1 gives a rigorous bound on the pair-resonance kernel; piece 2 shows the pair
kernel does not dominate `Var(W)`, so the Fourier-truncation route to `Var ≤ c·R2` is ruled out.
The genuinely rigorous path stays klein-S179b LEM-005 (phase-vector DISCREPANCY), whose uniform
explicit rate is the one remaining analytic mile. The overlap Fourier laws derived here are the
correct tool for THAT computation (they give the exact `far` integrand), not for a truncation.

## The TARGET is per-shape at HIGH R2, not a uniform c·R2 (klein-S183) — agrees the uniform route is dead

mac-mini's correction above rules out the *method* (Fourier truncation) for a uniform `Var(W) ≤
c·R2`. Independently, the *target* is ill-posed: `c = Var(W)/R2` is NON-uniform (sampled
5.0–8.7·10⁻⁵ over `R2 ≤ 614`, PEAKING at LOW R2 — the mid-spread `{0,2,5,9,14,20,27,35,44,54,65}`
has `R2 = 170, c = 8.7·10⁻⁵ > 6.45·10⁻⁵`). So NO uniform `c ≤ 6.45·10⁻⁵` exists; both the method
and the target for "uniform `c·R2`" are dead. But **brick (B) never needed it.** The per-shape PZ
condition is

> **`Var(W) ≤ (1−bar)/bar · E[W]² = 2.02·E[W]²`**  (i.e. `PZ = E[W]²/(Var+E[W]²) ≥ bar`),

BINDING at HIGH `R2` (near-block / far-point, `Var(W)` largest absolutely); LOW-`R2` high-`c` shapes
clear comfortably (small absolute `Var`, high `E[W]`, spread ⟹ lonely). The whole tail
(`prim-diam ≥ 16`, `R2 ≤ 614`) clears: far-point `PZ = 0.390` (`+0.059`), the binding tail shape;
block `PZ = 0.347` the global min (compact prim-diam 10, exhaustive-covered).

**Consequence.** This LOCALIZES the remaining analytic mile: the LEM-005 discrepancy estimate
(the surviving rigorous route, per mac-mini) is needed ONLY at HIGH `R2` (the far-point family),
not uniformly — and there the far-point clears by `+0.059`. The per-shape framing converts "prove
a uniform `c`" (impossible) into "the high-`R2` tail satisfies `Var ≤ 2.02·E[W]²`," a bounded
family where mac-mini's exact overlap-Fourier `far` integrand + the LEM-005 discrepancy apply.
Files: `lrc14_support_reduction_klein_S183.out`, `lrc14_brickB_reframe_klein_S183.out`.

## Pushing the discrepancy rate for far ≤ E[W]² (mac-mini-S57): structure + leading term

Target: an EXPLICIT uniform rate for `far → far_iid` (klein LEM-005's open mile). Via the master
law applied to `B = A₁ ∪ A₂` (measure 2ℓ): `q₂ = P(A₁∪A₂ empty) = Σ_{(n_i): Σ n_i e_i = 0} Π ĝ_B(n_i)`,
`ĝ_B(0) = 1−2ℓ = 5/7`. Integrating over disjoint `(y₁,y₂)`:

> **`far = far_iid + far_dev`**, `far_iid = (5/7)^{k+1}` (the `n=0` term), `far_dev = Σ_{(n)≠0: Σn·e=0}
> ∫∫_disjoint Π ĝ_B(n_i)` — the RESONANCE sum (nontrivial integer relations `Σn_i e_i = 0`).

**The leading deviation is PAIRWISE and difference-based.** The `|S|=2` part of `far_dev` is
`(5/7)^{k−2} Σ_{i<j} ∫∫_disjoint Cov_x(f_i, f_j)`, `f_i = (1−χ_i(A₁))(1−χ_i(A₂))`, and
`Cov_x(f_i,f_j)` is a sum of pair joint-window masses at difference `d_{ij} = e_i − e_j` — EXACTLY
klein's THM-638 / my THM-641 offset-pair law, which deviates from `ℓ²` by `O(1/d_{ij})` (Koksma).
So the leading rate is `|far_dev| ≲ (5/7)^{k−2} Σ_{i<j} |dev(d_{ij})|`, DIFFERENCE-based hence
shift-invariant, and `→ 0` as the differences grow (spread). This is the correct form of the
"discrepancy rate": it is the sum of THM-638/641 pair-mass deviations, damped by `(5/7)^{k−2}`.

**Correction (self-caught):** an element-based weight `Σ gcd(e_i,e_j)²/(e_i e_j)` I first tried is
WRONG — it is not shift-invariant and fails for `0`-containing families (e.g. `{0,7,17,…,115}`:
tiny weight `0.12` but `far_dev = 0.009`). The `0` phase is deterministically fixed, an artefact
removed by translation; the genuine rate is the DIFFERENCE-based pair-mass sum above.

**Status.** The exact resonance decomposition and the leading pairwise (difference-based,
Koksma-decaying, THM-638/641) rate are established — this is the correct scaffold for LEM-005's
explicit uniform rate. The remaining piece is the higher-support tail (`|S| ≥ 3` resonances
`Σ n_i e_i = 0`), which shares the non-truncation concern of the variance; but the LEADING rate is
now explicit and difference-based, reducing the mile to a tail bound over the higher resonances.

## Bounding the higher-support resonance tail (mac-mini-S57)

The tail is `Σ_{r≥3} T_r`, `T_r = Σ_{(n): |supp|=r, Σn_i e_i=0} ∫∫_disjoint Π ĝ_B(n_i)`. Two facts:

**(a) Per-term geometric decay (rigorous).** Replacing a zero entry (`ĝ_B(0)=5/7`) by a nonzero
one (`|ĝ_B(n)| ≤ 2/(π|n|) ≤ 2/π`) costs a factor `(2/π)/(5/7) = 14/(5π) = 0.8913 < 1`. So each term
of support `r` obeys `|·| ≤ (5/7)^{k+1}·(0.8913)^r/Π|n_i|` — support is geometrically penalized
(UNLIKE the variance's balanced sum). The generating function `Σ_{Σn·e=0} Π w(n_i) = ∫₀¹ Π_i
W(θe_i) dθ`, `W(φ)=1−1.783·ln(2 sin π|φ|)`, packages the whole tail in one integral.

**(b) The absolute-value bound loses the cancellation (the obstruction).** Because `T_r`'s terms
carry signs (`ĝ_B` complex, `y`-dependent), the geometric per-term decay does NOT give a summable
bound: the resonance COUNT grows, and `∫₀¹ Π_i |W(θe_i)| dθ` DIVERGES-ish (`|W|` has a `ln`
singularity at `φ=0`, and near `θ=0` all `θe_i→0` simultaneously, giving `∫(ln)^k`). So the tail,
like the variance, needs the SIGNED cancellation — the `T_r` are individually large and cancel.
The exact higher-support masses (the triple/quad `L̂_S` laws derived above, WITH their signs) are
what carry that cancellation; the absolute route is the same barely-covers wall.

**(c) Clean sub-result — dissociated (B_h) families.** If the speed set is `B_h` (no nontrivial
`Σ n_i e_i = 0` with `Σ|n_i| ≤ h`), then ALL resonances have `Σ|n_i| > h`, so every tail term has
`Π|n_i| ≥` (a growing bound) and the leading resonance is pushed to high `|n|`: `|far_dev| ≤
(5/7)^{k+1}·C_h` with `C_h → 0` as `h → ∞`. So for STRONGLY dissociated spread families the tail
(and `far_dev`) is provably small — a clean rate on a genuine (if restricted) spread class.

**Net.** The tail has rigorous geometric per-term decay (a) and is provably small on dissociated
families (c); the uniform bound over ALL spread families (b) shares the variance's signed-
cancellation wall — the exact higher-support overlap laws (derived above) are the tool, but the
uniform sign control is the same one remaining mile. The leading (pairwise, THM-638/641) rate plus
the dissociated tail bound close `far ≤ E[W]²` for strongly-dissociated spread families explicitly.

## Pushing dissociated → general spread: it IS the variance wall (mac-mini-S57, rigorous)

Trying to extend the dissociated `far ≤ E[W]²` bound to general spread yielded a clean rigorous
STRUCTURE theorem (and corrected two earlier claims):

**far_dev is supported ENTIRELY on DOUBLY-BALANCED resonances.** The disjoint-arc y-integral of a
resonance term is `f(P_T,Q_T)`, nonzero only when `P_T+Q_T = 0`, and `P_T+Q_T = Σ_i n_i` (the total
of the frequency vector). Hence a resonance `Σ n_i e_i = 0` contributes to `far_dev` ONLY IF also
`Σ n_i = 0` (numerically verified: `f(1,1)=f(2,-1)=0`, `f(t,-t)≠0`). So

> **`far_dev = Σ_{n≠0: Σ n_i = 0 AND Σ n_i e_i = 0} ∫∫_disjoint Π ĝ_B(n_i)`** — the DOUBLY-BALANCED
> resonances, which are EXACTLY the `Var(W) = Σ_m|Ŵ(m)|²` structure (`Ŵ(m)` at the resonant `m`).

So `far ≤ E[W]²` and `Var(W) ≤ near` are not just equivalent as scalars — they are built from the
SAME doubly-balanced resonance sum. **Pushing the dissociated bound to general spread is therefore
IDENTICAL to resolving the variance cancellation** — not a separable task. Dissociation is special
precisely because it kills all doubly-balanced resonances (no `Σn_i=0 ∧ Σn_i e_i=0` with small `n`).

**CORRECTIONS to my earlier claims (both self-caught here):**
- The "leading deviation is PAIRWISE (support-2, THM-638/641)" is WRONG: a support-2 resonance
  `n=(t,-t)` on a pair has `Σ n_i e_i = t(e_i−e_j) ≠ 0` (not a resonance) OR the frequency-pair
  `(te_j/g, −te_i/g)` has `Σ n_i ≠ 0` (unbalanced ⟹ y-integral 0). Either way support-2 contributes
  ZERO. The LEADING `far_dev` term is SUPPORT-3: the 3-term-AP pattern `n=(1,−2,1)` on triples
  `(e_i, e_j, e_l)` with `e_i − 2e_j + e_l = 0`. So `far_dev ≈ (5/7)^{k−3}·(weighted 3-AP count) +
  higher` — governed by the family's ADDITIVE (3-AP / additive-energy) structure, consistent with
  the `R2` link. Dissociated ⟺ few 3-APs ⟺ `far_dev` tiny.
- The element-based `Σ gcd²/(e_i e_j)` (earlier) and the difference-based pairwise rate are both
  superseded: the correct object is the doubly-balanced (3-AP-and-up) resonance sum.

**Net.** The general spread case is RIGOROUSLY the variance wall (doubly-balanced resonances,
leading = 3-APs); the dissociated bound is the resonance-free special case. This unifies the two
open miles into ONE (attack together) and pinpoints the leading obstruction as the 3-AP structure.
The overlap Fourier laws + the doubly-balanced characterization are the exact tools; the signed
cancellation of the 3-AP-and-higher balanced sum is the single remaining analytic mile (= the
variance).

## Attacking the doubly-balanced 3-AP cancellation (mac-mini-S57): it needs a 2D resummation

Directly attacking the doubly-balanced resonance sum (the unified far/variance target):

**Support-truncation FAILS here too (verified).** Computing the doubly-balanced far_dev by support
(analytic y-integral: contribution `= (5/7)^{k-r} Π(−ârc(n_i)) Σ_{T⊆supp} c_{P_T}`, `c_P = −sin(2πPℓ)/(πP)`,
`c_0 = 5/7`): the terms ALTERNATE in sign and GROW then decay — `(s3,s4,s5) = (−0.0025,+0.0047,−0.0047)`
spread, `(−0.014,+0.050,−0.033)` block; `|s4/s3| = 1.9–5.8 > 1`. The partial sum overshoots the tiny
true far_dev (`0.0002` spread), so support-4-and-up are essential and cancel. The doubly-balanced
localization does NOT make the sum truncate — the 3-AP (support-3) leading term is smaller than
support-4, and the whole thing only converges by inter-support cancellation. Same wall.

**The correct handle is a 2D resummation.** The doubly-balanced sum (two constraints `Σn_i=0`,
`Σn_i e_i=0`) is captured EXACTLY by a DOUBLE generating integral
`far + far_iid-type = ∫₀¹∫₀¹ Π_i W(θ + φ e_i) dθ dφ` (θ enforces `Σn_i=0`, φ enforces `Σn_i e_i=0`;
`W` the per-phase weight with the `c_P` y-kernel folded in). The `(θ,φ)` map fills a 2-DIM subset of
the k-torus, so `Π W` averages toward `far_iid` and the deviation `far_dev` is the 2D DISCREPANCY of
`(θ,φ) ↦ (θ + φ e_i)_i` — more equidistributed (2-parameter) than the 1-D `x`-average, hence a
genuinely cleaner discrepancy object. This double integral carries the cancellation the support
sum cannot.

**Net (attack outcome).** The doubly-balanced 3-AP cancellation does NOT yield to support truncation
(verified: alternating, non-truncating, support-4 dominant) — it is the same essential signed
cancellation. Its correct resummation is the 2D generating integral `∫∫ Π_i W(θ+φe_i)`, whose
deviation from `far_iid` is a 2-parameter discrepancy → 0 as the speed set spreads. This is the
precise, correct form of the one remaining mile (sharper than the 1-D `x` discrepancy of klein's
LEM-005): a 2D equidistribution estimate for the lattice map `(θ,φ) ↦ (θ + φ e_i)`. NOT closed;
fully mapped.

## CORRECTION: the "2D discrepancy resummation" is CIRCULAR (mac-mini-S57)

Attempting to PROVE the proposed 2D discrepancy bound exposed that it does not exist as a
simplification. Introducing duals `θ` (for `Σn_i=0`) and `φ` (for `Σn_i e_i = m`) and resumming
gives `Var(W) + E[W]² = ∫₀¹∫₀¹∫₀¹ Π_i V(θ₁+φe_i) V(θ₂−φe_i) dθ₁ dθ₂ dφ`, `V = 1 − 1_arc`. But the
`θ₁, θ₂` integrals each collapse to `g(±φ) = W(∓φ)`, and `∫₀¹ W(φ) W(−φ) dφ = ∫₀¹ W² = E[W²]`
(verified: `0.040428 = 0.040428` exactly). So the "3D integral" is just **Parseval** — it returns
`E[W²]` with NO new tractable object. The `(θ,φ)`-map "2D discrepancy" of the prior section is
therefore NOT a genuine reduction; that claim is RETRACTED.

**Honest final state of this thread.** The doubly-balanced 3-AP cancellation is the barely-covers
wall, and it resists every reformulation tried: Bonferroni (diverges), support truncation of the
variance (grows), support truncation of the doubly-balanced sum (alternates, support-4 dominant,
non-truncating), absolute-value bounds (lose cancellation, diverge), and now the 2D resummation
(circular = Parseval). The cancellation is genuinely essential and does not admit a low-complexity
closed form. The one route that sidesteps it is klein-S179b LEM-005's DIRECT 1-parameter estimate
(`far → (5/7)^{k+1}` as the phase vector `{e_i x}` equidistributes), which gives the CROSSOVER
(`far ≤ E[W]²` for prim-diam ≥ 36 at k=11, numerically) but not a proven uniform rate. A full
uniform rate is a genuine open analytic problem (a k-fold exponential-sum/discrepancy estimate in
the barely-covers regime `k/7 > 1`), not reducible to the pairwise or 2D shortcuts.

What is RIGOROUS and survives from this thread: the equivalences (`far ≤ E[W]² ⟺ Var(W) ≤ near`;
`Var(W) = Σ|Ŵ(m)|²`); the exact S-arc overlap Fourier mass laws (verified); the proof that far_dev
is supported ONLY on doubly-balanced resonances (leading = 3-APs, support-2 contributes 0); the
geometric per-term decay `(2/π)/(5/7) = 0.891`; and the dissociated (`B_h`) explicit bound. The
uniform spread case remains the barely-covers wall — now mapped exhaustively and honestly.
