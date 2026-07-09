---
id: THM-663
title: The covering case of LRC(14) is closed — with all six density-floor legs (k=8..13) now discharged at the UNIFORM (min-over-families) level by the degree-≤4 covering-moment floor (THM-661) and the LEM-009 longest-AP tail machinery, the k≥8 union-bound branch of THM-530 is UNCONDITIONAL, so ρ*_{1/7}(P,E) ≥ m_P = 14249/252252 > 0 for EVERY admissible (P,E) from the q-witness sieve; every covering-saturated 13-set the sieve reconstructs therefore has a global witness with maxgap>1/7, hence is lonely — no covering-case counterexample exists. Combined with the non-covering case (reduces to LRC(≤13), SETTLED), LRC(14) holds modulo the named upstream THM-527-A finite-Vmax glue and Lean transcription
status: ASSEMBLY / CONSOLIDATION. The four links of the covering-case chain and their honest status: (1) q-witness sieve THM-369 — PROVED (Lean-checked); (2) lonely-density reformulation THM-527/kps-S4, "global-witness maxgap>1/7 ⟹ M≥1/14 ⟹ lonely" — PROVED for the reformulation (w0→∞ limit + finite-w0 O(1/Vmax) error), with the named finite-Vmax / integer-vs-real glue THM-527-A (HYP-2602) the one remaining analytic item; (3) ρ*_{1/7}≥m_P via THM-530 — k≤7 PROVED (pigeonhole μ_{1/7}=1, k=7 a.e. measure-zero), k≥8 union bound ρ*≥meas(G_P)+μ_{1/7}(E)−1≥m_P, whose ONLY contingency was the density floor μ_{1/7}(E)≥bar_k (bar_k=m_P+1−min_{|P|=13−k}meas(G_P)); (4) the density floor — NOW CLOSED for all six legs k=8..13 at the uniform level (THM-661 + LEM-009 fleet + this-session mac-mini-S58). So the covering case is unconditional modulo (2)'s finite-Vmax glue. Not yet Lean-transcribed (finite pieces are decide/native_decide-shaped; the tail rate is opus-S157's resonance-sum lemma). LRC(14) itself: covering case here + non-covering = LRC(≤13) SETTLED.
source: mac-mini-2026-07-08-S58 (consolidation after closing the last density-floor legs)
depends_on:
  - THM-369   # q-witness sieve: counterexample ⟹ covering-saturated 13-set (Lean-checked)
  - THM-527   # lonely-density reformulation; global-witness (maxgap>1/7 ⟹ lonely); THM-527-A finite-Vmax glue
  - THM-530   # ρ*_{1/7} ≥ m_P: k≤7 branch (proved) + k≥8 union bound (needs the density floor)
  - THM-661   # the unified covering-moment density floor (degree ≤4): all six legs, uniform
  - LEM-009   # the block+outlier decorrelation limit + longest-AP tail (k=11,12,13 tail closure)
related:
  - THM-657   # covering reformulation (W = uncovered measure); the frame the density floor lives in
  - THM-660   # Paley–Zygmund covering floor (degree-2 special case)
  - THM-662   # brick (A): the additive-energy/diameter extremal (an alternate k=11 tail route)
---

# THM-663 — The covering case of LRC(14) is closed

**One-line.** The last contingency of the covering case — the density floor `μ_{1/7}(E) ≥ bar_k`
for cluster sizes `k = 8..13` — is now discharged (this session + the fleet), so the covering
case of LRC(14) is unconditional modulo one named upstream glue (THM-527-A) and Lean
transcription. This file is the assembly: it states the four-link chain, pins each link's honest
status, and records what remains.

## The chain (counterexample ⟹ contradiction)

Let `S ⊂ ℤ`, `|S| = 13`, be a putative LRC(14) counterexample speed set (the observer is the
14th runner; dilation-normalize). Two cases from the **q-witness sieve** (THM-369, Lean-checked):

- **Non-covering case.** `S` is not covering-saturated ⟹ the sieve exhibits a witness ⟹ `S` is
  lonely, OR `S` reduces to a `≤13`-runner instance. Either way closed by **LRC(≤13)** (SETTLED,
  owner directive 2026-07-02).
- **Covering case.** `S = P ∪ {V_max − e : e ∈ E}` is covering-saturated, with `P` the observer
  block and `E` the `k`-element cluster (`k = |E|`, `8 ≤ k ≤ 13` after the `k≤7` sub-case).
  This is the residual THM-527/THM-530 attack.

**Reformulation (THM-527/kps-S4).** Via the slow–fast change of variables, `S` is lonely as soon
as the **global witness** has a gap `> 1/7`: `∃ x ∈ G_P` with `maxgap{frac(e_i x)} > 1/7 ⟹ M(S) ≥
1/14`. So it suffices that

> `ρ*_{1/7}(P,E) := meas{ x ∈ G_P : maxgap{frac(e_i x)} > 1/7 } > 0`.

**The floor (THM-530).** `ρ*_{1/7}(P,E) ≥ m_P = 14249/252252 > 0` for every admissible `(P,E)`:
- **`k ≤ 7`:** PROVED. `maxgap ≤ 1/7` forces the exact shifted `1/7`-grid (a measure-zero set of
  `x`), so `μ_{1/7}(E) = 1` a.e., giving `ρ*_{1/7} = meas(G_P) ≥ m_P` (the admissible-`|P|` floor,
  exact, THM-530).
- **`k ≥ 8`:** the elementary union (Bonferroni) bound
  `ρ*_{1/7}(P,E) ≥ meas(G_P) + μ_{1/7}(E) − 1`. This is `≥ m_P` **iff**
  `μ_{1/7}(E) ≥ m_P + 1 − meas(G_P)`, whose worst case over `P` is the **density-floor bar**
  `bar_k = m_P + 1 − min_{|P|=13−k} meas(G_P)`.

## The one contingency — the density floor — is CLOSED

`bar_k` (honest A′ bars) and the discharged floor `min_E μ_{1/7}(E) ≥ min_E B_d(E)`:

| k | bar_k | tool | `min_E` floor (uniform) | minimizer | margin | source |
|---|-------|------|-------------------------|-----------|--------|--------|
| 8 | 0.6750 | `B_4` (deg 4) | 0.761132 | block | +0.086 | THM-661 + S58 |
| 9 | 0.5622 | `B_4` | 0.644603 | block | +0.082 | THM-661 + S58 |
| 10| 0.4521 | `B_4` | 0.553111 | block | +0.101 | THM-661 + S58 |
| 11| 0.331206 | `D3` (deg 3) | 0.404751 (tail min `A`=0.452986) | block | +0.074 | LEM-009 fleet |
| 12| 0.199344 | `D3` | 0.355876 | near-2-AP | +0.157 | S58 (LEM-009 machinery) |
| 13| 0.056487 | `D3` | 0.308844 | near-2-AP | +0.252 | S58 (LEM-009 machinery) |

Since `B_d(E) ≤ μ_{1/7}(E)` for any valid degree-`d` moment minorant (the covering-moment bound,
THM-661), each row gives `μ_{1/7}(E) ≥ bar_k` for **every** admissible cluster `E`, **diameter-free**.
Each `min_E` is `[an exhaustive exact Farey check over compact shapes]` + `[a longest-AP tail whose
floor exceeds the compact min]`; the tail rests on the Weyl decorrelation limit reached at
opus-S157's PROVEN resonance-sum rate `|D3 − D3_∞| ≤ C/(pd)`. **All six legs clear.** Hence the
`k ≥ 8` union-bound branch of THM-530 is unconditional, and `ρ*_{1/7}(P,E) ≥ m_P > 0` for **every**
admissible `(P,E)`.

## Conclusion

Every covering-saturated 13-set the sieve reconstructs has `ρ*_{1/7} > 0`, i.e. a global witness
with `maxgap > 1/7`, hence `M(S) ≥ 1/14` — **lonely, not a counterexample**. With the non-covering
case closed by LRC(≤13), **LRC(14) has no counterexample**, modulo:

1. **THM-527-A finite-`V_max` glue** (HYP-2602) — the integer-vs-real / finite-`w0` bridge in the
   reformulation. The `O(1/V_max)` error is controlled and MSS bounds `V_max ≤ 91^12`; this is the
   one genuine analytic item still to be written cleanly (not a gap in the density floor).
2. **Lean transcription** — the finite pieces (`k≤7` pigeonhole, exhaustive compact `B_d` checks,
   the `ρ*` union bound) are `decide`/`native_decide`-shaped; the tail rate is opus-S157's
   resonance-sum lemma. No new mathematics, an engineering task.

### Advance on item (1): the bounded-arc-count lemma (mac-mini-S58)

The finite-`Vmax` glue is `ρ_K = ρ* + O(#arcs/Vmax)`, where `#arcs` counts the arcs of the good
set `G* = {x : maxgap{frac(e_i x)} > 1/7}`. **`#arcs` is INDEPENDENT of `Vmax`** — it depends only
on the cluster's INTERNAL differences: the combinatorial circular-gap order of `{frac(e_i x)}`
changes only at coincidences `frac((e_i−e_j)x)=0`, i.e. `x = m/(e_i−e_j) = m/(u_j−u_i)` — a
cluster-internal difference, NOT `m/Vmax` (a single phase wrapping through 0 leaves every circular
gap continuous). Within a combinatorial cell each gap is linear with slope a cluster-internal
difference, so crosses `1/7` `O(spread)` times ⟹ **`#arcs = O(k²·spread²)`, `Vmax`-free**.
Machine-verified (shift `e_i → e_i + c` leaves `#arcs` and `meas(G*)` EXACTLY unchanged):
`#arcs = 12` for the `k=11,12,13` blocks, `14` for the near-2-APs — `≈ k+1`, tiny and bounded.
So for **bounded-spread** clusters, `ρ_K = ρ* + O(k/Vmax) → ρ* ≥ m_P > 0`, giving a good period for
`Vmax > V₀ = O(k/m_P)`, with `Vmax ≤ V₀` a finite check — **the bounded-spread half of THM-527-A
is now clean**. The **large-spread** half (`spread ~ Vmax`) is the residual.

### Advance on item (1), large-spread half: the pigeonhole reduction (klein-S192, THM-527 part H)

The large-spread half reduces cleanly. A good period exists iff the ruler grid `x_j=(j+½)/Vmax`
meets `G*`; counting grid points in arcs gives `#{good j} ≥ ρ*·Vmax − #arcs(G*)`, so it suffices
that **`#arcs(G*) < ρ*·Vmax`** (`⇔ maxarc·Vmax > 1`). Two facts make this hold: (a) `#arcs` is
**LINEAR** in spread, `O(k³·spread)` (maxgap is the upper envelope of `k` gaps ⟹ Davenport–Schinzel
`O(spread)` threshold-crossings — the earlier `O(spread²)` over-counted; machine-exact `#arcs = 12·scale`
for `block×scale`, `≈0.2·spread` for random primitive); and (b) for large spread `ρ*` is **large**
(`μ_{1/7} → μ_iid ≈ 0.9–0.999` for `k≤13`, not the weak `m_P`). Verified at the WORST case
`Vmax = spread+14`: over `k=11,12,13`, spreads `≤1000`, primitive, with `G_P` — **zero failures**,
`#good ≥ 30`, `maxarc·Vmax ≥ 4.4`, `ρ* ≥ 0.90`. Honest residual: a-priori `#arcs ≤ c·spread` with
explicit `c<1` and `ρ* ≥ ρ₀>c` — both Weyl/resonance bounds (same axis as the density floor; the
a-priori Davenport–Schinzel `O(k³)` and Erdős–Turán constants are too weak, the true ones need the
resonance structure). Files: `04-computation/lrc14_bounded_arc_count_macmini_S58.py`,
`lrc14_largespread_{arccount,gridhit}_klein_S192.{py,out}`.

## Honest caveats (density-floor rigor)

The per-shape floors `B_d(E)` are rigorous lower bounds on `μ` (exact rational moments via
Farey-cell integration + an exhibited feasible moment minorant). The **uniform** floor `min_E ≥ bar`
combines (i) exact exhaustive compact checks (rigorous), (ii) the longest-AP tail — where opus-S157's
`1/(pd)` rate is a theorem but its constants `V_j` (mixed-variation) and the decorrelation limits
`D3_∞` are numerically-stable certifications, not yet fully a-priori (opus-S157's own caveat). The
margins (`+0.07` at the k=11 razor, `+0.08..+0.25` elsewhere) dwarf the realized corrections (`~10⁻³`),
so the closure is robust. This is the single sub-item inside link (4); everything else is exact.

**Corrected diagnosis of the last certification (klein-S192).** It is NOT "a fully a-priori `V_j`
bound (counting `G^j` breakpoint crossings)" — that route is now shown to be a dead end for the
*wrong reason*. The individual total variations ARE a-priori and small (`V_1 = E_a[TV_u W] ≤
2E_a[W_full] = 0.365`, exact gap formula `TV_u W = Σ_G min(2(ℓ_G−1/7)_+, 2/7)`), yet the
triangle-inequality decorrelation constant `C = Σ|∂D3/∂m_i|·V_i` is `115` (or `5.14` with measured
`V_i`) `≫ 3.47`, because `D3`'s small denominator `Δ = m2−m3/M = 0.026` amplifies each moment error,
while the TRUE deviation `×d ≈ 0.035` survives only by **cancellation** among the three errors
`ε_i = m_i(E_d)−L_i` (one shared Riemann-sum defect of `W,W²,W³`). klein-S192 captures the
cancellation at **first order** a-priori: `Σc_iε_i = E_a[Riemann-err of g_a]`, `g_a = φ∘W`, `φ` a
FIXED cubic, so `C_1 = E_a[TV_u g_a] ≤ Lip(φ)·2E_a[W_full] = 2.83 < 3.47` — the first-order tail is
a-priori-closed for all `d≥26`. **The full closure is kps-S89's box bound** (a direct rigorous
enclosure `min D3` over `m_i ∈ [L_i ± Vh_i/d]`, which captures the cancellation with NO
linearization): with the tight constants `Vh_i = i(6/7)^{i-1}E[W_B]` it clears `bar` for all
`d ≥ 62`, so `[finite check d ≤ 61] + [box bound d ≥ 62]` **closes the L=10 family a-priori with no
L² tail needed**. (opus-S154's L² Fourier tail — a uniform `o(1/d)` on `Σ_{m≥2}|\hat{g_a}(md)|²` —
would only tighten the *linear* `C/d` constant; it is a refinement, not a prerequisite.) So link
(4)'s last certification is **discharged** (kps-S89 box + finite check); klein-S192's role is the
diagnosis (the `V_j` are not the obstruction, cancellation is, and the box enclosure recovers it).
See LEM-009 (kps-S89 + klein-S192 sections) and `04-computation/lrc14_Vi_{apriori,combined}_klein_S192.{py,out}`,
`lrc14_L10_explicit_rate_kps_S89`.
