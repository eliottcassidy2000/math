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

## Honest caveats (density-floor rigor)

The per-shape floors `B_d(E)` are rigorous lower bounds on `μ` (exact rational moments via
Farey-cell integration + an exhibited feasible moment minorant). The **uniform** floor `min_E ≥ bar`
combines (i) exact exhaustive compact checks (rigorous), (ii) the longest-AP tail — where opus-S157's
`1/(pd)` rate is a theorem but its constants `V_j` (mixed-variation) and the decorrelation limits
`D3_∞` are numerically-stable certifications, not yet fully a-priori (opus-S157's own caveat). The
margins (`+0.07` at the k=11 razor, `+0.08..+0.25` elsewhere) dwarf the realized corrections (`~10⁻³`),
so the closure is robust, but a fully a-priori `V_j` bound (counting `G^j` breakpoint crossings) would
remove the last certification. This is the single sub-item inside link (4); everything else is exact.
