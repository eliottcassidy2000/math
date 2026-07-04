---
id: THM-614
title: Partial multi-tightener confinement (THM-612 open gap). For a primitive tight family with q*=28 (m=2) and f=2 tighteners, the EXTREMITY LEMMA holds — on the U-loose region R, exactly one odd tightener is ≤1/14 and the other ≥3/7 — forcing every component of R to be single-tightener and giving the MAGNITUDE BOUND w_i ≤ 12·u_max (both tighteners bounded by the even part's scale). General f: f ≥ 2 and f ≥ 7·meas(R). These CONSTRAIN the multi-tightener case but do NOT close it — the residual is the global tightness (M(S)≤1/14 off R) + primitivity, not captured by the covering-of-R shift argument. Confinement (primitive tight ⟹ q*=14) remains a CONJECTURE.
status: PARTIAL — the Extremity Lemma, the single-tightener component structure, the magnitude bound w_i ≤ 12 u_max (f=2), and f ≥ max(2, 7·meas R) are PROVED (elementary; verified). Full confinement is NOT proved. Independent search (938 structured even-block+odd-tightener families, e=10,11,12) finds 0 primitive tight with q*>14, consistent with THM-612.
source: opus-2026-07-03-S61
depends_on:
  - THM-612   # the tower + Lemma C (f=1); this extends the shift obstruction toward f≥2
  - THM-610   # tight ⟹ 14|q* (so q*=14m)
  - LRCUpTo13 # M(U) ≥ 1/(e+1): forces U loose and bounds the magnitude
related:
  - THM-613   # the same Lipschitz slope idea powers the magnitude bound here
  - HYP-4062  # kps: tight locus rigidity (confinement ⟹ the finite mod-14 problem)
  - HYP-2913  # g(14)≤3, the other open gap
results:
  - 04-computation/lrc14_confinement_setup_opus_S61.py
  - 05-knowledge/results/lrc14_confinement_setup_opus_S61.out
---

# THM-614 — partial multi-tightener confinement

**Setup (THM-612, m=2).** `S` primitive, tight (`M(S)=1/14`), `q*=28`. `E` = even runners `= 2U`,
`F` = odd tighteners, `f=|F|`, `e=|E|=13−f`. `R = {t : g_U(2t) > 1/14}`, `g_U(τ)=min_{u∈U}‖uτ‖`.
On `R` the `E`-runners are strictly safe, so tightness forces `R ⊆ ⋃_{w∈F} D_w`, `D_w={‖wt‖≤1/14}`.
Two facts: `R` is `(+1/2)`-invariant (`g_U(2(t+1/2))=g_U(2t+1)=g_U(2t)`); and for odd `w`,
`‖w(t+1/2)‖ = 1/2 − ‖wt‖` (so `D_w ∩ (D_w+1/2) = ∅`) — a single tightener covers ≤ one of each pair
`{t, t+1/2}`, whence **`f ≥ 2`** (reproving Lemma C) and, by the measure/union bound, **`f ≥ 7·meas(R)`**.

## Extremity Lemma (f=2, PROVED)
> For every `t ∈ R`: **exactly one of `‖w_1t‖, ‖w_2t‖` is `≤ 1/14`, and the other is `≥ 3/7`.**

**Proof.** `t` and `t+1/2` both lie in `R`, so both are covered. If `‖w_1t‖≤1/14`: then covering
`t+1/2` needs `‖w_1(t+1/2)‖≤1/14` or `‖w_2(t+1/2)‖≤1/14`; the first is `1/2−‖w_1t‖ ≥ 6/14 > 1/14`,
impossible, so `‖w_2(t+1/2)‖≤1/14`, i.e. `‖w_2t‖ = 1/2 − ‖w_2(t+1/2)‖ ≥ 3/7`. Symmetrically if
`‖w_1t‖>1/14` then `‖w_2t‖≤1/14` and `‖w_1t‖≥3/7`. ∎ *(verified on 3728 doubly-covered points.)*

## Corollary — single-tightener components and the magnitude bound (f=2, PROVED)
Since `‖w_1t‖` is continuous on `R` and never lands in `(1/14, 3/7)`, each connected component `I` of `R`
has `‖w_1t‖` entirely `≤1/14` (a "`w_1`-danger" component, `w_2` then `≥3/7`) or entirely `≥3/7`
(a "`w_2`-danger" component); the `(+1/2)`-shift swaps the two types. A single-type component lies inside
ONE arc of the corresponding danger set, of length `1/(7w_i)`.

Apply this to the component `I_0 ∋ t_g := τ_0/2` where `τ_0 = argmax g_U` (`g_U(τ_0)=M(U)`). Since
`g_U(2·)` is `2u_max`-Lipschitz and peaks at `t_g`, `I_0 ⊇ (t_g−ρ, t_g+ρ)` with
`ρ=(M(U)−1/14)/(2u_max)`, so `|I_0| ≥ (M(U)−1/14)/u_max`. If `I_0` is `w_1`-danger,
`|I_0| ≤ 1/(7w_1)`, giving `w_1 ≤ u_max/(7(M(U)−1/14))`. With `e=11` runners, `M(U) ≥ 1/12` (LRC≤13),
so `M(U)−1/14 ≥ 1/84` and
> **`w_1, w_2 ≤ 12·u_max`.**

So a `q*=28`, `f=2`, primitive tight family is **compact**: its odd tighteners are bounded by `12×` the
even part's largest speed. (General `e`: the tightener covering the largest `R`-component is
`≤ [2(e+1)/(13−e)]·u_max`.)

## What this does and does not do
- **Does (proved):** reduces a hypothetical `q*=28` `f=2` tight family to a *compact* one (tighteners
  `≤ 12 u_max`) with a rigid extremity/component structure; reproves `f≥2`; gives `f ≥ 7 meas(R)`.
  Independent search (938 structured families) finds none, matching THM-612.
- **Does NOT (honest):** close confinement. The shift argument only uses that `F` **covers `R`**; it does
  not use the full tightness `M(S) ≤ 1/14` *off* `R`, nor that the max is *attained* at denominator `28`,
  nor primitivity beyond `E` even. A single-loose-arc `U` makes covering-of-`R` by two tighteners
  structurally easy, so the real obstruction lives in the global/attained/primitivity conditions the shift
  argument cannot see. `m ≥ 3` is untouched (the clean `1/2−‖wt‖` reflection is `m=2`-specific).

**Confinement (`primitive tight ⟹ q*=14`) remains a CONJECTURE** (THM-612). This note contributes the
extremity structure, the compactness bound `w_i ≤ 12 u_max`, and independent search confirmation — a
proper subset of the gap, clearly scoped.
