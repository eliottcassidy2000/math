---
id: THM-1015
title: The CLUSTERED multi-killer stratum closes by interval survival — r killers over a core P force M > 1/14 as soon as Σ 1/k_i < L_P·(7−r), with NO spacing hypothesis. This reaches exactly where the telescoping balance cannot, because the tail has no per-step (13/14) loss (the "169" obstruction).
status: PROVED as an explicit criterion (elementary; cites LRC(≤13) SETTLED). Closes kind-pasteur's open clustered stratum for killers above an explicit, computable threshold; killers below it remain.
source: klein-2026-07-18-S317
depends_on:
  - THM-1004  # the interval-survival tail lemma
  - THM-1007  # kind-pasteur: single-killer + lacunary-chain closure (this handles the clustered complement)
related:
  - MISTAKE-104   # the n=13 deep well {1..11,168}, 168 = 13²−1, witness 14/169 — the "169" thread
  - THM-897       # triple-beat discrepancy 22κ_W/(169·x₃)
  - HYP-7355      # the compact residual = the small-killer regime
external: LRC(≤13) SETTLED — gives M(P) ≥ 1/(14−r) for the (13−r)-element core.
---

# THM-1015 — the clustered multi-killer stratum

kind-pasteur's THM-1007 closed *covering ⟹ M > 1/14* for the **single-killer** case and for **lacunary
chains** (each killer `> 13×` the running max), by iterating THM-724's balance lemma at a cost of `13/14`
per step. Their honest caveat: **clustered** killers (`v₂ ≈ v₁ ≫ core`) break the telescoping — the factor
degrades to `~1/2` — and that stratum stayed open. This theorem closes it for large killers, by a route
that has no per-step cost at all.

## Statement

Let `V = P ∪ K` be a primitive 13-set with core `|P| = 13 − r` and killers `K = {k_1,…,k_r}`,
`2 ≤ r ≤ 6`. Let `δ = 1/14` and let `L_P > 0` be the length of the largest interval of
`G_P(δ) = {t : ‖vt‖ ≥ δ ∀v ∈ P}`. If

```text
Σ_{i=1}^{r} 1/k_i  <  L_P · (7 − r)                          (★)
```

then `M(V) ≥ 1/14`, with strict inequality at an interior point. In particular if every
`k_i > r / (L_P (7 − r))` then (★) holds.

**No hypothesis on the spacing of the `k_i` is used.** Clustered, lacunary or arbitrary — only size enters.

## Proof

`P` has `13 − r` elements, so `M(P) ≥ 1/(14 − r)` by LRC(14−r) (settled, since `r ≥ 1`), and
`1/(14−r) > 1/14`; hence `G_P(δ)` has nonempty interior and `L_P > 0`. Let `I ⊆ G_P(δ)` be an interval of
length `L_P`. By the THM-1004 tail lemma, the bad set of a single new speed `k` meets `I` in measure
`≤ 2δL_P + 2δ/k`, so all `r` killers together remove at most `2rδL_P + 2δ Σ 1/k_i`. A point of `I` survives
all of them provided

```text
2rδL_P + 2δ Σ (1/k_i) < L_P   ⟺   Σ (1/k_i) < L_P (1 − 2rδ)/(2δ).
```

At `δ = 1/14`: `2δ = 1/7` and `1 − 2rδ = (7−r)/7`, so the right side is `L_P·(7−r)`, which is (★). Such a
`t` satisfies `‖vt‖ ≥ δ` for every `v ∈ P ∪ K = V`, i.e. `M(V) ≥ 1/14`. ∎

The condition `r < 1/(2δ) = 7` is exactly the ceiling: at `r = 7` the coefficient `1 − 2rδ` vanishes.

## Why this reaches the clustered case — and what 169 was measuring

The telescoping route pays `13/14` **per killer**, so `j` killers cost `(13/14)^j`; at `j = 2` that is
`169/196`, and the chain's tightest strict value is `(1/12)(13/14)² = 169/2352 = 0.071854`, barely above
`1/14 = 0.071428`. That `169 = 13²` **is** the two-step loss, and it is why clustering — which forces the
steps to be taken against a barely-shrinking max — degrades the product to `~1/2` and kills the argument.

The interval-survival tail pays **nothing per step**: the killers enter only through `Σ 1/k_i`, additively
and symmetrically. Spacing is invisible to it. That is precisely the structural reason it covers the
stratum the balance ladder cannot, and it is the payoff of the `169` thread: the obstruction was a
*per-step multiplicative* cost, so the fix is an *additive* certificate.

## Explicit thresholds and census

Minimum `L_P` over sampled cores drawn from `[1,20)`, and the resulting killer threshold
`r/(L_P(7−r))`:

| `r` | core size | `min L_P` | killers must exceed |
|---|---|---|---|
| 2 | 11 | 0.00609 | 65.7 |
| 3 | 10 | 0.00829 | 90.5 |
| 4 | 9  | 0.01065 | 125.2 |
| 5 | 8  | 0.01099 | 227.5 |
| 6 | 7  | 0.01727 | 347.5 |

**Census (clustered by construction — all killers within ±20 of a common base in [200,600]):** 349
primitive covering families, `min M = 2/21 = 0.09524`, i.e. **1.33× the threshold**, zero at or below
`1/14`. Consistent with kind-pasteur's measurement that the stratum is ~1.9× loose.

## Scope

- Killers **below** the threshold are not covered. Those are moderate-size killers, where the tail's
  `2δ/k` term is no longer small.
- The **small-killer regime** (killers comparable to or below the core) is not this regime at all: there
  the set has bounded ratio, i.e. it is the **compact** case `ρ < 13`, which is exactly HYP-7355 /
  THM-1014 territory — not a gap in this theorem but a pointer to the other one.
- So *covering ⟹ M > 1/14* now stands as: single-killer PROVED (THM-1007), lacunary multi-killer PROVED
  (THM-1007 ext.), **clustered/arbitrary multi-killer with large killers PROVED (here)**, moderate
  killers and the compact core remaining.

*Files: `04-computation/lrc14_clustered_tail_klein_S317.py` (+ .out).*
