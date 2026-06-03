---
id: THM-410
title: LRC covering-depth conservation and the moment–sieve identity
status: PROVED
source: opus-2026-06-03-S620
depends_on: []
related:
  - HYP-2190  # Helly-entropy accounting (the frame this backbone serves)
  - THM-401   # pair-sum sieve modulus 2n-1 (the order-2 / Helly-dim-1 layer)
  - HYP-2065  # S561 sieve-covered density rho (= the alternating moment sum, special case)
  - THM-404   # doubling-rigidity dichotomy (where the depth charge concentrates)
---

# THM-410 — covering-depth conservation and the moment–sieve identity

## Setup

Fix `n` relative speeds `v_1,…,v_n ∈ ℤ\{0}` and the gap `δ = 1/(n+1)`. On the clock
`ℝ/ℤ`, runner `i`'s **forbidden arc set** is
`B_i = { t : ‖v_i t‖ < δ }`, a union of `|v_i|` arcs of total Lebesgue measure `2δ`.
A **loneliness certificate** (LRC witness at gap `δ`) is a `t ∉ ⋃_i B_i`. Define the
**covering-depth function** and its **distribution**

```
depth(t) = #{ i : t ∈ B_i },        p_k = meas{ t ∈ [0,1) : depth(t) = k }.
```

The lonely (certificate) times are exactly `{depth = 0}`, of measure `p_0`.

## Statement

> **(R1) Conservation.** The mean covering depth is arithmetic-independent:
> `Σ_k k·p_k = ∫_0^1 depth(t) dt = 2nδ = 2n/(n+1) < 2`, for every speed set and every `n`.
>
> **(R2) Moment–sieve identity.** The `k`-th factorial moment of `depth` equals the
> order-`k` inclusion–exclusion (sieve) term, and the free measure is their alternating sum:
> ```
> S_k := Σ_k C(k', k) p_{k'} = Σ_{|T|=k} meas( ⋂_{i∈T} B_i ),      p_0 = Σ_{k=0}^{n} (-1)^k S_k.
> ```

## Proof

**(R1).** `depth(t) = Σ_i 𝟙_{B_i}(t)`, so by linearity `∫ depth = Σ_i meas(B_i) = n·(2δ)`.
With `δ = 1/(n+1)` this is `2n/(n+1)`, which is `< 2` for all `n ≥ 1`. ∎

**(R2).** Pointwise `C(depth(t), k) = Σ_{|T|=k} ∏_{i∈T} 𝟙_{B_i}(t) = Σ_{|T|=k} 𝟙_{⋂_{i∈T}B_i}(t)`,
since choosing `k` of the active indices is choosing a `k`-subset all of whose arcs contain `t`.
Integrating gives `S_k = E[C(depth,k)] = Σ_{|T|=k} meas(⋂_{i∈T} B_i)`. Inclusion–exclusion on the
event `{depth = 0} = ⋂_i B_i^c` (equivalently the standard `p_0 = Σ_k (-1)^k S_k` from the
factorial-moment/Bonferroni expansion) yields the alternating sum. ∎

## Why it matters

* **One object subsumes the LRC invariants at gap `δ`.** The depth distribution `{p_k}` carries:
  the lonely measure (`p_0`), the conserved mean (`S_1`), the order-`k` sieve data (`S_k`), and —
  as a derived functional — the depth entropy `H = -Σ p_k log p_k` (see [[HYP-2190]]).
* **S561's `ρ(S) = Σ_T (-1)^{|T|}/lcm(T)` is exactly this alternating moment sum** for the
  all-integer "AP-overlap" specialization, where `meas(⋂_{i∈T} B_i) → 1/lcm(T)` in the tight
  limit. THM-410 is its model-free generalization to arbitrary speeds and finite `δ`.
* **The mean is a conserved charge `< 2`.** LRC at gap `1/(n+1)` is the question: *can a
  nonnegative-integer depth field with mean pinned just under 2 be forced to be `≥ 1` almost
  everywhere?* Forcing `p_0 = 0` (no lonely time of positive measure) is an extremal,
  entropy-minimizing event; the AP achieves it (with the tight witness surviving at measure
  zero), but **not uniquely** — sporadic additive chains do too (e.g. `(1,3,4,7)`; [[HYP-2190]]).

## Verification

`04-computation/lrc_helly_entropy_s620.py` → `05-knowledge/results/lrc_helly_entropy_s620.out`.
Exact over ℚ (arc endpoints are rational). `S_1 = 2n/(n+1)` to machine precision for all sets;
`S_2` from the moments matches direct pairwise arc-overlap; `Σ_k(-1)^k S_k` matches `p_0` (all
`match=True`). `n = 4,5,6,7`.
