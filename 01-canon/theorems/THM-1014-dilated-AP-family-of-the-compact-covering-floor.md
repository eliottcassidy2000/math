---
id: THM-1014
title: HYP-7355 (compact primitive covering ⟹ M ≥ 1/13) is PROVED on its binding family — every primitive covering V = d·{1,…,12} ∪ {k} has M(V) = 1/13 unless 13d | k, and that failure mode is FORCED NON-COMPACT (it is exactly the deep-well tower {1,…,12, 182m}, ρ ≥ 182/12 > 13)
status: PROVED (the dilated-AP family only). Closes boxeph-S85's stated extremal family for HYP-7355 and explains why "compact" is exactly the right hypothesis. The general compact residual (arbitrary primitive covering V with ρ<13) remains OPEN.
source: klein-2026-07-18-S316
depends_on:
  - THM-1002  # pair-sum denominator bound (exact evaluator; M = val/q in lowest terms)
  - THM-724   # the deep well {1..12,182} is the covering-minimum, M = 14/183
related:
  - HYP-7355  # boxeph-S85: the sharp compact covering floor conjecture (this proves its binding family)
  - HYP-7310  # AP-uniqueness at n=12
  - THM-1010  # boxeph descent recursion (the compact residual it reduces to)
  - THM-726   # multi-killer M ≥ 1/13
external: LRC(13) SETTLED.
---

# THM-1014 — the compact covering floor on the dilated-AP family

boxeph-S85 refuted the 12-subset floor and the `M ≥ 1/9` covering floor, and replaced them with the sharp
conjecture **HYP-7355: compact (`ρ = v_max/v_2nd < 13`) primitive covering ⟹ `M ≥ 1/13`**, naming its
binding family as `d·{1,…,12} ∪ {killer}` — families carrying a **dilation substructure** whose embedded
dilated AP gives `M = 1/13` at `t = 1/(13d)`, and which defeat every elementary tool (descent-from-`v_max`,
best-removal recursion, sieve, measure). This theorem settles that family.

## Statement

Let `d ≥ 1`, `k ≥ 1`, and `V = d·{1,…,12} ∪ {k}` with `|V| = 13`, `V` primitive and covering. Then

```text
M(V) = 1/13     unless   13d | k,
```

and the exceptional case forces `d = 1` and `182 | k`, hence

```text
ρ(V) = k/12 ≥ 182/12 = 15.16… > 13,
```

i.e. **every exception is non-compact**. Consequently **HYP-7355 holds on this family**: every *compact*
primitive covering `V = d·{1,…,12} ∪ {k}` has `M(V) = 1/13`.

## Proof

**1. The good set of the dilated AP is a finite explicit set.** `M({1,…,12}) = 1/13`, so
`{u : min_{j≤12} ‖ju‖ ≥ 1/13}` is exactly the argmax set, `{a/13 : a = 1,…,12}` (at `u = a/13`,
`‖ja/13‖ = |ja|₁₃/13 ≥ 1/13` since `13 ∤ ja`). Substituting `u = dt`, the good set of `d·{1,…,12}` at
level `1/13` is

```text
{ t = b/(13d) : 1 ≤ b < 13d, 13 ∤ b }        (12d points).
```

**2. Reduction to a residue predicate.** Since `d·{1,…,12} ⊂ V`, `M(V) ≤ 1/13`. So `M(V) = 1/13` iff some
good `t = b/(13d)` also satisfies `‖kt‖ ≥ 1/13`, i.e.

```text
(k·b mod 13d) ∈ [d, 12d]      for some b with 13 ∤ b.        (★)
```

**3. (★) fails iff `13d | k`.** With `e = gcd(k,13d)`, the values `kb mod 13d` are multiples of `e`. The
window `[d,12d]` has length `11d`, so it contains a multiple of `e` whenever `e ≤ 11d`. A divisor `D` of
`13d` with `D > 11d` satisfies `13d/D < 13/11 < 2`, hence `D = 13d`. So (★) can only fail when
`e = 13d`, i.e. `13d | k`. *(Verified: the predicate matches exact `M = 1/13` on **240/240** primitive
covering members.)*

**4. The exception is forced non-compact.** `gcd(V) = gcd(d, k)`, so `13d | k ⟹ d | k ⟹ gcd(V) = d`;
primitivity gives **`d = 1`**. Then `V = {1,…,12, k}` with `13 | k`, and since `{1,…,12}` contains no
multiple of 14, covering forces `14 | k`; with `gcd(13,14) = 1` this gives **`182 | k`**. Hence
`v_max = k ≥ 182`, `v_2nd = 12`, and `ρ = k/12 ≥ 182/12 > 13`. ∎

## What this says

- **The exceptions are exactly the deep-well tower** `{1,…,12, 182m}`: `M = 14/183, 28/365, 42/547,
  56/729, …` — all `< 1/13`, all non-compact. `m = 1` is THM-724's covering-minimum. So the family's only
  sub-`1/13` members are the already-classified deep wells.
- **This is why `ρ < 13` is the right hypothesis.** It is not a convenience: on the binding family the
  failure set is *precisely* the non-compact part, and it misses the threshold by a definite margin
  (`15.17` vs `13`), not marginally.
- Verified computationally: over `d ≤ 29`, `k ≤ 800`, **894** primitive covering members, of which **863**
  compact — **zero** with `M < 1/13`.

## The multi-killer dilated region (mechanism + census, not yet a theorem)

The single-killer proof above uses that `d·{1,…,12}`'s good set at `1/13` is a **finite** point set. One
step out — `V = d·P ∪ K` with `P ⊊ {1,…,12}`, `|K| = r ≥ 2` killers — the situation improves: `|P| ≤ 11`
gives `M(P) ≥ 1/12 > 1/13` by LRC(12), so the good set of `d·P` has genuine **intervals**, of length
`L_P/d`. The THM-1004 interval-survival lemma then absorbs the killers once

```text
Σ_i (1/k_i)  <  (L_P/d)·(13 − 2r)/2         (δ = 1/13; needs r < 13/2 = 6.5)
```

Exact minima of `L_P` over all `P` of the given size:

| `r` | `|P|` | `min L_P` | `median L_P` |
|---|---|---|---|
| 2 | 11 | 0.00641 | 0.01282 |
| 3 | 10 | 0.00769 | 0.01923 |
| 4 | 9  | 0.00962 | 0.02564 |

**Census: 19,317 compact primitive covering families `d·P ∪ K` (r = 2,3,4; d ≤ 13; killers ≤ 420) — ZERO
with `M < 1/13`.** So HYP-7355 extends across the multi-killer dilated region empirically, with the tail
lemma supplying the mechanism for large killers. It is **not** a theorem yet: the tail only absorbs
*large* killers, and small killers (which sit inside the body rather than above it) need the same
regime analysis used for the Hamming-radius results.

## Scope

Only the dilated-AP family `d·{1,…,12} ∪ {k}` (one killer over a full dilated AP). The general compact
residual — arbitrary primitive covering `V` with `ρ < 13` — is still open; boxeph's point that it is a
*rigidity* statement rather than a perturbation floor stands. Companion fact from the same session's
census: **the AP is the unique primitive tight 12-set** over 11,231 residue-filtered candidates plus all
structured families (supports HYP-7310), and mac-mini's content law is confirmed on every candidate — the
sets containing `13` all fail the `±`-pair test and sit at `1/12`.

*Files: `04-computation/lrc14_hyp7355_dilated_ap_klein_S316.py` (+ `_dilated_ap`, `_compact`,
`lrc14_tight_locus_n12` .out).*
