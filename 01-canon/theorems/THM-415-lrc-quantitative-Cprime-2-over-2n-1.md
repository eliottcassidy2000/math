---
id: THM-415
title: Quantitative C′(n) — every multiple-of-n config has M ≥ 2/(2n-1), the s=2 Kravitz value
status: VERIFIED (exhaustive boxes n=4..8, range-stable); general n = CONJECTURE (= Kravitz second-gap)
source: claudebox-2026-06-03-S623
depends_on:
  - THM-398   # C′ reduction: multiple-of-n ⟹ loose
  - THM-401   # the 2n-1 shell (pair-sum-sieve modulus) — where the extremal witness lives
related:
  - THM-412   # the twisted-shell dodge (S622) — corrected/contextualised here
  - HYP-2220  # the corrected constructive route to LRC(14)
  - HYP-2055  # n=14 tight family / witnesses on (ℤ/m)*
---

# THM-415 — quantitative C′(n): the multiple-of-n loneliness gap is 2/(2n-1)

## Statement

Convention as in [[THM-398]]: `n` runners, `S` = `n-1` distinct positive integers (primitive WLOG by
scale-invariance `M(cS)=M(S)`), `M(S)=max_t min_{v∈S}‖vt‖`. A config is **multiple-of-n** if `n | v`
for some `v ∈ S` — exactly the class on which the `1/n`-clock witness fails (C′, the whole of LRC).

> **THM-415.** `min { M(S) : S primitive, multiple-of-n } = 2/(2n-1)`.
>
> Equivalently the *quantitative* C′(n): **`n | v` for some `v ∈ S` ⟹ M(S) ≥ 2/(2n-1)`** (sharp).

The extremal value `2/(2n-1)` has the minimal witness `t = a/(2n-1)` on the **`2n-1` shell**
(the pair-sum-sieve modulus, [[THM-401]]) at **band-distance 2**.

## Evidence

Exhaustive over primitive multiple-of-n configs in a box, exact over ℚ
(`lrc_cprime_gap_s623.py` → `.out`):

| n | min M | `2/(2n-1)` | an extremal config |
|---|---|---|---|
| 4 | 2/7  | 2/7  | (1,3,4) |
| 5 | 2/9  | 2/9  | (1,3,4,5) |
| 6 | 2/11 | 2/11 | (1,3,4,5,18) |
| 7 | 2/13 | 2/13 | (1,2,5,6,7,8) |
| 8 | 2/15 | 2/15 | (1,4,5,6,7,11,16) |

Range-stable (n=4,5 checked to box `6n`). The minimizers vary in shape (the AP-with-2-deleted
`(1,3,4,…,n)` is extremal only for n=4,5), but the **value `2/(2n-1)` is robust** and forces the
witness denominator to be `2n-1`.

## The Kravitz normalization (why the user keeps pointing here)

In Kravitz's "barely/very lonely runners" normalization (`n_K = n-1` speeds, `ML ≥ 1/(n_K+1) = 1/n`),
the conjectured loneliness ladder is `s/(n_K s + 1)`. The tight value is the `s=1` rung
`1/(n_K+1) = 1/n`; the next rung is `s=2`:
```
2/(2 n_K + 1) = 2/(2(n-1)+1) = 2/(2n-1).
```
THM-415 says the multiple-of-n class (the configs that cannot be tight) lands **exactly on the
`s=2` rung**, and no lower. So it **confirms Kravitz's second-gap for this class** and identifies the
mechanism: the gap is realised on the `2n-1` shell at band-distance 2.

## Relation to the twisted-shell dodge (THM-412, corrected)

[[THM-412]] (S622) proposed escaping every multiple-of-n config via a twisted clock on a shell
`m ≤ 2n-1`. The honest statement (small-n exhaustion, S623): the `m ≤ 2n-1` ceiling holds for the
**extremal / hardest** configs — they sit *exactly* at `2/(2n-1)` on the `2n-1` shell — but **not**
for all loose configs: configs with margin (e.g. `(1,4,5,9)` at n=5, `M = 3/13`) realise their
optimum on higher shells. They are all still `≥ 2/(2n-1)`. So `2n-1` is the **extremal shell**, not a
universal dodge ceiling (S622's "0 residual" was a sampling artifact). THM-415 is the correct,
shell-anchored, quantitative form.

## Consequence for n=14

`multiple-of-14 ⟹ M ≥ 2/27` (the `s=2` rung, witness on the `27 = 3³` shell). This is the sharp
constructive target for C′(14) ⟹ LRC(14) — see [[HYP-2220]].

## Verification
`04-computation/lrc_cprime_small_s623.py`, `lrc_cprime_gap_s623.py` → `05-knowledge/results/*_s623.out`.
