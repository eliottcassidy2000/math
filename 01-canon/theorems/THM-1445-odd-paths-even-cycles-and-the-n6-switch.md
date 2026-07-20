---
id: THM-1445
title: "ODD PATHS, EVEN CYCLES: the parity duality is FALSE, and the reason is that only one side has a distinguished empty object — but an exact low-n formula falls out. (0) THE DUALITY PROPOSED AND REFUTED: Rédei says hp is always ODD; Hamiltonian CYCLES are literally even subgraphs (every degree 2) where Hamiltonian PATHS are not (two degree-1 vertices), so 'hc is always EVEN' is the natural dual. It is FALSE. Exhaustively at n = 5, strong tournaments split 144 even / 400 odd; at n = 6, 7920 / 14400. The even half is trivial and uninformative: hc = 0 exactly on the non-strong tournaments (Camion), and 0 is even. (I) WHY IT FAILS, structurally: hp = I(Ω,2) = 1 + 2α₁ + 4α₂ + … is odd because the EMPTY independent set contributes the 1. Cycles are counted by no such polynomial and have no distinguished empty object, so there is nothing to force a parity. The even-graph connection is real but is a DEGREE statement, not a counting-parity one. (II) THE POSITIVE THAT FELL OUT, and it is proved rather than observed: two vertex-disjoint odd cycles need ≥ 3+3 = 6 vertices, so α_k = 0 for k ≥ 2 whenever n ≤ 5, whence hp(T) = 1 + 2·c_odd(T) EXACTLY for every tournament on at most 5 vertices — verified on all 8 + 64 + 1024 of them — and the identity breaks at exactly n = 6, witnessed by a tournament with disjoint odd cycles {0,3,4} and {1,2,5}, residual 4"
status: >
  (0) REFUTED-BY-EXHAUSTION over all 2^C(n,2) tournaments, n ≤ 6, with the strong /
  non-strong split computed separately so the trivial half is visible; corroborated by
  sampling at n = 7, 8.
  (I) Structural explanation, argued not computed, and labelled as such.
  (II) PROVED — the vanishing α_k = 0 for k ≥ 2 at n ≤ 5 is a vertex count, not an
  observation; the exhaustive check over n ≤ 5 confirms it and the n = 6 break is
  exhibited with an explicit disjoint pair.  This is the reverse of the usual order:
  the theorem comes from the counting bound and the computation only confirms the
  break lands where predicted.
  Nothing here advances a named open problem.  One proposed duality is killed and one
  exact small-n formula is pinned.
source: kind-pasteur-2026-07-20-S128c114 (owner: see the relation between odd-valued functions and tournament-adjacent ideas, and how both relate to even concepts like even graphs and even functions)
depends_on: []
related: [THM-165, THM-159]
script: 04-computation/odd_paths_even_cycles_kps_S128c114.py, ocf_low_n_exact_kps_S128c114.py (+ .out)
---

# THM-1445 — the duality is false, and here is what is true instead

## 0. The proposal, and its refutation

Rédei: `hp(T)` is always **odd**. And there is a genuine sense in which cycles are the
*even* object:

> a Hamiltonian **path** has two vertices of degree 1 — **not** an even subgraph;
> a Hamiltonian **cycle** has every vertex of degree 2 — it **is** an even subgraph,
> i.e. a point of the cycle space, the repo's `E_n`.

So "hc is always even" is the natural dual. **It is false.**

| n | strong: hc even / odd | non-strong: even / odd |
|---|---|---|
| 3 | 0 / 2 | 6 / 0 |
| 4 | 0 / 24 | 40 / 0 |
| **5** | **144 / 400** | 480 / 0 |
| **6** | **7920 / 14400** | 10448 / 0 |

The law survives to `n = 4` and dies at `n = 5`. And the even column is worthless
where it holds: **`hc = 0` exactly on the non-strong tournaments** (Camion — a
tournament has a Hamiltonian cycle iff it is strong), and `0` is even. So the apparent
"evenness" at small `n` is one trivial fact plus a coincidence that expires.

## I. Why it fails — only one side has a distinguished empty object

`hp = I(Ω(T), 2) = 1 + 2α₁ + 4α₂ + 8α₃ + ⋯`, where `α_k` counts independent sets of
size `k` in the odd-cycle conflict graph. **The oddness is the empty independent
set.** It contributes the `1`; every other term carries a factor of 2. Rédei's theorem
is that observation.

Hamiltonian cycles are counted by no such polynomial, and there is no distinguished
empty cycle to supply an unpaired `1`. Nothing forces a parity, and nothing does.

The even-graph connection is real but it is a statement about **degrees**, not about
counting parity. Conflating the two is what made the duality look plausible. A useful
way to hold it: *"even" as a property of the objects being counted does not propagate
to "even" as a property of how many there are.*

THM-165's coefficient — changing `hc` by `δ` changes `hp` by exactly `2δ` — is the
correct and much weaker true statement in this direction: cycles enter `hp` through the
even part, without their own count being even.

## II. What is true: an exact formula below the threshold

Two vertex-disjoint odd cycles need at least `3 + 3 = 6` vertices. Hence

> **`α_k = 0` for all `k ≥ 2` whenever `n ≤ 5`**, and therefore
> **`hp(T) = 1 + 2·c_odd(T)` exactly**, for every tournament on `n ≤ 5` vertices,

with `c_odd` the total number of directed odd cycles. This is a counting bound, not an
observation. Verified on all `8 + 64 + 1024` tournaments at `n = 3,4,5`: **zero
violations**.

And it must break at exactly `n = 6`, which it does. Explicit witness:

```
hp = 37,  c_odd = 16,  1 + 2·c_odd = 33,  residual = 4
disjoint odd cycles  {0,3,4}  and  {1,2,5}
```

`3 + 3 = 6` is the first room for two disjoint odd cycles, so `α₂ ≥ 1` and the residual
`4α₂ + 8α₃ + ⋯` becomes visible — divisible by 4, as the grading requires.

**So the OCF's higher terms switch on at precisely `n = 6`, and the reason is a vertex
count rather than an arithmetic accident.** That also explains, in one line, why so many
of this repo's tournament patterns hold to `n = 5` and break at `n = 6`: below the
threshold `hp` is an affine function of a single cycle count, and above it the
independence structure of `Ω` starts to matter.

## Named next

- The same threshold argument gives the next one for free: `α_3 = 0` for `n ≤ 8`, since
  three disjoint odd cycles need `9` vertices. So `hp = 1 + 2α₁ + 4α₂` exactly for
  `n ≤ 8` — worth confirming, and it would extend the exact regime by three.
- The general statement is `α_k = 0` for `n < 3k`, i.e. the OCF truncates at
  `⌊n/3⌋` terms. That is a clean grading bound and I have not seen it stated.
