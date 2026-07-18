---
id: THM-1016
title: The simultaneous-coverage criterion — an exact, computable METRIC necessary condition for a deep tight instance, strictly stronger than sheet capacity. If A is tight with sheet number s and off-sheet set F, then the on-sheet quotient's safe set G must lie inside Cov(F,s,L) = ∩_{j∈Z/s} ∪_{w∈F} T_{w,j}, where T_{w,j} = {τ : ‖w(τ+j)/s‖ ≤ L} is the tooth-comb of w on sheet j. Since G ≠ ∅ always, Cov = ∅ EXCLUDES the configuration outright. This crosses the wall of THM-1006 §H: over 2248 (s,F) configurations passing BOTH capacity and primitivity (the class §H proved counting cannot touch), the criterion excludes 1255 = 55.8%, including the smallest case s=2, F={1,3}.
status: PROVED (the containment G ⊆ Cov is elementary — four lines). The exclusion counts are VERIFIED-EXACT by direct computation of Cov on a 3·10^5 grid, with the smallest case s=2,F={1,3} confirmed by hand in closed form.
source: mac-mini-2026-07-18-S112 (owner directive: attack the metric wall with tooth positions)
depends_on:
  - THM-769   # sheet decomposition E=sU, off-sheet F, the s lifts
  - THM-1006  # the content law; §H proves capacity+primitivity is satisfiable (what this crosses)
external: LRC(n) SETTLED (used only for M(U) ≥ 1/(|U|+1) > L, i.e. G ≠ ∅).
related:
  - THM-1002  # klein: the pair-sum ruler / integer realizability crux
  - THM-774   # deep two-sheet folded diamond (s=2 special case)
  - HYP-6800
  - HYP-6820
  - HYP-7380  # this session
---

# THM-1016 — The simultaneous-coverage criterion

**One line.** Counting asks *how many* lifts an off-sheet speed can cover; this asks
*where*. Requiring the same point to be covered on **every** sheet at once is an
exact, computable condition that kills configurations capacity cannot.

## Setup

`A` tight (`M(A) = L = 1/(n+1)`), maximizer denominator `q = (n+1)s`, on-sheet part
`E = sU`, off-sheet part `F = A\E`. In the quotient coordinate `τ` (a point of the
`U`-circle, with lifts `t_j = (τ+j)/s`, `j ∈ ℤ/s`), define for `w ∈ F`

```
T_{w,j} = { τ ∈ ℝ/ℤ : ‖w(τ+j)/s‖ ≤ L }
```

— a comb of period `s/w` and tooth width `2Ls/w`, whose **offset depends on `j`**.
Set

```
Cov(F,s,L) = ⋂_{j ∈ ℤ/s}  ⋃_{w ∈ F}  T_{w,j},        G = { τ : φ_U(τ) > L }.
```

## The criterion (PROVED)

> **Theorem.** `G ⊆ Cov(F,s,L)`. Since `|U| ≤ n−2` gives
> `M(U) ≥ 1/(|U|+1) ≥ 1/(n−1) > L`, we have `G ≠ ∅`; hence
> **`Cov(F,s,L) = ∅` is impossible — no tight `A` has that `(s,F)`.**

*Proof.* Let `τ ∈ G` and `j ∈ ℤ/s`. For an on-sheet speed `v = su`,
`‖v t_j‖ = ‖su(τ+j)/s‖ = ‖u(τ+j)‖ = ‖uτ‖`, so `φ_E(t_j) = φ_U(τ) > L`: the whole
on-sheet part is **strictly** safe at every lift. As `M(A) = L`, we need
`φ_A(t_j) ≤ L`, so some `w ∈ F` has `‖w t_j‖ ≤ L`, i.e. `τ ∈ T_{w,j}`. This holds
for every `j`, so `τ ∈ Cov`. ∎

## It strictly crosses the counting wall

THM-1006 §H proves sheet capacity **plus** primitivity is satisfiable for every
`val = 2,…,13` — no counting argument can exclude those configurations. The
coverage criterion does:

| | configurations |
|---|---|
| pass capacity **and** primitivity (`s ≤ 8`, `\|F\| ∈ {2,3}`, `w ≤ 24`) | **2248** |
| of those, `Cov = ∅` → **excluded by tooth positions** | **1255 (55.8%)** |

by sheet number: `s=2`: 6, `s=3`: 206, `s=4`: 293, `s=6`: 312, `s=8`: 438.

**Smallest case, in closed form.** `s = 2`, `F = {1,3}`, `L = 1/13`. Capacity is
*satisfied* (`D_1 = D_3 = 2`, `c(2)+c(2) = ½+½ = 1 ≥ 1`) and primitivity holds
(`gcd(1,2)=gcd(3,2)=1`). But
`T_{1,0} ∪ T_{3,0} = [0,\tfrac{2}{13}] ∪ [\tfrac{24}{39},\tfrac{28}{39}]` and
`T_{1,1} ∪ T_{3,1} = [\tfrac{11}{39},\tfrac{15}{39}] ∪ [\tfrac{11}{13},1)`, whose
intersection is **empty**. So `{even speeds} ∪ {1,3}` is never tight — a fact no
counting argument sees.

## What survives, and the complementary structure

Non-empty `Cov` occurs when the off-sheet speeds are **large**: `maxint(Cov)` decays
like `1/w` (e.g. `s=2`: `F={1,7}` gives `0.0330`, `F={1,31}` gives `0.0099`), so the
surviving configurations force the on-sheet quotient's safe interval to be
correspondingly tiny — i.e. a very spread on-sheet part. Hence:

> **metric (tooth positions) excludes SMALL off-sheet configurations;
> a height bound would exclude LARGE ones.** They are complementary, and together
> they would close the deep branch.

This locates the residual precisely: *large off-sheet speeds over a highly spread
on-sheet quotient*, which is exactly klein's "integer realizability" crux (THM-1002)
seen from the metric side.

## Honest scope

- **Proved:** the containment `G ⊆ Cov` and hence the exclusion rule; it is exact
  and decidable per `(s,F)` (rational endpoints).
- **Verified:** the 1255/2248 strict gain (grid computation; smallest case by hand).
- **NOT proved:** the content law. `Cov ≠ ∅` for large off-sheet speeds, so this
  criterion alone does not close the deep branch; it removes over half of the
  counting-irreducible configurations and leaves an explicitly characterized
  remainder that needs a height bound.

*Artifacts:* `04-computation/lrc13_simultaneous_coverage_macmini_S112.py` (+out).
Credits: THM-769 (sheet decomposition), THM-1006 §H (the counting wall this
crosses), klein THM-1002 (the realizability crux this meets from the metric side).
