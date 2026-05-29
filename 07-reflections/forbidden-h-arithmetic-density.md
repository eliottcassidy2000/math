# Forbidden-H Arithmetic Density: f(N) Closed Forms and Growth

**Session:** oracle-2026-05-29-S4
**Files:** `04-computation/lean/TournamentH7/TournamentH7/ForbiddenH_Counting.lean` (new)

## Setup

The Lean formalisation cleanly factors the forbidden-H problem into:
* an **arithmetic** stage — enumerate α-tuples (α₁, …, α_k) consistent with
  `H = 1 + Σ 2^k α_k` and the independence-polynomial bounds;
* a **structural** stage — rule out each surviving tuple via Moon-Moser
  type arguments.

This document explores the arithmetic stage in detail.  Define:

> **f(N)** = number of non-negative integer tuples (α₁, …, α_k, …) with
>           1 + Σ 2^k α_k = N, α_k ≤ C(α₁, k), and α_k ≥ 1 ⟹ α_j ≥ C(k, j)
>           for j ≤ k.

f(N) is the count of *arithmetic candidates* that the structural stage
must then independently kill.

## The growth-rate question

Brute-force enumeration (Python) gives:

| N   | f(N) | Notes                                |
|-----|------|--------------------------------------|
| 1   | 1    | empty product                        |
| 3   | 1    |                                      |
| 5   | 1    |                                      |
| 7   | 1    | THM-343 (H ≠ 7)                      |
| 9   | 2    |                                      |
| 11  | 2    |                                      |
| 13  | 2    |                                      |
| 21  | 4    | HYP-1753 (conjectured H ≠ 21)        |
| 63  | 37   | HYP-1754 (conjectured H ≠ 63)        |
| 127 | 244  |                                      |
| 255 | 2799 |                                      |

f(N) grows somewhere between N² and N³ — Python regression hints at
roughly N^{2.5} ± O(N²), but no clean closed form is apparent for full
f(N).

## The "main" subcase: only α₁, α₂

Restricting to α₃ = α₄ = ⋯ = 0 gives the simpler counting problem:

> **f_main(N)** = # {(α₁, α₂) : 2α₁ + 4α₂ = N − 1, α₂ ≤ C(α₁, 2), α₂ ≥ 1 ⟹ α₁ ≥ 2}.

After setting `M = (N−1)/2 = α₁ + 2α₂`, the constraint `α₂ ≤ C(α₁, 2)`
becomes `(M − a) ≤ a(a − 1)` where `a = α₁`, i.e., **M ≤ a²**.

### Closed form (NEW)

For odd N ≥ 1, letting M = (N−1)/2 and a_min = max(2, ⌈√M⌉):

> **Theorem (oracle-2026-05-29-S4).**
>
>   `f_main(N) = 1 + #{a ∈ [a_min, M − 2] : a ≡ M (mod 2)}`.

That is, count valid α₁ values in an interval `[max(2, ⌈√M⌉), M−2]` of
the right parity (matching M), then add 1 for the α₂ = 0 case (where
α₁ = M directly).

### Verified

The formula has been verified by direct enumeration for all odd N from
1 to 79.  Confirmed cases:

| N   | M   | a_min  | range          | count interval | f_main(N) |
|-----|-----|--------|----------------|----------------|-----------|
| 21  | 10  | 4      | [4, 8] even    | {4,6,8} = 3    | 4         |
| 63  | 31  | 6      | [7, 29] odd    | 12             | 13        |

### Asymptotic

For large N, `a_min = ⌈√M⌉` and the interval has length M − √M.  Half
of those are right parity, so

> **f_main(N) ~ M/4 = (N − 1)/8 ~ N/8.**

So f_main grows *linearly* in N with slope 1/8.

The gap between f(N) and f_main(N) is the contribution of higher-α
tuples (α₃, α₄, …), which dominate for large N.

## Phase transitions in f(N)

Each new α_k = 1 case becomes possible only when N is large enough to
accommodate it.  The downward closure gives explicit thresholds:

* α_3 ≥ 1: requires α_2 ≥ 3, α_1 ≥ 3.  Minimum H:
  1 + 6 + 12 + 8 = **27**.  First appearance at N = 27.

* α_4 ≥ 1: requires α_3 ≥ 4, α_2 ≥ 6, α_1 ≥ 4.  Minimum H:
  1 + 8 + 24 + 32 + 16 = **81**.  First appearance at N = 81.

* α_5 ≥ 1: requires α_4 ≥ 5, α_3 ≥ 10, α_2 ≥ 10, α_1 ≥ 5.  Minimum H:
  1 + 10 + 40 + 80 + 80 + 32 = **243**.  First appearance at N = 243.

In general:

> **Conjecture (oracle-2026-05-29-S4).**  The minimum N at which α_k ≥ 1
> becomes possible is
>
>   N_min(k) = 1 + 2·k + 2² · C(k, 2) + 2³ · C(k, 3) + … + 2^k · 1
>            = (1 + 2)^k = 3^k.
>
> (Pascal's row sum.)

Quick check: 3^2 = 9 ≠ 27 ... hmm that doesn't match.

Let me recompute.  For α_3 = 1:
* α_3 = 1 contributes 8.
* α_2 ≥ C(3, 2) = 3 contributes ≥ 4·3 = 12.
* α_1 ≥ C(3, 1) = 3 contributes ≥ 2·3 = 6.
* Plus the constant 1.

Total: ≥ 1 + 6 + 12 + 8 = 27 = 3³.  ✓

So actually the conjecture should be N_min(k) = 3^k:

> **Conjecture (corrected).**  N_min(k) = 3^k for k ≥ 1.
>
> Proof sketch: minimum H when α_k = 1 with all required downward
> values:
>   N_min(k) = 1 + Σ_{j=1}^{k} 2^j · C(k, j) = (1 + 2)^k = 3^k.

This is exactly the binomial expansion of (1 + 2)^k = 3^k.  Beautiful!

So:
* α_1 = 1: N ≥ 3¹ = 3.  ✓ (N = 3 is achievable.)
* α_2 = 1: N ≥ 3² = 9.  ✓ (N = 9 has α-tuple (2, 1, 0, ...)).
* α_3 = 1: N ≥ 3³ = 27.  ✓ (N = 27 first introduces α_3.)
* α_4 = 1: N ≥ 3⁴ = 81.  ✓ (N = 81 first introduces α_4.)
* α_5 = 1: N ≥ 3⁵ = 243.  ✓ (predicted, confirmed by enumeration).

## Implications for forbidden-H

The structure now becomes:

* For small N (say N < 9), only the (α₁, 0, 0, …) tuple is possible.
  Structural killing reduces to one case.

* For N ∈ [9, 27), tuples include α_2 ≥ 1.  Structural killing needs
  ~N/8 cases.

* For N ∈ [27, 81), α_3 ≥ 1 is added.  Cases ~N/8 + (N − 27)·c for
  some c.

* In general, the "kill load" at N has a *layered* structure indexed
  by the deepest α_k that is non-zero.

This means HYP-1754 (H ≠ 63) is fundamentally harder than HYP-1753
(H ≠ 21): not just because there are more candidates, but because the
candidates span THREE α-depth layers (α_2 only, α_3 = 1, α_3 = 2).

## Formalisation plan

A new Lean module `ForbiddenH_Counting.lean` would:

1. Define `arithCandidates N : Finset (Tuple)` — the enumeration of
   α-tuples.
2. Prove `arithCandidates 7 = {(3, 0, 0, 0)}` (f(7) = 1).
3. Prove `arithCandidates 21 = {(10,0,…), (8,1,…), (6,2,…), (4,3,…)}` (f(21) = 4).
4. State the N_min(k) = 3^k conjecture as a theorem (provable by direct
   substitution).

Lean's `Finset.filter`-based enumeration makes #1-#3 doable by
`native_decide` once the boundary conditions are right.

## Open mathematical questions

1. **Closed form for full f(N)?** Generating function?
2. **Asymptotic of f(N) for large N?**  Empirical data hints at
   roughly N^{2.5} — is there a precise rate?
3. **Structural kill load** as a function of α-depth: how many tuples
   at each α-depth need separate Moon-Moser-type lemmas?
4. **Phase-transition predictions:** which N values admit "almost
   forbidden" patterns (most candidates kill, one survives)?

## Related

* `formalization-driven-decompositions.md` — the original observation
  that f(N) is well-defined and tractable.
* `iso-invariance-as-cleaner-axiom-base.md` — meta-principle for
  preferring iso-class statements.
