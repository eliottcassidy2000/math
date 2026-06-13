---
source: monad-explorer-2026-06-06-S702
status: LOWER BOUND PROVED (elementary, explicit witness, all a<n/2);
        UPPER BOUND VERIFIED exact (n=4..48, 551/551; unified profile 777/777,
        n=4..40) + proved in the principal sub-case. CORRECTS the S701/HYP-2260
        "2/(odd) pinch ladder" claim. LRC(14) NOT proved.
tags: [LRC, AP_n, deletion-criticality, lower-half, binder, guard, coprime,
       next-coprime-denominator, three-gap, doubling, profile, reframe]
---

# The lower-half deletion ladder of AP_n is "2 / next denominator coprime to a"

## What S701 left open, and the wrong guess it left

S701 (HYP-2260) proved the **upper half** of the single-deletion profile of
`AP_n = {1,…,n-1}` exactly:

> `M(AP_n \ {a}) = 1/a` for every `n/2 ≤ a ≤ n-1`  (witness `t = 1/a`).

For the **lower half** `a < n/2` it only *verified* values and described them as a
"`2/(odd)` pinch ladder (n=14: `a=1→2/15, 2,3→2/17, 4,5→2/19, 6→2/23`)". That
description is **wrong as stated** — its own census has `n=12, a=3 → 1/8 = 2/16`
(even denominator). The lower half is not a 2/odd ladder. This session found what
it actually is.

## The law (S702)

> **For `1 ≤ a < n/2`:**
> `M(AP_n \ {a}) = 2 / D*(n,a)`, where **`D*(n,a) = min{ D ≥ n+a : gcd(D,a) = 1 }`**
> — the smallest integer `≥ n+a` that is **coprime to the deleted runner `a`**.

Special case `gcd(a,n)=1`: then `gcd(a, n+a)=gcd(a,n)=1`, so `D* = n+a` and
**`M = 2/(n+a)`** — a perfectly clean coprime law. The "2/odd" appearance at
`n=14` was an artifact: the *non*-coprime deletions there are `a∈{2,4,6}` (even),
and the next denominator coprime to an even `a` is forced odd — hence "2/odd"
*looked* like the rule when it is really "2 / next-coprime-to-`a`". At `n=12`,
`a=3` (coprime to nothing-mod-3) the next denominator coprime to `3` past `15`
is `16` — even — giving `1/8`, exactly the counterexample.

**Verified exactly** (`Fraction`/integer search) `551/551` lower-half cases for
`n = 4..48`, and the **unified profile** (both halves) `777/777` for `n = 4..40`.
Zero failures.

## Why it is true — the binder/guard picture, sharpened

Every optimum here is a **`c=2` two-binder optimum**: at the optimal `t = p/D*`
the gap is `2/D*`, realized by **two runners `b₁ < b₂` straddling the origin at
residues `±2`**, with `b₁ + b₂ = D*` (verified for the whole lower half). So the
denominator is literally the **sum of the two binding runners**, and maximizing
the gap `= 2/(b₁+b₂)` means **minimizing the binder sum** subject to feasibility.

The deleted runner `a` is exactly the obstruction that, if present, would sit at
residue `−1` (distance `1/D*`) and cap the gap at `1/D*`. Deleting it lets the gap
open to the next ring, `2/D*`.

### Lower bound — PROVED (elementary, explicit, all `a<n/2`)

Let `D = D*(n,a)` and choose `p ≡ a^{-1} (mod D)` (exists because `gcd(a,D)=1` by
definition of `D*`). Take `t = p/D`. For any runner `v ∈ {1,…,n-1}\{a}`:

- `vp ≡ 0 (mod D)`: impossible, since `gcd(p,D)=1` and `1 ≤ v ≤ n-1 < D`.
- `vp ≡ +1 (mod D) ⟺ v ≡ p^{-1} ≡ a`: that is the **deleted** runner. Excluded.
- `vp ≡ −1 (mod D) ⟺ v ≡ −a ≡ D-a`: but `D-a ∈ [n, D-1]` (since `D ≥ n+a`), so
  `v = D-a` is **out of range** `[1,n-1]`. Excluded.

Hence every surviving runner has `vp \bmod D ∉ {0, 1, D-1}`, i.e. residue
distance `≥ 2`, so `‖vt‖ ≥ 2/D`. Therefore `M(AP_n\{a}) ≥ 2/D*`. ∎

This *is* the whole reason coprimality enters: the construction needs `a` to be a
**unit mod D**. At `gcd(a,n)=1` the unit appears already at `D=n+a`; otherwise you
climb to the next `D` coprime to `a`. The lower-bound witness realizes `2/D*` in
all `551/551` tested cases.

### Upper bound — VERIFIED (`n ≤ 48`), proved in the principal case

The matching upper bound `M ≤ 2/D*` has **no small-subset certificate** (a search
over `|W|≤4` subsets returns none): it is a global statement, true but requiring
an analytic argument. Two pieces are clean:

- **Doubling case (coprime, `‖aθ‖ ≤ 1/(n+a)`).** Apply `M(AP_{n+a}) = 1/(n+a)`:
  some `w ∈ [1,n+a-1]` has `‖wθ‖ ≤ 1/(n+a)`. If `w ∈ R` we are done; if `w=a`
  then `2a ∈ R` (as `a<n/2`) has `‖2aθ‖ ≤ 2‖aθ‖ ≤ 2/(n+a)`. Done.
- **Residual.** The only remaining regime is `w ∈ {n,…,n+a-1}` *and*
  `‖aθ‖ ∈ (1/(n+a), 1/n]`. Continued-fraction analysis shows here `a` is the
  convergent of `θ` just below `n+j=w`, with `‖aθ‖ < 1/(n+j) ≤ 2/(n+a)`; the bound
  is confirmed by dense rational sampling (no violations) but the matching
  in-range witness is a three-gap statement left for a follow-up.

So: **lower bound proved, upper bound verified `n≤48` and proved in the principal
sub-case.** Honest status — not yet a closed THM.

## The complete single-deletion profile of AP_n (both halves)

Combining S701 (upper) and S702 (lower):

> `M(AP_n \ {a}) = ` &nbsp; `1/a` &nbsp; if `n/2 ≤ a ≤ n-1`
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `2/D*(n,a)` if `1 ≤ a < n/2`, &nbsp; `D*=min{D≥n+a: gcd(D,a)=1}`.

The two branches are genuinely different mechanisms: the upper half is a **`c=1`,
denominator-`a`** optimum (`t=1/a`, residues mod `a`); the lower half is a **`c=2`,
denominator-`(b₁+b₂)`** optimum.

### Corollaries (verified n=6..24)

- **Even `n`:** the unique maximal deletion is the **guard `n/2`**, `M=2/n`
  (recovers S701 — and now seen as the *peak* where the upper-half `1/a` law tops
  out, the lower-half `2/D*` ladder being strictly below it).
- **Odd `n`:** there is **no guard**; the maximal deletion is `M = 2/(n+1)`,
  attained as a **tie** at `a = 1` (lower, `D*=n+1`) **and** `a = (n+1)/2` (upper,
  `1/a`). The even/odd split of "which deletion helps most" is a clean parity
  fingerprint, complementary to S700's two-fixed-point picture.
- **All `n`:** the **minimal** deletion is the **fastest runner `a = n-1`**,
  `M = 1/(n-1)` — the *most load-bearing* single runner (its removal raises the gap
  the least). So deletion-criticality spans the band
  `[ 1/(n-1) , (2/n or 2/(n+1)) ]`: fastest runner at the floor, central
  runner/`a=1` at the ceiling.

## Connections

- **Sharpens the binder/guard dichotomy (S700/S701).** S701 showed grid-binders =
  units `(Z/n)^×`. Here the role of `gcd(a,n)` is made *quantitative*: the deletion
  gap is `2/(n+a)` exactly when `a` is a unit mod `n` (so a unit mod `n+a`), and
  degrades to `2/D*` (climb to the next unit modulus) precisely by the amount `a`
  fails to be a unit. Coprimality is the *engine* of the construction, not just a
  classification.
- **Doubling theme.** The lower-half optimum is a `c=2` (two-sided) structure, the
  upper half `c=1`; the guard `n/2` is the crossover. The `b₁+b₂=D*` law is a
  literal "two runners pin the origin" picture — the additive face (THM-401 pair
  sums) seen from the deletion side.
- **Res_C machinery.** `D*` is an *additive* modulus tied to `n+a`, transverse to
  the multiplicative shell `C=2n-1`; the `gcd(D,a)=1` gate is exactly the
  unit-condition the carry/owner machinery keeps re-encountering.

## Status / handoff

- **PROVED:** lower bound `M ≥ 2/D*` (explicit witness, all `a<n/2`); coprime
  upper-bound doubling case.
- **VERIFIED (exact):** the law `M = 2/D*` for `n=4..48` (551/551); unified profile
  `n=4..40` (777/777); `c=2`, `b₁+b₂=D*`; the corollaries.
- **OPEN:** the residual upper-bound three-gap statement (close it → full THM,
  reserve THM-411). Then: does the same "next-coprime-denominator" law govern
  lower-half deletions of the **other** tight families (`V*`, `2·AP`) — S701 showed
  the *guard* generalizes; does the whole ladder? And: **double** deletions
  `M(AP_n\{a,b})` — is there a two-variable `D*`?

Artifacts: `04-computation/lrc_lower_half_ladder_s702.py`,
`04-computation/lrc_lower_half_law_s702.py`,
`04-computation/lrc_full_deletion_profile_s702.py`,
`04-computation/lrc_upper_bound_probe_s702.py` (+ `.out`s).
Builds on S701 (HYP-2260), S700 (HYP-2259), THM-401, the S555 even-fold.
