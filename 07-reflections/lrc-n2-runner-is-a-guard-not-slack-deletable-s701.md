---
source: monad-explorer-2026-06-06-S701
status: PROVED (elementary) headline + VERIFIED census (exact, even n=4..20)
        Corrects the S700 (HYP-2259) handoff: the self-conjugate runner n/2 is
        NOT slack-deletable; it is the unique maximal GUARD. LRC(14) NOT proved.
tags: [LRC, n14, self-conjugate-runner, negation-fixed-point, guard, binder,
       deletion-criticality, doubling, half-system, even-fold, units, reframe]
---

# The second apex n/2 is a GUARD, not a slack runner: deleting it doubles the gap

## One-line correction of my own S700 handoff

S700 (HYP-2259) observed that the self-conjugate runner `n/2` is **slack**
(distance `1/2`, never binding) at every tight division-grid time of `AP_n`, and
conjectured: *slack off-grid too ⟹ `n/2` is deletable ⟹ even-n LRC reduces to the
odd-n single-apex residual.*

**That reduction fails, and S701 says exactly why.** Slack at the achieved
optimum does **not** imply deletable. Deleting `n/2` from `AP_n` does not preserve
the gap `1/n` — it **doubles** it:

> **Theorem (S701, proved, elementary).** For even `n`, `M(AP_n) = 1/n` but
> `M(AP_n \ {n/2}) = 2/n = 1/(n/2)`. Among all single-runner deletions, `n/2`
> gives the **unique maximal** gap.

So `n/2` is **load-bearing of the maximal degree**, even though it is the
**most-slack** runner at the optimum. It is a *guard*, not a *binder*.

## The proof (two lines)

Write `n = 2m`, so the runner is `m` and `R = AP_n \ {m} = {1,…,2m-1}\{m}`.

- **Lower bound.** At `t = 1/m`: for any `a ∈ R`, `a` is not a multiple of `m`
  (the only multiple of `m` in `1..2m-1` is `m` itself, which is deleted), so
  `‖a/m‖ ≥ 1/m = 2/n`. Hence `M(R) ≥ 2/n`.
- **Upper bound.** The even runners `{2,4,…,2m-2} = 2·AP_m ⊆ R` form a sub-config,
  and adding runners can only lower the gap, so `M(R) ≤ M(2·AP_m) = M(AP_m) =
  1/m = 2/n` (scaling speeds is gap-invariant).

Equality: `M(R) = 2/n`. ∎ Verified exactly (`Fraction`) for `n = 4,6,…,20`.

## Why "slack" and "load-bearing" are not in conflict — the guard mechanism

The optima of `AP_n` live on the **fine grid** `t = j/n` (`gcd(j,n)=1`, `j` odd).
There `‖(n/2)·j/n‖ = ‖j/2‖ = 1/2`: maximally slack, never a binder.

The optimum of `R = AP_n \ {n/2}` jumps to the **coarse grid** `t = j/(n/2)`.
There `‖(n/2)·j/(n/2)‖ = ‖j‖ = 0`: the runner `n/2` sits **on the origin**.

So `n/2` is the one runner that, while sitting harmlessly at `1/2` on the fine
grid, **single-handedly occupies the origin across the entire coarse grid** — it
*guards the gate* to the doubled optimum. Remove the guard and the gap collapses
upward from `1/n` to `2/n`. Slack at the achieved optimum, lethal at the optimum
it forbids: that is the definition of a guard.

## The full deletion census reveals a clean dichotomy

Computing `M(AP_n \ {a})` for **every** `a` (script `lrc_deletion_census_s701b`):

1. **`AP_n` is single-deletion-critical.** *Every* single deletion strictly
   raises `M` above `1/n` (verified n=4..16). No runner is redundant.
2. **`n/2` is the unique maximal guard.** It gives `M = 2/n`; every other deletion
   gives strictly less. (`a=[n/2]` is the unique argmax at every tested even n.)
3. **Binder/guard dichotomy on the grid.** Over the fine-grid optima, the runners
   that bind at some optimum are **exactly the units `(Z/n)^×`** (the binders at
   `t=j/n` are the negation pair `±j^{-1}`, both units — S700). The **non-units**
   (the even speeds) **never bind at any grid optimum** — they are pure off-grid
   guards. `n/2` is the maximal-slack non-unit (`1/2` at every optimum) and the
   maximal guard.

So `AP_n` splits its `n-1` runners into a **multiplicative binding skeleton** (the
`φ(n)` units, which carry the `1/n` floor on the fine grid) and a **guard set**
(the non-units, essential only off-grid). The self-conjugate `n/2` is the apex of
the guard set.

### The deletion profile is exactly known on the upper half — and it pins the peak at n/2

> **Lemma (proved). For every `a` with `n/2 ≤ a ≤ n-1`,  `M(AP_n \ {a}) = 1/a`,**
> realized at `t = 1/a`.

*Proof.* **(≤)** `{1,…,a-1} = AP_a ⊆ AP_n\{a}`, and `M(AP_a) = 1/a`; a superset of
runners has no larger gap, so `M(AP_n\{a}) ≤ 1/a`. **(≥)** At `t = 1/a`: for
`v ∈ {1,…,a-1}`, `‖v/a‖ = v/a ≥ 1/a`; for `v ∈ {a+1,…,n-1}` we have
`v-a ∈ {1,…,n-1-a} ⊆ {1,…,a-1}` (since `a ≥ n/2 ⇒ n-1-a ≤ a-1`), so
`‖v/a‖ = ‖(v-a)/a‖ ≥ 1/a`. Hence every remaining runner is `≥ 1/a` at `t=1/a`. ∎

So on the upper half the deletion profile is the **decreasing** function `1/a`,
maximized at its left endpoint `a = n/2`, value `1/(n/2) = 2/n` — recovering the
main theorem as the *peak* of a clean monotone law. The lower half `a < n/2`
stays strictly below: `M(AP_n\{a}) < 2/n` (verified n=4..18; the values form a
`2/(odd)` off-grid pinch ladder, n=14: `a=1→2/15, 2,3→2/17, 4,5→2/19, 6→2/23`).
Together: **`n/2` is the unique argmax of the single-deletion profile** — proved
on the upper half, verified on the lower. The self-conjugate runner is exactly the
crossover where the `1/a` law tops out.

## The genuine even→odd link (what survives of the S700 dream)

The reduction even→odd does **not** go by deleting `n/2`. But an even→odd link
*is* present, hiding in the upper-bound step: the **even sublattice**
`{2,4,…,n-2} = 2·AP_{n/2}` is a gap-invariant scaled copy of `AP_{n/2}`. For
`n=14` that is `2·AP_7` — the **odd prime** case `AP_7` (`M=1/7`). This is exactly
the "even-fold" `M_n(S) ≤ M(fold S)` of S555 (`fold` = halve the evens), now with
a precise role for the center: `n/2` is the unique runner *not* in the even
sublattice's preimage skeleton that nonetheless guards the fine grid against
collapsing onto it. The honest statement of the reduction:

> Even-n `AP` tightness `= ` (fine-grid floor `1/n` carried by the unit binders)
> `+ ` (the guard `n/2` blocking the coarse `2/n` optimum of the even half-system
> `2·AP_{n/2}`).

The half-system `AP_{n/2}` is where the odd case enters — but as a **ceiling the
guard defends**, not as a residual you reach by deletion.

## Connections

- **Project doubling theme.** The gap doubling `1/n → 2/n` on removing the unique
  self-conjugate runner is a clean instance of the recurring `n → n/2`
  "half-system / Mode-B" descent (Cayley–Dickson tower, `2·AP_{n/2}` skeleton,
  doubling law S612). The center `n/2` is the obstruction that keeps the full
  system off its half.
- **Sieve/apex picture (S679/S700).** The `0`-apex is the external sieve
  obstruction; the `n/2`-apex is now sharply characterized as the **internal
  guard** of the even half-grid. Two fixed points, two distinct mechanisms
  (external sieve vs internal guard) — the cleanest separation yet.
- **`Res_C` machinery (HYP-2164..2253).** `C = 2n-1` is odd, blind to `n/2`
  (S700). The guard role explains *why* the Lebesgue/owner machinery keeps finding
  `n/2` as a measure-zero wall (HYP-2252/2253): a guard contributes a
  codimension-1 origin-locus on the coarse grid, exactly a measure-zero wall.

## Status

- **PROVED (elementary):** `M(AP_n \ {n/2}) = 2/n` for all even `n`; the upper-half
  deletion law `M(AP_n\{a}) = 1/a` for `n/2 ≤ a ≤ n-1` (hence `n/2` is the unique
  argmax of the single-deletion profile, peak `2/n`); binders = units (S700).
- **VERIFIED (exact, n=4..18):** full deletion census; single-deletion criticality
  of `AP_n`; lower-half `M(AP_n\{a}) < 2/n` (the `2/(odd)` pinch ladder); the
  binder/guard dichotomy (non-units never bind on the grid).
- **CORRECTED:** the S700 handoff "slack ⟹ deletable ⟹ even→odd by deletion" is
  false; replaced by the guard picture.
- **NOT proved:** LRC(14). The off-grid recovery for *arbitrary* configs is
  untouched; this only dissects the canonical `AP_n` tight example.

Artifacts: `04-computation/lrc_n2_runner_loadbearing_s701.py`,
`04-computation/lrc_deletion_census_s701b.py` (+ `.out`s), HYP-2260, T754.
Builds on S700 (HYP-2259), S679 (HYP-2255), HYP-2231, S555 even-fold, S612
doubling law, THM-401 (C=2n-1).
