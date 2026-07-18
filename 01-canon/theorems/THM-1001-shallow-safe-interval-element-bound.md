---
id: THM-1001
title: The safe-interval element bound and uniform single-coordinate shallow winding exclusion
status: PROVED (elementary danger-tooth covering; the single-coordinate corollary is a proved uniform-in-height exclusion, delta-bound + a per-residue at-most-one exact value)
source: mac-mini-2026-07-17-S109 (shallow-branch winding, uniform in lift height)
depends_on:
  - THM-769   # shallow branch = complete nonzero residue system mod 13 (used, NOT re-inferred from one point)
  - THM-759   # ratio bound (this refines it via the exact safe interval)
  - LRCUpTo13 # M(C) >= 1/12 for the 11-element complement
related:
  - THM-770   # finite full-residue classification through lift height 12 (this is the uniform-in-height complement on the single-coordinate slice)
  - THM-763   # global finite height for primitive tight instances
  - HYP-6800  # the n=12 sporadic branch
  - HYP-6820  # the uniformity audit
  - HYP-7330  # this session's assembly
verification: 04-computation/lrc13_shallow_safe_interval_bound_macmini_S109.py
  (+ 05-knowledge/results/lrc13_shallow_safe_interval_bound_macmini_S109.out)
---

# THM-1001 — The safe-interval element bound (shallow branch)

**One line.** In a tight full-residue 12-set, no single speed can be much larger
than the "safe geometry" of the other eleven allows: `w ≤ 2/(13·δ(A\{w}))`, where
`δ(C)` is the widest arc on which the complement stays `1/13`-lonely. A larger `w`
would have a danger-comb too fine to plug that arc, opening a `>1/13` witness.
Consequence: **winding a single residue to any height never produces a tight set**
— a uniform-in-height exclusion where THM-770 is finite (height ≤ 12).

## Setup

`δ := 1/13`, `φ_C(t) = min_{v∈C} ‖vt‖`. For a finite speed set `C` put

```
δ(C) = length of the largest open arc J ⊆ ℝ/ℤ on which φ_C(t) > 1/13.
```

By THM-769 §2 a **shallow** tight 12-set `A` (`M(A)=1/13`, no speed divisible by
13) is exactly a **complete nonzero residue system mod 13**; every `a/13` is a
global maximizer. (We use THM-769's global characterization — completeness is
**not** inferred from a single maximizing point; cf. the S107 correction in
HYP-6820.)

## (A) The safe-interval element bound (PROVED)

> **Lemma.** Let `A` be a shallow tight 12-set and `w ∈ A`, `C = A\{w}`. Since
> `|C| = 11`, `M(C) ≥ 1/12 > 1/13` (LRC(12)), so `δ(C) > 0`. Then
> **`w ≤ 2/(13·δ(C))`.**

*Proof.* Suppose `w > 2/(13 δ(C))`, i.e. `δ(C) > 2/(13w)`. Let `J` be a largest
arc of `{φ_C > 1/13}`, `|J| = δ(C)`. The `1/13`-danger set of `w`,
`{t : ‖wt‖ ≤ 1/13}`, is a union of closed teeth each of width `2/(13w)`, spaced
`1/w`. If `J` contained no `w`-safe point it would lie inside the danger set, and
being connected inside a union of teeth separated by safe gaps it would lie inside
one tooth, forcing `δ(C) = |J| ≤ 2/(13w)` — contradicting `δ(C) > 2/(13w)`. So `J`
has a point `t*` with `‖w t*‖ > 1/13`; and `t* ∈ J` gives `φ_C(t*) > 1/13`. Hence
`φ_A(t*) = min(φ_C(t*), ‖w t*‖) > 1/13`, so `M(A) > 1/13`, contradicting tightness.
Therefore `w ≤ 2/(13 δ(C))`. ∎

**This refines THM-759.** The interval around `C`'s maximizer has width
`≥ 2(M(C) − 1/13)/max(C) ≥ 2(1/12 − 1/13)/max(C) = 1/(78·max C)`, so
`δ(C) ≥ 1/(78 max C)` and (A) yields `w ≤ 2/(13)·78·max(C) = 12·max(C)` — exactly
the ratio bound `a_12 ≤ 12 a_11`. The **exact** `δ(C)` is usually far larger than
the crude `1/(78 max C)`, giving a strictly sharper, configuration-aware bound.
Verified: `{1,…,12}` satisfies (A) at every `w` (bounds `4,8,12,16,20,24` for
`w=1..6`; `w` comfortably below each).

## (B) Uniform single-coordinate winding exclusion (PROVED, all heights)

> **Corollary.** For every residue `j ∈ {1,…,12}` and every integer `k ≥ 1`,
> `M(\{1,…,12\}\setminus\{j\} ∪ \{j+13k\}) > 1/13`. Hence **no** tight full-residue
> 12-set is obtained by winding a single coordinate — for **all** heights `k`.

*Proof.* The set is full-residue (`residues = {1,…,12}`); if it were tight, (A)
applied to `w = j+13k` with fixed complement `C = {1,…,12}\{j}` would force
`j+13k ≤ 2/(13 δ_j)`, `δ_j := δ(C)`. The exact `δ_j` (computed) gives thresholds
`2/(13δ_j) ∈ \{4,8,12,16,20,24,8,11,10,14,18,22\}` for `j=1..12`. For every `j`
except `j=5,6` already the smallest wound value `j+13` exceeds the threshold, so
**all** `k ≥ 1` are excluded by (A) alone; for `j=5,6` the single value below
threshold (`w=18, 19`) is checked exactly (`M > 1/13`). ∎

This is the **uniform-in-height** complement to THM-770's finite full-residue
classification (which certifies lift heights `≤ 12`): the single-coordinate slice
is closed for **every** height, including the `> 12` heights outside THM-770's box
and unreduced by its gcd-descent (descent lands on the primitive quotient, not a
lower height).

## (C) Reduction of the shallow residual

For any bound `H`, the 11-element complements `C` with all lift heights `≤ H` form
a **finite** family, so `δ_min(H) := min_C δ(C) > 0`. By (A), a shallow tight set
in which eleven coordinates have height `≤ H` has its twelfth coordinate bounded by
`2/(13 δ_min(H))` — hence lies in a finite, exactly-checkable set. Therefore:

> **Any shallow tight 12-set outside THM-770's height-12 box must have at least
> two coordinates wound above height 12.** The single-free-coordinate direction is
> closed uniformly.

## Honest scope

- **Proved:** the element bound (A) (elementary), its refinement of THM-759, and
  the uniform single-coordinate exclusion (B) (all heights). (C) is an immediate
  finiteness reduction.
- **Not proved here:** the full shallow branch (two or more coordinates wound
  high simultaneously), and nothing about the **deep** branch (THM-769 §3,
  THM-774/775/776/836). The bound (A) reduces to THM-759 in the worst case, so it
  does not by itself bound the multi-coordinate shallow spread — that residual
  needs the exact `δ(C)` controlled uniformly, or the descent/CSP machinery.
- The whole `n=12` sporadic-branch emptiness remains OPEN (HYP-6800/6820); this
  closes one uniform-in-height slice of the shallow half and sharpens the element
  bound.

*Artifacts:* `04-computation/lrc13_shallow_safe_interval_bound_macmini_S109.py`
(+out). Credits: THM-769 (shallow = full residue), THM-759 (the ratio bound
refined), THM-770 (the finite-height classification complemented), THM-763 (global
finite height), the codex-S64 audit (quantifier discipline).
