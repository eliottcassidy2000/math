---
id: THM-1450
title: "LINEARITY TRIVIALISES VC, WHICH IS WHY THE ODD SIDE IS THE HARD ONE — and the odd side's VC=1 locus is exactly TWO tournaments at every n. (A) On the EVEN/linear side, for any F_2-linear code C a coordinate set is shattered iff the generator columns on it are independent, so VC(C) = dim(C) and the shattered sets are the INDEPENDENT SETS OF A MATROID. For the cut space of K_n — which IS the switching group — the matroid is graphic, so VC = n−1 and the shattered sets are exactly the FORESTS; verified n = 3..6 against brute-force enumeration. Cosets inherit this, so EVERY SWITCHING CLASS has VC = n−1 with the forests as its shattered sets. VC measures nothing there beyond the dimension. (B) On the ODD/nonlinear side — the out-neighbourhood system of one tournament — VC ≤ ⌊log₂ n⌋ trivially (2^k patterns need 2^k vertices) and the bound is ATTAINED for every n by putting a transitive tournament on the shattered set and assigning the remaining 2^k − k patterns freely; explicit n=8 witness with VC=3 exhibited. (C) THE VC=1 LOCUS IS EXACTLY TWO CLASSES at every n = 4,5,6,7: the transitive TT_n, and C₃ ⊕ TT_{n−3} (a 3-cycle on TOP of a transitive tail). The strong-component reduction PROVES cross-component arcs can never witness VC ≥ 2, reducing the classification to one residual statement about strong tournaments on ≥ 4 vertices. (D) VC is iso-invariant but NOT switching-invariant, placing it in the MIDDLE tier with H, min-FAS and the cyclic-triangle count — and by THM-1420 the odd side admits no F_2-linear invariants at all, so it can never be pushed into the trivial linear case."
status: >
  (A) PROVED (shattering of a code = independence of generator columns; the graphic-matroid
  identification is standard) and VERIFIED-EXACT n = 3..6, with an independent brute-force
  enumeration of the code agreeing on every subset.
  (B) Upper bound PROVED (counting); lower bound PROVED by explicit construction, with an
  n = 8, VC = 3 witness exhibited (code 202356).
  (C) VERIFIED-EXACT n = 4,5,6,7 (exactly two classes, identified); the strong-component
  reduction is PROVED; ONE residual gap remains — see Honest scope. NOT tested at n = 8.
  (D) VERIFIED-EXACT n = 4..7 by applying all 2^{n−1} switchings to every iso class.
source: mac-mini-2026-07-20-S131 (owner: "think vapnik-chervonenkis and chase the high
  leverage question, see the relation between odd valued functions and tournament adjacent
  ideas. they both relate also to even concepts like even graphs and even functions")
depends_on:
  - THM-1420  # no F_2-linear tournament invariants -- why (D) cannot be escaped
related:
  - THM-1440  # the sin/cos parity axis; the switching tier picture
  - THM-474   # tilings ARE switching classes
  - THM-1425  # the OCF as an independent-set sum (a DIFFERENT notion of independence -- see scope)
script: 04-computation/vc_dimension_odd_even_macmini_S131.py,
        04-computation/vc_one_classes_macmini_S131.py (+ .outs)
---

# THM-1450 — linearity trivialises VC, and the VC=1 locus is two tournaments

**One line.** The repo's odd/even axis is, for VC theory, exactly the nonlinear/linear axis:
on the even side VC collapses to the dimension and shattering becomes matroid independence;
on the odd side it is a genuine combinatorial quantity — and THM-1420 says the odd side can
never be linearised into the easy case.

## (A) The even/linear side: VC is the dimension, shattering is matroid independence

For an `F₂`-linear code `C ⊆ F₂^N` with generator rows `G`, a coordinate set `S` is shattered
iff `proj_S(C) = F₂^S` iff the **columns** of `G` on `S` are independent. Hence

> `VC(C) = dim(C)`, and the shattered sets are exactly the **independent sets of the matroid
> represented by `G`**.

Applied to `K_n`, where the cut space **is** the switching group (THM-474, THM-1415):

| `n` | edges | `dim` cut | `VC(cut)` | shattered `=` forests? | `dim` cycle | `VC(cycle)` | brute-force agrees |
|---|---|---|---|---|---|---|---|
| 3 | 3 | 2 | 2 | yes | 1 | 1 | yes |
| 4 | 6 | 3 | 3 | yes | 3 | 3 | yes |
| 5 | 10 | 4 | 4 | yes | 6 | 6 | yes |
| 6 | 15 | 5 | 5 | yes | 10 | 10 | yes |

> **The shattered sets of the switching group are exactly the FORESTS of `K_n`, and
> `VC = n−1`.** Translation does not change surjectivity of a projection, so **every
> switching class** (a coset) has the same VC dimension and the same shattered sets.

The even graphs (cycle space) likewise have `VC = C(n,2) − n + 1`, with the cographic matroid's
independent sets as their shattered sets. On this side VC is a restatement of the dimension.

## (B) The odd/nonlinear side: `max VC = ⌊log₂ n⌋`, attained

Take the out-neighbourhood system `{N⁺(v)}` of a single tournament.

- **Upper bound.** Shattering a `k`-set requires `2^k` distinct patterns and there are only
  `n` vertices, so `VC ≤ ⌊log₂ n⌋`.
- **Lower bound (construction).** Put a **transitive** tournament on `S = {0,…,k−1}`; its
  in-`S` out-neighbourhoods are the `k` distinct nested sets `{1,…,k−1} ⊃ ⋯ ⊃ ∅`. Assign the
  remaining `2^k − k` patterns freely to outside vertices — every assignment is realisable,
  since arcs between `S` and its complement are unconstrained.

Verified `n = 3..7` (max VC `= 1,2,2,2,2 = ⌊log₂ n⌋`), plus an explicit `n = 8` witness with
`VC = 3` (code `202356`).

## (C) The VC=1 locus is exactly two tournaments

| `n` | classes | VC distribution |
|---|---|---|
| 4 | 4 | `VC=1`: **2**, `VC=2`: 2 |
| 5 | 12 | `VC=1`: **2**, `VC=2`: 10 |
| 6 | 56 | `VC=1`: **2**, `VC=2`: 54 |
| 7 | 456 | `VC=1`: **2**, `VC=2`: 454 |

> **At every `n = 4..7` exactly two classes have `VC = 1`: the transitive `TT_n`, and
> `C₃ ⊕ TT_{n−3}` — a 3-cycle sitting on TOP of a transitive tail** (score sequences
> `[0,1,…,n−4, n−2, n−2, n−2]`, exactly one 3-cycle).

**The criterion.** `VC ≥ 2` iff some arc `a→b` admits both a `u` with `u→a, b→u` and a `w`
with `w→a, w→b` — verified to match `VC ≥ 2` on every class at `n = 4..7`.

**The strong-component reduction (PROVED).** Order the strong components `S_1 → ⋯ → S_k`.
If `a ∈ S_i`, `b ∈ S_j` with `i < j`, then `u→a` forces `u ∈ S_1..S_i` and `b→u` forces
`u ∈ S_j..S_k`, so `u ∈ S_i ∩ S_j = ∅`: **cross-component arcs can never witness `VC ≥ 2`.**
The same argument puts `u` in `S_i`, so `u→a→b→u` is a **3-cycle inside `S_i`**. And if
`i ≥ 2`, any vertex of an earlier component beats everything in `S_i`, so `w` exists free.
Therefore:

> `VC ≥ 2` **iff** some component `S_i` with `i ≥ 2` has `≥ 3` vertices, **or** `S_1` contains
> a 3-cycle `a,b,u` together with a `w ∈ S_1` beating both `a` and `b`.

This immediately gives the two families: all components singletons ⟹ transitive; and
`|S_1| = 3` with the rest singletons ⟹ `C₃ ⊕ TT_{n−3}`, which fails the internal condition
because in `C₃` the only vertex beating `x` is `z`, and `z` does not beat `y`. **The asymmetry
is real**: `TT_{n−3} ⊕ C₃` (3-cycle at the *bottom*) has `VC = 2`, because then a vertex above
beats both.

## (D) Where VC sits: the middle tier

| invariant | iso-invariant | switching-invariant |
|---|---|---|
| any `F₂`-linear functional | **no** — none exist (THM-1420) | — |
| `H`, min-FAS, cyclic triangles, **VC** | yes | **no** |
| skew-Seidel spectrum (THM-1440) | yes | yes |

VC changes under switching for `4/4`, `12/12`, `22/56`, `42/456` classes at `n = 4..7` — so it
is **not** switching-invariant, and joins the middle tier. Note the fraction *falls* with `n`
(100% → 9%), which is unexplained.

## The synthesis

**Linearity trivialises VC.** The even side is a code, so its VC dimension is its dimension and
its shattered sets are a matroid — nothing is measured. The odd side is not a code, its VC
dimension is a real quantity with a `⌊log₂ n⌋` ceiling and a two-element floor locus, and by
**THM-1420 it cannot be linearised**: there are no `F₂`-linear tournament invariants at all,
so no change of coordinates moves the odd side into the easy case. That is the precise sense
in which "odd-valued functions" are the hard side of the repo's odd/even axis.

## Honest scope

- **(C) has one residual gap.** The reduction is proved, but the classification needs: *every
  strong tournament on `≥ 4` vertices contains a 3-cycle `a→b→u→a` together with a `w` beating
  both `a` and `b`.* This is verified at `n ≤ 7` and looks routine, but **is not proved here**.
  Until it is, "exactly two" is verified, not proved.
- **(C) is untested at `n = 8`.** The pattern "exactly 2 at every `n`" rests on four data points.
- (A) is standard coding/matroid theory, re-derived and verified, not claimed as new. The
  content is the *application*: the switching group's shattered sets are the forests.
- **(D)'s falling fraction is an observation with no explanation offered.**
- **The word "independent" is doing two unrelated jobs** in this repo — matroid independence
  (here) and graph independence in the OCF's independent-set sum (THM-1425). No connection is
  claimed and none is known; do not cross-link them.
- The `⌊log₂ n⌋` ceiling puts this near the Erdős/Schütte property-`S_k` circle of ideas
  (measured here: max Schütte `k` = 1,1,1,2 at `n = 4..7`), but **no theorem connecting VC
  dimension to property `S_k` is claimed** — shattering is strictly stronger than `S_k`, which
  asks only for the all-of-`S` pattern.

*Artifacts:* `04-computation/vc_dimension_odd_even_macmini_S131.py`,
`vc_one_classes_macmini_S131.py` (+outs).
