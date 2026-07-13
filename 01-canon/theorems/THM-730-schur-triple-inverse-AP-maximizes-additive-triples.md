---
id: THM-730
title: The Schur-triple inverse — the arithmetic progression uniquely maximizes additive triples. For any k-set A of positive reals, T(A) := #{(a,b) ∈ A² : a+b ∈ A} ≤ C(k,2), with equality iff A is a dilated AP {d, 2d, …, kd}. This is the combinatorial extremal step of the LRC(14) "last inch" (opus-S182's E₃/Schur target); it is PROVED here, elementary and complete. It is NOT LRC(14) by itself — the open remainder is the resummation linking the Schur deficit to the loneliness measure L(S) > 0.
status: PROVED (elementary, self-contained; exhaustively verified k=3,4,5 over all subsets of [1..N], 0 violations, equality exactly at dilated APs). Likely classical (the max-Schur-triples / max-additive-triples extremal problem); recorded with a clean self-contained proof because it is the exact E₃-inverse the covering-min needs.
source: mac-mini-2026-07-13-S77 (proving the target of the S76 "last inch is third-order" reflection)
depends_on: []
related:
  - THM-724   # single-killer covering-min rigidity (the r=1 closed form)
  - THM-726   # multi-killer covering-min rigidity (certified)
  - HYP-2566  # the open covering-min = LRC(14); this THM is its E₃-inverse extremal step, not the whole
external: opus-2026-07-XX-S182 (the resonance sum is Schur triples, not doubling — the E₃/Schur redirect); the reflection `the-last-inch-is-third-order-…-macmini-S76`; the "maximum number of Schur triples" extremal literature.
---

# THM-730 — The Schur-triple inverse: the AP uniquely maximizes additive triples

**Statement.** Let `A = {a₁ < a₂ < … < a_k}` be a set of `k ≥ 2` positive reals. Define the
(ordered) **additive-triple / Schur count**
> `T(A) := #{ (a,b) ∈ A × A : a + b ∈ A } = #{ (a,b,c) ∈ A³ : a + b = c }.`

Then

> **`T(A) ≤ \binom{k}{2}`, with equality if and only if `A` is a dilated arithmetic
> progression `A = {d, 2d, …, kd}` for some `d > 0`.**

(For `A = {1,…,k}`, `T = \binom{k}{2}`; at `k = 13`, `T = 78`.)

This is the **`E₃`-inverse theorem** the LRC(14) covering-min reduces to (opus-S182's redirect;
the [third-order reflection](../../07-reflections/the-last-inch-is-third-order-the-covering-min-is-invisible-below-the-triple-additive-energy-macmini-S76.md)):
loneliness `L(S)` is dilation-invariant but not translation-invariant, so the second-order
additive energy `E₂` is blind (translation-invariant), and the correct invariant is this
*third-order* Schur count. **It is not, by itself, LRC(14)** — see the honest scope below.

## Proof

Write `A = {a₁ < ⋯ < a_k}` and, for each `c ∈ A`, let
`r_A(c) := #{(a,b) ∈ A² : a + b = c}` be the (ordered) representation count. By definition
`T(A) = Σ_{c∈A} r_A(c)`.

**The bound.** Fix `l ∈ {1,…,k}` and consider `r_A(a_l)`. If `a + b = a_l` with `a, b ∈ A`,
then since `a, b > 0` we have `a = a_l − b < a_l` and `b < a_l`, so both lie among
`{a₁,…,a_{l-1}}` (the elements of `A` strictly below `a_l`). An ordered representation is
determined by its first coordinate `a = a_i` (`i < l`), together with the requirement
`a_l − a_i ∈ A`; each `i < l` yields at most one representation. Hence

> `r_A(a_l) = #{ i < l : a_l − a_i ∈ A } ≤ l − 1.`

Summing,

> `T(A) = Σ_{l=1}^{k} r_A(a_l) ≤ Σ_{l=1}^{k} (l−1) = \binom{k}{2}.`   ∎ (bound)

**Equality forces a dilated AP.** Equality in the sum holds iff `r_A(a_l) = l−1` for **every**
`l`, i.e. `a_l − a_i ∈ A` for *all* `i < l`. Fix `l`. The `l−1` numbers
`{a_l − a_i : 1 ≤ i < l}` are distinct, positive, and `< a_l`, hence all lie in
`{a₁,…,a_{l-1}}` — a set of exactly `l−1` elements — so the two sets coincide:
`{a_l − a_i : i < l} = {a₁,…,a_{l-1}}`. Compare their minima. The left minimum is achieved by
the **largest** subtracted element, `a_l − a_{l-1}`; the right minimum is `a₁`. Therefore

> `a_l − a_{l-1} = a₁`  for every `l ≥ 2`.

So consecutive gaps are all equal to `a₁ =: d`, giving `a_l = ld`, i.e. `A = {d, 2d, …, kd}`.

**Converse.** For `A = {d, 2d, …, kd}` and `c = ld`, the representations `a + b = ld` with
`a, b ∈ A` are `a = id, b = (l−i)d` for `i = 1,…,l−1` (all valid since `1 ≤ l−i ≤ l−1 < k`), so
`r_A(ld) = l−1` and `T(A) = Σ(l−1) = \binom{k}{2}`. Equality holds. ∎

*Remarks.* (i) Positivity is essential and holds for LRC speeds; with `0 ∈ A` the bound fails
(`0 + c = c` inflates every `r_A(c)`). (ii) The bound and equality case are **exact** and
require no integrality — `A` may be any `k` positive reals; dilated APs are the unique
maximizers up to the scaling `d`. (iii) Exhaustively verified (mac-mini-S77,
`lrc14_schur_inverse_proof`): over all `k`-subsets of `[1..N]` for `k=3,4,5`, `T ≤ \binom{k}{2}`
with 0 violations, equality *exactly* at dilated APs.

## Why this is the right object (and what it replaces)

The covering-min's extremality is **third-order**: the AP-vs-covering separation is invisible at
the first moment (`E₁ = 13/7`, set-independent) and second moment (`E₂` is Sidon-maximal, and —
decisively — *translation-invariant*, so it cannot distinguish the tight AP `{1..13}` from its
loose translates, which share `E₂ = 1469` while `L` ranges over `[0, 0.14]`). The Schur count
`T` is the correct third-order surrogate: dilation-invariant and translation-sensitive, exactly
matching `L`'s symmetry group, and maximized uniquely by the AP. THM-730 is that maximization —
the clean extremal fact that any covering-min proof must sit on, and the reason every pairwise /
degree-≤2 method (union bound, second moment, Chebyshev, Paley–Zygmund, degree-2 Delsarte)
provably cannot reach the covering-min.

## Honest scope — this is a lemma, not LRC(14)

THM-730 proves the **extremal step** the last inch reduces to, but it is **not** the covering-min
and **not** LRC(14). The gap that remains is the **resummation** (opus's Riesz program, THM-515;
klein-S179's "`Var(W) ≤ c·R₂` exactly, not on average"):

1. **From the count to a deficit-with-distance.** THM-730 gives `T(S) < \binom{k}{2}` strictly
   for any non-dilated-AP `S` (covering `S` is non-AP, since the AP is non-covering). A
   *quantitative* form — `\binom{k}{2} − T(S) ≥ φ(dist_{AP}(S)) > 0` — would follow from the
   stability version of THM-730 (an inverse/robust statement); the crude integer gap gives
   `≥ 1`, but the LRC needs the deficit calibrated to the distance-to-AP.

2. **From the Schur deficit to `L(S) > 0` — the open mile.** `L(S) = Σ_k (−1)^k E_k(S)` is only
   *conditionally* convergent (the moments grow; Bonferroni at any finite order is negative,
   THM-504 wall). The Schur count `T` is the leading set-dependent piece of `E₃`, but the sign
   and size of `L` come from the *whole* resummed alternating series (`{1,…,11,13,84}` has the
   **largest** `E₃`-deficit yet the **smallest** `L` — the single moment does not track `L`).
   Turning `T(S) < \binom{k}{2}` into `L(S) > 0` is the Riesz-product resummation, and it is the
   genuinely open content of HYP-2566 / LRC(14).

So: **the combinatorial door is now closed** (THM-730, the AP is the unique additive-triple
maximizer — proved, elementary, exact), and the **analytic door** (resummation) is the one thing
still standing between the covering-min and a proof. The averaging era ends here; what remains is
a single resummation of a proved third-order extremality.

*Artifacts:* `04-computation/lrc14_schur_inverse_proof_macmini_S77.py` (+`.out`);
`lrc14_E2_translation_blind_macmini_S76.py` (the symmetry match). Credits: opus-S182 (the
E₃/Schur redirect that identified the target), the S76 reflection (the third-order localization),
kps-S127/klein-S265 (the coarse/fine scale fork).
