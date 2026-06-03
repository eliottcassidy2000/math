---
source: claude-2026-06-03-S606
status: deep study of twisted involutions + verified reduction mechanism; n=14 NOT closed (obstruction pinned)
tags: [lonely-runner, n14, twisted-involution, sign-reversing, nerve, helly, redei, observer-coupling, half-shift, negation, vitali-wall]
---

# Twisted involutions for LRC@14: which ones can certify, and why the obvious ones can't

Goal: a proof or improvement for LRC@14, via S599's resolution template — turn
the permanent-shaped `p_0 = Σ_S (−1)^{|S|} m_S` into a certifiable (determinant)
form by a **sign-reversing involution**, the way Rédei does over GF(2). The
session's real output is a *deep* understanding of which involutions can do this,
a verified reduction mechanism, and a precise statement of the remaining obstacle.

## 1. Geometric clock involutions cannot certify (they are measure-preserving)

There are exactly two natural involutions of the clock `ℝ/ℤ` that respect the
speed structure:

- **Negation `t ↦ −t`.** `‖v(−t)‖ = ‖vt‖`, so `depth(−t) = depth(t)` *exactly*
  [verified]. Negation **fixes** the lonely set `Λ` and is **measure-preserving**
  — hence **sign-preserving** on `p_0`. (This is the `κ`-even fact of S605 in
  geometric clothing.) It is the *only* symmetry of `Λ`.
- **Half-shift `τ: t ↦ t+½`** (even `n`). It fixes even-runner distances and
  reflects odd ones (`‖vt‖ ↦ ½−‖vt‖`, origin↔antipode). It is measure-preserving
  (a translation) but maps `Λ` to a *different* set — `{even-far} ∩
  {odd-antipode-far}` — so it is **not even a symmetry of `Λ`** [verified: lonely
  status differs on a positive fraction of points]. Useful for reformulation,
  useless for cancellation.

**Lesson.** No geometric involution is sign-reversing on the loneliness
*measure*. A measure cannot be cancelled by a measure-preserving map. So the
Rédei trick has no geometric analogue on `p_0` directly.

## 2. The correct twisted involution acts on *subsets*, not the clock

The signs live in the inclusion–exclusion sum itself:
`p_0 = Σ_{S}(−1)^{|S|} m_S`, `m_S = meas(⋂_{v∈S} F_v)`. The right involution flips
subset parity. The **pivot involution** `σ_p : S ↔ S △ {p}` is sign-reversing,
and a pair cancels iff `m_S = m_{S∪p}`, i.e. iff the pivot `p` is **redundant**
on `S` (`⋂_S F ⊆ F_p`). The survivors — subsets where the pivot is *not*
redundant — carry the whole value:

> [verified] `Σ_{survivors} (−1)^{|S|}(m_S − m_{S∪p}) = p_0` exactly,
> for `(1,2,3,4), (1,3,4,7), (1,3,4,5,9)` (≈ 6–12 survivor pairs vs `2^n` subsets).

Iterating pivots `σ_{p_1}, σ_{p_2}, …` drives `p_0` to the **irredundant core** —
the **nerve** of the cover `{F_v}` (nerve lemma), which for LRC@14 is exactly the
**two-block Helly** structure codex-S599 isolates. This is the LRC realization of
S599's permanent → determinant template: the subset involution is the
"Gaussian elimination" that collapses the `(★)` permanent toward a certificate.

## 3. Redundancy is the fuel — and the Vitali wall is still there

The involution only cancels where there is redundancy. [verified]

- **Dissociated `(1,2,4)`: zero redundancy.** The involution is trivial, nothing
  cancels, `p_0` stays at the positive independence value — *never tight*. This
  recovers the S604 dichotomy (relation-poor ⇒ not tight) *through the
  involution*: no redundancy ⇒ no cancellation ⇒ no collapse to `0`.
- **Structured sets have redundancy**, but the *amount* does **not** cleanly
  decide tightness: non-tight `(1,2,3,5)` has `5/28` redundant pairs, more than
  tight `(1,3,4,7)`'s `3/28`. So redundancy-count is not the invariant — the
  **Vitali wall** (S603) persists. The residual that the involution cannot
  cancel is the **quantitative observer-coupling** (S580o), not a combinatorial
  count.

## 4. Why Rédei is solved and LRC is not — stated exactly

Rédei's reversal `T ↦ T^op` is sign-reversing on the **GF(2) parity** of the
Ham-path sum (a *signed* object); it cancels non-self-converse classes, the
self-converse fixed points have forced parity, and `#HamPaths` is odd — proved.
The *same* reversal, applied to LRC's **unsigned measure** `p_0`, is the
sign-*preserving* negation of §1. **That is the precise reason the two problems
sit on opposite sides of the line:** Rédei's natural object is signed (parity),
LRC's is an unsigned measure, and only the subset/nerve involution restores the
signs — at the price of leaving the un-cancellable coupling residual.

## 5. The n=14 target, pinned

LRC@14's worry-set reduces to `2^6 = 64` self-converse round classes (S576o) —
*exactly Rédei's fixed-point set*. Closing it needs **a second, sign-reversing
involution on those 64** (a twist on the `2^6` cube) — equivalently a nerve
collapse of the two-block Helly core (codex-S599). The subset pivot involution
(§2) is that mechanism; what it cannot cancel is the observer-coupling defect
(S580o), active for `m ≥ 5`. So the concrete, well-posed target is:

> **Find an involution on the 64 self-converse round classes that pairs them
> sign-reversingly modulo the observer coupling** — or show the two-block Helly
> core is its own certificate (a singleton/pair wall per codex-S599) so that no
> residual survives.

That is a *finite* object (64 classes / bounded Helly witnesses), not an
unbounded-speed quantifier — the involution lens converts the astronomically
large speed bound into a finite cube-symmetry search.

## Honest status

Not a proof of LRC@14. What is new and verified: (i) the geometric involutions
are provably sign-preserving / non-symmetries, so they cannot work; (ii) the
subset pivot involution *is* a genuine sign-reversing reduction of `p_0` to the
Helly nerve; (iii) redundancy is its fuel, absent for dissociated sets, present
but not tightness-deciding otherwise. The remaining gap is exactly the
observer-coupling residual — and it now has a finite, cube-shaped target.

## Next

1. **Parametrize the 64 self-converse round classes (HYP-2094/2086) as a `2^6`
   cube** and search bit-flip involutions for one that is sign-reversing on the
   loneliness defect — a finite search.
2. **Compute the irredundant nerve core for the n=14 AP and `V*`** directly
   (iterate the pivot involution) and check it equals the two-block Helly
   witnesses — a bridge from this reduction to codex-S599's certificates.
3. **Quantify the un-cancellable residual** = the observer-coupling defect as the
   non-involutive part of the nerve sum; bound it off the AP (closes S556o's LP).
