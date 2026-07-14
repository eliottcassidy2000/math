# The functor does not exist: the runner is a level without a curve, and its curve is Eisenstein, not 14a

*mac-mini-2026-07-13-S94. The capstone of the S89→S94 arc. Owner asked, in sequence: how tournaments
relate to the last bit (S89); to build the odd-graph→cusp transfer (S90, provably blind); to explore
the metagraph V₄ (S91, complement irreducible under Sₙ); to find the arc-transitivity breaker (S92,
the circle Z/14 where complement factors); to prove the moduli fixed-point arithmetic (S93, proved,
"curving creates the class numbers"); and finally to **prove the functor runner/circle → X₀(14)**.
The honest answer to the last one is: **it does not exist**, and the obstruction is clean enough to be
the right place to end. This reflection records why, and what the whole arc actually established.*

---

## The obstruction, in one line

A point of `X₀(14)` is a **pair** `(E, C)` — an elliptic curve `E` together with a cyclic subgroup `C`
of order 14. The runner side canonically supplies `C` (the danger circle `Z/14 = Z/2 × Z/7`, S92) — and
**no `E`**. A moduli of curves-with-level has no image from an object that has the level but not the
curve. There is nothing to be functorial *into*.

## Sharper: the runner's own curve is the wrong one

One might hope the runner side implicitly determines a curve. It does — but the wrong one. The
covering-min value is `n/Φ₆(n)`, and `Φ₆(n) = n²−n+1 = N_{ℚ(ζ₆)/ℚ}(n − ζ₆)` is the **Eisenstein norm**
(verified: `Φ₆(6,7,14) = 31, 43, 183`). So the runner's intrinsic arithmetic lives on the **Eisenstein
integers `ℤ[ζ₆]`** — the field `ℚ(√−3)`, the `j = 0` CM curve. But `X₀(14)`'s curve is `14a`, which has
`j ≠ 0` (not CM at all), and whose own CM fixed points (S93) sit at `disc −56 → ℚ(√−14)` and
`disc −7, −28 → ℚ(√−7)`. The runner's natural moduli point is in `ℚ(√−3)`; `X₀(14)`'s arithmetic is in
`ℚ(√−7)` and `ℚ(√−14)`. **Different curves, different fields.** The functor would have to realize `14a`
from runner data that instead realizes the `j = 0` Eisenstein curve. It cannot.

## The "coincidence at 14," made precise

The whole thread has carried a warning — "coincidence at n=14," "the two twelves," "do not conflate
different 14's." Here it resolves into an exact statement. **Two distinct arithmetics meet at the
integer 14:**

- The **level** `14 = 2·7` — the danger threshold `1/14 = 1/(2·7)`, the 2-adic descent times the
  apex-7 — points to `X₀(14) = 14a`: the level structure, the cusp form `f₁₄`, the Atkin-Lehner V₄.
- The **value-denominator** `Φ₆(14) = 183 = n²−n+1` — the covering-min `n/Φ₆(n)` — points to the
  Eisenstein integers `ℤ[ζ₆]`: `ℚ(√−3)`, `j = 0`.

These share only the number 14. The LRC(14) covering-min is an **Eisenstein / `Φ₆` object wearing a
level-14 (`2·7`) coat**; `X₀(14)` is the coat, not the body. No map turns one into the other — the
severance is that they are different objects that happen to be indexed by the same integer.

## What the arc did establish (the real ledger)

The functor is the one plank that fails, but the arc built a real bridge up to it:

1. **The last bit is an AP-rigidity** (S89): covering ⟹ `conc < 7` uniformly ⟺ the covering-min, with
   the AP the unique tight configuration.
2. **The tournament cannot see it** (S90): the difference tournament is rotation-invariant, so the AP
   and its covering perturbations are identical tournaments with different loneliness — order forgets
   measure.
3. **The root is Sₙ-transitivity** (S91): the complement is the only tournament involution and is
   *irreducible* — the arcs form one Sₙ-orbit, so the `2·7` factorization has no combinatorial image.
4. **The cure is the level** (S92): the danger circle `Z/14 = Z/2 × Z/7` breaks arc-transitivity by
   CRT, and there complement **factors** into the Atkin-Lehner V₄ — `W₂·W₇ = W₁₄`.
5. **The moduli arithmetic is proved** (S93): `X₀(14)`'s fixed points are `(0, 4, 4)`, the Fricke count
   is `h(−56) = 4`, and "lifting" the circle to the curve is the **curving that creates the class
   numbers** — the 2-part lifts as a fixed-point-free translation, the reflections gain their CM points.
6. **The functor does not exist** (S94): the runner is a level without a curve, and its curve is
   Eisenstein, not `14a`.

The bridge is real from the combinatorics to the *level structure* and its *fixed-point arithmetic*;
it stops at the *curve*. That is not a small gap papered over — it is the exact and honest boundary:
`X₀(14)` is where the LRC(14) level and cusp form live, and the runner combinatorics reach the level
but not the curve, because the covering-min's body is Eisenstein.

## Coda

I was asked to prove a functor, and the right result is that there isn't one — for a reason precise
enough to be worth more than a forced construction would have been. The Atkin-Lehner V₄, the class
numbers, the cusp form: these are genuine features of `X₀(14) = 14a`, and the runner side genuinely
carries the level `Z/14` and its `2·7` V₄. But the elliptic curve `14a` is not a shadow of the runner
problem — the runner problem's own arithmetic shadow is the Eisenstein `j = 0`. The last inch of
LRC(14) is a metric fact (the covering-min, the Dedekind/`f₁₄` residual); the modular curve organizes
its *symmetries* (the level, the fixed points), not its *values*. The functor's non-existence is the
statement that symmetry and value, here, are carried by two different curves — and that is the truest
thing the arc found.

---

*Cross-links: the arc — S89 (HYP-6545), S90 (HYP-6555), S91 (HYP-6565), S92 (HYP-6575), S93 (HYP-6580),
S94 (HYP-6585, this). The two curves — `Φ₆` = Eisenstein/`ℚ(√−3)`/`j=0` (HYP-3768), `14a` = conductor
14/`ℚ(√−7)`,`ℚ(√−14)` CM (S93). The caution vindicated — kps "two-twelves / coincidence-at-14."
Proof/obstruction: `04-computation/functor_runner_to_x014_macmini_S94.py`.*
