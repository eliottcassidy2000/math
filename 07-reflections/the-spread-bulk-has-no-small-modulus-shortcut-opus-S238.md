---
source: opus-2026-07-11-S238
status: A clean pigeonhole reduction + an HONEST NEGATIVE result on the spread-bulk residual. Not-clearing at
  a prime modulus q∈{17,19,23} ⟺ the folded residues occupy all (q−1)/2 fold-classes. But the {17,19,23}
  shortcut FAILS: an adversarial spread divisor-complete family occupies all fold-classes mod all three at
  once, clearing only at q=26. So the anti-concentration is irreducible to a few primes — the full
  ~7-modulus window disjunction is genuinely needed.
tags:
  - lrc14
  - spread
  - anti-concentration
  - pigeonhole
  - negative-result
  - no-shortcut
---

# The spread bulk has no small-modulus shortcut

**opus-2026-07-11-S238.** Working the residual's hard core — the **spread** (longest-AP ≤ 7) divisor-complete
families (99% of the residual, S237). The open statement: every such family clears at a non-14 modulus
`q ∈ [15,31]`, giving `M > 1/14` (band-edge lemma). This is the S230/S231 anti-concentration. I probed whether
a **small fixed modulus set** suffices — which would make the proof target far more tractable.

## The pigeonhole reduction (clean)

For a **prime** `q ∈ {17,19,23}` the danger band is `{0, ±1}`, and `×p` permutes the `(q−1)/2` antipodal
fold-classes. So:

> **A family with no multiple of `q` does NOT clear at `q` ⟺ its 13 folded residues
> `{min(v_i mod q, q − v_i mod q)}` occupy *all* `(q−1)/2` fold-classes** (no empty class ⟹ every multiplier
> hits `{0,±1}`). Clearing ⟺ some fold-class is empty.

With 13 speeds and `(q−1)/2 ∈ {8, 9, 11}` classes, a random spread family usually leaves several classes
empty (occupancy ≈ 7–9 of 11 mod 23), so it clears. The natural hope: *some* prime in `{17,19,23}` always has
an empty fold-class.

## The honest negative result: the {17,19,23} shortcut fails

Adversarial hill-climb **finds a counterexample**:

`v = [42, 48, 60, 108, 125, 154, 195, 206, 210, 245, 252, 259, 294]`

is primitive, divisor-complete, longest-AP = 3 (spread), with no multiple of 17, 19, or 23 — and it occupies
**all** fold-classes mod 17 (8/8), mod 19 (9/9), **and** mod 23 (11/11) *simultaneously*. So it clears at
**none** of `{17,19,23}`. It clears only at `q = 26` in the full window (`M ≥ 2/26 = 1/13 > 1/14`).

A random sample gives **0** such families (occupying all fold-classes mod all three primes at once is rare),
so the earlier "0/425" reading was misleading — it is rare but **achievable**. The minimal covering set of
non-14 moduli for spread divisor-complete families is **~7 moduli** (e.g. `{15,16,19,23,25,29,31}`), with **no
2–3-element subset sufficing**.

## Consequence

The spread-bulk anti-concentration is **irreducible to a few primes**. Proving it requires the full window
**disjunction** — clearing at *some* `q ∈ {15,17,19,23,25,27,31}` — where each `q` is a genuine
anti-concentration condition, and no fixed small prime set is a shortcut. This **rules out** the tempting
"prove clearing at 2–3 fixed primes via fold-class occupancy" approach.

## Honest scope

This session is a **negative result plus a clean reduction**: the per-prime clearing criterion is the
fold-class occupancy (proved, elementary), but the *disjunction* over the window cannot be collapsed to a few
moduli — the adversarial family defeats any 3-prime subset. So the spread bulk remains the genuine
anti-concentration wall the fleet has long circled; what this adds is (a) the exact fold-class reduction for
the prime moduli, and (b) the firm knowledge that the window cannot be shrunk to a provable-by-hand subset.
The honest path forward is either the full-window disjunction (a real anti-concentration theorem) or the
finite census at the LEM-010 diameter bound — not a fixed-few-moduli shortcut.

→ opus-S237 (the spread bulk = 99% of the residual), opus-S235 (band-edge lemma — the margin), opus-S230/S231
(the anti-concentration — the wall), opus-S232 (the fold-class / summand-shell structure at prime moduli),
kps cont.36 (the spread hard core). Files: `lrc14_spread_bulk_no_shortcut_opus_S238.py` (+`.out`).
