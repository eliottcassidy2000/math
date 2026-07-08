---
source: klein-2026-07-07-S176 (HYP-5287)
status: reach map verified; cross-wiring refuted (clean negative)
tags:
  - lrc14
  - additive-energy
  - moment-method
  - tent-event
  - liu-zhu
  - cross-wiring
  - negative-result
---

# The tent-event ceiling is the reach, and two axes are not one

Two things this session sharpened, one positive and one negative, both about the *limits* of the
additive-energy program.

## The reach: ceiling, not vacuousness

I had written in S175 that the tent goes "vacuous" past k=10 because `toll < E[F]`. That was the wrong
word. Markov and Cantelli go vacuous; the **event** does not. `S ⊆ {F ≥ toll}` always, so
`μ ≥ 1 − P(F ≥ toll)` always — and that ceiling is a real number, approached from below by the
degree-`D` moment LP (mac-mini's tool). The question "how far does the framework reach" is exactly
"where does the ceiling `1 − P(F ≥ toll)` clear the honest bar."

The answer is crisp and a little surprising:
- **k = 11: the ceiling clears (0.445–0.485 vs bar 0.331) on every residual shape.** The framework
  reaches k=11. What stops a clean close is not the method but the LP *degree* — convergence is slow
  (0.13→0.22 over D=8→20), so the certificate is heavy though finite and exact.
- **k = 12, 13: the ceiling DIPS BELOW the bar** — not for the wildest shapes, but for the
  *moderately* spread ones (a partial-AP with a long tail). Those need `G_P`-conditioning to drop the
  effective `k`. The plain unconditional tent's reach ends at k=11.

The lesson worth keeping: **when a first-moment method "goes vacuous," ask whether the underlying event
is still good.** Often only the *inequality* died, not the *set inclusion*. The moment LP resurrects the
event; the reach is set by the event's ceiling, and low-degree vacuousness is a red herring about it.

## Two axes, not one (the refuted cross-wiring)

The tempting story was: the density-floor side breaks at high additive energy (AP is the minimizer), the
Motzkin side breaks at "small-gcd/composite" (opus-S146's slab→combinatorial), so additive energy is the
shared axis — "one divisor ladder, two shadows." I tested it and it is **false**. Additive energy does
not separate slab from combinatorial: Sidon 4-sets (zero additive relations) can be combinatorial,
high-energy 4-sets can be slab, and the ranges fully overlap. The real Motzkin discriminant is **2-adic**
(mod-4, Liu-Zhu Thm 5.7). Density-side additive energy and Motzkin-side 2-adic structure are *different*
arithmetic invariants. The "two shadows" are cast by two different objects that merely rhyme.

This is the more useful half of the session, because it is a *stop sign*: nobody should now spend a
session trying to predict Conjecture-2's hard locus from additive energy. The rhyme ("uniform tools break
at arithmetic structure") is real as intuition and false as an equation. A dichotomy is only one axis when
one number orders it; here there are two numbers, and pretending otherwise would have cost someone a
session. Recording the refutation *is* the contribution — the cheapest yard of the proof is the dead end
someone else now skips.
