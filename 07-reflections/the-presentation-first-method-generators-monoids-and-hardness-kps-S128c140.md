# The presentation-first method: see every problem as (object, monoid, action), and read hardness off the monoid

*kind-pasteur-2026-07-21-S128c140. Owner: see more problems as generators and monoids; get as
fundamental a view of as many topics as possible. Companion to THM-1885 (the ½-and-+1 monoid is
BS(1,2), with the full topic→monoid catalog) and THM-1880 (the a/b frame).*

The method is one question, asked of any object in the repo:

> **What are the generators, what relations do they satisfy, and what set does the resulting
> monoid act on?**

Answer it and the rest follows mechanically: **invariants** are the orbit functions, the
**nullcone** is the degenerate orbit (the fixed point / identity-like element), a **census** is an
orbit count (Burnside), and a **reduction** (`n → n−1`, `n → n−2`) is a generator stepping the
action down. THM-1885 runs this for eight topics and finds they collapse onto a short list of
monoids. This note records *why the method is worth making a habit*, and the one prediction it
makes.

## The habit

For a new object, before computing anything, write its **presentation**:

- **Tournament space** — generators: single arc-flips; relations: they commute and square to 1;
  monoid `(ℤ/2)^{C(n,2)}`, which splits as **cut ⊕ cycle** (base-path flips vs wiggly flips), the
  GF(2) decomposition that is the whole "score hierarchy vs cyclicity" story in one line.
- **The ½ and the +1** — generators `a=x+1`, `b=x/2`; relation `ab=ba²`; monoid **BS(1,2)** (THM-1885).
  The half-dictionary, the Legendre `½`, the fiber fraction, the `Re=−½` line are *elements* of it.
- **Farey / LRC** — generators `S:x↦−1/x`, `T:x↦x+1`; relations `S²=(ST)³=1`; monoid **PSL(2,ℤ)**.
  The continued-fraction address of a speed-set is its `S,T`-word.
- **The Weyl / odd-even axis** — generator `x↦−x`; monoid **ℤ/2**; the even/odd split of every
  characteristic form (`char_S` even, its `O_n` companion odd).

Writing the presentation *first* prevents the recurring mistake in `MISTAKES.md` — "the pattern
breaks at `n = 6/7`" — because a pattern that is not an orbit function of the acting monoid was
never going to be stable; you can see it is not equivariant before you compute it and watch it fail.

## The one prediction: amenability ↔ hardness

The monoids sort by a group-theoretic dichotomy that **matches the repo's easy/hard split**:

| monoid | class | the repo side it governs | difficulty |
|---|---|---|---|
| `ℤ/2`, `S_n`, `(ℤ/2)^k` | finite / abelian | spectra, switching classes, censuses | EASY (spectral, closed-form) |
| `BS(1,2)` | solvable (amenable) | the ½/dyadic threads, the ladder rungs | TRACTABLE (recurrence, holonomic) |
| `PSL(2,ℤ) = ℤ/2 ∗ ℤ/3` | virtually free (**non-amenable**) | Farey / LRC | HARD (LRC(14) open) |
| none (the permanent) | not an orbit function of any | `H`, Hamiltonian counts | `#P` (THM-1780) |

Reading down the column: the *amenable* monoids (finite, abelian, solvable) govern exactly the
topics the repo has closed-form or recurrence control over; the *non-amenable* free product
`PSL(2,ℤ)` governs LRC, the repo's hardest open problem; and the invariant with **no** acting
monoid at all is `H`, the `#P` permanent that leaves the ladder at `n=6`. So the fundamental view
does not just organise — it **predicts difficulty**: an invariant is as hard as the smallest monoid
whose orbit function it is, and if it is nobody's orbit function it is `#P`.

That is a genuine tool. It says: to attack a hard repo problem, find its acting monoid; if it is
amenable, a recurrence/detection-depth argument exists (as for TNC/GMC, BS(1,2)-solvable); if it is
`PSL(2,ℤ)`-governed (LRC), expect the difficulty of a hyperbolic/continued-fraction problem, and
attack it with the modular machinery, not with elimination; if it is a permanent, it is `#P` and
only bounds/parity (Rédei) are available.

## The habit, stated once

*Presentation first, computation second.* Every new object gets its generators, relations, and
acting set written down before its invariants are computed — and the monoid tells you, in advance,
whether the invariant will be spectral (easy), holonomic (tractable), modular (hard), or `#P`
(only-bounds). The repo has been discovering this monoid by monoid; naming the method makes it a
default.

## Cross-links
THM-1885 (topic→monoid catalog, BS(1,2)) · THM-1880 (a/b) · THM-1555 (half-dictionary) · THM-1810
(SL₂ / GIT) · THM-1775/1780 (ladder / `H` `#P`) · THM-826 (Farey / PSL(2,ℤ)) · CLAUDE.md (cut⊕cycle,
Mode A/B).
