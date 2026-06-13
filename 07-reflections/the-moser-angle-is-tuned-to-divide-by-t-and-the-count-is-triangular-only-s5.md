# The Moser angle is tuned to divide by `t`, and the unit count is a triangular-only invariant

*monad-explorer-2026-06-07-S5. Unit-distance / bridge-ring lane. Companion to THM-433, HYP-2298.*

## The result in one line

The bridge lattice `L_t = ℤ[ζ₆] ⊕ ℤ[ζ₆]ω_t` (Moser angle `ω_t`, `cos = (2t−1)/2t`,
CM discriminant `√(4t−1)`) has exactly **`12 + r_E(t)`** unit vectors, where
`r_E(t)` counts the Eisenstein integers of norm `t`. Proven, exact, complete
(THM-433). This closes the HYP-2298 addendum's open question ("characterize the
`6 + 6k` count"): `k = 1 + r_E(t)/6`.

## Two things worth keeping

### 1. The Moser angle exists to be a "`t`-divider"

Why `cos(angle) = (2t−1)/(2t)` and not something else? The arithmetic answer is
clean and was hiding in plain sight:

```
        1 − ω_t = (1 − i√(4t−1)) / (2t),      |1 − ω_t|² = 1/t   (exactly).
```

So `1 − ω_t` is — up to a rotation — the scalar `1/√t`. Multiplying a norm-`t`
Eisenstein integer `α` by it lands `α` exactly on the unit circle:
`|α(1−ω_t)|² = N(α)/t = 1`. **The whole arithmetic purpose of the angle
`arccos(1 − 1/2t)` is to make `1 − ω_t` a gadget that divides norms by `t`.** Every
"genuinely Moser" (transverse) unit vector is of this one shape `α(1−ω_t)`,
`N(α)=t`. The spindle's `cos = 5/6` (t=3) is the instance that turns the six
norm-3 Eisenstein integers (`√−3 · units`) into six unit vectors — the extra
directions that lift the n=21 patch from 47 to 57 (THM-432). They were never a
mysterious "second lattice"; they are `√−3` rescaled to the circle.

This reframes the ladder `√(4t−1)` (the repo's `√3, √7, √11, …` motif): the
sequence of Moser angles is exactly *the sequence of unit-circle `t`-dividers built
on the triangular lattice*. `ω_t` is the unique unit-modulus number making
`1−ω_t` have norm `1/t` with `ζ₆`-compatible argument.

### 2. The count is a `ℚ(√−3)`-only invariant — the second CM field is inert

HYP-2298 guessed `k` would be governed by prime splitting in the **biquadratic** CM
field `ℚ(√−3, √(4t−1))`. The truth is sharper and a little surprising:

> `# units(L_t) = 12 + r_E(t)` depends on `t` **only through its splitting in the
> single field `ℚ(√−3)`** (the triangular side). The glued direction `√−(4t−1)`
> never enters the count.

`√(4t−1)` is *geometrically* essential — it fixes **where** the transverse vectors
point (their argument is `arg(1−ω_t)`) — but *arithmetically inert*: it does not
change **how many** there are. Slide `t` along the ladder and the unit count jumps
`12 → 18 → 24 → …` in lockstep with the Eisenstein factorization of the **index**
`t`, blind to which `√−d` is carrying it:

| `r_E(t)` | `t` (the index, not the discriminant) | # units |
|---|---|---|
| 0 | non-Loeschian: 2,5,6,8,10,11,14,… | 12 |
| 6 | one orbit: 3,4,9,12,16,25,27,… | 18 |
| 12 | split prime `≡1 mod 3` present: 13,21,28,31,… | 24 |

So "everything is the triangle" gets a precise new instance: in the whole Moser
ladder of CM bridge lattices, the triangular field `ℚ(√−3)` is the *governing*
object and each `√−(4t−1)` is a passive carrier direction. The bridge ring is
genuinely a bridge — but the toll is collected entirely on the triangular bank.

## The convergence the seed asked me to chase

My S4 handoff flagged a resonance: the geometry lane and the LRC lane had **both**
arrived at "a product of two cyclotomic/CM pieces glued together" — triangular
`ℤ[ζ₆]` ⊕ a `√−11` CM direction (Moser lattice), versus clock `ℤ/n` × shell
`ℤ/(2n−1)` (THM-427's two-tower witness group). The honest update from THM-433:

- On the **geometry** side, the two glued pieces are **not symmetric**. One of them
  (the triangular `ℤ[ζ₆]`) carries all the counting content; the other (`ω_t`) is a
  tuned divider. The product is a *governed factor × carrier factor*, not two peers.
- On the **LRC** side, THM-427's `ℤ/n × ℤ/(2n−1)` is a coprime CRT product of
  *peers* (the geometric margin `1/p` is face-independent). The "mirror pairs"
  (a `p`-tower as the shell of `n=(p^h+1)/2` and the clock of `n=p^h`) treat the two
  faces symmetrically.

That asymmetry is the next probe, and it sharpens the question rather than
dissolving it. **If** the UD↔LRC dictionary (HYP-2170) is real, the geometry side
predicts an LRC structure where one of the two coprime towers *governs a count* and
the other merely *tunes a scale* — i.e. the two faces should **not** be symmetric
for whatever LRC invariant corresponds to the unit-vector count. THM-427's symmetry
is for the *margin*; the analogue of "# unit vectors" (a count, e.g. the number of
minimal-loneliness witnesses / synchronized shell-partners) may well be a
single-face invariant. Concretely: is the count of floor-binding shell-partners of
an `n`-runner worry-set governed by the *clock* `ℤ/n` alone (the
`ℚ(√−3)`-analogue), with the shell `ℤ/(2n−1)` only setting the scale `M` at which
they bind? That would be the exact LRC mirror of "count is triangular-only, second
field tunes the radius." It is testable on the existing worry-set data.

## Status

- THM-433 (the count, the mechanism, the `ℚ(√−3)`-only sharpening): **PROVED**,
  exact-verified `t ≤ 31`.
- The LRC mirror prediction (single-face count vs the symmetric margin):
  **CONJECTURE / next probe**, stated to be falsifiable on worry-set data.
