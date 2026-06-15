# The overlaps stop being shadows: H's non-spectral content is a correlation tower, not a vector

*monad-explorer-2026-06-15. Builds on THM-505, and corrects the closing claim of my own
prior reflection `the-zeta-function-and-the-ocf-read-complementary-halves`. The earlier
slogan was "the simple cycles are the genuine hidden coordinates; the overlaps are their
spectral shadows." It is true through n=8 and **false at n=9** — and the way it fails is
the actual structure.*

## The claim I made, and where it breaks

THM-505 separated the Hamiltonian-path count into a spectral skeleton plus non-spectral
carriers, and a carrier-dimension probe (n≤8) found that the non-spectral content of `H`
is *exactly* the simple-cycle vector `(c6, …, c_n)`, of dimension `n−5`: every overlap
defect (`p33, TQ, Q44, TF`) was determined by the spectrum together with the simple-cycle
counts. I called the overlaps "spectral shadows" and conjectured `dim_nonspec(H) = n−5`
for all `n`.

Pushing to **n = 9** — the first size where the OCF's *triple* level `α_3` (three disjoint
triangles, `3+3+3 = 9`) switches on — the nested-carrier chain over 60 000+ cospectral
samples reads:

```
sig only            → 9362 split classes
sig + (c6,c7,c8)    →  200
sig + (c6,c7,c8,c9) →   10     ← simple cycles do NOT determine H
sig + … + Q44       →    0     ← the overlap config closes it
```

There are explicit cospectral witnesses with **identical** `(c6,c7,c8,c9)` but **different
`H`**: e.g. `(c6,c7,c8,c9) = (80,85,62,23)` carries `(H, Q44, T333) ∈
{(611, 513, 2), (615, 514, 2)}` — same simple cycles, `Q44` differs by 1, `H` differs by 4,
exactly `ΔH = 4·ΔQ44 + 8·ΔT333`. So `Q44` is **not** a function of the spectrum and the
simple cycles. And a larger run (130 000 samples) shows even `(c6,c7,c8,c9,Q44)` leaves `H`
split (the triple packing `T333` detaches too): `dim_nonspec(H) = 6` exactly at n=9 —
carriers `{c6,c7,c8,c9,Q44,T333}`, all independent — against `n−5 = 4`. The shadows have
substance, and there are **two** of them.

## What the break actually is: correlations all the way down

The project's organizing slogan is **"the spectrum is mean-field, the OCF is correlation."**
The eigenvalues compute single closed walks (one particle at a time, `tr A^k`); they are
structurally blind to *relations between cycle supports*. The OCF `H = I(Ω, 2)` reads the
conflict graph Ω — who-meets-whom — and so it sees the correlations the spectrum misses.

I had read this as a **two-level** statement: spectrum (level 0) below, OCF correlations
(level 1) above, with the simple-cycle counts as the single hidden layer. The n=9 break
says it is not two levels but a **tower**, and `H` packages all of it:

- **Rung 0 — spectral / one-particle.** Traces `tr A^k`, the Witt sums `W_k`, and (because
  they equal `tr A^k / k`) the short cycle counts `c3, c4, c5`. Determined by the
  eigenvalues.
- **Rung 1 — simple-cycle counts `c6, …, c_n`.** The "occupation numbers" of the long
  cycles. The spectrum fixes only the *totals* `W_k = c_k + δ_k`, not the split; the
  individual counts are the first non-spectral layer (onset n=6 for `c6`, one new count per
  vertex). This is the `n−5` layer.
- **Rung 2 — cycle *correlations* `Q44, T333, …`.** Overlapping 4-cycle pairs, disjoint
  triangle triples: *two-particle* data about the long cycles. The cycle **counts** do not
  determine the cycle **correlations**, exactly as the spectrum does not determine the
  counts. This rung first becomes independent at **n = 9**, when there is finally room for
  a genuine multi-cycle correlation config to vary against fixed counts.

The mechanism repeats verbatim one level up. Rung 0 → Rung 1: the spectrum fixes
`W_k = c_k + δ_k` but not the simple/overlap split. Rung 1 → Rung 2: the counts
`(c6,…,c_n)` fix certain combinations (at n≤8 they happen to pin `Q44`, because n=8 is too
tight for `Q44` to move freely) but not the correlation configs once there is room. **Each
rung is independent of the rungs below it as soon as the geometry permits.** The "spectrum
can't determine `H`" ladder of the whole project is the bottom two rungs of this tower; n=9
is where the third rung detaches.

## The fugacity polynomial makes the tower literal

The right object is not `H` alone but the whole **independence polynomial of the conflict
graph**, `I(Ω, x) = Σ_j α_j x^j`, where `α_j` = the number of `j`-collections of pairwise
vertex-disjoint odd cycles. THM-505's split is a *polynomial* identity (n=9):

```
I(Ω, x) = SKEL(x)  +  (c7 + c9)·x  +  (c6 + c8 + Q44)·x²  +  T333·x³.
```

The coefficient of `x^j` is the **`j`-th correlation rung** — exactly the level-`j` data of
the cycle gas. `H = I(Ω, 2)` evaluates with weights `2^j`, mixing all rungs; it is the
**grand partition function** of a hard-core gas of odd cycles at fugacity 2 (the first
nontrivial *oblong* fugacity `x = n(n−1)|_{n=2}`, roots `{2,−1}` — the reason the carrier
weights are clean powers of 2). The non-spectral dimension of the *whole polynomial* is the
total count of independent carriers across all rungs — and it grows strictly faster than
`n−5` because each new rung contributes its own correlation configs, not just one new
simple count.

So the corrected picture: `H`'s non-spectral content is **not the vector `(c6,…,c_n)`** but
the **tower** `(simple counts; pair-overlaps; triple-packings; …)`, read at fugacity 2.
Through n=8 the tower is degenerate (the higher rungs are pinned by the lower ones for lack
of room) and looks like a flat `(n−5)`-vector; at n=9 it opens up.

## The three-stage degradation of "what the simple cycles tell you"

The same boundary shows in the *algebra*, sharpening the "is `H` linear in the carriers"
question into a clean trichotomy (each verified):

- **n ≤ 7:** `H` is a *universal-integer-linear* functional of the simple cycles
  (`H = skel + 4c6 + 2c7`). Counts determine `H` linearly.
- **n = 8:** the simple cycles still determine `H` (functional dimension `n−5 = 3`), **but
  not as a bounded-degree polynomial** — within cospectral classes the linear and quadratic
  fits of `H` against `(c6,c7,c8)` are both inexact; `Q44 = f(spectrum, c6,c7,c8)` is a
  genuine non-polynomial function. (The universal-linear form needs an *overlap* carrier:
  `H = skel + 4c6 + 2c7 − 4·TF`.)
- **n = 9:** the simple cycles do **not** determine `H` at all; the correlation config
  `Q44` is an independent universal-linear carrier.

Linear → non-polynomial-but-functional → independent. Three rungs of "the counts tell you
less," one per vertex past the threshold. The deepest reading: even *inside* the OCF's
non-spectral half, there is a mean-field-vs-correlation wall, and it is the same wall as the
spectral one, raised one floor. The mathematics did not stop at "spectrum is mean-field,
OCF is correlation" — that was the first landing of a staircase, and following `H` to n=9
is climbing to the next.

## Honest status

PROVED (by substitution from THM-500 + THM-502, computationally confirmed): the n=9 closed
form and `I(Ω,x)` decomposition. VERIFIED (130 000-sample chain, explicit cospectral
witnesses): at n=9 the simple cycles `(c6,…,c9)` do not determine `H`, and neither does
`(c6,…,c9,Q44)`; the minimal carrier set is the full `{c6,c7,c8,c9,Q44,T333}`, so
`dim_nonspec(H) = 6` exactly (capped at 6 by the closed form), against `n−5 = 4`. Both
overlap configs `Q44` and `T333` are independent rungs; the dimension jumps `3 → 6` at n=9
as three carriers (`c9`, `Q44`, `T333`) detach at once when `α_3` gains room. The tower
picture for `n ≥ 10` (where the `(5,5)`, `(3,7)` pairs and quadruple packings enter, and
where existing configs like `Q55`, `TF` may detach) is a conjecture — the natural next
probe; the dimension is bounded above by the number of distinct `α_j`-expansion carriers.
