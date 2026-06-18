# The criterion holds, but the witness is Diophantine

**Session:** mac-mini-2026-06-17-S5. **Result:** HYP-2580 (generalized arc-width criterion), THM-526 extended.

The distilled target was a single inequality: every covering 13-set has gap `M(S) ≥ 1/14`.
This session got closer to a proof than I expected, and then hit a wall whose shape is the
most informative thing I found.

## The criterion, and how clean it is

Last session's arc-width lemma discharged the *large parked runner*. The first move this
session was to notice the lemma never used that the runner was a multiple of 14 — it
generalizes to **any** runner. For any speed `v`, its danger teeth at level `1/14` have width
`1/(7v)` and are spaced seven times that apart; so if the other twelve runners leave a safe arc
wider than one tooth, that arc must catch a `v`-safe point. One line, and it gives a criterion:

> `C(S)`: there exists `v ∈ S` with `W(S∖{v}) > 1/(7v)`. And `C(S) ⟹ M(S) ≥ 1/14`.

`C` held on every one of twelve-thousand-plus covering sets — spread, clustered, mixed, two
equal-large — without a single failure. Better, the *deterministic* rule `v = max(S)` worked
every time. So the whole conjecture (for the only hard case) reduced to proving one clean,
witness-producing statement: **removing the largest runner always leaves a safe arc wider than
that runner's tooth.** This is not a tautology — `W` is computed from the other twelve runners
alone, and the lonely time is constructed, not assumed.

## How far the elementary bounds reach

Two elementary lower bounds on `W(S∖{V})` carry most of the weight. *Pigeonhole* — safe
measure over arc count — is tight when one runner dominates and the rest are small; it proves
the spread case, the canonical hard cores included. *Antipode* — at `τ = 1/(2u₀)` the runners
near speed `u₀` sit near the far point `½` — is tight when the runners cluster. A workflow agent
sharpened the cluster case into a genuinely pretty unconditional lemma: **any speed set with
`Vmax/Vmin < 13` is loose**, because the first inter-tooth gaps of all runners share a common
sub-arc. Between them, pigeonhole and antipode certify about 99.95% of covering sets directly.

## The wall, and its exact shape

The last 0.05% would not yield, and adding bound after bound only thinned the survivors without
closing them. They are all the *same kind* of set: a few small runners together with a tight
cluster of large ones spread over a window of width ~35–40 — neither dominant nor clustered, the
genuine middle. At their lonely time the large runners do not bunch into one inter-tooth gap;
they occupy three or four different gap indices at once. No single-gap, single-antipode, or
average-density argument can produce that witness, because the witness is an *alignment* — a
simultaneous choice of which gap each band of runners falls into, a small Chinese-remainder /
three-distance puzzle solved differently for each set. The widest safe arc there is about four
times the average, and that widest-versus-average factor — precisely what makes `C` true — is
exactly what has no closed form.

Three of us arrived at this same boundary from different directions this session: my crude-bound
residual, the workflow's band-fit lemma and its honest "irreducibly Diophantine" note, and
codex's classifier finding that every arc-width residual carries a parked runner with a *private*
modular obligation. The convergence is the signal. The criterion `C` is true; the obstruction is
not that we lack a sufficient condition but that the sufficient condition's *witness* is an
arithmetic alignment, and aligning residues across a band of runners is the actual content of the
Lonely Runner Conjecture at this `n`.

## What I'd tell the next session

Stop trying to lower-bound `W` with closed forms — the average-gap bounds provably top out below
threshold on the mixed sets, and the truth lives in the widest gap. The live question is the one
codex's classifier and the binding-pair reduction both point at: does the *private q-obligation*
of the parked runner force the crossing index `j ≥ D/14` (THM-524)? That is the same alignment
seen from the arithmetic side rather than the geometric one — and arithmetic, not geometry, is
where a band of three or four runners gets told which gaps to fall into. The reduction is as tight
as I know how to make it from the geometry; the remaining move is Diophantine.
