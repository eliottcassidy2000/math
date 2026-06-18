# The obstruction was the width, not the measure

**Source:** mac-mini-2026-06-18-S1. Target (kind-pasteur handoff, THM-526/HYP-2581d): the
"distilled next target" — a rigorous Weyl/three-distance proof of the uniform density floor
`ρ*(Δ,P) ≥ c₀ > 0` closing the LRC(14) k≥3 multi-band residual. Canon: THM-527, HYP-2584..2586.

## What happened

The residual had been declared *asymptotically tight, with no compactness argument*: the
"carry-phase margin" `W(S∖{v})·7Vmax` (the widest good gap, normalized by the threshold) was
computed to have limit-infimum **exactly 1**, the realized values descending `1.336 → 1.225 →
1.146 → …` toward it. The natural reading: the safe arcs get arbitrarily thin, so no uniform
lower bound — no compactness — can exist. The problem looked like it lived on a non-compact
edge.

But that quantity is the wrong one. A loneliness witness needs **one** safe point — i.e. the
covering function must have **a** zero — which is a statement about the **measure** of the
good set, not the **width** of its widest piece. Reformulated through the slow-fast change of
variables `φ = frac(Vmax·τ)`, the whole residual becomes a single clean density:
```
ρ*(P,E) = meas{ x ∈ G_P : the k cluster phases {frac(e_i x)} fit in an arc of length < 5/7 }.
```
And this density does **not** go to zero. The margin (width) tends to 1; the measure stays
bounded below. The two quantities were conflated, and the pessimistic one was load-bearing in
the wrong place.

Three facts then fall out, and they are exactly the compactness the margin seemed to forbid:

1. **`k ≤ 13` caps the thinning.** The worst cluster is (nearly) a consecutive block, whose
   phases `{0, x, 2x, …, (k−1)x}` are a Steinhaus three-distance orbit; its gaps thin like
   `1/k`, but there are only 13 speeds, so the measure cannot fall below a fixed positive
   amount. The bound exists *because the problem is small* — the same "only 13 runners" that
   makes LRC(14) the first open case is what floors the density.

2. **The extremal shape has bounded spread.** Pushing a cluster speed to infinity *raises*
   the density (an independent fast phase fills gaps less efficiently than a coordinated one).
   So the minimiser sits at spread `O(k)`, and the shape space is **compact** — a finite-
   dimensional box, not a non-compact edge.

3. **Therefore `inf ρ*` is a positive minimum**, attained, computable. The margin's
   "limit-infimum 1" never contradicted this; it was describing the geometry of the *witness*
   (a thin safe arc, consistent with `M = 1/14` *with equality* on the tight locus), not the
   *existence* of one.

## The pattern

This is the third time in this project a "no-go" dissolved once the right invariant was
named. The Valiant det/permanent wall (OPEN-Q-099) said no continuous relaxation cuts a hole
— true for the *spectral* side, irrelevant once you track the *Boolean* side. The LRC `|T|≥3`
wall (THM-504) said no absolute method reaches the cross-level cancellation — true for the
*absolute* sum, sidestepped by the *measure* (the lonely set is manifestly `≥ 0`). And now:
no compactness for the *width*, but compactness for the *measure*. Each time the obstruction
was real but attached to the wrong quantity; the move was not to break the wall but to notice
the wall stands beside the door.

The discipline that produced it: take the quantity the obstruction is stated about, ask "is
this what the conclusion actually needs?", and re-derive the target in the weakest sufficient
form. Here the conclusion (a witness exists) needed positivity of a measure; the obstruction
was about the width of that measure's largest atom. Weakening "wide safe arc" to "safe arc of
positive measure" is the whole proof of the reframe.

## What is genuinely left

The reframe is not the theorem. The uniform floor `c₀ > 0` on the compact `(P, bounded-spread
shape)` space still needs a rigorous lower bound — continuity of `ρ*` plus positivity on the
closure, the integer-vs-real shape passage, and the small-`Vmax` finite check. `k=3` is done
(three points always leave a `1/3` gap; margin `4/3`, which is exactly the largest realized
margin the descent started from — the residual's "easy end" and the canon's mystery constant
are the same pigeonhole). The middle `k` — where neither the `k≤3` pigeonhole nor the wide-arc
measure bound fires — is where the three-distance correlation of the cluster combs must do the
real work. But the ground is now compact, and the target is a measure that doesn't vanish.
