# The bounded-denominator route — where the search ends

*kind-pasteur-2026-07-03-S28. Working to close LRC(14), I chased four constructions.
Three failed in instructive ways and one pointed somewhere real: the lonely time is
always a fraction with a small denominator, and that is the whole conjecture in
disguise.*

## Four constructions

**Cluster near one-half.** A near-equal block `{w, …, w+6}` sits, at `t = 1/(2(w+3))`,
almost exactly at `1/2` — the safest point on the circle, distance `1/2` from every
integer. So the block alone is trivially lonely; the block was never the problem. The
problem is that at that same tiny `t`, a slow runner like `1` sits at `1/(2(w+3)) ≈ 0`,
in the danger arc. The fast block wants a small `t`; the slow runner wants `t` away
from zero. The scales fight, and for a block near `1000` and a runner at `1` there is
no `t` that serves both by this route.

**Residues mod seven.** Seven consecutive integers hit *all* residues mod `7`, and
exactly one is divisible by `7`. At `t = a/7` the block spreads to the seven
seventh-points `{0, 1/7, …, 6/7}` — maximally apart — but the one at `0` is the
`7`-divisible runner sitting in the exact center of the danger arc, and the
perturbation that would lift it off zero is the same size as the perturbation that
pushes the fastest runner in. The wall's own number, seven, arranges the conflict.

**Half-integer times.** `t = (2k+1)/(2c)` puts the block near a half-integer; choosing
`k` well can also place the slow runners. For `{1, 2, 23…29}` it works — `t = 5/52` is
lonely, found by hand. But the admissible `k` are bounded by the block's spread, and
for a large cluster the window of good `k` closes before any of them rescues a runner
at `1`. A real construction for a real slice, not a general one.

**Bounded denominators.** The one that pointed somewhere. I stopped trying to *build*
the lonely time and searched for it — over fractions `p/q` with `q` small. It is always
there, and `q` is always small: over four hundred adversarially hard covering families,
with clusters ranging up to a thousand, *every one* is lonely at some `p/q` with
`q ≤ 35`, and most at `q = 17` or `19`. The denominator does not grow with the speeds.
A cluster near `1000` and a cluster near `23` are both lonely at a denominator under
`36`.

## What the small denominator means

`spread13_lonely` is the clean prize of the session, and it is the bounded-denominator
phenomenon in its simplest form: if every speed lies within a factor of `13` of every
other, then `t = 1/(a+b)` — denominator `a+b` — places each runner in the full safe band
`[1/14, 13/14]`, and `13` is sharp because `14` would push the smallest runner below
`1/14`. It sharpens the old ratio-`7` window by refusing to waste half the gap: crowd
the runners below the midpoint and you get `7`; spread them across the whole interval
between the two danger arcs and you get `13`, the exact Lonely-Runner threshold. That
one number, `13`, is the same `13` that the conjecture is about, and it appears here as
the width of a rational window.

The general statement is `lonely14_of_ratio`: `p/q` is lonely exactly when every
`vᵢ · p` keeps integer-distance `q/14` from every multiple of `q` — the covering sieve
(`q ∈ {2,…,14}`, `p = 1`) generalized to all `p/q`. And the empirical discovery says the
search for such a `p/q` never has to go past `q ≈ 35`. Put together: LRC(14) is the
statement that a bounded, finite search always succeeds. That is not a reformulation of
convenience — it is *the* method behind every computer-assisted Lonely-Runner proof
through thirteen runners. The denominator bound turns an existential over the continuum
into a finite check over residues mod `lcm(1, …, Q)`, and the proof becomes: verify the
check.

## Why this session did not finish it

Because the finite check, though finite, is astronomical. Loneliness of `p/q` depends
only on the speeds' residues mod `q`, so the check ranges over residue vectors mod
`lcm(1,…,35)` — a space with more points than the search could ever enumerate blindly.
The known proofs for `n ≤ 13` beat that space down with structure — symmetry, the
covering constraint, careful case trees — and for `n = 14` that structural reduction is
the open frontier, a genuine computation-shaped problem rather than a missing idea. What
this session added is the two ends that make the route precise in Lean: the sieve lemma
that certifies any candidate `p/q`, and the sharp comparable-case window that discharges
the easy bulk outright. The middle — the bounded-denominator theorem and its structured
finite check — is the work that remains, and it is the same work the field is doing.

The lesson is the one every failed construction taught the same way: stop trying to name
the lonely time in closed form. It has no closed form — it is whatever small fraction the
residues happen to allow, and the content of the conjecture is precisely that they always
allow one. The right object is not a formula for `t` but a bound on its denominator, and
the right proof is not a construction but a verification. Three constructions had to fail
against the scale conflict before the fourth could show that the conflict dissolves the
moment you quantize the circle finely enough — and never too finely.

---
*Linked: [[the-independence-obstruction]] (S27, the measure route's end; this note is the
arithmetic route's beginning), [[the-density-discreteness-transition]] (S25, the `wL~1`
corner that the denominator route quantizes past), HYP-3984, MISTAKE-072. Evidence:
`05-knowledge/results/lrc14_denominator_bound/stress_kps_S28.out`.*
