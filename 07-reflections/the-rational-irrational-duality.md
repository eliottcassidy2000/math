# The rational/irrational duality — why the base floor is a density argument

*kind-pasteur-2026-07-03-S31. The user asked me to work the density bridge and step-5
against the old dualities: irrational/rational, odd/even, positive/negative,
addition/multiplication. The density bridge turned out to BE the rational/irrational
duality made mechanical, and formalizing it (`base_floor_of_cite`, then
`base_floor_of_cite_gen`/`base_floor_list_of_cite`) is the clearest statement of what
that duality does in the LRC(14) proof.*

## The two-sided architecture, restated as a duality

`lrc14_of_magnitude_split` (S29) cuts every hard family at a magnitude `M`:

- **Below `M`** — a bounded-denominator census. The lonely time is a fraction `p/q`
  with `q` small (spread13, the sieve, the denominator route). This is the RATIONAL
  side: loneliness is a finite check over residues mod `lcm(1,…,Q)`.
- **Above `M`** — the analytic / renormalization sweep. The lonely time is chosen by a
  free construction `t* = (⌈N(t₀ − δ/V) − 1/14⌉ + 1/14)/N`, generically IRRATIONAL.
  This is the IRRATIONAL side: loneliness comes from a measure floor on an interval, not
  from a residue.

The two sides look like different mathematics — number theory versus measure theory,
finite check versus continuous sweep. The density bridge is the hinge that shows they are
the same statement seen from two directions.

## The bridge: an open set has a rational point

`base_floor_of_cite` proves `0 < length (goodRegion2 base (1/14))` from the LRC(≤13)
citation. The mechanism is exactly the duality:

1. **Irrational side (what the citation gives).** LRC for 13 runners hands us a REAL
   `t₀`, possibly irrational, that is lonely for the 12 base speeds at the generous
   margin `1/13`. This is a point in the continuum, with no denominator.

2. **The open-set fact.** The good set — the `t` at which every base runner keeps
   distance `> 1/14` from the integers — is OPEN. It is a finite intersection of open
   arcs. So around the real witness `t₀` there is a whole interval of good points.

3. **Rational side (what we extract).** By density of ℚ in ℝ (`exists_rat_btwn`), that
   interval contains a rational `x`. The margin slack `1/13 − 1/14 = 1/182` is exactly
   the room needed: a rational within `1/(182·ΣB)` of `t₀` is still strictly `1/14`-far
   from every base multiple. Its fractional part is a rational strict-good point in
   `[0,1)`, which is a *member* of `goodRegion2`, which forces positive length.

So: **a real (irrational) lonely point forces a rational one, because the good set is
open.** The irrational side of the architecture is never actually needed to *name* the
time — its only role is to certify that the rational side is non-empty. The floor is a
density argument, not a construction. That is the whole content of "bounded magnitude =
rational census, large magnitude = irrational sweep" (opus-S50): the sweep proves the
census non-vacuous, the census is what you verify.

## Why this is the right shape for the endgame

The peels consume this floor and nothing else. `far_peel_lonely` (one far runner) and
`goodRegion2_simul_peel` (opus, ≤6 far runners) both reduce loneliness of the full
13-family to (i) `0 < length (goodRegion2 B (1/14))` for the bounded base `B`, and (ii) a
size threshold on the far speeds. The density bridge closes (i) for every base of ≤ 12
positive speeds — which is exactly the arity the peels leave behind (`13 − |far|`). Once
the floor is a citation corollary, the peel is pure measure arithmetic: the base has
room, and a far-enough runner's dense comb cannot eat all of it.

This is why generalizing the bridge from the fixed 12-arity (`base_floor_of_cite`) to all
arities `≤ 12` (`base_floor_of_cite_gen`, `base_floor_list_of_cite`) mattered this
session: the single-far peel leaves 12, but the ≤6-far simultaneous peel leaves anywhere
from 7 to 12, and each of those bases needs its floor discharged the same way. The margin
`1/(j+1) ≥ 1/13` holds for every `j ≤ 12`, so the same `1/182` slack works at every
arity. One density argument, every base.

## The other dualities, briefly

- **Positive/negative.** Already mechanical: `lonely_abs_iff` (LRC14CertRoute) folds
  every sign choice onto `|v|`, because `|(−v)·t − m| = |v·t − (−m)|` and `−m` ranges
  over all integers. This is why the whole endgame works with positive speeds. It is the
  cheapest of the dualities and the one the reductions lean on first.
- **Addition/multiplication.** The covering constraint is multiplicative (`q ∣ vᵢ` for
  each `q ∈ {2,…,14}`); the sweep is additive (subtract danger combs from `[0,1)`). A
  family that is NOT covering misses some `q` and is lonely at the purely multiplicative
  `t = 1/q` (the sieve). A covering family has closed the multiplicative escape and must
  be beaten additively, by the peel. The magnitude split is the seam between the two
  operations: below `M`, multiply (residues mod `q`); above `M`, add (comb subtraction).
- **Odd/even.** The danger radius is `1/14 = 1/(2·7)`. The `7` is the odd prime that
  organizes the hardest near-equal obstructions (LRCSevenGap, the residues-mod-7
  construction). The simultaneous peel's floor is positive exactly for `|far| < 1/(2h) =
  7` — the same `7`, now as the number of far runners the union bound survives. The
  even factor `2` is the two-sided `dangerPair` (a tooth on each side of every integer).

## Where step-5 actually sits

Working the bridge against step-5, the honest boundary came out sharper than "22 < w ≤
threshold." The far-peel threshold is `≈ V²·const` (piece count `~V` times the `1/V`
floor), so the far-peel closes only genuinely huge runners; the certified deep-well
family `[1,…,12,182]` has `182 ≪ V²` and is a step-5 family, closed by census, not peel.
The clean decomposition of the covering-far surface is by OUTLIER COUNT `r`:

- `r ≤ 6` far runners, all large → the simultaneous peel closes it, and as of this
  session its base floor is a citation corollary. This slice is done modulo the fee check.
- `r ≤ 6`, some far runner bounded (`22 < w ≤ fee threshold`) → bounded-magnitude census.
- `r ≥ 7` → the union floor `1 − 2h·r` goes non-positive; the simple peel cannot survive
  seven outliers at once. This is the residual tower — the renormalization / iterated
  peel (mac-mini's THM-608 scale separation), the genuinely open middle.

The rational/irrational duality closes the first slice cleanly and explains why the third
is hard: seven is both the wall's denominator and the number of outliers the additive
union bound can absorb. Past seven you cannot subtract all the combs at once, and you are
forced back onto the multiplicative structure — the same finite residue check that the
field is grinding down for `n = 14`.

---
*Linked: [[the-bounded-denominator-route]] (S28, the rational side's beginning, corrected
S29), [[the-independence-obstruction]] (S27, the measure side's end). Corpus:
`base_floor_of_cite`, `base_floor_of_cite_gen`, `base_floor_list_of_cite` in
`LRCFarPeelGood.lean` (kernel-pure). HYP-4044, MISTAKE-096.*
