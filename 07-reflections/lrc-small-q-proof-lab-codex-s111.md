# LRC Small-q Proof Lab

This session applied the LRC14 proof machine to smaller even denominators
`q=8,10,12` and then back to `q=14`.

The main outcome is a route discriminator.  The Bonferroni witness floor,
bounded-bank `nu` minimization, bounded p0 cap, and AP-facing
difference-profile all survive exactly.  Scalar additive energy does not.

## Exact Results

The audit is in:

- `04-computation/lrc_small_q_proof_lab_codex_s111.py`
- `05-knowledge/results/lrc_small_q_proof_lab_codex_s111.out`

For each even `q`, it uses threshold `1/q` and cluster gap `2/q`.  It scans
primitive anchored shapes `E` with `max(E)<=q` in the hard cluster range
`k=q/2+1,...,q-1`.

The smallest Bonferroni floor in the lab is

```text
(q,k)=(12,7):  nu(AP_k)+cap_q,k-1 = 11/36 ~= 0.305556.
```

The smallest bounded p0-cap margin is

```text
(q,k)=(10,6): cap_q,k - max p0_q(E) = 1/25.
```

For q=14, the tight bounded margin is the known k=8 row:

```text
cap_8 - p0(AP_8) = 319/5880 ~= 0.054252.
```

Consecutive/AP is always the bounded-bank `nu` minimizer and dense-set `D`
maximizer.  AP difference-profile majorization has zero failures.  This is the
same aggregate-low-denominator signal as HYP-2866.

But p0 is not literally AP-maximal in every bounded bank: there are eight
`p0_AP_not_unique_or_not_max` cases.  They remain cap-safe.  This matters
because it separates the true theorem from the tempting but false theorem:

```text
false target: AP maximizes p0 everywhere.
true target: every p0 leader stays below cap, with AP/dense prefixes as the
             controlling boundary model.
```

Scalar additive energy is worse.  The lab reports

```text
scalar p0 inversions = 12706
scalar D inversions  = 12139.
```

So a scalar energy monotonicity proof is already dead in smaller q.  The q=14
worst p0 inversion at k=9 is the familiar HYP-2890 bridge row

```text
(0,2,4,6,7,8,10,12,14).
```

That is a labelled even-AP-plus-midpoint residual, not a random exception.

## Synthesis

The small cases do not reveal a shortcut that bypasses Part A or the p0/floor
work.  They do something more useful: they show which objects are stable as q
varies.

Stable:

- Bonferroni floor.
- Consecutive/AP `nu`/dense-set extremality.
- Bounded p0-cap margin.
- AP-facing difference-profile/Fejer majorization.
- Labelled residual leak after the positive same-frequency packet.

Unstable:

- Scalar additive-energy monotonicity.
- Local compression monotonicity.
- Pointwise p0 AP maximality.

The likely proof skeleton remains:

```text
GOOD floor via Bonferroni
  + cap-safe p0 bound
  + finite Part-A/arc-count conversion
  + labelled residual-leak inequality
  + wide tail/decorrelation.
```

The smaller q lab gives confidence that this is not overfit to 14.  The same
proof carriers already govern the easier denominators; only the margins change.

Concurrent mainline work makes the division of labor clearer.  HYP-2895/S108
uses smaller `N` to classify AP/Goddyn-Wong/apex-denominator exact-tiler
boundary phenomena and compact covering margins.  HYP-2896/S109 closes the
one-tail zeta `-1/12` branch.  This lab is the analytic companion: after those
boundary and one-tail branches are removed, scalar additive energy still fails,
while the labelled Fejer/residual route remains stable.

## Assumption Challenge

Candidate vertex sets considered: runners, speed gaps, q-sector labels,
small-part sets `P`, cover arcs, dense-set components, difference profiles,
additive-energy levels, Fourier modes, and proof obligations.

Chosen vertices for the tournament analysis are proof carriers:

```text
Bonferroni_floor,
bounded_nu_minimizer,
bounded_p0_cap,
AP_difference_profile,
scalar_additive_energy.
```

This preserves the LRC predicates that actually feed the proof:
positive witness measure, p0-cap safety, and labelled Fejer residuals.  It
destroys runner identity and scalar energy rank.  That destruction is deliberate:
the exact data show scalar energy rank does not preserve the needed inequalities.
