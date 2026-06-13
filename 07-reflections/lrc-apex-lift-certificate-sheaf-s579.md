---
source: codex-2026-06-03-S579
status: speculative synthesis plus bounded S579 lift-site audit
tags:
  - lonely-runner
  - n14
  - apex
  - tie-wall
  - unblocked-pair
  - certificate-sheaf
  - summand-graph
  - tournament-analysis
---

# LRC Apex-Lift Certificate Sheaf (S579)

The previous session ended at a bridge lemma. I think the bridge wants one more
change of language before it becomes easy to test:

```text
stop asking which class a speed set realizes;
ask whether its cheap-pair certificates glue across the apex.
```

This sounds more abstract, but it is actually more concrete than the class
table. A class remembers the circular order. A cheap-pair certificate remembers
the thing that proves loneliness: a small pair, a pinch time, and the exact
reason blockers fail.

## The Base Object

At the tight time, n=14 sees a 14-gon, not a generic 13-runner round
tournament. The six nonzero mod-7 lanes come in antipodal pairs, and residue
0 is the self-antipodal apex.

That makes the base:

```text
six paired lanes: +/-1, ..., +/-6 mod 7
one seam:         0 mod 7
```

The seam is where all the annoying threads meet: V*, zero-divisor behavior,
multiple-of-7/apex coupling, perturbation-direction ambiguity, and the failure
of a clean "the optimal class is self-converse" sentence.

## What A Section Is

A section is not a runner. It is a local proof certificate:

```text
(cheap pair, witness time, shield blockers, endpoint anchors, owner labels)
```

For HYP-2095, the important bit is whether the cheap pair is unblocked. For
HYP-2096 and HYP-2099, the important bit is which unit-spine and slack owners
carry the D/U/N and endpoint debt. For HYP-2083, the important bit is whether
the summand shell is visible to multiplication by a unit.

Those are different projections of the same local object.

## Why This Feels Like The Right Extension

Opus S570 says all `8191` gcd-1 transversals mod `27` are lonely and all have
an unblocked small pair. The incoming HYP-2100 unit-spine exchange sieve says
the same cheap-pair-first principle survives one-unit lifts: `13169/13169`
full D/U/N covers had an unblocked small pair. S578 also says the canonical
composite slack fibre behaves the same way through slack `42`: AP and V* are
the only floor rows, and both use `(1,13)` at `1/14`; the only block-all
controls are positive measure.
HYP-2102/Opus S571 adds a complementary positive-measure signal: bounded rows with a
multiple of `n` were all loose for `n=4,6,8,10,12,14`, but the crude lower bound
failed, leaving a thin-arc-cover residual. That is exactly the kind of failed
gluing branch that should use endpoint labels rather than a scalar bound.

So transversals and V* should not be separate proof chapters. They look like:

```text
transversals:      cheap section glues without crossing a nonunit seam;
unit-lift covers:  cheap section fires before exchange is needed;
mult-of-n rows:    apex seam goes positive-measure rather than tight;
V*:                cheap section survives the first nonunit apex lift;
block-all:         no section glues, so an endpoint/cover obstruction opens measure.
```

That is a much cleaner shape.

## What The Bounded Audit Added

The companion script `04-computation/lrc_n14_apex_lift_sheaf_s579.py` made this
finite over the HYP-2100 unit-spine lift site.  Objects are unit-shell lift
subsets, restrictions lower one unit representative, and sections are:

```text
ledger_failure
cheap_d14
cheap_side
positive_measure
residual
```

The all-slack one-lift scan has no residual section:

```text
full covers=13169
cheap_d14=7943
cheap_side=5226
restriction_residual=0
```

That is the key correction.  Denominator `14=2*7` is the big apex chart, but it
is not the whole sheaf.  Side denominators are not noise; they are required
charts.

The named two-lift stress rows show the gluing behavior:

```text
AP:       {cheap,cheap}:150, {cheap,ledger}:67, {ledger,ledger}:1
V*:       {cheap,cheap}:250, {cheap,ledger}:57, {ledger,ledger}:1
open-gap: {cheap,cheap}:150, {cheap,ledger}:67, {ledger,ledger}:1
```

The single `{ledger,ledger}` row in each layer is the warning label: a full
cover can appear only at the union object.  So a certificate sheaf must include
ledger-failure local sections, not just already-full local witnesses.

## Creative Extensions

### 1. Apex monodromy

Move a cheap-pair certificate around the six lanes and back through the apex.
If it returns to itself, the section glues. If it returns changed, the defect
is an obstruction cycle. The conjecture is that every nontrivial obstruction
cycle carries a private endpoint owner and therefore positive measure.

This is the "cohomology" version, but the actual data are small: pair labels,
shield labels, anchor labels, and D/U/N pivots.

### 2. Nonunit descent as certificate transport

The old language says a speed moved into gcd-3/gcd-9 mod `27` becomes invisible
to unit clocks. The new language says the unit-clock section must be transported
to an n-clock or pair-sum section.

V* is the example: the mod-27 unit witness cannot see the missed nonunit pair,
but the n-clock cheap pair `(1,13)` still proves the floor.

### 3. Positive measure as failed gluing

Instead of making block-all rows a residual bad set, declare them failed
gluing attempts. If every cheap local section is shielded or anchored, then the
shield/anchor data should form a cover-circuit. HYP-2095 and THM-397 say those
circuits have endpoint owners. S574 says private D/U/N pivots are abundant.

The proof dream:

```text
failed gluing -> cover-circuit -> private endpoint/pivot -> peel -> positive measure.
```

### 4. The staircase warning

The new all-0 staircase values are enormous by k=10. That is a warning not to
confuse path complexity with proof complexity. Near-regular class-like objects
can have explosive Hamiltonian-path counts while a small local law still
governs the feature we care about.

For LRC n=14, that local law is not H. It is whether a cheap section survives
the apex lift.

### 5. A minimal S580 script

The next computation should be modest. For AP, V*, transversal representatives,
and the smallest gcd-3/gcd-9 lifts, record:

```text
best cheap pair,
reduced sum,
mod-7 lane ownership,
mod-27 shell,
shield blockers,
endpoint anchors,
D/U/N private pivots,
exact M and safe measure flag.
```

Then compare rows by restriction maps, not by speed-set inclusion. The output
should say things like:

```text
certificate (1,13) survives AP -> V*;
certificate X dies at apex but creates endpoint owner Y;
certificate Z changes reduced sum but preserves unblocked status.
```

That would be a proof-relevant table rather than another census.

## Tournament Analysis

The vertices I would use are certificate germs and failed transports.

The observable should ask:

```text
which germ preserves the cheap pair under more lifts?
which failure creates a lower-denominator endpoint owner?
which obstruction has fewer private pivots left to peel?
```

Orient toward the harder object: no cheap section first, then lower exact M,
then fewer private pivots, then deeper nonunit lift. A directed cycle would be
interesting only if it is a loop of mutually incompatible certificate
transports around the apex. Otherwise the tournament is expected to be a
transitive proof ledger.

## Bottom Line

The most creative extension is also the most conservative one:

```text
Keep the tie-wall.
Keep the apex.
Keep the cheap-pair labels.
Forget class identity whenever it stops preserving those labels.
```

The proof of n=14 may be a gluing theorem:

```text
every tight-boundary speed set carries a global cheap-pair section,
or the obstruction to gluing creates positive measure.
```

That statement is wild enough to be useful and local enough to be attacked.
