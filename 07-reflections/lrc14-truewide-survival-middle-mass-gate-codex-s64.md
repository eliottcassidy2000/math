# True-wide cap is a survival middle-mass problem

**codex-2026-06-20-S64.**  The useful move in HYP-2701 is not a new
inequality; it is the coordinate change that makes the old cap-floor target
look like the right kind of inequality.

THM-556 says

```text
U4 = p0 + p5 + 5p6.
```

Since the missed-sector probabilities sum to `1`,

```text
U4 <= (k-6)/7
```

is exactly

```text
p1+p2+p3+p4 - 4p6 >= (13-k)/7.
```

So the true-wide branch is not asking for a raw Bonferroni bound.  It is asking
whether true-wide geometry forces enough middle missed-sector mass to pay for
the fully-missed tail.  The cap comparison has the same left side; only the
right side changes to `1-cap_k`.  That is why the k=8 dividend is visible as
currency, rather than as an exception that feels external to the argument.

The exact scan makes the next target much sharper.  In the audited boxes there
are no cap failures, no k>=9 floor failures, and only three k=8 floor failures.
All three are cap-safe.  The worst one,

```text
(0,3,6,9,12,14,15,18),
```

has floor debt `107/8820` and cap slack `295/3528`.  More important than the
numbers is the profile: it is barely true-wide, with exactly two far speeds.
Every tight leader in the scan has `far_count=2`; the `far_count>=3` rows are
looser in every audited layer.

That suggests a proof split:

```text
two-far true-wide rows     -> survival-currency lemma
three-or-more far rows     -> separate margin lemma
k=8 floor failures         -> finite dividend templates
```

The two-far addendum makes the first line more precise.  If a bounded core `B`
is followed by two decorrelated far runners, the missed-count profile evolves
by a literal death chain: each far runner picks one of seven sector colors, and
the residual count drops only when a currently missed inner sector is hit.  The
resulting boundary currency is already safely above the floor for every bounded
core in the audited universe.  The smallest margins are still positive:

```text
k=9: 569/3430
k=10: 5717/36015
k=11: 5317/24010
k=12: 35543/123480
```

So the two-far proof is not "find the missing middle mass from scratch."  The
boundary value already has it.  The real task is to bound how much a resonant
far pair can subtract from that boundary value:

```text
actual currency = boundary currency + signed two-far deviation.
```

The tight rows all have small relation-distance far pairs, usually consecutive
or near-consecutive.  That aligns perfectly with THM-548/HYP-2679: off
resonance, signed Abel/BV should control the deviation; on resonance,
Freiman/scale reduction should hand the pair to a finite atlas.  This is a more
workable theorem target than a global scalar survival bound.

This also clarifies how to use Tournament Analysis.  The vertices should not
default to runners, arcs, gaps, colors, or far elements.  Those are all possible
coordinates, but the quotient that preserves the current LRC predicate is a
risk order on proof obligations and row families.  The pairwise observable is
floor slack; the switch is how much cap-floor currency a family spends.  In
this quotient the tournament is transitive, and its Hamiltonian path starts
with the three k=8 dividend rows before the k=9 and k=10 two-far leaders.

The challenged assumption is that a scalar like `p0`, `U4`, additive energy, or
a binary tournament statistic should be the final true-wide observable.  HYP-2701
keeps the scalar only after choosing the survival basis.  Near equality must
still be lifted back to sector labels, transfer-tax packets, residual-profile
coordinates, and relation-lattice phase.

The incoming HYP-2702 sparse-tail automaton gives the same lesson in another
coordinate system.  Raw sparse residual coordinates are too strong, but
miss-zeta generated-context words look tractable.  Here raw Bonferroni language
was too opaque, but survival middle mass exposes the proof currency.  In both
branches, the work is less about discovering another scalar and more about
choosing the coordinate system before scalarizing.
