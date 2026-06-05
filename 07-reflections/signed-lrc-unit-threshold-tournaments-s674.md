# Signed LRC And Unit-Threshold Tournaments

The negative-speed idea is sneakier than it first looks.  If a runner changes
from speed `v` to speed `-v`, the stationary observer sees exactly the same
distance at every time.  So negative speed does not create a new LRC geometry.
It creates a gauge cover of the same geometry.

That sounds like a dead end until the pair clocks are inspected.  For two
signed runners the relative clock is

```text
s_i v_i - s_j v_j.
```

Same signs give differences; opposite signs give sums.  So the sign gauge
turns the repo's `C=2n-1` pair-sum ledger into ordinary relative motion.  This
is the real trick: observer-blind, pair-visible.

In the n=14 floor rows, the finite scout gives a crisp signal.  AP, Vstar, and
2AP all keep the same observer-safe mask under the tested sign patterns, as
they must.  But Vstar has a sign lift with a zero pair clock because

```text
3 + 24 = 27.
```

AP and 2AP do not have such a signed zero clock in the same audit.  That makes
Vstar's old nonunit seam less mysterious: in the signed double cover it is not
only a missing/nonunit shell fact, it is a literal pair reset.  The observer
does not see the reset directly, but the owner/carry proof layer should.

This suggests a proof route for LRC14.  Search sign gauges as if they were
Ising variables on the complete runner graph.  Maximize zero clocks, gcd-9
clocks, or gcd-3 clocks modulo `C=27`.  Then ask whether the resulting signed
pair-clock pressure descends to a legal owner/carry fiber.  If it does, it
should land in AP, Vstar, 2AP, or a scalar lift.  If it does not, HYP-2241's
owner-private deletion bit should be exactly the address that constructs a
strict-tax witness.

So the sign gauge is not a replacement for the existing pair-sum support
calculus.  It is a way to animate it.  Pair sums become clocks; clocks have
zeros, gcd shells, cycle structure, and local update rules; those can be fed
into tournament and automaton tools already living in the repo.

The unit-distance half of the session has the same moral.  The old binary
mapping, "flip unit edges" or "flip nonunit edges", silently assumed that
lengths only mattered as unit versus nonunit.  The more honest threshold
tournament has three states:

```text
d < 1    close / pressure
d = 1    unit / wall / trienerment tie
d > 1    far / slack
```

For normalized unit-distance constructions where minimum distance is `1`, the
`d<1` state vanishes.  Then the three-state model collapses back to the binary
unit/nonunit flip, which explains why the earlier lattice and Moser unit-spine
work was not actually missing a third state.  It was looking at a face where
the third state is empty.

The scout's triangular-spine rows through `n=10` make this visible.  There are
no close pairs, the unit graph is traceable, and both choices of unit sign have
a directed all-unit Hamiltonian path.  The sign only chooses whether the spine
is read forward or backward.  The compressed-line and square-with-center toys,
where close pairs exist, behave differently and lose the all-unit directed
path.  That is the regime where the threshold tournament has real information.

This sharpens the unit-spine question: a genuine flop is not caused by choosing
unit-positive versus unit-negative.  It should be caused by nontraceability of
the extremal unit graph, or by a deletion/perturbation stage where close
compression edges enter before normalization.

The shared abstraction is a hidden side channel that is invisible to the
primary scalar.  Signed LRC keeps observer loneliness fixed while changing
pair clocks.  Unit thresholds keep point count and unit count fixed while
changing close/wall/slack pressure.  In both cases, the right tournament
vertices are not necessarily runners or points; they may be gauges, pair
states, spine steps, owner deletions, or proof obligations.

## Concrete Leads

1. Extend S674 over all HYP-2164 `Res_27` survivors and record the sign pattern
   that maximizes `(zero, gcd9, gcd3, opposite)` pair-clock pressure.
2. Add owner-private deletion and carry-vector labels to that sign table, then
   test whether HYP-2241's no-leak channel separates all signed pressure
   fibers.
3. Treat sign choices as a MaxCut/Ising optimization problem on residues, with
   edge weights given by `gcd(shell(s_i r_i - s_j r_j), 27)`.
4. Build a unit-distance threshold audit for Moser `P_2^-` (`n=21`) and
   `P_2^+` (`n=22`), first normalized, then under one-point deletion and
   small perturbation before rescaling.
5. Identify the threshold states with the repo's three-state automaton:
   close = left/pressure, unit = middle/wall, far = right/slack.
