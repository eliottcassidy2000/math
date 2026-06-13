---
source: codex-2026-06-03-S576
status: proof-program synthesis + finite audit
tags:
  - lonely-runner
  - parity-ladder
  - simplex-polygon
  - proof-program
  - nonunit-descent
  - tie-discharge
  - converse-quotient
  - endpoint-owners
---

# LRC Parity-Ladder Proof Program

HYP-2091 changes the shape of the proof search.

Before it, a lot of the repo was asking:

```text
Which clocks matter?
Which quotient forgets too much?
Where do the 2n-1 antipodal witnesses fit?
Why do odd/even and addition/multiplication keep reappearing?
```

HYP-2091 answers the geometry part:

```text
even LRC n -> clean polygon tournament
odd LRC n  -> wall/mesh tie-resolution tournament
```

That is not a proof, but it is a sorting machine. Once the geometry is sorted,
the remaining proof burden is no longer mysterious. It is the clock burden:
unit versus nonunit summand shells, D/U/N covers, pair-sum endpoint owners, and
private pivots.

## Four Lanes

Use `C=2n-1`.

```text
clean unit lane
  n even, C prime or all shells unit-visible.
  No antipodal tie burden and no nonunit shell burden.
  Target: D/U/N exchange on the converse-merged round seam.

clean composite lane
  n even, but C composite.
  Geometry is clean; arithmetic is not.
  Target: D/U/N exchange plus gcd-stratum descent.

wall unit lane
  n odd, C prime or all shells unit-visible.
  Arithmetic is clean; geometry is not.
  Target: discharge antipodal tie resolutions.

wall composite lane
  n odd and C composite.
  Both burdens are present.
  Target: tie-discharge and gcd descent must commute.
```

This four-lane split is more useful than the old single ladder. It says which
lemma has to fire before we spend time on a row.

## What n=14 Becomes

The important line of the S576 table is:

```text
n=14
m=13
C=27=3^3
lane=clean_polygon
unit/nonunit shells=9/4
D/U/N obligations=34
round/converse seam=190 merged nodes
```

So fourteen is not hard because the regular polygon has antipodal ties. It
does not. The runner tournament has odd size `m=13`, so HYP-2091 puts it on
the clean polygon ladder.

Fourteen is hard because the clean geometry is sitting over composite
arithmetic:

```text
C=27
gcd 1 shells: 9
gcd 3 shells: 3
gcd 9 shells: 1
```

This makes the n=14 theorem target wonderfully narrow:

```text
attach D/U/N labels to the 190 converse-merged round nodes;
track endpoint owners for pair-sum blockers;
prove every gcd-3/gcd-9 nonunit defect descends or produces a second clock;
then apply private-pivot exchange to force a floor or pair-sum witness.
```

That feels like progress because it names what is still alive after all the
quotients. The bare round class, the binary time word, and the runner
tournament each forget one of these labels.

## Which Clocks Matter

The clocks worth preserving are the ones that already have theorem-level
escape meanings.

```text
pair-sum clocks
  The exact maximin of a 1D lower envelope occurs at pair-sum or antipodal
  active-set times in the S570 audits. THM-397 says failed small pair-sum
  covers expose endpoint blockers.

n-clock
  The equality/floor branch. This is the AP wall, the observer-source escape
  clock, and the place where evenness re-enters as a midpoint/apex.

2n-1 unit clock
  The summand graph clock. Addition creates antipodal shells; multiplication
  by units makes missed shells visible.

D clocks
  Small-denominator divisibility. Failing one gives an immediate rational
  witness.

labelled event clocks
  Runner owner, endpoint owner, pair-sum denominator, wall side, and private
  obligation labels. These are the clocks that survive quotient loss.
```

Clocks to ignore until they carry labels:

```text
primitive reset length
  In normalized integer rows it is period 1. It sorts commensurable from
  dense torus cases, but it does not distinguish hard integer rows.

binary cyclic lonely-word stabilizer
  S571/S573 found it usually trivial. Reflection symmetry is real, but binary
  words forget owners.

unlabelled round class
  It is the correct open body, but not the proof object.

runner-vertex tournament alone
  It is often a lossy shadow of endpoint and obligation incidence.
```

## Wild Remodels

Here are the more aggressive frames that now look mathematically disciplined
rather than decorative.

### Sheaf View

Base space:

```text
converse-merged round nodes
```

Fibre:

```text
D/U/N owner labels, pair-sum endpoint owners, wall tie resolutions
```

The proof asks whether a counterexample can be a global section of this fibre
bundle without ever exposing a source gap. The likely obstruction is not a
cohomology group in any literal current sense, but the shape is right: every
quotient forgets local labels, and the proof must show the forgotten labels
cannot be glued consistently around the boundary seam.

### Automaton View

State:

```text
(parity lane, C gcd strata, D/U/N cover, private labels, endpoint owners)
```

Transitions:

```text
n -> n+2
speed exchange
gcd descent
tie resolution
pair-sum pinch
```

Accepting states:

```text
floor witness or pair-sum witness
```

Bad states:

```text
full cover, no private exchange, no endpoint discharge, no source gap
```

The conjecture is simply that bad states are empty.

### Matroid/Circuit View

The D/U/N hypergraph should be treated like an incidence matroid shadow:

```text
minimal full cover = circuit-like object
private obligation = leaf/pivot
exchange           = basis move
```

The AP row is the lower clean circuit. Open-gap lifted rows are not failures;
they are different circuits whose refined pair-sum witness sits above the
floor. The proof target is not "all covers are AP"; it is "all cover circuits
have a witness."

### Tie-Discharge View

Odd LRC `n` has even runner size `m`, hence antipodal tie diagonals. Instead
of enumerating `2^(m/2)` mesh resolutions, the wall proof should try to
discharge each tie pair:

```text
choose tie orientation
-> either observer source gap appears
-> or the oriented tie becomes an owner-labelled obligation
-> or the row reduces to a clean neighboring lane
```

This is where HYP-2091 makes a real promise: the odd ladder's exponential tie
choices should not be globally enumerated if each tie has a local discharge.

## The Immediate Computation We Need

The next computation should be the n=14 fibre table:

```text
190 converse-merged round node
  -> multiset of D/U/N owner labels
  -> pair-sum endpoint owners from THM-397
  -> active floor/pair-sum witness type
  -> private-pivot profile
```

This table would test the clean-lane theorem directly. It would also tell us
whether AP/V* are isolated as special labelled fibres or merely two visible
points in a larger exchange-connected component.

The next theorem target should be:

```text
clean composite lane theorem for n=14:
every labelled fibre over the 190-node round/converse seam either has an
n-clock witness, has a pair-sum witness, or descends along a gcd-3/gcd-9
private pivot.
```

That is still hard. But it is much smaller than the original continuum.

## Tournament Analysis

For S576:

```text
vertices: LRC n-ladder rows
observable: (tie pairs, nonunit shells, D/U/N total, merged round nodes, D_ext)
switch: larger remaining proof burden wins
tie path: increasing n
```

The fingerprint is transitive:

```text
directed_3_cycles=0
sccs all singletons
```

That is expected. A routing table should be transitive. The next nontrivial
tournament should live inside a fixed lane, probably with vertices:

```text
converse-merged round nodes
```

or

```text
D/U/N owner-labelled obligations
```

The pair observable should then be private-pivot exchange reachability or
endpoint-owner pressure. That is where cycles, SCCs, and edge flips can become
proof information rather than decoration.

## Short Version

HYP-2091 says:

```text
split geometry first.
```

S576 says:

```text
then split clock burden.
```

For `n=14`, the geometry is already clean. The proof should stop worrying
about antipodal tie walls there and focus on the composite `C=27` nonunit
descent plus owner-labelled D/U/N exchange over the 190-node converse-merged
round seam.
