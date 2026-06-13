# Signed LRC And Unit-Distance Trienerment, S674b

The useful thing about negative speeds is almost perversely simple:

```text
||v t|| = ||-v t||.
```

So independent direction choices do not change the LRC observer predicate.
That means they are a gauge symmetry of loneliness, not a new version of the
problem.

But the pairwise runner movie is not gauge-invariant.  If runner `i` has sign
`eps_i`, then the relative clock between `i` and `j` is:

```text
eps_j v_j - eps_i v_i.
```

Same signs give differences.  Opposite signs give sums.  That is the trick.
The observer sees the same clearance sequence; the pair tournament sees a
different additive address.

## The Clean AP13 Witness

For total `n=14`, the shell modulus is:

```text
C = 2n - 1 = 27.
```

S674b samples AP13 under all-plus and alternating directions.  The observer
clearance sequence is identical.  The pair tournament is not:

```text
all plus:         100 unique tournaments, 18 fingerprints
alternating plus: 355 unique tournaments, 29 fingerprints
```

The exact pair-clock split is the gem:

```text
same-sign differences = 2,4,6,8,10,12
opposite-sign sums    = 3,5,7,9,11,13,15,17,19,21,23,25
```

So alternating direction inserts the odd interior pair-sum sieve directly into
the runner-runner dynamics over `C=27`, while the observer shadow is unchanged.

This is a beautiful quotient test.  If a proposed LRC14 summary remembers the
pair-sum shell, it should split all-plus AP13 from alternating AP13.  If it
does not, then it is still too observer-scalar.

## What This Says About Pair-Sum Modulus

The old slogan was that the pair-sum sieve modulus wants to be literally
`2n-1`, or its odd part.  S674b does not prove the theorem.  It does something
more local and maybe more useful:

```text
signed runner-runner pair clocks expose sums modulo C=2n-1
without changing the LRC observer predicate.
```

That gives a route to theorem form.  The additive face is not a separate
resonance energy.  It is a side-channel over the same observer projection.
Alternating signs are the laboratory switch that makes the side-channel
visible.

This also explains why the address-coordinate language keeps winning.  The
observer quotient is true but incomplete.  The pair-clock address is invisible
to loneliness at a fixed row, but it can be proof-relevant when classifying
which rows are coherent floor atoms and which are strict perturbations.

## Unit Distance Needs Three States

The unit-distance prompt asks for:

```text
distance = 1     as trienerment/tie edges
distance < 1     one orientation
distance > 1     the opposite orientation
```

That is better than a binary flip.  The binary map forces equality into one
side too early.  But the equality layer is exactly where the geometry lives.

S674b uses:

```text
S = short, follows base order
U = unit, tie layer, resolved positive or negative
L = long, reverses base order
```

The square-center toy already shows the impairment.  It has `S=4,U=4,L=2`,
but the center is isolated in the unit graph, so no unit Hamiltonian spine can
exist.  Any Hamiltonian path must pay a short or long impurity.

The Eisenstein rows show the opposite regime.  Hex-7, hex-19, and a simple
Eisenstein `n=21` fringe toy all have a unit spine in the equality layer.
For the `n=21` toy:

```text
pair states: U=44, L=166
base path states: U=20
unit-positive tournament: c3=157, strongly connected
unit-negative tournament: c3=0, reverse transitive
```

So the sign of the unit tie layer is not cosmetic.  With no short edges, making
unit and long share the negative side collapses the whole object to a
transitive order.  Making unit positive leaves a strong cyclic carrier.

## The Flop Should Be An Impurity Profile

HYP-2202 already says the fixed canonical tiling order can flop before the
geometric unit graph loses a Hamiltonian path.  S674b sharpens the diagnostic.

Instead of asking only:

```text
is there a unit Hamiltonian path?
```

ask:

```text
what is the S/U/L word along the path?
```

The first true flop, if it exists, should be measured by the minimum impurity
profile:

```text
0 impurities: pure unit spine
1 short impurity: one close-pair bridge is forced
1 long impurity: one nonunit jump is forced
mixed impurities: geometry has left the Eisenstein-spine regime
```

That gives a smaller-size program.  Build toys where a single short edge
isolates a point; then build toys where the unit graph remains traceable but
every natural tournament insertion path pays long impurities; then compare
against the actual `n=21/22` candidates.

## Assumption Challenge

For LRC, vertices should not default to runners.  Signed runners are useful,
but the predicate-preserving object is:

```text
observer clearance sequence + hidden pair-clock address.
```

Candidate vertices remain runners, pair clocks, residues, cover arcs, fixed
sections, wall-crossing events, Fourier modes, and proof obligations.  S674b
uses proof lenses and pair clocks because they preserve the distinction the
prompt cares about.

For unit distance, vertices should not default to points either.  The likely
proof vertices are:

```text
distance states,
unit-spine obligations,
Hamiltonian path edges,
deletion owners,
construction moves.
```

Points are where distances live.  Path obligations are where the theorem will
probably live.

## Next Moves

For LRC14, use sign gauges as a quotient stress test:

```text
same observer shadow, different pair-sum address.
```

Any proposed no-leak invariant should split those when it claims to remember
additive shell data.

For unit distance, write a spine-impurity scanner:

```text
input: exact point set
classify all pairs as S/U/L
find min impurity Hamiltonian path word
record unit-spine/deletion-owner profiles
```

Then run it on the known Moser/Eisenstein `n=21/22` lanes and any non-lattice
candidate sets we can import.

The deeper common pattern is now clear:

```text
gauge symmetry preserves target predicate;
side-channel changes proof address;
trienerment/equality layer is where the address enters.
```

That is a good shape.  It feels theorem-friendly.
