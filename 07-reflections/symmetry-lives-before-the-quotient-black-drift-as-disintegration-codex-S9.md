# Symmetry lives before the quotient

The most useful reframe of the merged metagraph is to stop asking whether the
picture is symmetric and ask **at which level the symmetry lives**.

Concurrent THM-787 found the same score-energy axis from the opposite
direction.  Its variable is exactly eight times the reversed `C3` coordinate;
its unmerged-class census and this session's converse-merged-node census agree
once that quotient level is kept explicit.
HYP-6855's centered-shift proof and half-tiling count give a third description
of the blue fixed axis.  Agreement among energy, endpoint flux, and reflection-
orbit counts is evidence that this coordinate is structural rather than a
drawing convention.

At the labelled tiling level, complement is brutally exact.  With

```text
a=d0-n/2,  epsilon=d0+dlast-(n-1),
```

it sends `(a,epsilon)` to `(-a,-epsilon)`.  Blue tilings lie on the axis
`epsilon=0`.  Black tilings fill transverse layers.  The flux toward lower
score variance is the oblique functional `2a-epsilon`.  So a useful mental
image is not a one-dimensional chain but a centrally symmetric two-dimensional
cloud, observed through an oblique projection.

The iso-class quotient then folds points with unequal multiplicities into
fibres, merges converse classes, and forgets the sign of the endpoint defect
unless the line instance is retained.  Pure-blue, mixed, and pure-black are
not three intervals in the original cloud; they are three types of quotient
fibre.  Conditioning the symmetric measure on those fibre types need not be
symmetric.  The `n=7` reverse black drift is therefore a disintegration
phenomenon:

```text
symmetric labelled measure
  -> condition on an asymmetric family of quotient fibres
  -> orient by the oblique C3 functional
  -> asymmetric categorical current.
```

This viewpoint connects several old threads that otherwise look unrelated.

- Score majorization supplies the horizontal potential: `C3` is exactly
  inverse quadratic score energy.
- The half-tiling reflection thread identifies the blue fixed axis.
- The line-metagraph and simultaneous-isomorphism threads warn that endpoint
  node pairs do not determine a line orbit or its multiplicity.
- The Hamiltonian-path inverse theorem says each node fibre is an
  automorphism quotient of observer cuts; stabilizers are therefore plausible
  sources of conditional bias.
- Continued fractions and centered mechanical words show how an apparently
  symmetric static fibre acquires directed transport once an endpoint clock is
  retained.
- The prime-seven monodromy thread gives the analogous LRC warning: a static
  cyclic node can be symmetric while labelled owner transport has nontrivial
  carry.
- Concurrent THM-784 sharpens that warning: a fixed slow-owner rainbow chamber
  can contain arbitrarily many fast-owner walls.  Static node complexity and
  raw event count both miss metric refinement; the transport stalk must retain
  extent and owner incidence.

The right next object is consequently a **flux kernel**, not merely another
node statistic.  Its rows should be simultaneous-isomorphism line orbits and
its columns should retain

```text
(source node, target node, C3 step, signed epsilon,
 endpoint path orbits, stabilizers, multiplicity, colour).
```

Successively forget fields and measure the first stage at which the oriented
current changes.  This turns “what information must be preserved?” into an
objective minimality experiment.  It may also explain the visual left/right
effect without fitting it after the fact.

For the four-coordinate LRC object, this is a model rather than a solution.
The tournament node is a chamber label; the tiling fibre is an observer-cut
fibre; endpoint mechanical words give directed wall transport; owner and sheet
coordinates give monodromy; the metric threshold remains a separate stalk.
The lesson is that symmetry of the base never licenses deletion of transport
coordinates above it.

The sharp conjectural question is now: does the exact symmetric endpoint-bit
measure, disintegrated over `HP(T)/Aut(T)` fibres, force the sign of the black
categorical current?  If yes, the imbalance is structural.  If no, the first
size where it changes sign will identify the missing orbit datum more cleanly
than any drawing can.
