# LRC Restricted Tournament Quotient Stack

codex-2026-06-01 S535

The user's request was not just "map LRC to tournaments again."  The sharper
instruction was: keep the form of the question as a statement about which
tournament isomorphism classes can be exhibited, but make the class set more
meaningfully restricted.

That means the ambient object should not be:

```text
all tournament classes A000568(n)
```

It should be:

```text
the image of arithmetic clock geodesics inside a stack of marked, thresholded,
kinetic, and pressure-labelled tournament quotients.
```

The proof target becomes more abstract but also more honest:

```text
Every primitive speed system exhibits at least one class in the pure source
substack, or in its compactified wall boundary.
```

## The New Mappings

### 1. Raw Circular Body

Vertices are moving runners.  Edges are the half-turn circular comparator.
This is the ordinary A000568-style body.  It is the wrong proof object by
itself because good and bad states share classes.

S535 repeats the defect:

```text
n=4: 2 raw phase classes, both mixed
n=5: 4 raw phase classes, all mixed
n=6: 11 raw phase classes, all mixed
```

The raw body is still useful as the base manifold.  It carries H, SCCs,
round/circular structure, and chord-channel geometry.  But it does not know
the observer threshold.

### 2. Source-Deleted Target Menu

Restrict to states where the observer is a source, then delete the observer.
The remaining runner tournament asks:

```text
Which runner-subtournament classes can be exhibited at lonely/source times?
```

This is the cleanest answer to the user's class-exhibition framing.  It is a
small target set inside A000568(n-1), not a classification of every state.
HYP-1999 identifies the open part as the Ferrers interval menu; S535's sampled
version includes walls too, so it sees compactified boundary classes.

### 3. Observer-Source Marked Quotient

Keep the observer as a colored vertex.  Orient observer incident edges by the
actual LRC threshold:

```text
observer -> runner iff ||v_i t|| >= 1/n.
```

This is THM-381/THM-385 in quotient language.  Goodness becomes class-local:
the observer is a source.  The price is a label tax: many more colored classes
than raw phase classes.

### 4. Gap-Threshold Necklace

Vertices are circular gaps, not runners.  Orient by gap length; color the two
observer-adjacent gaps and whether each gap is at least `1/n`.

This is the outside clasp view.  THM-384 says the LRC witness is exactly:

```text
both observer-adjacent gap vertices are long.
```

The class set is no longer just A000568 of runner phases.  It is a tournament
on the cyclic gap necklace with two marked clasp positions.

### 5. Kinetic Gap-Flow Necklace

Use the same gap vertices, but let the edge switch compare:

```text
(long?, growth sign, gap length)
```

This is the instantaneous THM-387 flavor: the source gap is not merely a static
object; it has directed flow.  In S535 this begins to be a stricter quotient at
`n=7`, where the static gap map sees `68` classes and the kinetic map sees
`50`.

That difference is an abstract metric:

```text
kinetic_torsion = image(static gap quotient) - image(kinetic flow quotient).
```

When kinetic_torsion is positive, monotonicity has collapsed several static
possibilities into fewer directed-flow possibilities.

### 6. Blocker-Deficit Shadow

Vertices are moving runners.  Edges rank endpoint deficit:

```text
max(0, 1/n - distance_to_observer).
```

Colors record whether a blocker is on the left or right side of the observer.
The good class is tiny: all colors safe.  This quotient is close to
tautological, but it is valuable as a repair layer: it measures the local
neighborhood around source, almost-source, and blocker handoff fibers.

### 7. Apex-Boundary Runner Body

Vertices are moving runners with the raw phase tournament, but the two runners
bracketing the observer are colored by side and by whether the observer-adjacent
gap is long.

This is HYP-2008 as a quotient map.  The observer-source gap is the
source-sink/apex arc.  Instead of marking every safe runner, this map marks the
two boundary actors that determine the source-sink opening.

It has higher label tax than blocker-deficit, but it preserves the actual
runner tournament body and the apex identity at the same time.

## Computed Shape

The S535 audit over primitive clocks gives:

```text
raw phase:
  mixed fibers persist, certification mostly zero.

source-deleted:
  small target menu, all sampled systems certified.

observer-source / gap-threshold / kinetic-gap / blocker / apex:
  zero mixed fibers in all audited windows.
```

The map meta-tournament is transitive:

```text
source_deleted_phase
> blocker_deficit_shadow
> gap_kinetic_flow
> observer_source_marked
> gap_threshold_necklace
> apex_boundary_runner
> phase_runner.
```

This ordering should not be read as "source_deleted is the only good map."
It says the maps occupy different levels:

```text
source-deleted = target language
blocker-deficit = local repair layer
gap-kinetic = flow monotonicity layer
observer-source = exact THM-381 layer
gap-threshold = outside clasp geometry
apex-boundary = source-sink actor identity
raw phase = base manifold
```

## More Abstract Metrics

I want future sessions to stop asking for one number.  The object is now a
bundle; it deserves bundle metrics.

```text
source_codimension_bits
  log2(A000568(m) / target_image_size).

label_tax_bits
  Colored-image expansion needed to make goodness class-local.

purity_curvature
  mixed_fiber_count times entropy.  Measures how badly a quotient folds good
  and bad chambers through the same class.

kinetic_torsion
  Static threshold classes minus kinetic-flow classes.

blocker_viscosity
  Entropy and depth of non-source endpoint-deficit layers.

apex_label_tax
  Extra classes caused by remembering which two runners bracket the observer.

compactification_index
  Wall target classes minus open target classes.

pressure_survival_rank
  Size of the THM-380 owner-compatible endpoint core after projection.

diameter_channel_state_entropy
  HYP-2015's independent-pair state entropy inside the hidden diameter channel.

monodromy_defect
  Classes gained by completing the full arithmetic period compared with one
  local chamber pass.

fiber_hysteresis
  Difference between entry classes and exit classes for the same source-gap
  wall, especially around THM-387 LS -> LL -> SL flow.

quotient_lens_resolvent
  The smallest subcollection of maps whose class data separates every sampled
  good state from every sampled bad state.
```

These metrics are deliberately abstract.  The point is not to crown them all.
The point is to make the project ask sharper questions:

```text
Which quotient kills mixed fibers?
Which quotient keeps the target menu small?
Which quotient preserves pressure debt?
Which quotient has a monotone transition law?
Which quotient exposes the AP wall as a single compactified class?
```

## Whackier Extensions

### CRT Aperture Tournament

For composite `n`, make vertices the CRT residue fibers.  Orient fiber `a` to
fiber `b` if `a` has more unit-wall aperture than `b`, breaking ties by
endpoint debt.  For `n=145=5*29`, this would measure the zero-residue embryo
against the unit-wall aperture from HYP-1989.

### Pressure-Sheaf Tournament

Vertices are quotient classes themselves.  Edge `C -> D` if every
source-avoiding lift of `C` exports endpoint debt into a lift of `D`.  LRC
becomes acyclicity or sink-peeling of this sheaf tournament.

### Rapidity Defect Tournament

Use the rapidity coordinate from HYP-1992:

```text
r(g) = atanh(1 - 2g).
```

Vertices are gaps or blockers.  Orient by rapidity surplus over the threshold
rapidity `0.5 log(n-1)`.  This might linearize gap competition enough for a
convexity/Farkas certificate.

### Boundary-Word Tournament

Encode a full clock period as a word in wall crossings:

```text
endpoint wall, phase wall, adjacent swap, reset wall.
```

Vertices are word factors.  Orient by forced precedence.  Isomorphism classes
become factor-order tournaments, and LRC asks whether every primitive word
contains an LL/source factor from the allowed menu.

### Diameter-Channel Fiber Product

Take the HYP-2015 diameter independent-pair state and fiber it over the
apex-boundary quotient.  This may be the finite classifier needed for:

```text
wall-only at n=14 => AP.
```

## Conclusion

The right sentence is no longer:

```text
LRC is a statement about A000568 tournament classes.
```

The better sentence is:

```text
LRC is a forced-exhibition theorem inside a restricted quotient stack over
A000568, where raw circular classes are the base and threshold/kinetic/pressure
labels carve out a small pure source substack.
```

That is more complicated, but it is also more constrained.  A counterexample
can no longer be "a loop in tournament space."  It must be a loop that avoids
source-deleted Ferrers classes, avoids THM-384 long-long gap classes, survives
THM-387 kinetic monotonicity, keeps blocker-deficit nonzero, preserves
apex-boundary closure, and still carries a THM-380 endpoint-pressure core.

That is a much smaller monster to hunt.
