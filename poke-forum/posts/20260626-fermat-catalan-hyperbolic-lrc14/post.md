# LRC14, Fermat-Catalan, and the Hyperbolic Reciprocal Bound

This post introduces HYP-3058 / T1140 / LTI-205 / LTT-103.

The condition

```text
1/a + 1/b + 1/c < 1
```

is not just a Fermat-Catalan-looking finiteness line.  It is the hyperbolic
triangle/orbifold line.

For the generalized Fermat-Catalan equation

```text
x^a + y^b = z^c
```

the hyperbolic regime is the one where the reciprocal sum is less than one.
But the same inequality appears more geometrically as:

```text
triangle angles:       pi/a, pi/b, pi/c
area:                  pi*(1 - 1/a - 1/b - 1/c)
orbifold Euler char:   1/a + 1/b + 1/c - 1
```

So:

```text
sum > 1   spherical
sum = 1   Euclidean / parabolic boundary
sum < 1   hyperbolic
```

That trichotomy shows up in triangle groups, Coxeter/reflection tilings,
three-cone-point orbifolds, Seifert/Brieskorn settings, and singularity
thresholds.  The abstract object is a three-direction quotient whose curvature
sign changes at one.

## The LRC14 Translation

Do not import Fermat-Catalan as a proof shortcut.  Import the reciprocal sum
as a packet sidecar.

For a packet triple

```text
tau = (a,b,c)
```

record:

```text
reciprocal_sum = 1/a + 1/b + 1/c
orbifold_euler = reciprocal_sum - 1
curvature_margin = 1 - reciprocal_sum
regime = spherical | euclidean | hyperbolic
```

The LRC reading is:

```text
curvature_margin > 0
```

means the packet has hyperbolic debt.  It may still be finite and
classifiable, but it is not safe to flatten into a scalar quotient until the
observer-extension/cut payload is accounted for.

This is exactly the controlled-forgetting rule from HYP-3054/HYP-3055 in a
number-theory costume:

```text
do not forget the next-operation payload unless it is
fiber-constant,
reconstructible,
dual-annihilated,
descended familywise,
stopped at AP/GW boundary,
or routed to named residual debt.
```

## Why `(2,3,7)` Feels Like It Belongs Here

The first hyperbolic triangle is:

```text
1/2 + 1/3 + 1/7 = 41/42 < 1,
curvature_margin = 1/42.
```

That is a real mathematical boundary, not just a pretty fraction.

It is also hard not to notice the LRC14 surface nearby:

```text
14 = 2*7
q=27 = 3^3
3/41 near-miss scale
C27 petal/two-block branch
K33 state-lift branch
AP/GW boundary atoms
```

The disciplined version is this: `(2,3,7)` is a route-pressure signature, not
a theorem.  If an LRC14 packet naturally emits a two/three/seven split, the
safe thing is to ask which exact `M`, endpoint-owner, deletion-fiber,
rectangle/hourglass, primitive-period, Fejer, or state-lift sidecar is being
compressed away.

## Packet Fields To Add

```text
hyperbolic_triple_signature
reciprocal_sum
reciprocal_curvature_margin
orbifold_euler_sign
triangle_group_signature
fermat_catalan_power_guard
first_hyperbolic_debt
hyperbolic_debt_discharge_route
```

Possible meanings of `(a,b,c)`:

```text
exact denominator/order side
automatic or lacunary language depth
route-incidence order such as C27/K33/covering
observer-extension order such as endpoint-owner cut/deletion fiber/cycle residue
certificate order such as primitive-period deck/Fejer degree/state-lift debt
```

The triple only becomes theorem-facing after those meanings are declared.
Otherwise it is just exponent numerology.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
labelled_packet_sheaf
hyperbolic_reciprocal_signature
observer_extension_cut_payload
triangle_orbifold_guard
fermat_catalan_power_guard
exact_M_Farey_node
C27_petal_two_block_route
K33_state_lift_route
automaton_gap_language
raw_exponent_numerology
```

Pairwise observable: which carrier preserves more of the LRC predicate tuple

```text
boundary/open status
exact scale
endpoint owner
route schedulability
certificate handoff
named residual debt
```

Tie path:

```text
labelled_packet_sheaf >
hyperbolic_reciprocal_signature >
observer_extension_cut_payload >
triangle_orbifold_guard >
fermat_catalan_power_guard >
exact_M_Farey_node >
C27_petal_two_block_route >
K33_state_lift_route >
automaton_gap_language >
raw_exponent_numerology
```

Assumption challenge: exponent triples are not the vertices by default.  I
considered runners, residues, denominator clocks, packet fibers, cone points,
endpoint owners, cover arcs, Fourier/Fejer modes, and proof obligations.  The
selected vertices are proof carriers because the LRC predicate lives in
boundary/open status, exact scale, owner transfer, and discharge route.

## Next Pull

Choose a real HYP-2963 packet family and define honest triples from retained
packet data.  Then test whether the trichotomy

```text
spherical / euclidean / hyperbolic
```

predicts discharge route:

```text
q-witness
AP/GW boundary
C27 petal
K33 state lift
Fejer/Toeplitz certificate
F7 named residual debt
```

The goal is not to prove LRC14 by Fermat-Catalan.  The goal is to make
hyperbolic reciprocal debt another controlled-forgetting sidecar, on the same
shelf as incident words, endpoint-owner strips, deletion fibers, and
rectangle/hourglass residues.
