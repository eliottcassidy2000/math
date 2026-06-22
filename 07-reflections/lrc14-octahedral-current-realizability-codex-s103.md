# Octahedral Current Realizability

The S101 packet graph had one more piece of structure hiding in plain sight.
The three affine-zero pairs are a perfect matching on the six residue vertices.
Deleting that matching from `K6` gives the octahedron graph, equivalently
`L(K4)`.

That changes the realizability question.  We should not ask for an arbitrary
tournament or arbitrary signed graph current.  A repeated-packet lift chooses a
height cochain on the six edges of a tetrahedron, and its nonzero packet
couplings are exactly the adjacent-edge pairs.  The affine-zero lane is the
opposite-edge relation.

The octahedron has `12` edges and `6` vertices, so its cycle rank is `7`.  It
also has `8` triangular faces with one global boundary dependence.  That is a
better `7` than a slogan: the support-six repeated-packet residual has a
seven-dimensional face-curl module before the analytic reciprocal lift.

The layer scan makes the point numerically.  At `H=10`, the start-aligned
constant gauge has L1 divergence `0.0219283`, the raised constant gauge has
`0.00754806`, and the best scanned cochain `210210` has `0.00225361`.  Wall
max and wall incidence correlate positively with divergence, while mixed
positive/negative local signs correlate with smaller divergence.  This matches
the old "opposite bounded signs" intuition in a realizable finite carrier.

The proof target should now be an octahedral Hodge lemma:

```text
finite Kirchhoff balance on L(K4)
  + low-height wall deletion
  + face-curl ledger on the 8 triangular faces
  + HYP-2636 additive-channel Abel summation
  => lifted packet divergence bound.
```

The key topological advantage is that the octahedron is a sphere.  Once vertex
divergence and triangular face curl are controlled, there is no independent
harmonic one-current left to hide a global obstruction.  That may be the
realizability principle the LRC route was missing.
