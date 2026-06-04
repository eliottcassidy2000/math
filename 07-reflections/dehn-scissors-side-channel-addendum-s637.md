# Dehn and Scissors Side Channels: S637

Incoming HYP-2211 already made the two-shadow correction precise: it is not
known that exactly one of `e+pi` and `e*pi` is rational.  What is known is
stronger in one direction and weaker in another: they cannot both be
algebraic.  If both elementary symmetric
coordinates

```text
s = e + pi
p = e*pi
```

were algebraic, then `e` and `pi` would be roots of `x^2-s*x+p`, so both would
be algebraic over the algebraic numbers.  That contradicts the
transcendence of both constants.  Thus at least one of `s,p` is transcendental;
if one of them were rational, the other would have to be transcendental.  The
"exactly one rational" version is not a theorem.

That correction does not kill the analogy.  It sharpens it.  Rational versus
irrational is not parity.  Odd/even is a real quotient `Z -> Z/2`; the
irrationals are not closed under `+` or `*`.  The better object is the side
channel that says which field the data lives over and what algebraic dependence
is being retained.  Addition and multiplication are not rival colors; they are
the two elementary symmetric handles on a hidden unordered pair.

This is the same move as cuboids versus simplices.  Volume is the scalar.  Dehn
data is the retained channel.  A cuboid has zero Dehn invariant; a regular
tetrahedron has nonzero Dehn data because its dihedral angle has
`cos(theta)=1/3`, so `theta/pi` is irrational by Niven's theorem.  A box can be
triangulated into simplices whose Dehn pieces cancel, and a simplex can be
packed inside a box, but equal volume alone is not equidecomposability.  Packing
is a weaker witness than scissors congruence.

That gives a clean bridge back into the repo's current spine:

- S633 says a hard sequence should be surrounded by shadows: fixed, merged,
  nonfixed, bisection, q-shadow, skinny quotient, and transporter quotient.
- S634 says a useful quotient must name the perspective/carrier it preserves:
  LRC source/sink target, rooted converse action, or unit-distance spine/bulk.
- HYP-2186/HYP-2187 say equinumerosity is count-level while equidecomposability
  is retained-fiber data.
- S630/S634 say the unit-distance `n=21` row is not a raw `H=21` echo.  It is
  `57 = 20 + 37`, with a unit Hamiltonian spine and centered-hex bulk.

So the working rule becomes: before collapsing to one number, ask what plays
the Dehn role.

For unit distance, the next creative impairments should be small and
side-channel targeted.  Fix an edge count and damage traceability.  Fix a
traceable graph and damage the direction-support mask.  Fix a `57`-edge
`n=21` carrier and vary the endpoint-compatible ears.  Build examples that
have the same scalar but different spine count, bulk shell, direction price, or
embedding obstruction.  These are not merely negative examples; they are the
geometry analogue of the `e+pi/e*pi` trap.  They show exactly which hidden
coordinate a proof method needs.

Tournament Analysis in S637 used proof lenses as vertices rather than numbers:
algebraic-dependence channel, Dehn invariant, symmetric-polynomial trap,
carrier compression, sequence-shadow packet, log bridge, parity analogy,
volume-only packing, and the raw exact-one claim.  The majority tournament was
transitive with `H=1`, ranking Dehn and carrier compression at the top and the
raw rationality claim at the bottom.  That is a good sign: the ranking agrees
with the corrected theorem and with the repo's recent carrier discipline.

The deepest connection is therefore not that rational/irrational, odd/even,
addition/multiplication, cuboid/simplex, LRC/unit-distance, and tournament
`H` are the same object.  The connection is that each problem becomes clearer
when a tempting scalar is treated as a projection and the missing coordinate is
made explicit.  The new technique is to search for scalar twins with different
side channels at the smallest sizes first, then let those controlled failures
tell us which carrier to preserve at scale.
