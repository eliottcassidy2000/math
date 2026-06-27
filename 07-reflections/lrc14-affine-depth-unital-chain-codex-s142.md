# LRC14 Affine-Depth Unital Chain

The repeated triangular/perfect-number prompt now has a more concrete LRC14
hook than HYP-2941 alone.  HYP-2941 said the equality is false as a scalar
identity but useful as an affine-depth warning.  HYP-2942 then supplied exactly
the kind of labels the warning needed: calibrated q=3 unital blocks for GW,
near/K33, and petal components.

The new packet grammar assigns a component depth

```text
1 + depth_gcd(hole) + depth_gcd(double)
```

to each marked C27 transfer.  This gives depths `3` for GW `12->24`, `4` for
near-miss `12->36`, and `1` for the unit petal `10->20`.  Composing those
components as affine words `b a^d` is noncommutative.  The S140 linked order

```text
GW -> near/K33 -> petal10
```

has suffix depths `[8,5,1]` and sum `14`.  The other five permutations of the
same depths miss `14`.  This is the useful signal: LRC14 appears as a calibrated
path-order invariant, not as a triangular-number slogan.

The proof route remains conservative.  Exact `M`/Farey branch still comes
first; C27 marked transfer and unital pair-completion still carry the local
packet; K33/octahedral/Clebsch still provide the state-lift address.  The affine
depth-14 signature is a test for whether a residual has entered the linked
nonunit path rather than a unit-only splice.  Unit-only splices should be killed
by C27 petal/two-swap rigidity.  Depth-14 linked packets should be routed to
HYP-2908 / THM-572.
