# Summand/Multiplicand Farey Basis Merge

The current additive-basis thread had one piece still too implicit.  HYP-2998
said the old Goldbach, Fermat polygonal, Zeckendorf, Fibonacci, and Farey
threads are representation economies.  HYP-2999 named the Pascal-slope packet
fields.  HYP-3000 fixed the user's row as a path-rank identity.  The missing
piece is that the two familiar graph shadows are not metaphors:

```text
p+q is a summand-graph node.
p*q is a multiplicand-graph node.
```

For a Farey packet `p/q`, the additive lane and product lane are two projections
of one labelled object.  The summand graph gives a dense antidiagonal fiber
`x+y=p+q`, exactly the old pinch-denominator world.  The multiplicand graph
gives a sparse hyperbola `xy=p*q`, exactly the factor/Kpq incidence world.  The
log map explains why these feel parallel, but the proof data they preserve is
different enough that quotienting one into the other is dangerous.

The Fibonacci row

```text
1, 1, 1+1, 1+2, 1+3+1, 1+4+3, 1+5+6+1, ...
```

is now best viewed as the finite additive fiber before canonicalization.  It is
the row `C(n-1-k,k)`, not merely the scalar `F_n`.  Zeckendorf then adds the
carry rule that turns this finite packet into a unique normal form.  That puts
Goldbach and Zeckendorf at opposite ends: Goldbach wants enough representations;
Zeckendorf kills all but one by local confluence.

The LRC14 branch split gets cleaner:

```text
p=1: q-parent/AP-GW boundary side
p=2: C27 petal / two-block summand side
p>=3: product-incidence / Kpq-K33 side
```

So `2/27` and `3/41` are not just nearby Farey children.  They are different
operation-graph failures: the first is still a second-gap summand/C27 packet;
the second is the first product-incidence wall.

Assumption challenge: the tournament vertices here are proof currencies and
representation fibers, not runners, arcs, or raw sequence entries.  The
quotient preserves the LRC packet only when exact Farey scale, additive fiber,
factor fiber, carry law, and proof economy have all been declared.  It destroys
which operation graph supplied the witness if those labels are collapsed to one
scalar.

Next use: add `summand_fiber_id`, `multiplicand_factor_fiber`,
`zeckendorf_carry_width`, and `farey_shadow_lane` to any packet classifier that
uses sequence shadows.  The user prompt is pushing toward a clean rule: sequence
data is useful only after we state whether it is abundance, bounded invoice,
carry normal form, additive antidiagonal, product hyperbola, or stress test.
