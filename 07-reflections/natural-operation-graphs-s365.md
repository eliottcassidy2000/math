# Natural Operation Graphs and Product-Sum Collisions

codex-2026-05-31 S365

## Starting Point

The user prompt asks for the old incomplete-tournament idea on natural-number
modes:

```text
x -> z and y -> z iff x + y = z
```

and for the sparser multiplicative analogue

```text
x -> z and y -> z iff x * y = z.
```

The repo already had several pieces of this picture:

- `04-computation/summand_graph_investigation.py` treats the distinct-summand
  graph and finds the `{2,3}` dependency module with exceptions `{1,4,6}`.
- `07-reflections/summand-graph-fermat-zeckendorf.md` connects that graph to
  Fermat polygonal decomposition, Zeckendorf structure, and triangular phase
  transitions.
- `04-computation/natural_numbers_s116i.py` says multiplication is linear in
  rapidity space while addition remains external.
- `04-computation/sum_product_tournament.py` studies the finite-field
  sum-product tension inside circulant tournament connection sets.

The new S365 computation is in
`04-computation/natural_operation_graphs_s365.py`, with output stored in
`05-knowledge/results/natural_operation_graphs_s365.out`.

## The First Correction: This Is A Cospan Hypergraph

The additive rule is not fundamentally a unary graph rule.  It is a labeled
two-input operation fiber:

```text
(x,y) --+--> z
       |
       +--> x -> z and y -> z
```

For addition, forgetting the label `y=z-x` destroys nearly all information.
On positive naturals:

```text
x -> z in the additive simple shadow
iff exists y >= 1 with x+y=z
iff x < z.
```

So the unlabeled additive graph on `{1,...,N}` is exactly the complete
transitive tournament/order.  Its transitive reduction is only the successor
chain.

For multiplication, after suppressing the unit edge convention:

```text
x -> z in the multiplicative simple shadow
iff exists y >= 2 with x*y=z
iff x is a proper nonunit divisor of z.
```

So the multiplicative simple shadow is the divisor DAG, and its Hasse skeleton
is

```text
x -> x*p,  p prime.
```

This is the conceptual asymmetry:

```text
addition forgets to total order;
multiplication forgets to prime factorization.
```

The additive graph only becomes interesting when we keep its fibers.  Those
fibers are partitions of `z`, Goldbach-style representation counts, and the
old distinct-summand graph.  The multiplicative graph remains interesting even
after forgetting labels because divisibility survives as a sparse partial
order.

S365 measured this collapse.  At `N=80`, the additive simple graph has all
`3160` transitive-order edges.  The multiplicative nonunit graph has only
`209`, about `6.6%` as many.

## Product-Sum Equations Are Operation Collisions

The family

```text
x + y = x*y
x + y + z = x*y*z
x_1 + ... + x_k = x_1 * ... * x_k
```

is exactly the locus where the additive k-input operation and multiplicative
k-input operation land at the same natural-number mode.

Examples:

```text
(2,2)             -> 4
(1,2,3)           -> 6
(1,1,2,4)         -> 8
(1,1,1,2,5)       -> 10
(1,1,1,3,3)       -> 9
(1,1,2,2,2)       -> 8
```

The collision target `z` is the common output:

```text
sum(tuple) = product(tuple) = z.
```

This suggests a new invariant of a natural mode:

```text
collision spectrum of z = all arities k and cores c for which
                         sum(c with ones) = product(c) = z.
```

The mode `8`, for instance, is hit at arity `4` by `(1,1,2,4)` and at arity
`5` by `(1,1,2,2,2)`.  That is a small resonance between a two-factor core
and a three-factor core.

Grouped by target through arity `20`, the richest small mode is already the
highly factorable number `24`:

```text
z=24:
  k=12 core=(2,12)
  k=15 core=(3,8)
  k=16 core=(4,6)
  k=17 core=(2,2,6)
  k=18 core=(2,3,4)
  k=19 core=(2,2,2,3)
```

Then come `16`, `12`, `18`, `20`, and `28`.  This is weak evidence for a
collision-spectrum invariant that favors numbers with many multiplicative
partitions, but the defect filter makes it subtler than divisor count alone.

## Defect-Core Theorem

Delete all 1s from a positive solution.  Call the remaining nonunit core

```text
c = (c_1,...,c_m),  c_i >= 2.
```

Let

```text
P(c) = product_i c_i
C(c) = sum_i c_i
D(c) = P(c) - C(c).
```

Then:

```text
c padded with q ones solves sum=product
iff q = D(c) >= 0.
```

Proof:

Adding 1s changes the sum and leaves the product fixed.  A tuple with `q`
ones and core `c` has

```text
sum    = q + C(c)
product = P(c).
```

Equality is exactly `q=P(c)-C(c)`.  Conversely, if that integer is
nonnegative, padding by that many 1s gives a solution.

So 1 is not a harmless identity here.  It is a defect absorber.  Multiplicative
identity becomes additive slack.

This is the cleanest bridge I see to the summand graph.  The distinct-summand
setting removes the binary collision `(2,2)`, and the first genuinely mixed
collision becomes `(1,2,3)`.  The old exceptional nodes `{1,4,6}` now look less
accidental:

- `1` is the defect absorber.
- `4` is the binary diagonal collision `(2,2)`.
- `6` is the first ternary collision `(1,2,3)` and also the first triangular
  packet `1+2+3`.

## Finite Prefix Gate

For sorted positive solutions

```text
1 <= a_1 <= ... <= a_k
sum_i a_i = product_i a_i,
```

we have

```text
a_1 * ... * a_{k-1} <= k.
```

Proof:

Let `P=a_1*...*a_{k-1}`.  Then

```text
P*a_k = sum_i a_i <= k*a_k.
```

Since `a_k > 0`, divide by `a_k`.

This gives a very small exact search.  Enumerate sorted prefixes of length
`k-1` with product at most `k`; once a prefix has product `P` and sum `S`, the
last entry is forced by

```text
x*(P-1) = S.
```

This is a practical speedup over a naive k-dimensional search.  The S365
script uses this gate to enumerate all sorted positive solutions through
arity `20`.

## Shifted Divisor Law For Two-Element Cores

For a two-element nonunit core `(a,b)`, the defect is

```text
D(a,b) = ab-a-b.
```

The full tuple has length

```text
k = 2 + D(a,b) = ab-a-b+2 = (a-1)(b-1)+1.
```

Therefore:

```text
(a,b) is a two-core at arity k
iff (a-1)(b-1)=k-1.
```

This is a tiny but satisfying bridge between the product-sum collision stack
and the multiplicative divisor graph.  The arity `k` sees the divisor pairs of
`k-1`, shifted by one in each coordinate.

Examples from S365:

```text
k= 5: (2,5)->10, (3,3)->9       since 4 = 1*4 = 2*2
k=13: (2,13)->26, (3,7)->21, (4,5)->20
k=17: (2,17)->34, (3,9)->27, (5,5)->25
```

The universal family

```text
(1,1,...,1,2,k)
```

comes from the trivial divisor `1*(k-1)`.

## Representation-Theoretic Reading

There are two monoidal structures on the same basis of natural modes:

```text
e_x tensor_+ e_y -> e_{x+y}
e_x tensor_* e_y -> e_{xy}.
```

The additive shadow is the induced reachability preorder of `tensor_+`.  It is
too coarse: every lower mode points to every higher mode.  The multiplicative
shadow is the induced reachability preorder of `tensor_*`.  It is still
structured because the primes remain as irreducible generators.

The product-sum collision stack is where these two monoidal products agree
after arity matching:

```text
tensor_+^k(e_{a_1},...,e_{a_k})
=
tensor_*^k(e_{a_1},...,e_{a_k}).
```

In rapidity coordinates, multiplication is linear:

```text
log(product) = sum logs.
```

Addition is external to that linearization.  Product-sum collisions are
therefore places where an external additive fiber lands exactly on a rapidity
lattice point produced by the multiplicative monoid.

The operator "adjoin a 1" acts like a creation operator for additive slack:

```text
core c with defect D -> D applications of the 1-creator -> collision.
```

This is crude, but it feels like the right abstraction: a Fock-like space of
nonunit multiplicative cores, graded by defect, whose grade determines how
many unit particles are needed to collide with addition.

## Why This May Be Useful

1. The old incomplete-tournament language should be upgraded to labeled
   operation cospans.  Otherwise the additive graph collapses to the order and
   looks falsely empty.

2. Multiplication is the sparse surviving skeleton.  The divisor DAG is the
   arithmetic incomplete tournament, with prime Hasse edges as generators.

3. Sum-product equations are not side curiosities.  They are the collision
   strata between the dense additive completion and the sparse multiplicative
   skeleton.

4. The defect-core theorem turns the whole equation family into a core
   enumeration problem:

   ```text
   choose c_i >= 2 with product(c)-sum(c) = k-len(c).
   ```

5. The two-core layer is exactly the shifted divisor graph of `k-1`.  Higher
   cores should similarly encode higher multiplicative partitions with a
   defect constraint.

## Next Questions

- Extend the target collision spectra past arity `20`.  Do highly composite
  numbers have unusually rich spectra after normalizing for size, or does the
  defect filter select a different sequence?

- Study the graph whose vertices are nonunit cores and whose edges append or
  split a factor while tracking defect.  Does it have Zeckendorf-like normal
  forms or triangular phase changes?

- Complete the multiplicative divisor DAG by additive complement edges and ask
  which cycles appear first.  The additive order supplies every missing
  comparison; multiplication supplies the sparse irreducible skeleton.

- In finite fields, compare this natural-number collision stack with
  `sum_product_tournament.py`: QR sets are multiplicatively closed, intervals
  are additively coherent, and product-sum collisions may be the integer
  shadow of that tension.

- Try a formalization in the semiring category: natural modes as objects,
  addition and multiplication as two symmetric monoidal products, and
  product-sum tuples as equalizer fibers between the two arity-k functors.
