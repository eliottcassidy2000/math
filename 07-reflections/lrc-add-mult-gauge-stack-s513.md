# LRC add/multiply gauge stack after S513

This pass makes the user's `x+2` and `x*2` picture precise enough to use in
Tournament Analysis.

The natural-number address is:

```text
N = 2^h * odd_core.
```

The `x+2` direction is horizontal motion through odd cores.  The `x*2`
direction is vertical dyadic motion.  The LRC hard rows `14` and `18` sit next
to each other on the first-even row:

```text
14 = 2 * 7,
18 = 2 * 9.
```

The pure dyadic test row is:

```text
16 = 2^4.
```

That already clarifies the recursive structure: row-parent and gate ladders
are not just bigger numbers.  They are vertical translations in a grid whose
odd core remembers which horizontal branch produced the obstruction.

Addition and multiplication play opposite roles.  Addition flattens.  As a
plain operation shadow, `x+y=z` only says `x<z`, so it collapses to a
transitive order.  The useful additive information lives in the fiber:

```text
even N: Goldbach pairs p + q = N,
odd N:  Levy/Lemoine pairs p + 2q = N.
```

Multiplication branches.  Divisors, largest proper divisors, dyadic height,
and odd core make a sparse skeleton inside the additive order.  Product-sum
equations are where the two operations collide:

```text
P = product(F)
S = sum(F)
D = P - S
D ones repair the equation: D + S = P.
```

For two nonunit factors at arity `k`, this becomes:

```text
(a-1)(b-1) = k-1.
```

That is why product-sum data belong in the LRC operation stack.  They are not
loneliness by themselves.  They are labels for places where multiplication
creates exactly the surplus that addition repairs.

The computational punchline is clean.  On the selected denominators, every
single scalar operation gauge becomes a transitive tournament:

```text
H=1, c3=0.
```

That includes additive abundance, additive scarcity, multiplicative branch
depth, product-sum collision, endpoint debt, A000568 odd survival, and a
scalar grid-stack score.  These are ledgers.  They order denominators, but
they do not yet contain the shape of loneliness.

The shape appears when the gauges disagree pairwise.  The majority tournament
combining additive abundance, divisor branching, dyadic height, product-sum
critical pairs, endpoint debt, and A000568 odd survival has:

```text
H=242407, c3=42, SCCs=(13,1).
```

The pressure-style majority gauge has:

```text
H=509, c3=15, SCCs=(9,1,1,1,1,1).
```

That is a useful distinction.  The loneliness metric is not "the Goldbach
count" or "the dyadic height" or "phi(N)".  It is the tournament shape created
when those coordinates vote against each other.

The A000568 analogy fits as the quotient layer.  Ordinary tournament
isomorphism classes forget labels through an odd-cycle Burnside filter.  LRC
cannot live only there, because the stationary observer and endpoint labels
matter.  But the odd-core/dyadic grid is the natural-number counterpart of
that same fixed-point survival: odd structure survives quotienting, while
dyadic motion tells how fibers lift and split.

So the practical metric object is:

```text
runner-level dynamic gauges
+ pair-cell danger and operation labels
+ denominator-level add/multiply/A000568 labels.
```

`H` should be read inside each layer, not as a universal scalar.  Score
histograms, directed triangles, SCCs, and edge flips tell whether a layer is
only a transitive ledger or a real switchboard.

The counterexample-shaped target becomes much narrower.  A bad row would need
low additive smoothing, high endpoint debt, product-sum/gate trouble, marked
A000568 projection trouble, and a nonpeeling pressure SCC at the same time.
That conjunction is the next object to try to make impossible.
