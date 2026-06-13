# LRC operation arc zoo after S511c

This pass was useful because it did not try to collapse everything into one
new arc rule.

The runner-level criteria split into three families.

First, scalar rankers are ledgers.  `origin_clearance_rank`,
`endpoint_deficit_rank`, and `two_neighbor_guard` usually complete to
transitive tournaments.  That makes their `H` uninteresting, but the score
sequence and marked observer score still say something concrete about where
the stationary runner sits and how the local two-neighbor bracket is arranged.

Second, phase-overridden arithmetic gates are perturbations.  `product_sum_gate`
and `goldbach_lemoine_gate` only change the half-turn tournament when the
runner-pair residues hit an operation relation.  On the hard `n=14` best
witness times they often agree with phase, which is itself information: the
current lonely corridor is not being cut by those prime/product gates.  On the
initial row they do perturb H and score histograms, so they still label the
arithmetic chamber.

Third, the bundled vote is a marked fiber coordinate.  It mixes phase,
endpoint, two-neighbor, dyadic, odd-core, Goldbach/Lemoine, and product-sum
votes.  Its H is far smaller than phase H on `n=14`, but that does not make it
worse.  It is measuring the agreement of many switchboards, not circular
spread alone.  The relevant readings are strict SCC size, score width, observer
score, and directed triangles.

The arithmetic report clarified the user-supplied picture:

```text
even N:  Goldbach gate       N = p + q
odd N:   Lemoine gate        N = p + 2q
all N:   multiplicative row  N = 2^h * odd_core
all N:   product-sum gate    (a-1)(b-1) = N-1
```

Addition is horizontal movement: fixed sums, adjacent odd rows, and `x+2`.
Multiplication is vertical movement: dyadic height, divisibility, and `x*2`.
The product-sum equation is the interface between them.

That makes the A000568 analogy cleaner.  A000568 counts unmarked complete
binary relations.  Tournament Analysis is the act of choosing the binary
relation from a richer pairwise metric.  LRC is not the unmarked quotient alone;
it is a marked switchboard fiber above it.  The observer label, endpoint
threshold, close-pair mass, two-neighbor data, and operation grid coordinates
are all part of the same object.

The next useful experiment is not another static operation label.  It should
combine this zoo with HYP-1980's danger weighting and HYP-1981's observer-source
target: for each exact endpoint chamber, emit phase H, source status, operation
bundle score histogram, close-pair tie count, and the product-sum/prime gate
activations.  Then ask whether mixed A000568 fibers split by those labels.
