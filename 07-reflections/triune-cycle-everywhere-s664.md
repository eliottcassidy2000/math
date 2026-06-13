# Triune Cycle Everywhere (S664)

The useful version of the cycle is not that every domain has a sum, a product,
and a fraction.  That is too cheap.  The useful version is a projection-repair
grammar:

```text
sum -> product -> fraction -> sum
```

where:

```text
sum      = additive aggregation: packets, traces, walls, counts, moments
product  = local factorization: gcd/norm shells, components, sieves
fraction = boundary memory: carry, owner, branch, lift, deletion, recursion
```

The three arrows mean:

```text
sum -> product
  Additive packets expose or localize the factor/shell data.

product -> fraction
  Local factors are still a quotient.  They need a lift, branch, owner,
  deleted-card mark, generic, or recursion boundary.

fraction -> sum
  The boundary state re-expands into the next additive packet ledger.
```

That last arrow is the important one.  It prevents the trinity from becoming a
static checklist.  The fraction face is the generative seam: it tells us which
next sum is being summed.

## Computation

`04-computation/triune_cycle_everywhere_s664.py` scans the repo and builds a
curated domain atlas.  The scan is intentionally broad and noisy, but it gives
a useful sanity check: thousands of research files mention all three faces, and
the top files are exactly the navigation/index structures plus recent carrier
scripts.

The curated face-dependency tournament has vertices:

```text
sum, product, fraction, raw_scalar.
```

Pairwise observable:

```text
domain-weighted dependency support.
```

Switch:

```text
weighted majority, with tie path sum -> product -> fraction -> raw_scalar.
```

It recovers the desired nontransitive face core:

```text
score_hist={0: 1, 2: 3}
directed_3cycles=1 [('sum', 'product', 'fraction')]
sccs=[['fraction', 'product', 'sum'], ['raw_scalar']]
hamiltonian_paths=3
```

The strongest weights are:

```text
sum -> product       support 54
product -> fraction  support 52
fraction -> sum      support 51
```

All three beat `raw_scalar`.

## Strong Instances

### LRC14

S663 is the exact core.  AP, `Vstar`, and `2AP` share the same additive/product
shadow at `C=27`, but local `+27` carry perturbations become strict.  Grouping
by only the sum/product shadow leaves:

```text
3 mixed groups = 1 floor + 91 strict each.
```

Adding the carry-continuant fraction face leaves:

```text
0 mixed groups.
```

This is the cleanest current proof target:

```text
fixed odd-wall packets
+ fixed C=27 gcd/product shell
+ fixed carry-owner continuant
=> AP, Vstar, 2AP, or strict looseness.
```

### Tournament Decks

S660 is the exact sibling.  Full decks collide through `n=6`; global repairs
like `deck+H`, `deck+score`, and `deck+scalar` do not split the collisions.
The repair is a fraction face:

```text
(card isomorphism type, deleted vertex outdegree).
```

That is the same owner rule as LRC carry state.  The product-like deck quotient
must remember which missing vertex owned which card.

### OCF / H(T)

OCF is the algebraic hinge.  Additive odd-cycle packet coefficients already
become the product-like partition function:

```text
H(T) = I(Omega(T), 2).
```

But `I(Omega,2)` is still a quotient.  Deletion/contraction/substitution
boundaries should be carried as a continuant-like fraction state.  This is the
best next pure tournament experiment: encode a small macro-word as a continuant
and compare it with Hamiltonian-path DP.

### Unit Distance

The triune ledger should be:

```text
sum      = unit-edge spine / edge-packet counts / tile flips
product  = direction support, Eisenstein norm shells, unit-vector factors
fraction = point-deletion frontier owner or ear-attachment state
```

This suggests a new impairment test for `n=21/22`: deliberately keep the same
edge-count and direction/norm product shadow, then vary deletion owner state.
If the best witnesses split only after the owner is kept, unit distance has the
same projection-repair structure as LRC14.

## Number-Theory Side Labs

The S664 script includes three small toys.

For Goldbach/Lemoine through the chosen window, the same prime pair gives both
additive columns:

```text
E = p + q
O = p + 2q
```

and the fraction/branch reconstruction is:

```text
q = O - E
p = 2E - O.
```

The diagonal examples are the fixed-pair shadows:

```text
(p, 2p, 3p) = (3,6,9), (5,10,15), (7,14,21), ...
```

For the finite-field Vieta lab over `F_17`:

```text
sum_only_max_fiber=17
product_only_max_fiber=33
sum_product_max_fiber=2
sum_product_branch_max_fiber=1
```

So trace/norm reduces the fiber to the two root orders, and the branch sheet
splits it completely.  This is the finite pi/e trace-norm model.

For aliquot/perfect numbers through `300`:

```text
sigma product formula mismatches = []
perfect fixed points = [6, 28]
aliquot two-cycle = (220, 284)
```

The divisor-sum carrier has an Euler product face, but the aliquot orbit is the
fraction face.  Perfect numbers are fixed points of that orbit state.

## Assumption Challenge

Candidate Tournament Analysis vertices considered:

```text
faces, domains, files, runners, residues, wall pairs, gcd shells, carry words,
deck cards, deleted vertices, unit-distance points, prime pairs, models,
generics, proof obligations.
```

S664 uses faces for the cycle tournament and domains for the route tournament.
The preserved predicate is not raw truth of each problem; it is:

```text
which carrier repairs a projection collision?
```

The repo scan destroys local labels.  The curated entries restore the known
owner state from S660/S663 and mark weaker domains as finite-lab or transfer
targets rather than theorems.

## New Rule

When a scalar quotient almost works, ask:

```text
What is its fraction face?
```

In this repo that usually means:

```text
who owned the deleted object?
which lift of the residue shadow was chosen?
which branch of the trace/norm polynomial is active?
which generic/model boundary was retained?
which recursive route produced the visible packet?
```

That is the common move across LRC, decks, OCF, unit distance, Goldbach/Lemoine,
pi/e, perfect numbers, CH, cauldrons, and Heegner-style polynomial primes.
