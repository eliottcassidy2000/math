# HYP-2226: Heegner Prime Polynomials Use a Tournament-Style `n-2` Witness Horizon

**Status:** OPEN structural bridge with exact finite scout and a classical
Rabinowitsch/Heegner side channel.  Refines HYP-2225 by adding the THM-410
long-edge interval-reversal model and the pair-witness split.

## Claim

The Euler/Rabinowitsch prime-generating primes

```text
{2,3,5,11,17,41}
```

and the Heegner class-number-one list

```text
{1,2,3,7,11,19,43,67,163}
```

should be connected to the tournament `n-2` witness layer through the typed
quadratic

```text
Q_p(x) = x^2 + x + p.
```

For the six Euler primes, the map

```text
d = 4p - 1
```

gives the Heegner tail

```text
{7,11,19,43,67,163}.
```

The base Heegner values `{1,2,3}` are not produced by this tail map; they are
the unit/Gaussian/Eisenstein base fields and should remain typed separately.

The useful tournament-shaped statement is:

```text
Q_p(0)       = p      boundary prime
Q_p(1..p-2) are prime interior witnesses
Q_p(p-1)    = p^2    forced boundary composite
```

Thus the `p-2` positive inputs are not an accident of indexing.  They are the
interior witnesses between two fixed boundary positions.

## Tournament Analogue

THM-410 gives the exact tournament model.  Start with the transitive
`p`-vertex tournament on `0<1<...<p-1` and reverse only the long edge
`0 -> p-1`.

Then every interior vertex `v` with

```text
0 < v < p-1
```

creates exactly one cyclic triangle

```text
0 -> v -> p-1 -> 0.
```

Therefore

```text
#C3 = p - 2.
```

Equivalently, for the reversed boundary edge `p-1 -> 0`, the pair-witness
decomposition is

```text
sigma = 0
lambda = p - 2
delta = 0
```

all witnesses are cyclic.  This matches the earlier pair-fiber identity
`sigma + lambda + delta = n-2`.

## Evidence From S650

`04-computation/heegner_prime_horizon_tournament_s650.py` checks primes
`p <= 500` and finds:

```text
Prime p <= 500 with first positive composite at x=p-1:
[2, 3, 5, 11, 17, 41]
```

For those six primes:

```text
p   d=4p-1   positive_horizon   p-2   first_bad   Q(first_bad)
2       7            0            0        1              4
3      11            1            1        2              9
5      19            3            3        4             25
11     43            9            9       10            121
17     67           15           15       16            289
41    163           39           39       40           1681
```

The same rows have interval-reversal triangle count `c3_long_edge=p-2`.

Non-Heegner controls break before the boundary.  For example:

```text
p=7:  Q_7(1)=9,  horizon=0 < p-2=5
p=13: Q_13(1)=15, horizon=0 < p-2=11
p=29: Q_29(2)=35, horizon=1 < p-2=27
```

## Proof Use

This should be used as a role-mismatch guardrail.

The Euler polynomial run length is not a Hamiltonian-path count.  It is a
boundary/interior witness horizon.  The tournament object it resembles is
THM-410's interval-reversal triangle count and the pair-witness fiber
`n-2`, not the global value `H(T)`.

The proof-route question is whether class-number-one factorization can be
modeled as a "no early interior collapse" condition, analogous to a tournament
interval whose only cyclic defects are the forced witnesses of one marked long
edge.

## Assumption Challenge

Do not choose tournament vertices only as numbers or primes.  Possible vertex
sets include:

- polynomial input slots `x`
- boundary positions `0,p-1`
- discriminants `4p-1`
- class group/factorization classes
- reversed intervals
- pair-witness roles `(sigma,lambda,delta)`
- proof obligations in the Rabinowitsch criterion

The quotient preserved here is the `p-2` interior-horizon predicate plus the
forced boundary collapse `Q_p(p-1)=p^2`.  It destroys the actual residue-level
primality mechanism unless the class-number-one side channel is retained.

## Tournament Analysis

S650 uses proof lenses as vertices.  The pairwise observable is:

```text
which lens preserves more typed information about the bridge?
```

The switch/gauge orients `A -> B` when `A` explains `B` without scalar
collapse.  The tournament is transitive with one Hamiltonian path:

```text
class_number_one_factorization
> euler_prime_horizon
> heegner_discriminant_map
> interval_reversal_long_edge
> pair_witness_n_minus_2
> first_composite_boundary
> raw_number_list
```

Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles = 0
scc_count = 7
hamiltonian_paths = 1
```

## See

`04-computation/heegner_prime_horizon_tournament_s650.py`;
`05-knowledge/results/heegner_prime_horizon_tournament_s650.out`;
`07-reflections/heegner-prime-horizon-tournament-s650.md`;
THM-410, HYP-2225, HYP-2224, HYP-2223, HYP-2222, HYP-1220, HYP-1202.
