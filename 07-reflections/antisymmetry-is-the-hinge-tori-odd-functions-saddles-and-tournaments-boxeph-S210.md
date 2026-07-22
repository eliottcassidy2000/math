# Antisymmetry is one coordinate, not one theorem

*Scope-repaired audit of boxeph-2026-07-21-S210 / HYP-8835.*

> **Correction:** MISTAKE-225 retracts the claims “pure saddle iff
> transitive,” “intransitivity is toroidal recurrence,” and “bagel deficit is
> reduced Euler characteristic.” It also separates three unrelated sign
> actions. The useful survivor is an equilibrium-support/positive-kernel lens.

## The odd-support theorem uses two coordinates

For a tournament adjacency matrix `A`, the zero-sum payoff

```text
M=A-A^T
```

is real skew-symmetric. Hence its rank is even and every odd-order principal
block is singular. That observation alone does **not** prove that an optimal
support has odd size: even-order skew matrices can also be singular.

The missing coordinate is tournament-specific. Modulo two, every principal
tournament payoff block of size `m` is `I+J`: it has zero diagonal and one in
every off-diagonal position. Its determinant is `1+m` modulo two. For even
`m` this is one, so the block is nonsingular. An optimal support block must be
singular; therefore its size is odd. This is the parity mechanism behind the
Fisher--Ryan tournament-game result, not “antisymmetry implies odd support.”

The script exhausts every tournament through `n=5`, asserts that exactly one
candidate optimal support is found, and sees support sizes `1,3,5` only. This
is a finite control, not a replacement for the published theorem.

## Pure means Condorcet, not transitive

A pure optimal strategy exists exactly when some vertex beats every other
vertex—a Condorcet winner. The rest of the tournament may still contain
cycles. The smallest decisive witness has four vertices:

```text
0 beats 1,2,3;       1 -> 2 -> 3 -> 1.
```

It is intransitive, but `e_0` is its unique optimal strategy. Under replicator
dynamics, every interior state satisfies

```text
x_0' = x_0(1-x_0),
```

so it flows to `e_0`. This single example refutes all three implications

```text
intransitive => no pure optimum,
intransitive => support at least 3,
intransitive => recurrent/toroidal dynamics.
```

“Game saddle,” “Morse saddle,” and “hyperbolic dynamical saddle” are different
definitions. A shared word is not a map.

## The correct dynamical hinge

For replicator dynamics of a skew payoff,

```text
x_i' = x_i (Mx)_i,
```

let `p` be a kernel vector on an invariant support face. Then

```text
d/dt sum_i p_i log(x_i) = p^T Mx = -(Mp)^T x = 0.
```

When `p>0`, this is the monomial first integral

```text
product_i x_i^(p_i).
```

For the regular three-cycle, `p=(1/3,1/3,1/3)` and `x_0x_1x_2` is conserved.
The positive level sets in the two-dimensional simplex are periodic circles
`S^1`. Their union is an annular region, not a closed surface `T^2`, and no
general intransitive tournament is thereby integrable.

Why did the original derivative happen to vanish? The regular three-cycle has
zero column sums. A generic skew matrix does not. The transferable object is a
positive kernel/equilibrium support, not antisymmetry or intransitivity alone.

This suggests a sound tournament program: stratify by support faces carrying
positive kernel vectors, record which faces survive deletion/join/SCC
operations, and ask which monomial integrals persist. That is more informative
than assigning every cycle a torus.

## Solid torus versus boundary torus

The cutting “bagel” counted by A003600 is a **solid torus** `D^2 x S^1`. Its
boundary is `T^2`; these have different first Betti numbers:

```text
solid torus:  b=(1,1),
boundary T^2: b=(1,2,1).
```

A perfect Morse function on the closed boundary `T^2` has one minimum, two
index-one saddles, and one maximum. Importing that `1,2,1` count into a solid-
region cutting problem needs a cell map that S210 did not supply.

The exact sequence identity

```text
bagel(n)-cake(n)=T_n-1
```

also varies with `n`, whereas `chi(T^2)=0` and the reduced Euler characteristic
of a connected `chi=0` space is `-1`. MISTAKE-222 already explains why the
matching minus one is not a geometric bridge. The identity remains valuable
as a binomial operation law, not as an Euler theorem.

## Three involutions, three proof obligations

The original synthesis silently identified:

1. tournament complement/transpose, sending `M` to `-M`;
2. torus inversion, sending `theta` to `-theta` and fixing `2^n` torsion points;
3. an odd coordinate permutation, changing the sign of the Vandermonde.

These are abstractly similar `Z/2` actions but live on different objects.
Calling them “one oddness” requires an equivariant map and a preserved target
predicate. None is currently known. Moreover the LRC far-set Fourier weights
are **even** in frequency. The periodic Bernoulli sawtooth is anti-invariant
only with its value at integer discontinuities specified consistently.

THM-473 also concerns an ensemble-average characteristic polynomial related
to Hermite polynomials; it does not say each tournament matrix has Hermite
characteristic polynomial. Individual real skew matrices merely have the
usual imaginary-pair spectral symmetry.

## What remains worth pursuing

The corrected cross-domain questions are now precise:

- Can positive-kernel support faces be made into an operation-compatible
  tournament invariant finer than SCC data?
- Which deletion or join operations preserve the monomial first integral?
- Is there an explicit equivariant map from an LRC signed wall complex to a
  tournament sign representation, and which LRC predicate would it preserve?
- Can a genuine cell valuation explain the solid-torus cutting correction,
  or does the boundary/solid distinction provide a no-go theorem?

Primary tournament-game reference: [Fisher--Ryan, *Tournament games and
positive tournaments*](https://doi.org/10.1002/jgt.3190190208). The solid
object is explicit in [OEIS A003600](https://oeis.org/A003600).

Artifacts: corrected HYP-8835; MISTAKE-225;
`04-computation/tori_odd_saddles_tournaments_boxeph_S210.py` and its frozen
output.
