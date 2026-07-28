---
title: "The same S3 plane in rainbow matchings and quartic Kummer theory"
date: 2026-07-28
status: CROSS-FRONTIER PROOF INTERFACE; no LRC or Keller realization follows
source: root-2026-07-28-c2-c3-v4-synthesis
---

# Three matchings and the `V4` character plane are one representation

The repaired THM-2648 matched-wall atlas contains a literal common object for
the prompt's binary/ternary and quartic-resolvent clues.  On its active
`K_(3,3)` block there are three one-factors

```text
M={m_0,m_+,m_-},
```

cyclically permuted by `C3`; reflection fixes `m_0` and swaps `m_+,m_-`.
Thus the local symmetry is `C3 semidirect C2=S3`.

Over `F2`, the permutation module has the canonical splitting

```text
F2^M = F2*1 direct-sum W,
W={x in F2^M: sum_m x_m=0}.                               (1)
```

The three nonzero elements of the two-plane `W` are
`110,101,011`, and `S3` permutes them faithfully.  Hence `W` is the standard
irreducible `F2[S3]` plane.  It is exactly the module

```text
Hom(V4,C2)                                                (2)
```

appearing in the proved quartic Kummer gate THM-2655.  This is a genuine
identification of representations, not merely a cardinality analogy.

## The two incarnations keep different data

In the rainbow problem, a vector of `W` is a parity difference between chart
**occurrences**.  THM-2648 proves that after retaining the occurrence bit,
both nontrivial `C2` chart parity and every nonzero `C13` carry character have
positive energy.  Carry projection alone forgets which edge occurrence
belonged to which chart; even equal hole sets and reflection symmetry do not
recover it.

In the quartic problem, `W` is the character plane of the semiregular `V4`
kernel.  THM-2655 places it in

```text
0 -> A^*/A^(*)2 -> H1_et(U,mu2) -> Cl(A)[2] -> 0.         (3)
```

Because `W` is projective over `F2[C3]` and `F2[S3]`, a copy in either end of
(3) is exactly sufficient for an **abstract** equivariant connected `V4`
Kummer carrier.  It does not make that carrier the graph quartic of a Keller
map.

So the common representation exposes a shared question but not a transfer:

```text
Does the standard S3 plane lift through the realization functor?           (4)
```

For LRC the functor is `physical same-base packet -> occurrence-labelled
matching kernel`.  For JC it is `polynomial Keller map -> graph-quartic
Kummer torsor on the resolvent normalization`.

## The exact binary/ternary fork

THM-2648 now separates two optima.

- A binary alternating four-cycle between two nonlinear rainbows gives the
  unrestricted sharp union of thirteen edges.  It has no distinguished
  affine chart.
- A ternary alternating six-cycle is the sharp repair after the affine chart
  is fixed.  Its union has fourteen edges, and adding the inverse repair
  completes the three one-factors carrying (1).

Thus `2` is cheaper for pure thinning, while `3` is the first move that keeps
the affine reference and exposes the `S3/V4` module.  This is a precise form
of “three old edges for three new edges repair a two-hole defect.”

The local action is still only `S3`.  It factors through

```text
C2*C3 -> S3
```

and imposes the extra Coxeter relation.  THM-2646 shows what is missing from
the finite quotient: the integral central/full-twist coordinate of `B3`.
Therefore the local matching frame does not by itself generate the Farey or
modular tower.  An iterated path/clock sidecar is required.

## Typed transfer table

```text
source:       three labelled matching occurrences / quartic V4 torsor
common map:   take the F2 augmentation/character plane W
preserved:    S3 action, its three nonzero vectors, C2/C3 symmetry
destroyed:    edge occurrence and positivity / polynomial affine realization
needed LRC:   lawful same-base function-valued kernel with endpoint labels
needed JC:    graph-quartic/Jelonek realization of the Kummer carrier
cheap test:   compute the W-valued packet before carry/discriminant projection
```

The warning is as useful as the bridge: equality of representation types is
not equality of carriers.  The next computation on either frontier should
test the lift in (4), not recompute an `S3` character table already known on
both sides.

