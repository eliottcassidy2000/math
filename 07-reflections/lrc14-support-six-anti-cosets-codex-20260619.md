# LRC14 support-six anti-cosets

Codex 2026-06-19.

The S10 frontier asked for "the Minkowski count of the support-6 relation."
Executing the first finite version changed the shape of the missing lemma.  The
object is not a free product tail.  It is a deleted anti-coset theta sum:

```text
Lambda(E) = { n in Z^(k-1) : sum_i n_i e_i = 0 }
```

with all coordinate hyperplanes and all nonzero 7-cosets deleted, then with the
THM-538 support floor `support(n) >= 6`.  This quotient preserves the exact
predicate `K(n) != 0`; it destroys the x-space component geometry, which is why
single-stranger Weyl and cluster-collapse remain parallel proof routes.

The first count is in
`04-computation/lrc14_support6_minkowski_count_codex_20260619.py`, with output
`05-knowledge/results/lrc14_support6_minkowski_count_codex_20260619.out`.

## What moved

The naive span-only Minkowski hope is false in the strongest possible finite
sense.  Wide shapes can have height-1 large-involving support-six relations.

```text
k=8:  E=(0,1,2,3,4,5,6,21)
      h_typeII=1 because 21 = 1+2+3+4+5+6.

k=10: E=(0,1,2,3,4,5,6,7,8,22)
      h_typeII=1 via a signed six/seven-term identity.

shifted cluster: E=(0,50,51,52,53,54,55,56)
      h_typeII=1 everywhere; span decay is the wrong language.
```

But dissociation behaves exactly as the contraction story predicts:

```text
E=(0,1,2,3,4,5,6,97): first large-involving shell jumps to h=5.
E=(0,1,2,3,4,5,6,7,68): first large-involving shell is h=3.
```

The unsigned `c1` envelope is also dramatically too loose.  For the k=10 worst
verifier at box `H=4`, the shell envelope is already `45863`, while the exact
`|K|` through height 3 is `0.61553` and the signed shell sum through height 3 is
near zero (`-0.002338 + 0.002030 + 0.000047`).  A proof that only counts
`64*c1^s/prod |n_i|` is paying for a ghost.  The sign and the deleted 7-coset
structure are doing real work.

## New split

The support-six tail should be split into three regimes.

1. Low-height anti-coset resonances.  These are finite additive-energy
   identities such as `N = sum six core offsets`.  They are not controlled by
   span, but they are enumerable and empirically safe.

2. Dissociated deleted-coset theta tail.  Once low-height anti-cosets are
   excluded, the first large-involving height rises quickly.  This is the true
   successive-minima/Minkowski-second-theorem target.

3. No-scale clusters.  These have many small relations at every height, but the
   exact sector measure is small.  They should be controlled by a collapse
   quotient/common-translation argument, not by span decay.

## Tournament note

I used proof obligations as tournament vertices rather than runners or arcs.
The pairwise observable was proof-closure strength across five traits:
uniformity, computability, signedness, multi-large coverage, and compatibility
with the finite bounded-spread check.  The resulting tournament was transitive:

```text
anti_coset_shell_count
-> bounded_finite_check
-> successive_minima_theta
-> cluster_collapse
-> single_stranger_weyl
-> free_coordinate_envelope
```

This is a useful warning.  Speed-vertex tournaments hide the support floor.  The
meaningful tournament for this subproblem lives one level up, on proof carriers
and anti-coset classes.

## Status

LRC(14) is still not proved.  The session did not execute the infinite
Minkowski bound.  It made the missing bound more precise and ruled out the
over-simple version: "large span forces large coefficients" is false without a
low-height anti-coset ledger.
