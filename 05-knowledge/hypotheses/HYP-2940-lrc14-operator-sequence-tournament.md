---
id: HYP-2940
title: LRC14 operator sequence tournament
status: PROOF-INTERFACE / synthesis warning; not a proof
source: codex-2026-06-23-S137
related:
  - HYP-2939
  - HYP-2938
  - HYP-2937
  - HYP-2936
  - HYP-2935
  - HYP-2934
  - HYP-2933
  - HYP-2932
  - HYP-2931
  - HYP-2930
  - HYP-2908
  - THM-572
---

# HYP-2940: LRC14 operator sequence tournament

S137 is a synthesis pass over S130-S136.  It keeps the prompt's four Farey
mutations

```text
p+q, p*q, q^p, p^q
```

but treats them as labelled operator sequences before using them as proof
carriers.  The live chain is the LRC14 unit-excess chain

```text
p/q = p/(14p-1),      M - 1/14 = 1/(14q).
```

## Computation

The script
`04-computation/lrc14_operator_sequence_tournament_codex_s137.py` stores output
at
`05-knowledge/results/lrc14_operator_sequence_tournament_codex_s137.out`.

It reuses the exact `M(S)` engine and graph routines, avoiding optional
network dependencies.  The default row bank is the S130 `749`-row bank through
replacement `70`; HYP-2937 remains the larger `140` frontier audit.

## Sequence readout

Along `p/(14p-1)`:

```text
q        = 14p - 1       linear, Delta=14
p+q      = 15p - 1       linear, Delta=15
p*q      = 14p^2 - p     quadratic, Delta2=28
```

So the `n+2` recursion lane is additive and theorem-facing, while the `n*2`
or product lane is a two-dimensional coimage/area side channel.  The power
payloads are deliberately non-polynomial magnitude amplifiers and remain
stress tests for magnitude-blind quotients, not proof denominators.

The first three unit-excess rows keep their established routing:

```text
1/13  -> star / q-parent
2/27  -> C=27 two-block shell branch
3/41  -> first K33 wall
```

## Row-bank readout

On the S130 row bank through replacement `70`:

```text
nonunit-excess           690
q-parent-star             54
tight-floor                2
C27-petal-two-block        2
K33-unit-excess            1
```

The low frontier is unchanged:

```text
M <= 3/41: AP, GW, near-miss 12->36
M <= 2/27: AP, GW, near-miss 12->36, swap 10->20, swap 13->26
```

The C27 transfers are:

```text
AP:       perfect transversal
GW:       H[12:g3] -> D[3:g3]
10->20:  H[10:g1] -> D[7:g1]
13->26:  H[13:g1] -> D[1:g1]
12->36:  H[12:g3] -> D[9:g9]
```

Transfer labels recur in safely loose rows, so exact `M`/Farey branch remains
the first coordinate.  The C27 quotient is a router, not a classifier.

## Graph and PZ carriers

S137 confirms the S132 carrier roles:

```text
octahedron L(K4): v=6, e=12, cycle_rank=7
Clebsch folded Q5: v=16, e=40, triangle-free, cycle_rank=25
halved 5-cube: v=16, e=80, complement(Clebsch), cycle_rank=65
```

The toy six-sector Paley-Zygmund lower bound loses about `0.162` at `k=8`,
`0.144` at `k=10`, and `0.114` at `k=12` against the exact independent union.
Thus PZ remains an existence gateway; tight cap comparison still needs labelled
degree-4/factorial moment data.

## Tournament Analysis

Vertices are operator/carrier roles:

```text
q, p+q, p*q, C27 shell, Kpq/K33 incidence,
octahedron L(K4), Clebsch/half-cube, Paley-Zygmund, power payloads.
```

Two gauges were tested.

The conservative proof gauge is transitive:

```text
q binding scale
> C=27 shell
> p+q additive lane
> Kpq/K33 incidence
> p*q product lane
> octahedron L(K4)
> Clebsch/half-cube
> Paley-Zygmund
> power payloads
```

Fingerprint:

```text
c3=0, SCCs all singleton, Hamiltonian paths=1
```

The majority-of-criteria gauge uses:

```text
theorem scale, sequence simplicity, branch separation,
state-lift fit, graph-packet fit, anti-scalar guard,
magnitude stress.
```

It is not transitive:

```text
directed 3-cycles = 2
nontrivial SCC = {p+q, p*q, octahedron, Clebsch/half-cube}
Hamiltonian paths = 5
```

The two directed cycles are:

```text
p+q -> p*q -> Clebsch/half-cube -> p+q
p+q -> octahedron L(K4) -> Clebsch/half-cube -> p+q
```

This is the useful warning.  Once sequence simplicity, product/coimage growth,
graph-packet fit, and magnitude stress are all allowed to vote, the additive
and packet layers form a real SCC.  No single scalar carrier should be trusted
inside that SCC without its side-channel labels.

## Proof target

The current best interface becomes:

```text
1. Keep exact M/Farey branch first.
2. Use C=27 shell transfer for p=2.
3. Use Kpq/K33 incidence for p>=3.
4. Attach octahedral/Clebsch packet data if a support-six or folded-mask
   residual remains.
5. Use Paley-Zygmund only as an existence gateway.
6. End, if possible, in the HYP-2908 / THM-572 forbidden-H state lift.
```

Candidate lemma shape:

```text
After the standard LRC14 finite-core reductions, every low-gap non-AP/GW atom
has either a unit-visible C27 hole discharged by petal/two-block rigidity, or a
sign-visible K33/octahedral/Clebsch packet whose connected OCF state-lift would
have activity-two value 7.
```

This is not a proof of LRC14.  It narrows the packet that a proof must
construct before invoking the forbidden-H endpoint.
