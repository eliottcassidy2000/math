---
id: THM-380
name: lrc-pressure-certificate-trilemma
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S503
depends_on:
  - THM-357
  - THM-359
  - THM-366
  - THM-369
  - THM-379
---

# THM-380: LRC counterexamples must pass sieve, core, and pressure-cycle tests

## Statement

Let `V={v_1,...,v_{n-1}}` be a Lonely Runner speed set at threshold `1/n`.
Let `F(V)` be the open forbidden union from THM-357, and let `(E_*,I_*)` be
the terminal endpoint/interval core from THM-359.

Let `G` be any directed pressure graph on the speeds in `V`.  Suppose that, on
the active owners of `(E_*,I_*)`, the graph `G` realizes the core protection
incidences in the sense of THM-379: for each core endpoint owner, at least one
strict core protector gives an edge

```text
protector speed -> endpoint-owner speed
```

in `G`, with no self-protection.

If `V` is a full-open-cover counterexample, then all three conditions hold:

1. `V` is small-denominator sieve-complete:

```text
for every m in {2,...,n}, some v in V is divisible by m;
```

2. the terminal protection core `(E_*,I_*)` is nonempty;
3. the active pressure graph `G` contains a directed cycle, equivalently a
   nontrivial strongly connected component.

Equivalently, any one of the following is a certificate that `V` is not a
full-open-cover counterexample:

```text
missing small denominator;
empty terminal endpoint core;
pressure-realized terminal core + pressure DAG.
```

## Proof

Assume first that `V` misses a denominator `m` with `2 <= m <= n`: no speed is
divisible by `m`.  By THM-366, and by the Lean-formalized master sieve
THM-369, every unit point `a/m` is a closed lonely witness; if `m<n`, it is an
open lonely witness.  Therefore `F(V)` cannot be the whole circle.  Thus a
full-open-cover counterexample must be sieve-complete.

Assume next that the terminal protection core is empty.  If `F(V)` were a
full open cover, then by THM-357 every forbidden endpoint would be strictly
protected.  The whole endpoint/interval system would then be a nonempty
protection core: every endpoint has a protector, and every interval boundary
endpoint lies in the endpoint set.  This contradicts THM-359, which says the
terminal core contains every core.  Hence a full-open-cover counterexample
must have nonempty terminal core.

Finally assume the terminal core is nonempty and pressure-realized by `G`.
THM-379 applies to the active owner set of the terminal core.  It says the
owner-protection graph contains a directed cycle, and any pressure graph
containing those realized edges contains a directed cycle as well.  Therefore
`G` has a nontrivial strongly connected component.

The three contrapositive certificates are exactly the negations of these
necessary conditions.

## Interpretation

This theorem packages the current LRC proof workflow into a finite certificate
ladder:

```text
denominator sieve -> endpoint core -> pressure SCC
```

The first two layers are arithmetic/topological.  The third layer is
Tournament Analysis: if the pressure lift is a DAG and it captures the
endpoint-core protection incidences, the supposed counterexample core cannot
exist.

The theorem does not claim that the current `k1`, `k2`, or deficit pressure
lifts always realize every endpoint-core incidence.  That realization is the
next proof obligation.  What is now formal is the payoff:

```text
realized core + DAG => no counterexample core.
```

## Verification Record

`04-computation/lrc_endpoint_pressure_formal_s503.py` prints the trilemma and
checks the finite graph engine behind the pressure-cycle layer.  It is not a
new LRC search; it is a sanity audit for the formal certificate shape.

Stored output:

```text
05-knowledge/results/lrc_endpoint_pressure_formal_s503.out
```

## Related

- THM-357: Lonely Runner endpoint-protection trichotomy.
- THM-359: endpoint/interval protection core peeling.
- THM-366: LRC small-denominator divisibility sieve.
- THM-369: Lean formalization of the denominator sieve.
- THM-377: selected n=14/n=18 pressure acyclicity.
- THM-379: owner-realized endpoint cores force pressure cycles.
- HYP-1960 and HYP-1961.
