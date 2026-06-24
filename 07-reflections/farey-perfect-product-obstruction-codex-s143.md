# Farey Perfect Products and the Obstruction Guardrail

Session: codex-2026-06-24-S143

The useful surprise is that the perfect-number material does not sit outside
the current Farey thread.  It lands exactly on the `n=2` unit-excess chain:

```text
2^(p-1)/(2^p-1),      product = 2^(p-1)(2^p-1).
```

So `2/3` in `F_3` is not only the first non-star product `6`; it is the first
even perfect product, with graph shadow `K_{2,3}`.  Then `F_4` adds `3/4`,
whose product `12` is not perfect but whose graph shadow `K_{3,4}` is the first
reduced complete-bipartite Farey term containing a `K_{3,3}` minor.  That makes
the `F_3/F_4` seam genuinely useful:

```text
F_3: perfect product / planar two-block.
F_4: first reduced K33 incidence wall.
```

The formula that keeps the analogy honest is:

```text
sigma(a(na-1))/(a(na-1)) = n(2a-1)/(na-1)
```

when `a=2^k` and `na-1` is prime.  For `n=2`, this is exactly `2`.  For
`n=14`, it is `2 - 12/(14a-1)`.  Thus LRC14 is not "explained by perfect
numbers"; it is a deficient parallel chain whose `n=2` sibling gives a clean
fixed-point control.

The graph side needed the Kuratowski/Wagner correction.  Edge count alone is
too lossy: `K5` has `10` edges, but `2/5 -> K_{2,5}` is planar with the same
edge count; `K3,3` has `9` edges, but `3/3` is not reduced.  Disjoint unions
and density mediants are bookkeeping, not new forbidden minors.  The structural
operation is minor/subdivision closure.

The Tournament Analysis result is the practical warning.  The carrier-role
majority tournament has a nontrivial SCC:

```text
{K33_incidence, farey_level, product_edges, unit_excess_chain}.
```

That is the layer where the proof can go wrong if one scalar is allowed to
stand for the whole packet.  The exact `M`/Farey branch still comes first; the
product ledger is useful only after the incidence/minor label and the
unit-excess address stay attached.

Next useful move: annotate the existing LRC14 row bank with a small
`unit_excess_abundancy_shadow` column for `p/(14p-1)`, then compare it to the
C27/K33 split from HYP-2937/HYP-2940.  A null result would still be useful
because it would certify the perfect-number bridge as a control metaphor rather
than a proof engine.
