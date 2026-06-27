# LRC14 Farey-Perfect Kuratowski Carrier Reflection

S144 separates three motifs that are easy to collapse:

```text
perfect-number product,
Farey mediant,
forbidden graph minor.
```

The exact split is useful.  `F_3` already has the perfect product `2*3=6`
through `2/3 -> K_{2,3}`, but that graph is planar.  `F_4` first has the
complete-bipartite Kuratowski wall through `3/4 -> K_{3,4}`, but its product is
`12`, not perfect.  So product perfection and nonplanarity are not the same
signal.

The even-perfect lane is still worth keeping because it gives a clean edge-load
family:

```text
2^(r-1)/(2^r-1) -> K_{2^(r-1), 2^r-1}.
```

But after the `2/3` seed, every such graph is nonplanar only by inherited
`K_{3,3}` containment.  That makes the lane a stress test for product-ledger
scalarization, not a source of new obstruction types.

The most important proof guardrail is Kuratowski/Wagner minimality.  `K5` and
`K3,3` are the two cores.  Their disjoint union is not a third core, and a
mediant is not an average graph.  The proof-relevant iteration is minor or
subdivision transitivity.

For LRC14, this reinforces the existing route:

```text
1/14 -> 2/27 -> 3/41
```

with `2/27` as the planar C27 two-block/petal branch and `3/41` as the first
unit-excess `K_{3,3}` incidence packet.  The next useful test is to attach the
mediant cross terms `ad` and `bc` to the C27 shell-transfer labels and ask
whether every remaining low-gap non-AP/GW atom either stays in the `p=2` petal
branch or crosses the `p>=3` minor-transitive packet needed by HYP-2908/THM-572.
