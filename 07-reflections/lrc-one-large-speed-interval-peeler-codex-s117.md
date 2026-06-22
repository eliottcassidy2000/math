# LRC one-large-speed peeling is an interval-width theorem

*codex-2026-06-22-S117*

The scale-separated branch does not need full equidistribution when there is
only one genuinely large speed.  It only needs one connected safe interval from
the smaller seed.

Let `B` be a seed with max speed `m`.  If some `tau` satisfies

```text
||b tau|| >= alpha > 1/n
```

for every `b in B`, then the seed remains threshold-`1/n` safe on an interval
of length

```text
2(alpha-1/n)/m.
```

The added speed `v` has danger teeth of length `2/(n v)`.  If the seed-safe
interval is longer than one tooth, it cannot be fully swallowed by the danger
comb.  Hence `B union {v}` is safe whenever

```text
v > m / (n*(alpha-1/n)).
```

Using `alpha=1/(n-1)` from the smaller LRC theorem gives the clean induction
gate

```text
v > (n-1)m.
```

For LRC14 this says: if the largest speed exceeds `13` times the second
largest, then the row is safe by LRC13.  Any counterexample must be
top-balanced:

```text
v_max <= 13 v_second.
```

The AP-core seed `{1,...,11,13}` is better because `tau=1/12` gives
`alpha=1/12`; the local gate is then `v>78`, so the committed lcm/radical
speeds are certified by a one-arc existence argument.  This is sharper than
the HYP-2904 global component-budget threshold `v>=768`, because HYP-2904
proves a positive measure lower bound while this lemma proves only existence.

Connection to the current proof split:

- HYP-2906 closes the one-large, high-ratio Node-3 branch before p0-cap or
  Bonferroni obligations enter.
- HYP-2904 remains the right object for density floors, multiple large speeds,
  and finite-Part-A budgets, because those need component or arc-count control.
- HYP-2903's missing-depth parity guard is now cleanly localized to the
  balanced/binding cap leg; a locally peelable far speed should be removed
  before attempting to prove a third-order Bonferroni tail sign.

Tournament Analysis uses proof carriers as vertices:
`connected_seed_safe_arc`, `single_comb_tooth_width`, `smaller_LRC_margin`,
`scale_ratio_gate`, `global_component_budget`, `runner_count_only_induction`,
and `raw_runner_vertices`.  The challenged assumption is that the scale
reduction must retain the whole safe set.  For existence, a single interval is
the right object.
