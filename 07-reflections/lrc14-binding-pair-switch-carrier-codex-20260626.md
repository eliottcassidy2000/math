# LRC14 Binding-Pair Switch Carrier

THM-524 is stronger than a max-finding trick but weaker than a scalar proof
route.  It says the maximum lives at a pair crossing.  It does not say that the
pair crossing alone carries the LRC predicate.

The exact audit in `lrc14_binding_pair_switch_carrier_codex_20260626.py` makes
the failure mode visible.  AP13 and GW are zero-open boundary atoms with the
literal complement optimum pairs `(1,13)`, `(5,9)`, `(3,11)`.  Petal and K33
rows move to positive-open active pairs such as `(7,20)`, `(1,26)`, `(5,36)`,
and `(1,36)`.  The covering row `{1..11,13,84}` is the off-grid-only witness
case, with `M=7/89` and no denominator-14 grid witness.

The raw-pair quotient is badly unsafe.  Large-tail rows have hundreds or
thousands of pair-good decoy times: the pair can sit at distance `1/2` while a
third runner collapses the true minimum.  So the actual proof carrier is:

```text
pair crossing + denominator lane + active blockers + all other-runner clearance
```

This suggests a practical next move for the LRC14 packet classifier.  Add
binding-switch fields to each HYP-2963 packet, then test whether they are
constant on families or whether they are the exact missing coordinate behind a
covering, petal, K33/state-lift, or F7 route.
