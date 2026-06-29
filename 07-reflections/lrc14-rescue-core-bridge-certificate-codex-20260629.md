# LRC14 Rescue-Core Bridge Certificate Reflection

HYP-3439 was useful because it changed the meaning of the HYP-3437 rank-`6`
overlap core.  Before the bridge, rank `6` looked like possible high-rank
noncanonical debt.  After joining HYP-3437 with HYP-3438 and HYP-3450/HYP-3451,
the rank-`6` signal on the AP/84m spine is exactly the canonical `m=1`
corridor-fence row, duplicated by name as `covering_AP_with_84` and
`ap_omit_12_tail_84x01`.

The AP tail now has a cleaner finite target.  Rows `ap_omit_12_tail_84x02`
through `x12` still have negative one-branch naive slack, but the minimum
overlap rescue core is always `(5,7,9,11,13)` of rank `5`, and the
component-cover ledger still leaves low-rank two-colour escapes.  So the next
proof should not chase raw rescue rank.  It should prove:

```text
canonical m=1 rank-6 base case
rank-5 AP-tail descent for m>=2
component-cover escape obstruction for arbitrary primitive rows
```

After HYP-3452/HYP-3454/HYP-3456/HYP-3457, the AP-tail clause is sharper than
the original HYP-3439 bridge.  HYP-3457 closes `m=1..4` as finite mixed
transients, HYP-3454 gives the rank-one `E:84m/E:84m` endpoint component for
`m>=5`, paired dead-cover rank is already `<=2` from `m>=3`, and HYP-3456
derives the mod-35 escape count.  Those sidecars should be used as the AP-tail
proof skeleton rather than re-scanning the bridge rows.

The nonnegative control `multi_far_84_154` also matters.  It shows the bridge
can route a row with no one-branch deficit while the component-cover ledger
still emits many low-rank escapes.  That makes raw component count and raw
dead-fraction unsafe as theorem carriers unless paired with the exact
one-branch slack state.

Assumption challenge: runner vertices, odd blockers, even-half gates, raw gaps,
fixed sections, wall crossings, residues, cover arcs, Fourier modes, matroid
circuits, survivor gates, component graph nodes, endpoint walls, and proof
obligations were all plausible tournament vertices.  The chosen carrier is the
bridge obligation because it preserves the actual LRC predicate being routed:
negative one-branch slack plus legal two-colour escape.  It destroys raw runner
order and most interval geometry, so endpoint labels, survivor-gate words, and
component addresses must remain sidecars.

HYP-3462 now formalizes the AP84 carrier/splice: HYP-3431 supplies the
branch-union carrier, `m=1` is the rank-`6` base, `m>=2` is the checked rank-`5`
AP-tail descent, and HYP-3454/HYP-3456/HYP-3457 supply the endpoint, floor, and
finite packets.  The next pull is broadening the bridge to non-AP primitive
covering rows through HYP-3453/HYP-3451 and the HYP-3455 gluing clause.
