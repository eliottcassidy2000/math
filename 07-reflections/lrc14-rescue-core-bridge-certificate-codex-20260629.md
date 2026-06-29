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

After rebasing over HYP-3452, the AP-tail clause is sharper than the original
HYP-3439 bridge.  HYP-3452 turns the rank-5 descent into a phase theorem:
`m=1..4` are finite mixed transients, `m>=5` has the rank-one `E:84m/E:84m`
endpoint component, paired dead-cover rank is already `<=2` from `m>=3`, and
escape counts follow a mod-35 Beatty correction.  That result should be used
as the AP-tail proof skeleton rather than re-scanning the bridge rows.

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

Next pull: formalize the canonical rank-`6` base case by combining HYP-3431's
corridor fence with HYP-3450/HYP-3451's four low-rank component escapes; then
use HYP-3452 to prove the rank-`5` AP-tail descent through finite transients
`m=1..4`, the rank-one `m>=5` phase, and the mod-35 escape clock.  Only after
that should the bridge be broadened to arbitrary primitive covering rows.
