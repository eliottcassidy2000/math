# LRC14 AP84 Finite Transients

HYP-3457 turns the last finite AP-tail instruction into a concrete packet.
For `m=1..4`, there are exactly four low-rank survivor windows, all with
endpoint rank `2`, and the formulas are mirror-paired across the two HYP-3431
low corridors.

The useful sign is

```text
(98m+13)/(588m) - 6/35 = (455-98m)/(2940m).
```

It is positive for `m=1..4`, so the moving high gap still crosses the fixed
`B1:5` wall and the best escape is mixed `E:84m/B1:5`.  At `m=5` the sign
switches, which is exactly the HYP-3454 pure `E:84m/E:84m` phase.

This narrows the AP-tail handoff.  HYP-3454 handles the infinite endpoint
interval, HYP-3456 handles the mod-`35` escape count, and HYP-3457 handles the
finite transient packet.  The remaining work is no longer "check m=1..4"; it
is the carrier/splice step: import HYP-3431 as the complete low branch-union
carrier and feed the three AP-tail sidecars into HYP-3439.

Assumption challenge: I considered organizing the transient proof by runners,
`m`, residues, high-grid gaps, four survivor windows, fixed corridor walls,
dead-cover graph summaries, and proof obligations.  The survivor-window
quotient is the right one here: it preserves the AP-tail escape predicate and
endpoint labels, while making clear that non-AP component geometry is outside
this sidecar.
