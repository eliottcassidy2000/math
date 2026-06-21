# LRC14 Period-Max / Boolean Bridge, S76

The useful outcome is a correction, not a new period-max proof.

Incoming `THM-563` and S6/S7 make the single-far branch a finite
endpoint-period problem.  HYP-2791 then tempted a simple bridge: reuse the
three-term Boolean/type cut

```text
Phi_low = 21*T1 + 57*T2sep + 2*T2adj
```

as a bounded-base slack certificate.  The S76 overlay says that is the wrong
projection.  On the one-far base ledger the k=8 frontier has negative
`Phi_low-AP` witnesses, with minimum `-4153/3080`.  The HYP-2791 cut lives on
final-row Boolean laws; after deleting the far runner and replacing the far
effect by a period-max numerator, the three low-depth atom weights are no
longer monotone.

What survives is more useful:

```text
period-max direct scan
  > AP/dilation orbit filter
  > q0 cover-atom slack on bases
  > skipped-period audit
  > Phi_low as final-row or size-shifted coordinate
```

All non-AP frontier bases in the overlay have positive `q0` gap, with global
minimum `71/5880`.  The S6 broad-scan worst checked k=8 row
`(0,4,6,8,10,12,14)` has `Phi_low` gap `493/294` and `q0` gap `5/49`, so the
worst observed pressure is not a counterexample to the Boolean story.  The
counterexamples are lower-pressure base profiles where `Phi_low` moves in the
wrong direction.

Assumption challenge: the vertices here are not runners, arcs, or final rows.
They are proof obligations: exact period scan, AP/dilation filter, q0 base
slack, skipped-period audit, and final-row Boolean slack.  The quotient
preserves missed-sector cyclic shape and period-frontier address, but destroys
the far runner's wall ownership.  That destruction is exactly why a final-row
cut cannot be moved blindly onto the base.

Next proof target: formalize the base-level `q0` slack statement on the S6/S7
frontier and test whether the negative k=8 `Phi_low` witnesses become positive
again after adjoining the actual far runner at the period-max witness.
