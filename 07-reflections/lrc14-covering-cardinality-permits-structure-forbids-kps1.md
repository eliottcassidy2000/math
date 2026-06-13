# LRC14 covering: cardinality permits, structure forbids

THM-497 is a useful correction to the way we talk about denominator shells.

The band criterion can be read as a cover problem: each runner contributes a
dilated symmetric band in `(Z/q)*`, and a strict witness is just an uncovered
unit.  That turns the LRC14 residual into a finite cyclic covering question
with a clean deficit function `D(q,S)`.

The first surprise is negative but clarifying.  Counting cannot prove the
deficit is positive: in band `k`, thirteen runners have enough total interval
mass to cover the unit group.  So the obstruction is not cardinality.  It is
alignment.  The hard rows are hard because their dilated bands over-correlate.

This fits the Q27/Q31 work neatly.  Plain shell ceilings are scalar shadows and
can be false globally; the retained fiber and valuation layers are where the
proof information survives.  HYP-2470 is a bounded eight-core theorem, HYP-2471
is the fiber repair for the same local exceptions, and THM-497 supplies the
global covering language in which the next discrepancy theorem should live.

---

## The second lesson (kps1): the tool-domain boundary

The dispatch asked to apply the repo's *novel* contributions to LRC 14. The
honest verdict from an adversarial pass: the celebrated machinery — the sum-free
grading dichotomy (THM-469), the pentagonal winding-positivity (THM-488), the
η-discriminant (THM-489), the OCF `H = I(Ω,2)` Paley character sums, mod_rank —
gives **essentially no transfer**. Every one lives on the *additive* structure of
`Z/q` or on tournament/code generating functions; the LRC deficit `D(q,S)` is a
*multiplicative*-character sum over the unit group (a dilate `v^{-1}B_q`). Even
the `14 = 2·7` CRT fibration does not cleanly separate the obstruction.

This is worth recording precisely *because* the repo's recurring miracle is that
its objects reduce to the same few atoms (φ, the pentagonal product, the dyadic
valuation). It is tempting to assume the atoms are universal. They are not: the
tournament/code toolkit is additive-and-generating-function; LRC is
multiplicative-character/covering, and the on-target instruments are
incomplete-character/Weil bounds plus the fiber-ladder descent (the Q27/Q31
route). A stated tool-domain boundary is as useful as a transfer — it stops the
cluster forcing the wrong key into the lock. (Contrast
[[the-three-twos-kp0611]], where the dyadic atom genuinely *did* transfer; the
difference between the two cases is exactly additive-vs-multiplicative.)

And the corrective coda: the naive `D(q,S) = q(6/7)^13 + O(√q · polylog)` bound
that looks like the obvious next step is FALSE for structured configs — they
over-correlate, the deviation grows faster than `√q`. The deficit theorem must be
proved in the over-correlated regime, not the independent one.
