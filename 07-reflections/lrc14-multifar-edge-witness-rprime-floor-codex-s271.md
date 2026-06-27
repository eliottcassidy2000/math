# LRC14 Multi-Far Edge-Witness Rprime Floor Reflection

Anchors: HYP-3125, HYP-3124, HYP-3121, HYP-3123, HYP-3122, T1199,
LTI-260, LTT-158, OPEN-Q-108.

The useful reframing in this pass is that the open `Rprime >= c` floor is not
just a scalar correlation estimate.  It is a directed proof edge:

```text
R-safe packet -> Q-safe packet
```

HYP-3121 says this diagonal intersection is the middle engine of LRC(14).
HYP-3124 says a directed edge is only a witness after both endpoint-deletion
children are retained.  Combining those statements makes the floor a two-ended
object: the tail recursion deletes or relaxes a small-part `R` constraint,
while the tip recursion deletes or conditions on a multiple-of-14 `Q`
constraint.  Losing either endpoint turns the floor back into the failed
Bonferroni scalar.

The S271 audit is deliberately modest, but it lines up with the incoming
picture.  On `{1..11,13,84}` and `{1..10,13,84,154}`, Bonferroni is negative
while `Rprime` is positive (`0.513784` and `0.925326` on the midpoint grid).
Tail deletion generally improves the ratio, while tip deletion exposes the
multiple-of-14 side as a sharper child.  Individual R/Q edge ratios stay near
`1`, so the obstruction is not one bad pair edge.  It is a packet-level
distribution question across many sector cuts.

That is where the wild analogies become disciplined:

- Elliott-Halberstam becomes only a level-of-distribution metaphor: prove an
  LRC-specific BV/L2 average over edge-sector residue errors outside a finite
  major-arc ledger.  Do not import a prime theorem.
- Gaussian functions become a positivity interface: prove the diagonal after
  heat/Selberg smoothing, then desmooth through finite-ruler and Cech
  component fields.
- Asano contraction becomes a legality test: contract tail and tip variables
  only when the two-ended sector polynomial carries a Lee-Yang zero-free
  sidecar, with HYP-3122's phi4 cumulant signal as a stabilizer.

The selected tournament path puts
`edge_witness_two_ended_RQ_packet` before the raw floor, followed by Gaussian
smoothing, Cech finite-ruler desmoothing, Asano contraction, chiral guard,
phi4 stabilization, H7 state-lift, EH proxy, and raw Bonferroni.  That ordering
is not saying the edge packet proves the floor by itself.  It says the floor
has to be carried in an edge packet before the analytic estimates are legal
proof data.

Post-rebase, HYP-3127 strengthens the same lane: Asano is no longer merely a
legality check after smoothing, but a candidate main contraction engine because
coverage is multi-affine.  That turns this S271 packet into the audit harness
for HYP-3127.  Its row fields should record the single-far zero-free polydisk,
the `SPEC`/signed-cancellation bound, and the recursion-monotonicity status
beside the tail/tip deletion-child ratios.  EH should stay as a
level-of-distribution/SPEC language unless a real LRC-specific bound is
proved; Gaussian should be read as the decoupled free-field limit that the
Asano contraction perturbs.

Next proof move: build an `edge_floor_packet` row for real `r=2..6` covering
packets, with `Rprime`, both deletion-child ratios, residue exception sets,
Gaussian smoothing width, finite-ruler threshold, Asano/Lee-Yang contraction
status, phi4 sign, normal-fan chamber, chiral guard word, and terminal exit.
The hoped-for theorem is a dichotomy: positive smoothed diagonal plus
finite-ruler desmoothing gives the uniform floor, or the exception is forced
into an Asano/Lee-Yang/phi4/Cech/H7 named debt class.
