# Focused review of the new graph-growth proofs

Date: 2026-11-11. This is a review of this packet's new statements, not
a fresh verification campaign on inherited repository claims. No numbered
canon promotion or literature-priority conclusion is made.

## Mathematical review

The summand, multiplicand, and paper-model lanes were developed separately.
Root integrated the prime-exponent chart and supplied the three-block
chase and shrinking-window formula. The summand lane proved the interval
family that completes saturation; the multiplicand lane proved the exact
independent-seed first-hole law. Other lanes cross-read the new interfaces.

* The paper-model reviewer checked the valuation chart, the block-chase
  bound, and the exact deficit window. In particular the tracked block
  forces the true second-largest gap into the retained window; its next
  lower edge lies above the old maximum, so no old vertex is silently lost.
* The multiplicand reviewer checked the additive observer and inverse,
  the synchronized scale formula, and both coefficient bounds. The
  root chart was corrected to say that extending the unmodified observer
  to the unit destroys its first finite-prefix diagonal; it does not
  globally annihilate the unit basis. The unit-present convention has a
  separate correction. The additive formal-series remainder was worded
  precisely at its first possible degree.
* The paper-model reviewer independently checked the entire stronger
  saturation proof, including overlapping input intervals and the fixed
  suffix threshold. The cardinality justification was tightened: the
  chosen T puts the maximum beyond K, and every permanent hole is below
  K. Early finite saturation alone would not imply that condition.
* Root checked the multiplicatively independent two-seed theorem, its
  single-mixed-branch lower bound, its numerical least-hole comparison,
  the first-failure criterion, and the dependent counterexample. The
  property is multiplicative independence, not simply being composite or
  lying on differently drawn graph branches.

The strongest hostile controls are `{2,3}` in addition (temporary
saturation can break), `{5,6,7}` (even stable top three do not by themselves
guarantee permanent saturation), `{2,3}` in multiplication (a forced serial
branch survives below a fast maximum), and `{2,4}` (dependence changes
the leaf counts and destroys that obstruction).

## Finite controls and their actual scope

* `summand/check.py`: exact closure memberships and residue updates;
  literal full-state versus moving-window comparison; additive matrix
  inversion and common-scale examples.
* `summand/top_suffix_check.py`: 832 interval sums, 38 seed sets, and
  exact bounded-window transients for consecutive seed triples through
  `{100,101,102}`. Large full states are not materialized for that last
  test. Both normal and optimized runs passed.
* `multiplicand/check.py`: 130,816 lower-triangular inverse equations
  through numeric cutoff 512, checked against the literal arc observer.
  Coefficient bounds are proved for all orders, not extrapolated from
  this matrix sample.
* `multiplicand/saturation_check.py`: 20 small initial states across
  independent examples and a dependent hostile. The all-round result is
  supplied by the proof, not this short sample.
* `papers/graph_norm_check.py`: actual permutation products, power/induction
  and norm identities, positive carries, paired coordinates, and finite
  filtration controls. The graph model does not satisfy every Raynaud
  relation; the explicit failure is retained in its proof.

The root reran all five checkers and checked relative Markdown links and
the repository documentation gate. Normal/optimized controls remain live
because the checkers raise explicit exceptions, rather than using Python
assert statements as the only validation.

The interactive node-birth view is separate response content. Its primary
slider interaction and displayed pairs were checked for every integer
2 through 40 against direct enumeration. Browser rendering was not
visually inspected. It illustrates numeric target births, not the
simultaneous seed-generation clock of the saturation theorem.
