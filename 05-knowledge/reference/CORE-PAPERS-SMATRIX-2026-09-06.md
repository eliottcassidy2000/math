# Li's S-matrix preprint: a scoped planar-JC method transfer

**Status: PRIMARY PREPRINT READ; local algebraic consumers independently
audited. The preprint's complete theorem and Lean development were not
independently certified in this session.** Checked September 6, 2026.

- Yinjie Li, *The S-matrix conjecture*,
  [arXiv:2608.29750v1](https://arxiv.org/abs/2608.29750v1), submitted
  August30,2026. [Primary full text](https://arxiv.org/html/2608.29750v1).
- The paper claims the remaining even-dimensional inverse-Frobenius bound,
  with a formalization conditional on Frankel--Urschel Lemma2.1. Its method
  combines a defect budget, binary rounding, and joint Gram measurements.
- The budget is under the contradiction hypotheses: even `n>2`, real
  nonsingular `B` with entries in `[0,1]`, and
  `||B^(-1)||_F<=2n/(n+1)`. Cube entries alone do not supply the bound.
- Sections2 and4 supply the shared-budget and projection mechanism;
  Section6 retains several column contributions in the exceptional order6.

The local consumers are the [complete response classification](../results/planar_jc_long_20260906_memory.md),
[sharp finite depth recognition](../results/planar_jc_long_20260906_depth.md),
[moving collision identity](../results/planar_jc_long_20260906_collision.md),
and [generic algebraic testers/integral Student response](../results/planar_jc_long_20260906_smatrix.md).
They use complete measurement systems and explicit discarded-data records.
Their proofs do not depend on the external S-matrix theorem.

The decisive transfer firewall is `(1,i)`: a nonzero complex vector can
have zero bilinear square. A Hermitian norm restores positivity but supplies
neither an algebraic complex-coefficient budget nor the paper's binary
rounding hypotheses. The local replacements are nonzero integer minors,
full polynomial coefficient ideals, and an actual moving-collision identity.
No S-matrix-to-JC implication or literature priority claim is made.
