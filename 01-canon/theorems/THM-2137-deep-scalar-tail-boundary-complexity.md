---
id: THM-2137
title: "Deep scalar tails force unit-comb boundary complexity"
status: >
  RESERVED / PROOF AUDIT IN PROGRESS. A fibre-discrepancy argument is being
  checked which would bound a common 13-power in the scalar tails by the
  number of guard/unit boundary components. The intended constants are
  191/6930 in the six-plus-one branch and 961/6930 in the five-plus-two
  branch. The almost-everywhere cover, endpoint conventions, and the use of
  THM-2091's distinct-charge spectrum still require a written audit here.
source: codex-2026-07-22-LRC-deep-fibre-boundary-complexity
depends_on:
  - THM-2091
  - THM-2133
  - THM-2135
---

# THM-2137 -- deep scalar tails force boundary complexity

This identifier is claimed for a prospective all-depth reduction of the two
scalar tails in THM-2135.  The proposed proof changes the quotient from one
thirteen-root needle to the complete root fibre of the largest common
thirteen-power in the deep coefficients.  A clean fibre outside the divided
deep danger combs must then miss the residual guard/unit set entirely.

What is presently known is the target statement

```text
B=H+sum_i q_i,
delta_6=191/6930,                 delta_5=961/6930,

six-plus-one:  B >= delta_6 13^nu13(13v),
five-plus-two: B >= delta_5 13^min(nu13(13v_1),nu13(13v_2)).
```

No claim is made until the null-set selection and shifted-grid discrepancy
argument have been recorded below.  Even after proof, these inequalities are
height-versus-depth reductions rather than a closure of LRC(14).
