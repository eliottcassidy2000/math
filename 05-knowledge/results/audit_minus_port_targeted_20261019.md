# Independent targeted audit: minus-port bounded-section atlas

**Verdict: PASS.** No proof-critical issue or required correction found.

Reviewed `square_sum_minus_bounded_section_contacts_20261019.md`, including
its newly added reciprocal-involution subsection, against the companion
plus-port note and the exact recurrence-bound theorem. This audit does not
repeat the broad finite census and makes no novelty assessment.

## Checks completed

- Independently reconstructed the cubic slope/intersection group law and
  checked the negative-sheet anchor, the displayed doubling abscissa,
  and both inverse coordinates of `3Q_-`, including the minus sign on z.
- Expanded the integer-polynomial square identity independently; it is
  identically zero after bringing both sides together. The resulting root
  has the required minus-input/plus-output square.
- Independently reproduced the shifted denominator, displacement, and
  endpoint polynomial coefficients and all three exact negative-tail sums.
  They prove `E>31`, `L=(N-tE)/(t-1)>63`, and `5E-4N>4` throughout
  `1<t<=17/16`. In particular, the bound of 63 concerns L, not v−t.
  These inequalities imply `1<t<v<5/4`. The derivative of f is positive
  there, and `f(5/4)=720/1681<1/2`; all claimed admissibility and
  nonvanishing facts follow, including at the closed cutoff endpoint.
- Confirmed the three exact polynomial reversals
  `t^11 N_+(1/t)=-E`, `t^11 E_+(1/t)=-N`, and
  `t^22 Z_+(1/t)=W`; hence `v_-(t)=1/u_+(1/t)`.
  Since `f(1/t)=-f(t)`, this transports the square identity correctly.
  The text properly limits this to algebraic duality: the plus-port
  admissibility interval does not contain `1/t` for `t>1`.
- Checked the fallback against the inherited theorem's precise anchor
  hypothesis. Returning to the positive-z anchor is necessary and is done
  correctly. For `t>17/16`, `f_0>35904/297025>1/10` and D=2 give
  `D/f_0<20<32`, so `J<=5`. The bound is `2 floor(900/7)=256` and the
  signed-even nonzero list contains exactly 256 candidates. The inherited
  theorem also supplies the strict changed-spread and ordinary-chart
  conditions, so no compactness assertion near the upper boundary is
  needed.
- Checked degrees `(11,10)`, coefficient norms `(446,760)`, and
  `446^2+760^2=776516`. Homogenization to degree 11 and primitive
  normalization therefore prove `C'<=776516(2C)^11` on the cusp chart.
  On the fallback, a bounded number of rational operations, including
  the finite identity/inverse/doubling branches, gives effective absolute
  polynomial numerical-height bounds and O(log C) operand/output/certificate
  bits. Integer gcd/reduction and comparisons preserve the conservative
  O((1+log C)^3) bit-cost claim. Taking maxima with the independently
  established plus-port atlas yields the stated both-port theorem.

## Reproducible controls

The independent checker `audit_minus_port_targeted_20261019.py` passes
28 exact symbolic/rational controls. Its normal and `python3 -O` outputs
are byte-identical. Results are saved in the matching `.out` and
`.optimized.out` files. The checker imports no companion construction code.

Only these uniquely named audit artifacts were written. No main research
files or repository files were changed, and no push was performed.
