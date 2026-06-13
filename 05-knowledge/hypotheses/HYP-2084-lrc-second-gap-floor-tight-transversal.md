---
id: HYP-2084
status: WEAKENED by S573; floor-only claim holds only in the original small boxes
source: codex-2026-06-03-S572
related:
  - HYP-2088
  - HYP-2083
  - HYP-2082
  - HYP-2153
  - HYP-2052
  - HYP-2055
  - HYP-2059
---

# HYP-2084: below `2/(2n-1)`, the bounded residual is already the floor-tight perfect-transversal branch

**Original claim, now weakened.** The right separator below the spectral-gap
edge `2/(2n-1)` is not coarse AP/sumset minimality.  In the original S572
bounded boxes, rows below `2/(2n-1)` were already `n`-clock-tight rows with
`M(S)=1/n`, and their residues mod `2n-1` formed perfect antipodal
transversals.

**S573 correction.** This does not globalize to larger lifted integer boxes.
`lrc_second_gap_bounds_s573.py` finds genuine open-gap rows:

```text
n=7: (1,5,6,11,16,17), M=5/33 in (1/7, 2/13)
n=8: (1,2,3,4,5,7,18), M=3/23 in (1/8, 2/15)
n=8: (1,3,4,5,7,13,18), M=3/23 in (1/8, 2/15)
```

The surviving corrected claim is local: S572's floor-tight statement is true in
its original small boxes, but the larger theorem target is HYP-2088's
three-clock blocker ledger, not a global floor-only separator.

**Exact bounded audit (`lrc_second_gap_transversal_audit_s572.py`).**
- Primitive boxes scanned:
  - `k=3, max_speed=20`
  - `k=4, max_speed=16`
  - `k=5, max_speed=13`
  - `k=6, max_speed=11`
- In every row with `M(S) < 2/(2n-1)`:
  - route is `n_clock`,
  - exact value is `M(S)=1/n`,
  - residues form a perfect antipodal transversal modulo `2n-1`.
- The full bounded menu below the edge is:

```text
n=4: AP only
n=5: AP and flip-set {2}
n=6: AP and flip-set {2}
n=7: AP only
```

So the bounded residual matches the known small-`n` floor-tight transversals
from S553/S553b.

**S602 additive-chain refinement.** `lrc_p0_collapse_additive_chains_s602.py`
rechecks the primitive canonical transversals for `n=4..8`.  The only non-AP
primitive floor transversals are exactly the flip-set `{2}` rows:

```text
n=5: (1,3,4,7)      = 1,3,1+3,3+(1+3)
n=6: (1,3,4,5,9)    = 1,3,1+3,1+(1+3),(1+3)+(1+(1+3))
```

Both are two-seed additive chains whose top is the sum of the previous two
chain terms.  Thus the local S572 floor-transversal residual is AP plus a
specific sparse additive-chain flip, not AP plus arbitrary small sumset.

**Important correction.** Sumset minimality is not the sharp separator. It is
sufficient but not necessary. Example:

```text
(1,3,4,5,9)
```

lies below `2/(2n-1)` while having positive sumset excess `3`. So "AP-like
small sumset" is weaker than "floor-tight perfect transversal."

**Interpretation.**
- Addition creates the odd summand shells `{a,2n-1-a}`.
- Multiplication by units makes unit shells visible as witness clocks.
- Odd `2n-1` removes the midpoint, so the gap story is shell coverage.
- Evenness re-enters only at the floor branch, through the `n`-clock midpoint /
  apex defect.

This suggests the second-gap proof should split:

```text
1. missed unit shell -> 2/(2n-1) witness;
2. residual perfect transversals -> floor-tight classification by n-clock;
3. composite nonunit holes -> second-clock / gcd-stratum / endpoint-core closure.
```

**Honest scope.** This is not a proof of the global spectral gap, and S573 shows
the global floor-only reading is false.  Treat this file as a useful small-box
observation and a stepping stone to HYP-2088.

**See:** `07-reflections/lrc-addition-multiplication-odd-even-second-gap-s572.md`,
`04-computation/lrc_second_gap_transversal_audit_s572.py` (+.out),
`04-computation/lrc_second_gap_bounds_s573.py` (+.out),
`04-computation/lrc_p0_collapse_additive_chains_s602.py` (+.out),
`07-reflections/lrc-2n-minus-1-summand-unit-bridge-s571.md`,
S553/S553b, HYP-2088, HYP-2083, HYP-2082, HYP-2153.
