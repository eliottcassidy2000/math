---
id: HYP-2084
status: SUPPORTED by exact bounded audit + prior small-n gap data
source: codex-2026-06-03-S572
related:
  - HYP-2083
  - HYP-2082
  - HYP-2052
  - HYP-2055
  - HYP-2059
---

# HYP-2084: below `2/(2n-1)`, the bounded residual is already the floor-tight perfect-transversal branch

**Claim.** The right separator below the spectral-gap edge `2/(2n-1)` is not
coarse AP/sumset minimality. It is floor-tight perfect-transversal structure:
rows below `2/(2n-1)` are already `n`-clock-tight rows with `M(S)=1/n`, and
their residues mod `2n-1` form perfect antipodal transversals.

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

**Honest scope.** This is not a proof of the global spectral gap. It is exact
bounded evidence that the hard branch below `2/(2n-1)` is already the
floor-tight branch, not a broad non-AP family.

**See:** `07-reflections/lrc-addition-multiplication-odd-even-second-gap-s572.md`,
`04-computation/lrc_second_gap_transversal_audit_s572.py` (+.out),
`07-reflections/lrc-2n-minus-1-summand-unit-bridge-s571.md`,
S553/S553b, HYP-2083, HYP-2082.
