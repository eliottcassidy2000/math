---
id: THM-2085
title: "Explicit height-57 rank-seven relation gate from Selberg box sandwiches"
status: >
  CLAIMED. This file reserves the explicit effective sharpening of THM-2083.
  The target is: G_Q subset E_h forces a nonzero relation
  a h+b q_i+c q_j=0 with max(|a|,|b|,|c|)<=57. The proof uses the classical
  degree-57 Selberg interval sandwich, a signed tensor box minorant, and the
  THM-2081 relative spanning-tree inequality. Exact constants and the
  distinction from the factorwise nonnegative-minorant no-go are under audit.
source: codex-2026-07-22-LRC14-Selberg-relative-Hunter
depends_on:
  - THM-2081
  - THM-2083
related:
  - THM-537
  - THM-2051
---

# THM-2085 -- explicit height-57 rank-seven Selberg gate

Reserved for the paper proof and exact rational referee.  The proposed
constants at degree `57`, with `epsilon=1/58`, are

```text
I_i >= 5363/164836,
w_ij >= 655135/66923416,
6 min(w_ij) - (2/7-7 min(I_i))
  =6435/8365427 >0.
```

The remaining checkpoint obligation is to spell out the pointwise signed box
minorant and verify every fraction independently before promoting the status.
