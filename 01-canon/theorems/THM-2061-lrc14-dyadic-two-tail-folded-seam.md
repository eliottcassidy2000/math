---
id: THM-2061
title: "LRC(14) dyadic two-tail folded seam and metric box"
status: >
  RESERVED WITH A PROOF CANDIDATE UNDER AUDIT. Sharp sheet capacity reduces a
  primitive thirteen-speed row with an imprimitive eleven-speed core to the
  seam 2C union {x,y}, x,y odd. Failure at 1/14 should be equivalent to a
  persistent opposite-parity cover over the full weak-safe set of C, and a
  maximizer interval should force x,y<=12 max(C) plus a two-centre separation
  inequality. This is the 1/14 analogue of THM-769's 1/13 folded criterion;
  it is not yet a proof that the seam is empty.
source: codex-2026-07-21-LRC14-dyadic-seam
depends_on:
  - THM-2060
related:
  - THM-761
  - THM-769
  - THM-775
  - HYP-8846
---

# THM-2061 -- LRC(14) dyadic two-tail folded seam

The reserved claim concerns

```text
S=2C union {x,y},       |C|=11,       x,y odd,
delta=1/14.
```

For `phi_C(theta)=min_(c in C)||c theta||`, let

```text
G_C={theta in R/Z: phi_C(theta)>=delta}.
```

The candidate exact obstruction is that `M(S)<delta` if and only if, for
every `theta in G_C`, both odd tails are `2delta`-eligible and their nearest
integer parities are opposite. The proposed metric consequence uses a
maximizer component of `G_C`; settled lower-dimensional LRC gives it radius
at least `1/(84 max(C))`, which should imply

```text
x,y<=12 max(C),
1/(xy)+2 rho <= 1/(7x)+1/(7y).
```

The strict/weak endpoint conventions, connected-component parity constancy,
and exact relation to THM-769/775 still require referee review.
