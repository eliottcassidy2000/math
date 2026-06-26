---
id: HYP-2948
title: LRC14 Euler pentagonal recurrence and Beurling-Selberg minorants are labelled cancellation gates
status: COMPUTATIONAL SCOUT / cancellation-gate synthesis; not a proof
source: codex-2026-06-24
tags: [lrc14, euler, pentagonal, sigma, beurling-selberg, minorants, tetrahedral]
related:
  - HYP-2947
  - HYP-2946
  - HYP-2944
  - HYP-2877
  - HYP-2506
  - HYP-2426
  - HYP-2221
  - HYP-2220
results:
  - 04-computation/lrc14_ph_j37_pentagonal_minorant_codex.py
  - 05-knowledge/results/lrc14_ph_j37_pentagonal_minorant_codex.out
external:
  - https://arxiv.org/abs/math/0411587
---

# HYP-2948: pentagonal recurrence and minorants are labelled cancellation gates

The Euler paper in arXiv:math/0411587 is directly relevant to the previous
Farey/perfect-number pass.  It uses the pentagonal number theorem to derive a
recurrence for the divisor-sum function `sigma(n)`.

The computation verifies the recurrence through `n <= 30` from the sparse
Euler product coefficients:

```text
prod_{m>=1}(1-q^m)
  = 1 + sum_{k>=1} (-1)^k
      (q^{k(3k-1)/2} + q^{k(3k+1)/2}).
```

The first generalized pentagonal supports are:

```text
1, 2, 5, 7, 12, 15, 22, 26, ...
```

The script reports:

```text
sigma recurrence failures through 30: []
```

## Minorant Reading

This is the same proof shape as a Beurling-Selberg minorant in the LRC analytic
thread:

```text
sparse support,
signed coefficients,
one-sided defect/cancellation label,
controlled tail.
```

The point is not that Euler's recurrence is a Beurling-Selberg function.  The
point is that both are useful only with their labels retained:

```text
Euler: generalized pentagonal support + signs.
BS minorant: bandlimit + low Fourier coefficients + one-sided defect + tail budget.
```

Dropping those labels turns both into raw scalar approximations, which is
exactly the mistake the current LRC14 proof route avoids.

## Pentagon Versus Tetrahedral

This also refines the earlier Pollock-versus-pentagonal split:

```text
pentagonal numbers:  degree 2, Euler/quadratic recurrence lane
tetrahedral numbers: degree 3, Pollock/cubic finite-exception lane
```

The script prints the two sequences side by side.  They cross at small values,
but they live in different proof categories.  The pentagonal side supplies
cancellation recurrences; the tetrahedral side supplies bounded-arity additive
basis and finite-exception machinery.

## LRC14 Use

The candidate theorem-shape is:

```text
Any analytic LRC14 certificate based on a Beurling-Selberg minorant must carry
the same kind of labelled sparse-support data that Euler's sigma recurrence
carries: support, signs, one-sided error, and tail budget.
```

In the current proof order, this sits below:

```text
exact M/Farey branch,
AP/Goddyn-Wong labels,
C27 shell transfer.
```

but above:

```text
raw product/perfect-number analogies.
```
