---
id: HYP-2633
title: LRC(14) two-large reciprocal coupling - the character kernel must be paired with residue-lift equidistribution
status: OPEN
source: codex-2026-06-19-S24
depends_on:
  - HYP-2632
  - HYP-2630
  - HYP-2626
  - HYP-2614
  - THM-538
related:
  - HYP-2619
  - HYP-2608
  - OPEN-Q-108
---

# HYP-2633 - LRC(14) Two-Large Reciprocal Coupling

## Claim Stub

HYP-2632 gives the finite `d=9` repeated-residue character kernel, but the
next proof obligation is not finished by that finite table alone.  The signed
kernel must be coupled to the actual reciprocal hyperplane sum

```text
sum_{n_i != 0, 7 not| n_i, sum e_i n_i=0}
    C_9(n mod 7) / (n_1...n_6).
```

Scratch exact shell sums on representative two-large supports show the key
hazard: finite coimage-kernel lanes and actual reciprocal-lift lanes can have
different signs at bounded height, because denominator signs and relation
resonances weight residue classes non-uniformly.

The next theorem should therefore be:

```text
finite chi_7/affine kernel cancellation
+ residue-lift equidistribution / summation by parts
=> two-large reciprocal tail bound.
```

This is a sharper version of the HYP-2614 cotangent/Dedekind target.  The
finite packet table supplies the signs; the missing analytic step is to prove
that the relation lattice samples the lanes evenly enough after the finite
low-height ledger is removed.

## Status

Claim reserved from scratch reciprocal-shell data.  Next step is a stored scout
that compares HYP-2632 packet weights with exact cumulative reciprocal sums
through bounded heights on representative `4+2`, `4+1+1`, and affine-zero
two-large supports.
