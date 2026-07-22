# HYP-2430 - Code support realization is a Tutte/matroid cancellation gate

**Status:** OPEN synthesis / next computational route.
**Source:** codex-2026-06-11-P2.
**Companions:** HYP-2415, HYP-2425, HYP-2429, OPEN-Q-063, Greene's theorem,
Tutte partition functions.

## Statement

The scalar weight enumerator of a binary linear code is a Tutte specialization of
the associated binary matroid. Thus the `[72,36,16]` problem should be attacked as
a matroid support-realization gate, not only as a polynomial/invariant-ring gate.

Gleason supplies the scalar target. A binary matroid/code support must realize
that target while satisfying:

- self-duality;
- doubly-evenness;
- minimum distance 16;
- the weight-16 `5-(72,16,78)` design layer;
- neighborhood and automorphism restrictions.

## THM-2069 sharpening: the first failed deletion radius

For a binary generator matrix, THM-2069 identifies the `k`-deletion wheel
with the initial code weight distribution, and identifies its first failure
radius with column-matroid cogirth. Thus the support-realization target can be
stated exactly:

> Find a rank-36 binary column matroid whose deletion wheel is full through
> radius 15 and whose first cocircuit shell, at radius 16, consists of exactly
> 249849 supports realizing the forced `5-(72,16,78)` design, while the full
> row space is self-dual and doubly even.

The equality

```text
78*C(72,5)/C(16,5)=249849
```

shows that the scalar `A_16` count and the design block count agree. It does
not supply the incidence structure. In particular, a weight enumerator gives
the number of bad projective directions at each radius but forgets which
coordinates support them; this is precisely why the scalar Gleason gate is
insufficient. The Paley `eQR(72)` code has cogirth 12, so its deletion wheel
fails four layers too early.

## Why This Extends the Eta Analogy

Eta asks whether a sparse signed denominator has a product structure that keeps
zeros out of the disk. The code asks whether a scalar weight enumerator has a
binary matroid behind it. In both cases, the scalar partition function is a
shadow; the hidden support object is the theorem.

## Proposed Computation

Build a "Tutte leakage" diagnostic for candidate binary matroids/support systems:

```text
leakage = first forbidden low dual weight or violated design incidence
```

Random supports with the same low layer should leak almost immediately. A true
Type II support gate should cancel all leakage up to the Gleason-forced threshold.

Tournament Analysis upgrade: vertices should be support-building moves, not just
sign laws. Pairwise observable should compare `(low-weight suppression,
design-compatibility)` so cycles reveal tradeoffs between cancellation and
realizability.
