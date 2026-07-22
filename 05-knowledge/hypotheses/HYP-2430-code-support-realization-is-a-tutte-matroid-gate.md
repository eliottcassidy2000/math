# HYP-2430 - Code support realization is a Tutte/matroid cancellation gate

**Status:** OPEN synthesis / next computational route.
**Source:** codex-2026-06-11-P2.
**Companions:** HYP-2415, HYP-2425, HYP-2429, HYP-2441, OPEN-Q-063,
Greene's theorem, Tutte partition functions.

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

## Exact linear non-identifiability control

Let `V` be the rational vector space on the `16`-subsets of a `72`-set and
let `W` be the rational vector space on its `5`-subsets.  The incidence map

```text
partial_5: V -> W,
partial_5(B)=sum_(T subset B, |T|=5) T
```

cannot be injective because `dim V=C(72,16)>C(72,5)=dim W`.  Clearing
denominators in a nonzero kernel vector gives an integral signed `5`-trade:
two formal block multisets with identical `5`-marginals.  This does not prove
the existence of two simple `5-(72,16,78)` designs, but it proves that the
linear five-incidence data alone cannot identify a formal support layer.

The small hostile control is the Pasch trade

```text
{123,145,246,356} <-> {124,135,236,456}.
```

Every vertex degree and pair codegree is preserved while every block changes.
It is not a parameter-`(72,16,5)` construction; it demonstrates exactly why
low marginals need a realizability sidecar.  For the code problem that sidecar
must at least retain binary linearity, self-orthogonality, simplicity, and the
column-matroid circuit structure.

This ambiguity survives the exact order-five symmetry of HYP-2441 and can be
made disjoint from every fixed block.  For an action of type `5^14 1^2`, a
nonidentity element fixes

```text
2*C(14,3)=728 sixteen-subsets,       14 five-subsets.
```

Hence Burnside gives `823261004433604` block orbits and `2798320`
five-subset orbits.  More sharply, the invariant subspace spanned only by
moving block-orbit sums has dimension

```text
(C(72,16)-728)/5 = 823261004432876 > 2798320.          (1)
```

Equivariance of `partial_5` therefore gives a nonzero `C_5`-invariant integral
signed `5`-trade supported entirely on moving block orbits, with kernel
dimension at least `823261001634556`.  It preserves every lower marginal as
well.  Thus HYP-2441's fourteen fixed minimum blocks may be held literally
fixed while its scalar and five-incidence ledger still has an enormous formal
kernel.  This is a non-identifiability theorem for formal block systems, not
the construction of another simple design or a binary code.

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
