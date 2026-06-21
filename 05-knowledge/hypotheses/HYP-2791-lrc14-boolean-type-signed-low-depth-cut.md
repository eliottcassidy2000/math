---
id: HYP-2791
title: LRC14 Boolean type quotient has a sharp low-depth signed aggregate cut
status: CONFIRMED finite k=8 bounded-bank certificate; OPEN as global proof route
source: codex-2026-06-21
depends_on:
  - HYP-2744
  - HYP-2747
  - HYP-2746
  - HYP-2738
related:
  - THM-534
  - HYP-2602
  - HYP-2726
  - HYP-2740
  - HYP-2750
---

# HYP-2791: Boolean Type Signed Low-Depth Cut

## Claim

In the k=8 bounded bank

```text
E = {0} + 7-subsets of [1,14],
```

after quotienting out the two-row AP type orbit

```text
(0,1,2,3,4,5,6,7),
(0,2,4,6,8,10,12,14),
```

the first nontrivial Boolean/type separator, beyond the already-known `q0`
cover atom, lives on exactly three low-depth dihedral missed-sector atoms:

```text
T1    = (1,(1))      one missed inner sector,
T2sep = (2,(1,1))    two separated missed inner sectors,
T2adj = (2,(2))      two adjacent missed inner sectors.
```

There is a positive low-depth functional

```text
Phi(E) = 21*T1(E) + 57*T2sep(E) + 2*T2adj(E)
```

such that AP strictly minimizes `Phi` among all non-AP-type rows in this bank.
Equivalently, the signed functional `-Phi` is a small signed aggregate cut on
which AP is the strict maximizer.

## Exact Certificate

Script:

```text
04-computation/lrc14_boolean_mobius_signed_cut_codex_20260621.py
```

Stored output:

```text
05-knowledge/results/lrc14_boolean_mobius_signed_cut_codex_20260621.out
```

The full 12-state dihedral type quotient first rediscovers the trivial `q0`
separator:

```text
q0-only exact min gap modulo AP orbit = 3/56.
```

Excluding `q0`, the LP support is exactly the three atoms above.  The sharp
L1-normalized positive separator has coefficients

```text
60601/231049, 164633/231049, 5815/231049
```

on `(T1,T2sep,T2adj)` and exact margin

```text
2324813/97040580.
```

The three active non-AP rows are

```text
(0,1,2,4,5,6,7,10)   diffs = (27/490, 11/980, 44/735)
(0,2,3,4,5,6,7,8)    diffs = (25/1176, 149/5880, 19/1470)
(0,4,5,6,8,9,10,14)  diffs = (157/8820, 23/980, 899/8820)
```

The compact integer certificate

```text
(21,57,2)
```

has exact minimum margin

```text
21*dT1 + 57*dT2sep + 2*dT2adj >= 8447/4410
```

before normalization, or

```text
8447/352800
```

after dividing by `21+57+2=80`.  Its tight witness is

```text
(0,4,5,6,8,9,10,14).
```

## Interpretation

This is a small signed aggregate cut in the 64-state Boolean lattice.  The
three features are sums of exact atom probabilities `q[M]` over all missed-set
masks with the indicated cyclic run type.  Thus the quotient preserves cyclic
adjacency of missed sectors, but destroys speed ownership, relation height,
and the `x`-location of the wall atoms.

The sign is important.  Most individual type atoms have many AP-beaters
(HYP-2744), so the useful statement is not monotonicity of a single positive
atom.  It is a three-term aggregate balance: AP has enough low-depth miss mass
in the separated and adjacent two-miss channels, relative to singleton misses,
to force a strict cut after the AP dilation orbit is removed.

This also clarifies the CJJ/Mobius hierarchy guardrail.  The generic SoS/theta
lift can collapse to the size quotient, and the per-subset Mobius certificate
can become circular for global extremality.  HYP-2751 is not claiming a generic
hierarchy proof.  It is a concrete finite quotient certificate in the generated
LRC sector laws.

Incoming `THM-563` and `HYP-2788` give the most promising lift target.  THM-563
turns the signed one-far deviation into a finite endpoint-period maximum, while
HYP-2788 routes wide near-cap rows to a single-perturbation bounded scaffold.
That is structurally parallel to HYP-2751: both replace an apparently analytic
or high-dimensional obstruction by a small signed finite ledger.  The natural
next check is whether the low-depth type functional above remains a stable
bounded-scaffold certificate on the finite non-consecutive base list that
survives the THM-563/HYP-2788 wide reduction.

## Formalization

Lean source:

```text
04-computation/lean/TournamentH7/TournamentH7/LRCBooleanTypeCut.lean
```

The Lean module records the exact arithmetic of the sharp coefficients, compact
coefficients, active rows, and margins.  It proves, by finite `native_decide`
arithmetic:

```text
optimalCoeff_sum = 231049
compactCoeff_sum = 80
all three sharp active rows saturate the sharp LP cut
the compact cut is valid on the sharp active rows
active2 saturates the compact cut
both margins are positive
```

The full 3430-row validation is carried by the exact Python script; the Lean
module is the arithmetic nucleus that should be extended if the finite bank is
later imported as data.

## Tournament Analysis

Vertices are proof lenses, not runners:

```text
full_64_boolean_law
compact_low_depth_cut
sharp_three_atom_LP
relation_level2_pin
q0_trivial_cut
generic_CJJ_level2
```

Pairwise observable:

```text
(separation, preserves_mask_shape, tractability,
 compatibility_with_incoming_hierarchy_work, formalizability)
```

The computed tournament is transitive:

```text
full_64_boolean_law
> compact_low_depth_cut
> sharp_three_atom_LP
> relation_level2_pin
> q0_trivial_cut
> generic_CJJ_level2.
```

## Assumption Challenge

The useful vertices here are not runners, tournament arcs, Tanner checks, or
Fourier modes.  They are low-depth Boolean atom types in the missed-sector
lattice.  This quotient preserves exactly the information needed for the cut
and forgets several things that are probably needed for the global proof:
height of relations, speed ownership of walls, and far-element phase data.

The next proof obligation is to lift this finite bounded-bank certificate into
one of the existing global ledgers: relation-level marginal matching,
state-word transport, or the THM-563 endpoint-period far-element ledger.
