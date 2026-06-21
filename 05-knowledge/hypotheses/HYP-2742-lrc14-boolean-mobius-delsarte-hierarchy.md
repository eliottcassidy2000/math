---
id: HYP-2742
title: LRC14 sector atoms admit a complete Boolean-Mobius Delsarte hierarchy
status: OPEN hierarchy target; exact k=8 bounded-bank scout
source: codex-2026-06-21
depends_on:
  - THM-534
  - HYP-2740
  - HYP-2738
  - HYP-2602
related:
  - HYP-2726
  - HYP-2739
  - HYP-2741
  - HYP-2735
  - HYP-2603
---

# HYP-2742: Boolean-Mobius Delsarte Hierarchy

## Claim

The higher-order Delsarte hierarchy for linear codes suggests the right LRC
analogue: stop keeping only the size of the missed-sector set.  For the six
inner sectors `U={1,...,6}`, define exact atom and containment variables

```text
q[M] = meas{x : the exact missed inner-sector set is M},
a[A] = meas{x : A subset M(x)} = sum_{M superset A} q[M].
```

Then Boolean Möbius inversion gives the exact atom law:

```text
q[M] = sum_{A superset M} (-1)^(|A|-|M|) a[A].
```

The THM-534 Delsarte depth law is just the size quotient

```text
p_t = sum_{|M|=t} q[M].
```

The full `64`-mask Boolean lattice is complete for the sector atom law, while
intermediate quotients such as dihedral run type give hierarchy levels between
the current `7`-state Delsarte LP and the full generated-sector object.

## External Inspiration

The code-hierarchy analogy is precise enough to be operational:

- Coregliano-Jeronimo-Jones introduce LP hierarchies generalizing Delsarte via
  higher-order Krawtchouk/MacWilliams constraints over interactions of `ell`
  words.
- The exact-completeness paper rewrites hierarchy variables as
  pseudoprobabilities and uses Möbius inversion on a subspace lattice to recover
  exact atoms.
- In LRC, the analogous lattice is not the subspace lattice but the Boolean
  lattice of missed sector subsets.  This is smaller (`2^6`) and fully exact.

## Exact Evidence

Script:

```text
04-computation/lrc14_boolean_mobius_hierarchy_codex_20260621.py
```

Stored output:

```text
05-knowledge/results/lrc14_boolean_mobius_hierarchy_codex_20260621.out
```

For consecutive/AP rows `k=8..13`, the script verifies exact
Möbius inversion and records the quotient sizes:

```text
k=8..13: mobius_inversion_ok=True
k=8..13: mask_support=32 of 64
k=8..13: dihedral type states=12
```

Thus the consecutive rows live on only half of the full Boolean atom lattice,
but their type quotient is still richer than the current depth profile.

The k=8 bounded bank

```text
E = {0} + 7-subsets of [1,14],  3432 rows
```

tests whether any single dihedral type atom supplies a monotone proof cut.
Result: no, except for the already-known cover/deep-tail extremes.

AP is maximal only for:

```text
(0, ())  = full cover depth 0,
(5,(5)) = five missed sectors in one run,
(6,(6)) = all six inner sectors missed.
```

For every other type, there are AP-beaters.  Examples:

```text
(1,(1)):      3354 AP-beaters
(2,(1,1)):    3427 AP-beaters
(2,(2)):      3361 AP-beaters
(3,(1,2)):    3339 AP-beaters
(4,(4)):        11 AP-beaters
```

## Consequence

The hierarchy cut is real, but it will not be a one-positive-atom monotonicity
proof.  This confirms HYP-2738's lesson from a stronger basis: the cap-closing
certificate must be a signed aggregate cut in the type/full Boolean-Möbius
basis, not a single positive monotone statistic.

The next target is:

```text
Find a small signed functional Phi(q[M]) such that
1. Phi is nonnegative on all generated LRC sector laws;
2. Phi separates the AP/consec extremal row with exact rational margin;
3. Phi has a compact certificate in the containment variables a[A].
```

This is the LRC version of cutting the Delsarte polytope: add exact
Möbius/nonnegativity constraints that are invisible after the size quotient but
still tractable enough to prove symbolically.

## Tournament Analysis

Vertices are hierarchy views:

```text
full_boolean_mobius
> dihedral_type_quotient
> size_delsarte_depth
> single_type_atom
> raw_tanner_carrier.
```

Pairwise observable:

```text
(completeness, preserves_generated_masks, tractability, proof_signal)
```

The computed tournament is transitive with `directed_3cycles=0` and one
Hamiltonian path.

## Assumption Challenge

The useful hierarchy vertices are not runners, arcs, or Tanner checks.  They
are missed-sector masks and containment events.  The full Boolean lift preserves
the generated sector law exactly and destroys none of the sector information,
but it hides speed ownership, relation height, and finite-denominator geometry.
Those must still be supplied by the HYP-2602/HYP-2603/HYP-2741 proof pipeline.
