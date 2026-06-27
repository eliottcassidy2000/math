---
id: HYP-2949
title: LRC14 real-factor recombination suggests branch-local C27 subset ledgers
status: COMPUTATIONAL SCOUT / recombination analogy with guardrails; not a proof
source: codex-2026-06-24
tags: [lrc14, recombination, subset-sum, c27, unital, farey, factorization]
related:
  - HYP-2948
  - HYP-2947
  - HYP-2946
  - HYP-2945
  - HYP-2942
  - HYP-2937
  - HYP-2934
results:
  - 04-computation/lrc14_ph_j37_pentagonal_minorant_codex.py
  - 05-knowledge/results/lrc14_ph_j37_pentagonal_minorant_codex.out
external:
  - https://arxiv.org/abs/2410.15880
---

# HYP-2949: real-factor recombination suggests branch-local C27 subset ledgers

arXiv:2410.15880 revisits integer polynomial factorization by first factoring
over the reals and then recombining linear/quadratic real factors.  The paper's
abstract says the recombination stage becomes an integer subset-sum problem.

The LRC14 analogy is not "polynomial factorization proves LRC."  The analogy is:

```text
candidate local packets are cheap,
the hard step is recombining them into globally legal integer/labelled objects.
```

For LRC14, the candidate packets are:

```text
Farey product factors,
C27 shell transfers,
q=3 unital branch-local blocks,
Kpq/K33 owner packets,
Beurling-Selberg low Fourier packets.
```

The recombination problem is to choose a compatible subset without violating:

```text
exact M/Farey branch,
AP/Goddyn-Wong marking,
C27 shell ownership,
unital lambda=1 pair uniqueness,
Paris-Harrington bad-child rank.
```

## Computation

The script records the scale difference between raw subset enumeration and
recombination-style savings:

```text
degree 20: 2^20 raw packets, 2^5 Horowitz-Sahni proxy
degree 40: 2^40 raw packets, 2^10 proxy
degree 80: 2^80 raw packets, 2^20 proxy
degree 100: 2^100 raw packets, 2^25 proxy
```

This is not a complexity proof for LRC14.  It is a proof-interface warning:
the productive object is a labelled subset ledger, not a raw product count.

## Relation To HYP-2942

HYP-2942 showed that one q=3 unital chart cannot globally carry both the tight
GW branch and the K33 near-miss branch under the raw residue-pair model,
because the unital has `lambda=1` pair uniqueness.

HYP-2949 turns that into a recombination rule:

```text
Do not put incompatible branch packets in the same recombination subset.
```

Use multiple branch-local charts, or split the repeated pair by an additional
branch coordinate.  That is the LRC version of preserving enough side data to
recombine real factors into honest integer factors.

## Candidate Test

Build a finite subset-sum ledger over the current low-gap packet library:

```text
GW H12->D3 block,
K33 H12->D9 block,
petal H10->D7 block,
petal H13->D1 block,
two-hole splice packets,
post-F4 product-excess leakage atoms.
```

Score a subset as legal only if it preserves:

```text
no repeated unital pair inside one chart,
exact Farey branch,
C27 gcd stratum,
PH-style bad-child rank decrease.
```

The hoped-for output is not a scalar optimum.  It is a short illegal-subset
certificate for every non-AP/GW residual recombination.
