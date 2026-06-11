# HYP-2412 - Sign-law Tournament Analysis ranks cancellation quotients

**Status:** METHOD CONFIRMED on the pilot set; mathematical content OPEN.
**Source:** codex-2026-06-11-P1.
**Artifacts:** `04-computation/pentagonal_lyapunov_code72_codex.py`,
`05-knowledge/results/pentagonal_lyapunov_code72_codex.out`.

## Setup

Tournament vertices are sign laws on generalized pentagonal support:

```text
Euler, all-plus, all-minus, period-3 paired, Thue-Morse paired,
Legendre-7 paired, and selected random paired samples.
```

Pairwise observable: compare finite-window reciprocal Lyapunov slopes. The lower
slope wins because it represents stronger cancellation. Ties use the listed order
as the Hamiltonian tie path.

## Pilot Fingerprint

For the nine-vertex pilot:

```text
score histogram: {0:1, 1:1, ..., 8:1}
directed 3-cycles: 0
SCC sizes: nine singleton SCCs
Hamiltonian path count: 1
```

The pilot tournament is transitive because the observable is a scalar. This is a
feature, not a theorem: it gives a clean ranking of cancellation quotients while
advertising that the quotient has forgotten root geometry and support design.

## Next Upgrade

Replace scalar slope by a two-coordinate observable:

```text
(nearest-zero radius, root-angle rigidity)
```

or by a code-facing pair:

```text
(low-weight suppression depth, support-design compatibility).
```

Those observables should create nontransitive cycles, identifying where "better
cancellation" and "more realizable support" diverge. That divergence is likely
where the `[72,36,16]` problem lives.
