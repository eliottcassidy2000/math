# HYP-2425 - The [72,36,16] Type II obstruction is support realization, not the scalar Gleason gate

**Status:** CONFIRMED for the scalar enumerator; OPEN for code existence.
**Source:** codex-2026-06-11-P1.
**Companions:** HYP-2415, HYP-2420, HYP-2421, OPEN-Q-061, T783.
**Artifacts:** `04-computation/pentagonal_lyapunov_code72_codex.py`,
`05-knowledge/results/pentagonal_lyapunov_code72_codex.out`,
`04-computation/cancellation_gate_atlas_codex.py`,
`05-knowledge/results/cancellation_gate_atlas_codex.out`.

## Statement

Building on HYP-2415's tournament-gauge formulation, the Type II Gleason
invariant ring completely determines the formal weight
enumerator of a putative extremal binary doubly-even self-dual `[72,36,16]`
code. With

```text
A = x^8 + 14 x^4 y^4 + y^8
B = x^4 y^4 (x^4-y^4)^4,
```

the unique degree-72 Type II formal enumerator with weights `4,8,12` suppressed is

```text
W = A^9 - 126 A^6 B + 3015 A^3 B^2 - 4398 B^3.
```

It has nonnegative integer coefficients, total `2^36`, and first nonzero
nontrivial coefficient

```text
A_16 = 249849.
```

Therefore the famous existence problem is not blocked at the scalar
modular-form level. The obstruction is the realization of this enumerator by an
actual
self-dual support system with the required orthogonality, doubly-evenness, and
design structure.

## Exact Data

The first coefficients of `W(1,y)` are:

```text
weight  0: 1
weight  4: 0
weight  8: 0
weight 12: 0
weight 16: 249849
weight 20: 18106704
weight 24: 462962955
weight 28: 4397342400
weight 32: 16602715899
```

The weight-16 words would form a `5-(72,16,78)` design; the derived lambda values
are:

```text
t=1: 55522
t=2: 11730
t=3: 2346
t=4: 442
t=5: 78
```

This is the finite analogue of the pentagonal cancellation gate: a small
invariant ring forces a forbidden prefix, but the hard problem is whether there
is a support object behind the scalar partition function. The atlas addendum
checks extremal Type II formal enumerators at lengths `24,48,...,240`; all are
integral, nonnegative, and sum to `2^(n/2)` in this range. That reinforces that
scalar feasibility is a weak gate; support realization is the real gate already
at length 72.

## Literature Check

As of this session's web check, the existence of an extremal `[72,36,16]`
self-dual code is still listed as unknown by Error Correction Zoo; Hannusch-Major
2023 gives newer necessary conditions via neighborhoods; Borello 2013 restricts
possible automorphism groups. This fits the support-realization reading: scalar
invariants allow the object, while local neighborhoods and automorphisms squeeze
possible supports.

## Next Tests

1. Build a local-support search for the `5-(72,16,78)` layer before trying to
   generate a full code.
2. Translate neighborhood constraints into a signed incidence matrix whose
   MacWilliams transform kills weights `4,8,12`.
3. Compare low-symmetry random support searches with the pentagonal random-sign
   model: generic supports should have positive "MacWilliams Lyapunov" unless
   self-duality forces a cancellation gate.
4. Treat the weight-enumerator as a matroid/Tutte-Greene partition function and
   ask which support matroids can pass the same low-weight suppression gate.
