# HYP-2429 - Type II extremal formal enumerators stay scalar-positive far beyond length 72

**Status:** COMPUTATION CONFIRMED through length 240 in the stored atlas; scratch
scan saw no negative coefficient through length 1200. Code existence remains a
separate support problem.
**Source:** codex-2026-06-11-P2.
**Companions:** HYP-2425, HYP-2430, HYP-2415, HYP-2420.
**Artifacts:** `04-computation/cancellation_gate_atlas_codex.py`,
`05-knowledge/results/cancellation_gate_atlas_codex.out`.

## Statement

For Type II formal weight enumerators at lengths `n=24m`, the scalar Gleason gate
continues to produce integral, nonnegative extremal enumerators well beyond the
famous length-72 case. Therefore scalar positivity is not a sharp proxy for code
existence.

Stored run:

```text
n= 24 d= 8  A_d=759
n= 48 d=12  A_d=17296
n= 72 d=16  A_d=249849
n= 96 d=20  A_d=3217056
n=120 d=24  A_d=39703755
n=144 d=28  A_d=481008528
n=168 d=32  A_d=5776211364
n=192 d=36  A_d=69065734464
n=216 d=40  A_d=824142912363
n=240 d=44  A_d=9826364199840
```

For every row in the stored run, `W(1,1)=2^(n/2)` and all coefficients are
nonnegative.

## Interpretation

Length 24 is realized by the extended Golay code. Length 72 remains open despite
its perfectly healthy scalar enumerator. The scalar ladder therefore separates
three questions:

1. invariant-ring feasibility;
2. support/design realization;
3. classification or construction of actual codes.

Only the first is automatic in this atlas.

## Next Tests

1. Store a longer positivity scan with checkpointed output if needed.
2. Compare the first place scalar positivity fails, if it fails, with known
   nonexistence or shadow-bound thresholds.
3. Build a support-realization invariant per length: minimum-word design
   parameters, neighborhood constraints, and automorphism restrictions.
