---
id: HYP-1797
status: EXPLORATORY
source: opus-2026-05-30-S4
related:
  - HYP-1795
  - HYP-1796
  - THM-002
  - THM-354
---

# HYP-1797: Completeness-Defect Principle

## Statement

The parabolic frustration law

```text
E[Delta c3 | score=s] = s(n-1-s)/2
```

is the zeroth-order checksum of tournament completeness.  For structures that
are almost tournaments but contain missing, tied, weighted, or noisy pairwise
comparisons, the first diagnostic should be a completeness-defect vector:

```text
observed frustration injection - tournament parabolic prediction.
```

The hypothesis is that the stability of tournament-native invariants such as
`H`, OCF-style residue rank, phase vectors, and endpoint incidence depends
first on this completeness defect.

## Why This Is Different From Residue Rank

Residue rank assumes a tournament is already present and asks what survives a
projection.  Completeness defect asks whether the tournament axioms themselves
are stable enough for residue and phase features to remain meaningful.

This is relevant for:

- trienerments and mixed graphs, where ties or absent arcs change parity;
- applied ranking data, where pairwise comparisons may be missing or noisy;
- soft tournament attention matrices, where arcs are weights rather than
  binary decisions;
- engineering uses of `tournament_tda.py`, where input data may only be a
  tournament after thresholding.

## Evidence

- The older "what each piece really represents" reflection identifies the
  parabolic law as completeness made quantitative.
- Good-cut/SCC, OCF parity, and endpoint-transfer parity all rely on complete
  binary pair decisions before quotienting or decomposition.
- Trienerment notes show that allowing ties can make forbidden tournament
  values such as `7` and `21` achievable, so the tournament obstruction is not
  merely arithmetic; it depends on completeness.

## Predictions

1. Small completeness defect should preserve low-order residue/phase features
   under thresholding more reliably than raw `H`.
2. Tied or incomplete structures that realize tournament-forbidden values
   should have large localized completeness defect around the forcing cycles.
3. In soft attention tournaments, training changes in tournament structure
   should be visible first as changes in completeness-defect features before
   hard `H` or homology features stabilize.

## Test Plan

1. Define a `completeness_defect` feature for partial, tied, or weighted
   pairwise data.
2. Test it on trienerment examples realizing `H=7` and `H=21`.
3. Test threshold stability: perturb a weighted matrix, threshold into
   tournaments, and compare variance of `H`, residue rank, and phase features
   against the completeness defect.

## Sources

- `07-reflections/what-each-piece-really-represents.md`
- `07-reflections/residue-phase-incidence-synthesis.md`
- trienerment messages and writeups in `06-writeups/trienerments.tex`.
