---
id: HYP-1783
status: PARTIALLY-TRUE
session: codex-2026-05-30
---

# HYP-1783: Quotient Gap/Residue Principle

## Statement

Many apparently unrelated tournament invariants are best understood as
support, parity, or residue data of quotient maps.

In this language:

- a **gap** is an empty quotient fiber;
- a **boundary** is a parity-visible collection of quotient fibers;
- a **projection defect** is residue left after a quotient forgets
  coordinates;
- a **transport feature** is where quotient mass moves between nonempty
  fibers.

## Confirmed Core

THM-355 proves the finite-set core in Lean. For any finite quotient
`q : alpha -> beta`, an empty source fiber gives an empty transport row, and an
empty target fiber gives an empty transport column. The formalization includes
both set-level emptiness and cardinal row/column-sum zero statements.

For the concrete good-cut quotient

```text
goodCutBucket : StTiling n -> Fin (n+1),
```

Lean proves that, for `n >= 3`, the image is exactly

```text
{0} union {2,...,n-1}.
```

So bucket `1` and bucket `n` are genuine quotient gaps.

## Open Reading

The broader claim is not yet a theorem: forbidden H-values, endpoint transfer
boundaries, projection defects, and good-cut height transport may all be
instances of one support/residue calculus for quotient towers.

The proposed feedback loop:

1. Identify the quotient map.
2. Compute its fiber support and parity boundary.
3. Use transport rows to distinguish empty, rigid, and leaking fibers.
4. Use projection defects to measure what survives after forgetting a layer.

## Evidence Threads

- THM-347: endpoint transfer sees odd child fibers but not necessarily full
  row rank.
- THM-349: good-cut polynomials count the fibers of a height quotient.
- THM-352/353: transport rows give a checksum for quotient movement.
- THM-355: empty fibers are zero source rows and zero target columns.
- Projection-defect notes: deletion and shadow maps measure surviving residue.

## Next Tests

- Treat forbidden H-values as gaps of the finite quotient `T -> H(T)` at fixed
  `n`, and compare zero transport columns with neighboring H-values.
- Track endpoint transfer and Hamming-mask transport with the same
  source/target gap vocabulary.
- Add `gap_count`, `boundary_parity`, and `residue_mass` feature blocks to
  `tournament_tda.py`.
