---
id: HYP-8803
title: "A3 cotangent pullback versus A2 nonfiltered quantum descent"
status: >
  RESERVED / OPEN. The externally replayed A3 certificate is exact but does
  not decide DC(2). The proposed experiment is to compare its inverse-Jacobian
  momentum formulas with the symplectic substitution behind THM-2044 and
  identify whether the HYP-8802 plus-six cascade is the Taylor shadow of a
  resumable nonfiltered descent or a genuine no-finite-descent obstruction.
source: codex-2026-07-21-DC2-filtered-pullback-wall
related:
  - THM-1300
  - THM-2044
  - THM-2046
  - HYP-8802
external:
  - https://github.com/techno-optimist/erdos-frontier-atlas/tree/main/certificates/dixmier-conjecture
---

# HYP-8803 -- descend the exact A3 certificate across the filtration wall

This namespace is reserved for a controlled comparison, not for a claimed
DC(2) counterexample.  The linked certificate has been replayed exactly and
constructs `Phi:A_3 -> A_3`; its own README and terminal verdict explicitly say
that `DC(1)` and `DC(2)` are untouched.

The missing experiment is to transport the exact inverse-transpose-Jacobian
momentum fields through the `ell=L+g` symplectic substitution of THM-2044 and
compare the resulting normal-ordering terms, weight by weight, with HYP-8802.
The key falsifiable alternatives are:

1. the weights `m -> m+6` are coefficients of a closed polynomial/nonfiltered
   conjugation and the cascade resums in `A_2`; or
2. every finite-order lift has a nonzero top weight, giving a rigorous
   no-finite-descent theorem for this witness.

The honest gap is the actual transport calculation and a proof covering all
allowed simultaneous corrections of `T,D,S`, rather than only the fixed-`T`
tangent-linear sector already audited in HYP-8802.
