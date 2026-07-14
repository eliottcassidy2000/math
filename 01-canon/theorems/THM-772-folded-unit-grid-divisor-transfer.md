---
id: THM-772
title: Folded unit-grid divisor transfer in the two- and three-sheet equality packets
status: CLAIMED (number reserved; exact unit-grid proof and independent audit in progress)
source: codex-2026-07-14-S7
depends_on:
  - THM-769
related:
  - THM-523
  - THM-593
  - THM-765
  - HYP-6820
---

# THM-772 — Folded unit-grid divisor transfer

This stub reserves the theorem number while the exact finite shell table is
independently audited.  The intended statement is:

1. If a primitive tight twelve-packet lies in THM-769's two-exception equality
   branch

   ```text
   A=2U union {x,y},       |U|=10,       x,y odd,
   ```

   then `U` contains a multiple of every `m=2,...,12`.  Moreover
   `x,y<=11 max(U)`.

2. In the three-sheet equality branch

   ```text
   A=3U union {x,y,z},     |U|=9,        3 does not divide xyz,
   ```

   the core `U` contains a multiple of every `m=2,...,11`.  If it has no
   multiple of 12, then, after permutation and independent signs,

   ```text
   {x,y,z} mod 36 = {+/-2,+/-10,+/-14},
   ```

   one residue from each signed pair.  Moreover
   `x,y,z<=10 max(U)`.

The proof carrier is the unit-grid obligation hypergraph.  If `U` has no
`m`-multiple for `2<=m<=12`, then every unit fraction `a/m` lies in its strict
`1/13`-loose set.  THM-769 therefore requires the exceptional runners to be
eligible and to own all sheet colours at **every** such fraction.  This is a
finite residue problem modulo `2m` or `3m`; the simultaneous all-unit condition,
not any single residue, supplies the contradiction.  The speed bounds come
from fitting the guaranteed core-safe interval at a maximizer into each
exception's eligibility tooth.

Until the table, strict endpoint convention, and the three-sheet modulus-6/12
splice have passed the independent audit, cite this file only as a claimed
reduction, not as a proved theorem.
