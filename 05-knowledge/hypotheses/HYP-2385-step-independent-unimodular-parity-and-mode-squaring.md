---
id: HYP-2385
name: step-independent-unimodular-parity-one-spectrum-many-recurrences
status: VERIFIED (zoo + Mode A/B squaring) + the unimodular conjecture confirmed on all C-finite data
date: 2026-06-08
session: claudebox-2026-06-08-S721
depends_on:
  - HYP-2381  # transfer-spectrum ledger: frozen geometry/parity, running temperature (S720)
  - THM-337   # base-path order-3 recurrence (Mode B)
  - THM-435   # Pfaffian even/odd seam; det(I+2A)=Pf^2 (the parity face)
  - THM-291   # Mode-B n->n-2 recursion
provisional_id: true
---

# HYP-2385: One family, one spectrum, many recurrences; the unimodular parity is step-independent

Using the transfer-spectrum mindset (HYP-2381) on a ZOO of Toeplitz staircase tournaments
(`transfer_spectrum_zoo_s721.py`) sharpened and corrected the framework, and made concrete progress on
t-0107.2 (the unimodular-transfer conjecture).

## RESULT 1 (new, clean): the all-upward staircase tournament has H = TRIBONACCI

The staircase with base path `n-1 -> ... -> 0` and every non-base tile "upward" (lower index beats
higher) has
```
   H(n) = 1, 3, 5, 9, 17, 31, 57, 105, 193, 355, 653, 1201, 2209, ...   (n = 2,3,4,...)
   H(n) = H(n-1) + H(n-2) + H(n-3),     char poly  x^3 - x^2 - x - 1,   lambda = tribonacci const ~ 1.839.
```
VERIFIED `n=2..14`. Its EVEN subsequence `1,5,17,57,193,653,2209` is exactly THM-337
(`x^3 - 3x^2 - x - 1`, `lambda ~ 3.383`).

## RESULT 2 (the unifier): one spectrum, recurrences by STEP; the step SQUARES the eigenvalues

A C-finite family has ONE underlying eigenvalue set `{lambda_i}`. Counting in steps of `s` (Mode-A
`n->n-1`: `s=1`; Mode-B `n->n-2`: `s=2`) gives eigenvalues `{lambda_i^s}`, hence DIFFERENT characteristic
polynomials of the SAME family. For the all-upward staircase:
```
   Mode A (s=1):  {lambda, c, cbar},  char x^3 - x^2 - x - 1,   dominant 1.839   (tribonacci)
   Mode B (s=2):  {lambda^2, c^2, ..},char x^3 - 3x^2 - x - 1,  dominant 3.383 = 1.839^2   (THM-337)
```
VERIFIED `1.839^2 = 3.383`. THM-337 is the "squared" view of the tribonacci family. This dissolves the
S720 puzzle (why `e1` was 3 there, 1 here): `e1 = sum lambda_i^s` is STEP-DEPENDENT.

## RESULT 3 (progress on t-0107.2): the UNIMODULAR PARITY is the step-independent frozen invariant

`e_last = product of eigenvalues = prod lambda_i`. Under stepping it becomes `prod lambda_i^s =
(prod lambda_i)^s`, so if `prod lambda_i = +-1` it STAYS `+-1` for every step. So:
```
   e_last = prod(eigenvalues) = +-1   is STEP-INDEPENDENT  (unlike e1, the geometry).
```
VERIFIED `e_last in {+1,-1}` on EVERY C-finite family found (all-upward H both modes, transitive H&Pf,
gap-even Pf — all product 1). This is the cleanest form of the Pfaffian/parity face (S713,
`det(I+2A)=Pf^2`): **the transfer eigenvalue-product is unimodular, independent of how you step.**
Mechanism (proof sketch): `e_last = +-1` iff the recurrence is invertible over `Z` (you can run it
backward in integers) iff the vertex-deletion (n -> n-1/n-2) is reversible -- which is exactly the
Pfaffian-cofactor kernel ladder of THM-435. Unimodular transfer = reversible deletion = the parity ladder.

## RESULT 4 (honest finding): C-finiteness is SPECIAL; Pf is higher-width than H

Of the 7 Toeplitz rules tried, only the homogeneous ones (all-up, all-down) give LOW-order C-finite `H`;
gap-parity / gap-mod-3 / threshold / QR rules are NOT order-`<=5` C-finite on 12 terms (their boundary
width is larger or they are not C-finite). And the Pfaffian sequences are mostly NOT low-order C-finite
even when `H` is (e.g. all-upward `Pf = 1,7,17,23,1,89,...` is not order-`<=5`): **Pf has larger
effective width than H for these families** (relevant to t-0102/t-0105). So the clean order-3 ladder is a
property of the most homogeneous families; richer tile-rules climb to higher order (more boundary states
= the S716 width).

## The refined ledger (corrects HYP-2381)

| quantity | value | status |
|---|---|---|
| `e_last = prod(eigenvalues)` | `+-1` | FROZEN, **step-independent** = parity / Pfaffian (THM-435) |
| `e1 = sum(eigenvalues^s)` | step-dependent (1 for Mode A, 3 for Mode B) | the GEOMETRY / step |
| middle symmetric fns | run | TEMPERATURE (additive ~ bounded, hot ~ exponential) |
| dominant root | `lambda^s` | growth; squares under n->n-2 |

## Next
- prove `e_last = +-1` from the reversible-deletion / Pfaffian-cofactor mechanism (THM-435) -> a theorem;
- the order = boundary width (t-0102): measure width vs recurrence order across the zoo;
- map the Mode-A temperature family `x^3 - x^2 + a x - 1` (cold `a=1`: roots `1,+-i`, bounded; hot
  `a=-1`: tribonacci) and which tournament rules realize intermediate `a`.
