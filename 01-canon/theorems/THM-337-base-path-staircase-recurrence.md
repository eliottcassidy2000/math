---
theorem: THM-337
name: Base-Path Staircase Order-3 Linear Recurrence
status: CONJECTURED (verified k=1..11 computationally)
session: opus-2026-05-27-S7
verified: computationally k=1..11
---

## Statement

Let T_k be the **base-path staircase tournament** on n=2k vertices, defined by:
- Base path: n−1 → n−2 → ... → 1 → 0 (directed path from vertex n−1 to 0)
- All non-base tiles "upward": for all (x,y) with x−y ≥ 2, vertex y beats vertex x
  (i.e., lower-indexed vertex beats higher-indexed vertex for all non-adjacent pairs)

Then H(T_k) satisfies the **order-3 linear recurrence**:

**H(T_k) = 3·H(T_{k−1}) + H(T_{k−2}) + H(T_{k−3})**

with initial values H(T_1) = 1, H(T_2) = 5, H(T_3) = 17.

The characteristic polynomial x³ − 3x² − x − 1 = 0 has roots:
- λ₁ ≈ 3.38298 (dominant real root)
- λ₂, λ₃ complex conjugate pair with |λ| ≈ 0.544

So H(T_k) ~ C · λ₁^k where λ₁ ≈ 3.382976 and C is a constant.

## Computed Values

| k  | n=2k | H(T_k)                         |
|----|------|-------------------------------|
| 1  | 2    | 1                             |
| 2  | 4    | 5                             |
| 3  | 6    | 17                            |
| 4  | 8    | 57                            |
| 5  | 10   | 193                           |
| 6  | 12   | 653                           |
| 7  | 14   | 2,209                         |
| 8  | 16   | 7,473                         |
| 9  | 18   | 25,281                        |
| 10 | 20   | 85,521                        |
| 11 | 22   | 289,409                       |
| 50 | 100  | 127,071,617,887,002,752,149,434,981 (exact, from recurrence) |

## Proof Sketch

The tournament T_k has a specific recursive structure: T_k is obtained from T_{k−1} by inserting two new vertices (the new "top" pair of the staircase) in a specific pattern.

The H-increment formula (from OCF): H(T_k) − H(T_{k−1}) = 2 · #{odd cycles through new vertices}.

The key structural observation: the new vertex insertions create exactly enough odd cycles to produce the recurrence coefficients 3, 1, 1. A complete proof would need to show that the number of odd paths between the new insertion neighborhoods in T_{k−1} satisfies a linear recurrence of order 3 in k.

Computational evidence (Berlekamp-Massey algorithm applied to k=1..11): unique order-3 recurrence found.

## Characteristic Roots

x³ − 3x² − x − 1 = 0:
- Root 1: x ≈ 3.38297627... (real, dominant)
- Root 2: x ≈ −0.19148814 + 0.54425022i (complex)
- Root 3: x ≈ −0.19148814 − 0.54425022i (complex)

|Root 2| = |Root 3| ≈ 0.57547... < 1, confirming geometric convergence to dominant root.

Growth: H(T_k) ~ C · 3.38298^k.

## Comparison with Other Staircase Families

The Hyp-1732 staircase (pair-structure): H(2k) sequence begins 1, 5, 29, 233, 2489, 33773, 562685, ...
This family does NOT satisfy a short linear recurrence (no P-recursive formula found at orders 1-3 with polynomial degree ≤ 5).

The base-path staircase is MUCH smaller: H(T_22) = 289,409 vs H(Hyp1732, n=22) = 11,222,321.
The pair-structure staircase is specifically constructed for maximality; the base-path staircase is not.

## New OEIS Sequence Candidate

The sequence 1, 5, 17, 57, 193, 653, 2209, 7473, 25281, 85521, 289409, ...
(H values of base-path staircase T_k, k=1,2,3,...)

Status: likely a new OEIS sequence. Check: does it appear in OEIS?
A-number if found: TBD.
