# THM-130: Complete Betti Numbers of Paley Tournaments

**Status:** VERIFIED (P_7, P_11) — CONJECTURED for general p

## Statement

For the Paley tournament P_p (p ≡ 3 mod 4 prime), with m = (p-1)/2:

The GLMY path homology has:
- β_d(P_p) = 0 for d ∉ {0, m, m+1}
- β_0 = 1
- **β_m = m(m-3)/2** (diagonals of regular m-gon)
- **β_{m+1} = m(m+1)/2 = C(m+1, 2)** (triangular number)

Consequences:
- χ(P_p) = 1 - m(m-3)/2 + m(m+1)/2 = 1 + 2m = **p**
- β_m + β_{m+1} = m(m-1)
- β_{m+1} - β_m = 2m = p - 1

## Per-eigenspace decomposition

Z_p acts on P_p by cyclic relabeling. The chain complex decomposes into p eigenspaces (k = 0, 1, ..., p-1):

**k = 0 eigenspace (Z_m-invariant):**
- β_m^{(0)} = m(m-3)/2 = m · β_m^{orb}
- β_{m+1}^{(0)} = m(m-3)/2 (from budget equality β_m^{(0)} = β_{m+1}^{(0)})
- All other β_d^{(0)} = 0 except β_0^{(0)} = 1

**k ≠ 0 eigenspaces (p-1 = 2m copies):**
- β_{m+1}^{(k)} = 1
- All other β_d^{(k)} = 0

**Orbit complex (Z_m = QR quotient of k=0):**
- β_m^{orb} = (m-3)/2 (diagonal types of m-gon)
- β_{m+1}^{orb} = (m-3)/2

## Mechanism (HYP-710 Rank Shift)

The face-0 boundary phase creates a rank-1 shift: R_d^{(k)} - R_d^{(0)} = (-1)^{d+1} for 1 ≤ d ≤ m.

At d = m+1: R_{m+1}^{(k)} - R_{m+1}^{(0)} = β_m^{(0)} - 1.
At d ≥ m+2: R_d^{(k)} = R_d^{(0)} (shift vanishes).

This automatically gives β_{m+1}^{(k)} = β_{m+1}^{(0)} - β_m^{(0)} + 1 = 1.

The ONLY free parameter is β_m^{(0)} = m(m-3)/2, which determines everything.

## Verification

| p | m | β_m | β_{m+1} | χ |
|---|---|-----|---------|---|
| 7 | 3 | 0 | 6 | 7 |
| 11 | 5 | 5 | 15 | 11 |
| 19 | 9 | 27(?) | 45(?) | 19(?) |

P_7 and P_11 fully verified computationally. P_19 predicted.

## Open question

WHY is β_m^{orb} = (m-3)/2? This is the fundamental parameter.
Combinatorial interpretation: diagonal types of regular m-gon.
Algebraic proof: unknown.

opus-2026-03-13-S71b
