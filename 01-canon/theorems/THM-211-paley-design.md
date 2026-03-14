# THM-211: Paley Tournament P₇ Three-Cycles Form a 2-(7,3,2) Design

**Status:** VERIFIED, PROOF NEEDED (likely follows from known QR properties)
**Found by:** opus-2026-03-14-S89b
**Verified in:** `04-computation/sum_h_proof_89b.py`

## Statement

The Paley tournament P₇ (where i→j iff (j-i) mod 7 is a quadratic residue, QR₇ = {1,2,4}) has exactly 14 directed 3-cycles. Their vertex sets form a **2-(7, 3, 2) balanced incomplete block design** (BIBD):

- 7 points = vertices {0,1,...,6}
- 14 blocks = 3-cycle vertex sets
- Each pair of points appears in exactly λ = 2 blocks
- This is precisely **twice the Fano plane** PG(2,2)

## The Resolution Structure

The 7 disjoint pairs of 3-cycles (d₃₃ = 7) form a **resolution** of this design:

| Pair | Triple A | Triple B | Leftover |
|------|----------|----------|----------|
| 1 | {0,1,3} | {2,4,5} | 6 |
| 2 | {0,1,5} | {3,4,6} | 2 |
| 3 | {0,2,3} | {1,5,6} | 4 |
| 4 | {0,2,6} | {1,3,4} | 5 |
| 5 | {0,4,5} | {1,2,6} | 3 |
| 6 | {0,4,6} | {2,3,5} | 1 |
| 7 | {1,2,4} | {3,5,6} | 0 |

Each vertex is the leftover vertex in exactly one pair. This is a **near-resolution** (each "parallel class" covers 6 of 7 points).

## IP Formula for P₇

H(P₇) = 1 + 2·(t₃ + t₅ + t₇) + 4·d₃₃
       = 1 + 2·(14 + 42 + 24) + 4·7
       = 1 + 160 + 28
       = **189**

P₇ is self-complementary: H(P₇) = H(P̄₇) = 189.

## Connection to Projective Geometry

- The 14 blocks as 2× Fano suggest that the **two directions** of each triangle (clockwise/counterclockwise = the two directed 3-cycles on each vertex set) produce the doubling.
- Since the Fano plane has 7 lines and 7 points, and the Paley tournament has 7 vertices, each triple of the Fano plane maps to exactly 2 directed 3-cycles (one for each orientation).
- This means the **Fano plane = undirected 3-cycle structure of P₇**.
