---
id: THM-340
name: SC-composition-formula
status: PROVED
date: 2026-05-28
session: opus-2026-05-28-S3b
---

# THM-340: SC Composition Formula — Q(d,k) = [xᵈ]B(x)ᵏ

## Statement

Define:
- B(x) = Σ_{a≥2} SC(a+1) · xᵃ = x² + 5x³ + 50x⁴ + 903x⁵ + ···
  where SC(n) = #{path-fixed strongly connected tilings on n vertices}

- Q(d,k) = number of labeled arrangements of k non-overlapping SC tiling blocks of total width d

Then:
  **Q(d,k) = [xᵈ] B(x)ᵏ**

Equivalently, let C(x) = B(x)/x² = 1 + 5x + 50x² + 903x³ + ··· = Σ_{n≥0} SC(n+3)xⁿ.
Then:
  **Q(2k+j, k) = [xʲ] C(x)ᵏ**

## Proof Sketch

Q(d,k) counts all ordered compositions d = a₁ + ··· + aₖ with each aᵢ ≥ 2,
weighted by SC(a₁+1)·SC(a₂+1)···SC(aₖ+1).

This is exactly the coefficient of xᵈ in (Σ_{a≥2} SC(a+1)xᵃ)ᵏ = B(x)ᵏ.

Interpretation: to place k SC blocks in a row with total width d, choose
lengths a₁+···+aₖ = d with each aᵢ ≥ 2 (minimum SC block has 2 edges),
and for each partition put SC(aᵢ+1) ordered SC tilings in slot i.

## Verification

Verified computationally for all Q(d,k) with d ≤ 20, k ≤ floor(d/2).

Known Q values (first few rows):
| d\k | 1 | 2 | 3 | 4 | 5 |
|-----|---|---|---|---|---|
| 2 | 1 | — | — | — | — |
| 3 | 5 | — | — | — | — |
| 4 | 50 | 1 | — | — | — |
| 5 | 903 | 10 | — | — | — |
| 6 | 30773 | 125 | 1 | — | — |
| 7 | 2032504 | 2306 | 15 | — | — |
| 8 | 264271477 | 73076 | 225 | 1 | — |
| 9 | 68184627441 | 4463038 | 4334 | 20 | — |
| 10 | 35047197032002 | 552760703 | 130659 | 350 | 1 |

## Diagonal Formulas (proved, see THM-342)

- Q(2k,k) = 1
- Q(2k+1,k) = 5k
- Q(2k+2,k) = 25k(k+3)/2
- Q(2k+3,k) = 903k + 250k(k-1) + 125·C(k,3)
- Q(2k+4,k) = 30773k + 4515k(k-1) + 2500·C(k,2) + 3750·C(k,3) + 625·C(k,4)

## Diagonal GF (proved corollary)

Σ_{k≥1} Q(2k+j, k) zᵏ = [xʲ] z·C(x) / (1 − z·C(x))

The full diagonal GF is 1/(1−z·C(x)) including the k=0 term (= 1).
