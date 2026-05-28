---
id: THM-342
name: Q-triangle-diagonals
status: PROVED
date: 2026-05-28
session: opus-2026-05-28-S3b
---

# THM-342: Q-Triangle Diagonal Formulas

## Setting

The Q-triangle has Q(d,k) = [xᵈ]B(x)ᵏ, the k-th power convolution of B(x) = Σ_{a≥2} SC(a+1)xᵃ.
The "j-th diagonal" consists of all Q(2k+j, k) for k = 1, 2, 3, …

## Main Theorem

Let C(x) = B(x)/x² = Σ_{n≥0} SC(n+3)xⁿ = 1 + 5x + 50x² + 903x³ + ···

Then: **Q(2k+j, k) = [xʲ] C(x)ᵏ**

All diagonal formulas follow from expanding C(x)ᵏ:

### Diagonal j=0 (main diagonal):
  Q(2k, k) = 1 for all k ≥ 1

  Proof: [x⁰]C(x)ᵏ = C(0)ᵏ = 1ᵏ = 1.
  Combinatorially: the only composition of 2k into k parts each = 2 is (2,2,…,2), with SC(3)ᵏ = 1.

### Diagonal j=1:
  Q(2k+1, k) = 5k for all k ≥ 1

  Proof: [x¹]C(x)ᵏ = k·[x¹]C(x)·C(0)^{k-1} = k·5·1 = 5k.
  Combinatorially: choose which of the k blocks gets length 3 (= SC(4)=5), rest get length 2.

### Diagonal j=2:
  Q(2k+2, k) = 25k(k+3)/2 for all k ≥ 1

  = 50k + 25·C(k,2)

  Proof: [x²]C(x)ᵏ = k·50 + C(k,2)·25. From SC(5)=50 and SC(4)²=25.

### Diagonal j=3:
  Q(2k+3, k) = 903k + 250k(k−1) + 125·C(k,3)

  Proof: Three types of compositions of 2k+3 into k parts ≥2 with total excess 3:
  - One part=5: k·SC(6) = 903k
  - One part=4, one part=3: k(k-1)·SC(5)·SC(4) = 250k(k-1)
  - Three parts=3: C(k,3)·SC(4)³ = 125·C(k,3)

### Diagonal j=4:
  Q(2k+4, k) = 30773k + 4515k(k−1) + 2500·C(k,2) + 3750·C(k,3) + 625·C(k,4)

  Using SC(7)=30773, SC(6)·SC(4)=4515, SC(5)²=2500, SC(5)·SC(4)²·3=3750, SC(4)⁴=625.

## General Formula

For any j, Q(2k+j, k) = Σ_{λ⊢j, parts≥1} multinomial(k; mult(λ)) · Π SC(λᵢ+3)

where:
- λ ranges over integer partitions of j with all parts ≥ 1
- mult(λ) = multiplicities of the distinct parts of λ
- multinomial(k; mult) = k! / (m₁! m₂! ··· (k−|λ|)!)  [distributing parts among k positions]
- SC(ℓ+3) is the SC tiling count for a block of width ℓ+2

This is equivalent to Q(2k+j,k) = [xʲ] C(x)ᵏ.

## Diagonal GF

Σ_{k≥0} Q(2k+j, k) zᵏ = [xʲ] (1/(1−z·C(x)))

## Verification

All formulas verified computationally for k = 1..8 (j=0), k = 1..8 (j=1),
k = 1..7 (j=2,3), k = 1..6 (j=4).
