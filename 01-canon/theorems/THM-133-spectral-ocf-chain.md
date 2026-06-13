# THM-133: Spectral-OCF Chain — Paley Maximizes C₅ via Jensen

**Status:** VERIFIED (computationally at p=7, 11)
**Author:** opus-2026-03-12-S60
**Date:** 2026-03-12

## Statement

For any circulant tournament on Z_p (p prime) with connection set S ⊂ Z_p\{0}, |S| = (p-1)/2:

1. **Universal real part:** All eigenvalues satisfy Re(λ_k) = -1/2 for k ≥ 1.

2. **Writing λ_k = -1/2 + iy_k**, the imaginary parts satisfy:
   - Σ_{k=1}^{p-1} y_k² = p(p-1)/4 (universal, from Parseval)
   - y_{p-k} = -y_k (conjugate pairing)

3. **C₅ formula:** The number of directed 5-cycles satisfies
   C₅ = f(p) - (1/2)Σy_k⁴
   where f(p) = [((p-1)/2)⁵ - (p-1)/32 + 5p(p-1)/16] / 5.

4. **Jensen inequality:** Σy_k⁴ ≥ [p(p-1)/4]² / (p-1) = p²(p-1)/16
   with equality iff all y_k² are equal (spectral flatness).

5. **Paley achieves equality:** For p ≡ 3 mod 4, the Paley tournament T_p has |λ_k|² = (p+1)/4 for all k ≥ 1 (from Gauss sums), giving y_k² = p/4, hence Σy_k⁴ = p²(p-1)/16.

**Corollary:** The Paley tournament maximizes C₅ among all circulant tournaments on Z_p (for any prime p).

## Proof

### Step 1: Re(λ_k) = -1/2

For any circulant tournament with connection set S, the complementary set is -S = {p-s : s ∈ S}, and S ∪ (-S) = Z_p\{0} (since S is a tournament). Therefore:

Re(λ_k) = Σ_{s∈S} cos(2πks/p) = (1/2) Σ_{d=1}^{p-1} cos(2πkd/p) = (1/2)(-1) = -1/2. □

### Step 2: Parseval constraint

Σ_{k=0}^{p-1} |λ_k|² = p·|S| = p(p-1)/2 (standard Parseval for circulants).
With |λ_0|² = ((p-1)/2)², we get Σ_{k≥1} |λ_k|² = (p²-1)/4.
Since |λ_k|² = 1/4 + y_k²: Σy_k² = (p²-1)/4 - (p-1)/4 = p(p-1)/4. □

### Step 3: C₅ = Tr(A⁵)/5

In a tournament, all closed walks of length 5 are simple directed cycles.
Proof: If a closed 5-walk revisits a vertex, it decomposes at that vertex into a cycle (length ≥ 3) and a walk (length ≥ 2). Since there are no 2-cycles: cycle length ≥ 3 and remaining ≥ 3, total ≥ 6 > 5. Contradiction. □

### Step 4: Tr(A⁵) expansion

λ_k⁵ = (-1/2 + iy_k)⁵. Taking real part:
Re[(-1/2+iy)⁵] = -1/32 + (5/4)y² - (5/2)y⁴

Summing: Σ_{k≥1} Re[λ_k⁵] = -(p-1)/32 + (5/4)·p(p-1)/4 - (5/2)·Σy_k⁴.

So: Tr(A⁵) = ((p-1)/2)⁵ - (p-1)/32 + 5p(p-1)/16 - (5/2)Σy_k⁴.
And: C₅ = Tr(A⁵)/5 = f(p) - (1/2)Σy_k⁴. □

### Step 5-6: Jensen + Paley

Standard Cauchy-Schwarz gives Σy_k⁴ ≥ (Σy_k²)²/(p-1). The Paley tournament achieves equality by Gauss sum theory. □

## Verification

| p | C₅(Paley) | C₅(flattest other) | Σy⁴(Paley) | Σy⁴(other) | Jensen min |
|---|-----------|-------------------|-------------|-------------|------------|
| 7 | 42 | 28 | 18.375 | 46.375 | 18.375 |
| 11 | 594 | 484 | 75.625 | 295.625 | 75.625 |

## Additional finding: H = 1587/8 - Σy⁴/2 at p=7

At p=7, H(T) itself is an exact linear function of Σy_k⁴:
H = 1587/8 - (1/2)Σy_k⁴

This arises because:
- H = 1 + 2α₁ + 4α₂ where α₁ = C₃ + C₅ + C₇ and α₂ = disjoint pairs
- C₃ = 14 (universal), C₅ = 51.1875 - Σy⁴/2
- 2C₇ + 4α₂ = 535/8 + Σy⁴/2 (empirically exact)
- Net: H = 198.1875 - Σy⁴/2

The partial cancellation (C₅ gives -Σy⁴, but C₇+α₂ gives +Σy⁴/2) reduces the
effective coefficient from -1 to -1/2.

## Limitations

- At p=11, H is determined by σ₂ = Σy⁴ but NOT linearly.
- At p=13, σ₂ alone does NOT determine H; σ₃ = Σy⁶ is also needed.
- The "spectral flatness → max H" argument works at p=7 because σ₂ determines H monotonically. At larger p, the relationship is more complex.
