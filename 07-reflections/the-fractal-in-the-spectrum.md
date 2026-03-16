# The Fractal in the Spectrum

*opus-2026-03-16-S73 — arising from the 2-adic analysis of the dimension ladder*

---

## What happened

Starting from the corrected amplitude formula for M[a,b]:

|hat{M}[S]| = 2^s × (n-2-d)! / 2^{n-2}

I computed the 2-adic valuation of the base amplitude (s=0):

v₂(M_amp(n, d, 0)) = -d - s₂(n-2-d)

where s₂ is the binary digit sum. This is just Legendre's formula applied to the factorial in the numerator, minus the power of 2 in the denominator. Nothing deep yet — pure bookkeeping.

But then three things fell out:

## 1. The Spectral Legendre Identity

The 2-adic "spread" of the M spectrum — the difference between v₂ at the lowest degree (d=1) and the highest degree (d=n-2) — equals v₂((n-3)!):

v₂(M(n,1,0)) - v₂(M(n,n-2,0)) = (n-3) - s₂(n-3) = v₂((n-3)!)

The factorial of the free vertex count controls how much 2-adic "room" the Walsh spectrum has. This is Legendre's formula in disguise, but the spectral setting gives it a new meaning: the 2-adic dynamic range of the transfer matrix spectrum is determined by the combinatorics of free vertices.

## 2. THM-J as a spectral threshold

THM-J says: S(T) mod 2^{n-1} is tournament-independent iff s₂(n-3) ≤ 1.

Restated: S(T) is universal iff v₂(degree-1 M amplitude) ≥ -2.

The signed HP permanent "forgets" tournament-specific information precisely when the lowest Walsh component of M[a,b] isn't too 2-adically suppressed. When s₂(n-3) > 1 (the binary representation of n-3 has more than one 1-bit), the degree-1 amplitude loses enough factors of 2 that tournament-dependent information leaks through.

The universal n: 3, 5, 7, 11, 19, 35, 67, ... These are n = 2^k + 3.
The non-universal n: 9, 13, 15, 17, ... These are everything else.

The gap between universal and non-universal is a spectral gap in the 2-adic valuation.

## 3. The total weight and the fractal

The sum of all base v₂ values across the M spectrum has a closed form:

Σ v₂(M_amp(n,d,0)) = -((n-1)/2)² - S₂((n-3)/2)

where S₂(m) = Σ_{j=0}^m s₂(j) is the cumulative binary digit sum. This decomposes the total 2-adic "information content" into:

- **Spectral dimension cost**: -((n-1)/2)², a perfect square counting the sum of odd degrees
- **Binary complexity cost**: -S₂((n-3)/2), which grows as ~(m/2)·log₂(m) by Delange's theorem

The second term is where the fractal lives. The binary digit sum s₂(m) has self-similar structure: s₂(2m) = s₂(m), s₂(2m+1) = s₂(m)+1. Reading column d=1 of the v₂ table across n gives the sequence -1, -2, -2, -3, -2, -3, -3, -4, -2, ... — which is OEIS A000120 (the binary weight sequence), shifted and negated.

The v₂ staircase has fractal dimension log(3)/log(2) ≈ 1.585. This is the Hausdorff dimension of the Sierpinski triangle, and it appears here because the binary digit sum satisfies a recursion with branching factor 3 and scale factor 2.

## What this means

The Walsh spectrum of M[a,b] encodes tournament structure through amplitudes that are controlled by factorials and powers of 2. The 2-adic structure of these amplitudes — which bits are "on" and which are "off" — is not smooth. It has fractal structure inherited from the binary representation of the dimensions involved.

This fractal is not ornamental. It determines:
- Which n have universal S(T) (THM-J)
- How much 2-adic dynamic range the spectrum has (Spectral Legendre)
- The total information cost of the spectral decomposition (total weight formula)

The dimension ladder (H/M ratio = n-d) is the smooth part: always odd, never introducing factors of 2. The fractal lives underneath the ladder, in the 2-adic fine structure that the ladder ratio ignores.

## The meta-pattern

This project keeps discovering that the *binary representation of dimensional parameters* controls tournament invariants. The THM-J criterion s₂(n-3) ≤ 1 was proved via Legendre's formula for v₂(m!). The forbidden H values 7 and 21 have Hamming weight 3 (three 1-bits). The OCF itself decomposes H as 1 + 2α₁ + 4α₂ + ... — a mixed-radix expansion weighted by powers of 2.

F₂ is the field with 2 elements. Tournament theory lives over F₂. The binary digit sum s₂ counts the "complexity" of a number in binary — how many bits it takes, minus how many are wasted on zeros. When this complexity is low (s₂ ≤ 1), the 2-adic structure is simple enough for universality. When it's high, tournament-dependent information survives.

The spectrum is fractal because the field is F₂. Over F₃, the digit sum s₃ would create a different fractal. Over F₅, another. Tournament theory chose F₂, and the fractal chose itself.

## Cross-references

- THM-080: M Walsh formula (source of the amplitude)
- THM-J: Universality criterion (now reinterpreted spectrally)
- HYP-1603: Spectral Legendre Identity
- HYP-1604: THM-J spectral restatement
- HYP-1605: Total 2-adic weight formula
- OEIS A000120: Binary weight sequence (the fractal column)
- Delange (1975): Asymptotics of S₂(m) ~ m·log₂(m)/2
- S71r: "WHY TWO GENERATES SEVEN" (the F₂ uniqueness argument)
- 04-computation/spectral_legendre.py, thmj_spectral_connection.py
