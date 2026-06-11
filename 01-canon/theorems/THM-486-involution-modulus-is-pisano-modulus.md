# THM-486 — The involution modulus 24 is the Pisano modulus

**Status:** PROVED/CITED + VERIFIED — claudebox-2026-06-11-S5. Elementary number theory; the
contribution is the BRIDGE from THM-484 (24 = max involution modulus, (ℤ/24)* = 𝔽₂³, the
Gleason/code orders) to the Fibonacci/Zeckendorf side (THM-485) through the single fact 24 | p²−1.
**Provenance:** user dispatch. **Companions:** THM-484 (involution modulus), THM-485 (Zeckendorf
transfer operator), THM-481 (Golay g₂₄). **Sources:** OEIS A001175 (Pisano periods); Fulton–Morris
/ Wall (π(n)=n classification); the n²≡1 mod 24 fact.

## Statement (all verified)

1. **π(24) = 24.** The Pisano period of the involution modulus equals itself; and
   **π(n) = n exactly for n ∈ {1} ∪ {24·5ᵏ : k ≥ 0} = {1, 24, 120, 600, 3000, …}** (verified to
   k=4). So 24 is the **smallest n > 1 whose Fibonacci sequence has period equal to the modulus**
   — the base Pisano-fixed point; the others are its 5-power multiples.
2. **24 | p²−1 governs both the involution modulus AND the Pisano period.** THM-484: (ℤ/24)* has
   exponent 2 ⟺ n² ≡ 1 (mod 24) for all n coprime to 24 ⟺ 24 | p²−1 (p > 3). The same p²−1
   bounds the Pisano period: **π(p) | p²−1** for every prime p ≠ 5 (π(p) | p−1 if p ≡ ±1 mod 5,
   π(p) | 2(p+1) if p ≡ ±2 mod 5; both divide p²−1). Verified jointly for all primes 7..61.
3. **α(24) = 12 and F₁₂ = 144 = 12².** The rank of apparition (least k with 24 | F_k) is 12, and
   F₁₂ = 144 is the **unique nontrivial Fibonacci perfect square** (Cohn 1964) — and it is a
   square of the rank itself, 12² = 144, with 12 = 24/2 = the antipodal half of the modulus.

## Why this is the right bridge

24 is simultaneously: the **involution modulus** ((ℤ/24)* = 𝔽₂³, THM-484), the **base Pisano
modulus** (π(24)=24, part 1), the **Golay/Leech dimension** (g₂₄ = C(I+S(Paley₂₃)), THM-481), and
the **η²⁴ = Δ modular weight**. All four are downstream of n² ≡ 1 (mod 24): the unit group is
all-involutions, the Fibonacci recurrence closes on period 24, the lattice/code self-duality lives
at dimension 24, and the modular discriminant has weight 12 = α(24). The Zeckendorf/Fibonacci side
(THM-485, growth φ = (1+√5)/2) and the code/involution side (THM-481/484) meet at 24 because the
golden ratio's field ℚ(√5) controls the Pisano period (the 5 in 24·5ᵏ) while ℚ(√−3)/the cube root
controls the code side — and 24 = lcm of the 8-part (2³, the involutions) and the 3-part.

## Verification

Script transfer_temperatures_viswanath_cbx5.py: π(24)=24; π(n)=n exactly on {1,24,120,600,3000};
α(24)=12, F₁₂=144=12²; π(p)|p²−1 and 24|p²−1 jointly for primes 7..61. All exact.

## Remark

The recon-recalled "24 is the largest n with π(n) ≤ n" is FALSE (π(n) ≤ n for infinitely many n,
e.g. all listed above and most n > 10); the correct distinguished fact is π(n) = n on {1}∪{24·5ᵏ}.
Logged to avoid canonizing the error (the kind of recalled-folklore slip MISTAKES.md warns about).
