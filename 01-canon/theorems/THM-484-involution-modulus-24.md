# THM-484 — The involution modulus 24: the arithmetic source of the 8/24 Gleason orders

**Status:** PROVED (claudebox-2026-06-11-S4) + VERIFIED. The arithmetic (parts 1-2) is
classical; the contribution is the BRIDGE: identifying 24 = the maximal involution modulus
as the common source of the lengths 8 and 24 of the two Gleason Type II generators, both of
which THM-480/481 realize as Paley tournament gauges, and the additive/multiplicative
grading of those two code families.
**Provenance:** claudebox-2026-06-11-S4, user dispatch (p²≡1 mod 24 + the happy-number
8-cycle puzzle). **Companions:** THM-480 (ê₈ = C(I+S(Paley₇)) = RM(1,3)), THM-481 (g₂₄ =
C(I+S(Paley₂₃)); eQR ladder), THM-482 (doubling thermalizes to d⁺), HYP-2415 ([72,36,16]),
THM-447 (skew-Walsh doubling). **Classical context (cited):** the 24 of η²⁴/Leech/Δ and of
the Gleason ring's second generator W_{g₂₄} (Conway–Sloane SPLAG; Gleason 1970).

## Part 1 — The maximal involution modulus (elementary, verified)

For p coprime to 24, p² ≡ 1 (mod 24): p is coprime to 8 ⟹ p² ≡ 1 (mod 8); p coprime to 3 ⟹
p² ≡ 1 (mod 3); CRT. Equivalently (ℤ/24)* has exponent 2.

**24 is the LARGEST n with (ℤ/n)* of exponent ≤ 2.** Proof: (ℤ/n)* has exponent ≤ 2 iff
every unit is its own inverse iff n | 24. (The 2-part: (ℤ/2^k)* is cyclic·(ℤ/2^{k-2}) for
k≥3, exponent 2 only for 2^k | 8; the p-part for odd p: (ℤ/p^k)* cyclic of order
p^{k-1}(p-1), exponent ≤2 only for p=3, k=1. So n | 8·3 = 24.) Verified: exponent-2 moduli =
{1,2,3,4,6,8,12,24}, max 24.

**Consequence:** φ(24) = 8 and (ℤ/24)* = {1,5,7,11,13,17,19,23} ≅ 𝔽₂³ — the elementary
abelian group of order 8, all eight units involutions.

## Part 2 — The two Gleason generators sit at φ(24) and 24

Gleason's theorem: every Type II (doubly-even self-dual) weight enumerator lies in
ℂ[W_{ê₈}, W_{g₂₄}], generators of degrees **8 and 24** (Molien 1/((1−t⁸)(1−t²⁴)),
Clifford–Weil group order 192 = 8·24). Both lengths are read off the involution modulus:
**8 = φ(24)** and **24** itself. THM-480/481 realize both generators as Paley tournament
gauges:
- ê₈ = C(I+S(Paley₇)) = RM(1,3): its 8 coordinates are the points of 𝔽₂³ ≅ (ℤ/24)* — the
  eight involutions index the first generator.
- g₂₄ = C(I+S(Paley₂₃)): length 24 = the modulus; defining prime **23 ≡ −1 (mod 24)**, the
  antipode unit (the σ:v↦−v involution of the perspective key, at the top of the modulus).

So the full Type II weight-enumerator ring is generated at the two scales {φ(24), 24} of the
maximal involution modulus, both by Paley tournaments whose primes are units of that modulus.

## Part 3 — Additive vs multiplicative grading; the two order-8 groups

The two tournament-gauge code families of THM-480/481 are graded by the two groups of order 8:
- **Doubling / d⁺ family** (THM-480/482): graded by the ELEMENTARY ABELIAN 𝔽₂^m (the skew-
  Walsh/Sylvester character group, THM-447); the n=8 grading group is 𝔽₂³ = (ℤ/24)*.
  ADDITIVE side.
- **Border / eQR family** (THM-481): graded by the CYCLIC ℤ/q (quadratic-residue/Paley
  structure); the n=8 grading group is ℤ/8. MULTIPLICATIVE side.

At length 8 the two families COINCIDE — ê₈ = RM(1,3) = eQR(8) = d₈ is the unique Type II
[8,4] code (verified: C(I+S(Paley₇)) ≅ RM(1,3) by explicit coordinate permutation). They
SPLIT upward: one doubling step sends the tower to d₁₆⁺ (THM-482) while the border family
continues as the arithmetic eQR codes. The cyclic-vs-elementary-abelian fork of the two
order-8 groups IS the repo's additive/multiplicative temperature axis (S720/729) at the
self-dual-code level: crystalline d⁺ (cold) vs arithmetic QR (hot).

## Part 4 — Squaring's two faces; the happy-number puzzle (with honest scope)

The shared operation across the dispatch is SQUARING. It has two faces:
- **Crystalline (mod 24):** the squaring map x ↦ x² collapses all of (ℤ/24)* to the
  identity (p² ≡ 1) — maximal self-duality, every element its own σ-partner. This is the
  arithmetic shadow of "doubly-even self-dual."
- **Chaotic (digit squares):** the base-10 sum-of-squares-of-digits map S(n) = Σ dᵢ² has a
  unique nontrivial attractor, the 8-cycle **4 → 16 → 37 → 58 → 89 → 145 → 42 → 20 → 4**
  (the "unhappy cycle"; the puzzle's unknown start node is **20**). As a directed object it
  is the cyclic 8-cycle = the ℤ/8 grading group of Part 3.

**HONEST SCOPE (flagged coincidence).** The cycle LENGTH 8 = φ(24) is NOT structural: a
base-by-base computation (script involution_modulus_happy_cbx4.py) shows the digit-square
map's nontrivial cycle lengths are 1 (base 2,4), {1,2} (base 3), {1,3} (base 5,11),
{1,2,3} (base 8,9), {1,4} (base 7), {1,8} ONLY at **bases 6 and 10**, {1,2,3,10} (base 12).
So a unique 8-cycle happens to occur at base 6 (= lcm(2,3), the repo's recurring modulus) and
base 10, but is not forced by 24. What IS structural is the cyclic-vs-elementary-abelian
duality of the two order-8 groups (Part 3); the puzzle is a memorable avatar of the cyclic
pole, the involution modulus of the elementary-abelian pole.

## Verification (script involution_modulus_happy_cbx4.py)

Exponent-2 moduli {1,2,3,4,6,8,12,24} (max 24); (ℤ/24)* = 𝔽₂³; ê₈ ≅ RM(1,3) by explicit
permutation; eQR ladder extremal d at 8,24,32,48 (THM-481/HYP-2415); digit-square 8-cycle and
its base-dependence; puzzle answer 20. All exact.

## Remarks

1. **Why 24, the deep version (cited, not proved here).** 24 is the involution modulus, the
   degree of the Gleason generator W_{g₂₄}, the Leech dimension, the exponent of η²⁴ = Δ, and
   the unique n>1 with 1²+…+n² a perfect square (n=24, cannonball 70²). The elementary
   involution-modulus fact (Part 1) is the most arithmetic of these and the one that pins the
   tournament-gauge orders. The trace identity Σμ_j² = n(n−1)/2 (THM-472) is the repo's
   sum-of-squares ledger; whether it meets the cannonball identity at any tournament order is
   open (Σμ² = C(n,2) is never a sum of consecutive squares 1²+…+k² for the flag spectra
   checked — a non-lead, recorded for honesty).
2. **The Paley primes are involution units.** Every prime q ≡ 7 (mod 8) (the eQR ladder)
   reduces to 7 or 23 (mod 24) — both involutions of (ℤ/24)*, with 23 = −1. The QR ladder
   lives inside the unit group of the involution modulus.
