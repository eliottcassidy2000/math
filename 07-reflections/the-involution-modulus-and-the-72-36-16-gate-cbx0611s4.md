# The involution modulus, squaring's two faces, and the [72,36,16] gate (claudebox-2026-06-11-S4)

Dispatch: continue the Gleason/d⁺ line, merge in p²≡1 (mod 24), and use a famous puzzle —
the 8-cycle 4,16,37,58,89,145,42,(?) returning to the start.

## The puzzle

It is the base-10 sum-of-squares-of-digits map (happy-number dynamics); the unique nontrivial
attractor is the 8-cycle 4 → 16 → 37 → 58 → 89 → 145 → 42 → **20** → 4, so the unknown start
node is **20**. The 8-cycle plus "p² ≡ 1 mod 24" both point at 8 and 24 — which is exactly the
pair of orders where THM-480/481 found the two Gleason Type II generators (ê₈ and g₂₄) as
Paley tournament gauges. That is not an accident at the level of 8 and 24; it IS one at the
level of the cycle's length.

## What is real (THM-484)

1. **24 is the maximal involution modulus**: (ℤ/n)* has exponent ≤ 2 iff n | 24, so the
   largest is 24, where (ℤ/24)* ≅ 𝔽₂³ (all eight units are involutions). Elementary, proved.
2. **The two Gleason generators live at φ(24)=8 and 24.** ê₈ = RM(1,3) is indexed by 𝔽₂³ =
   (ℤ/24)* (the eight involutions are the eight coordinates); g₂₄ sits at length 24 with
   defining Paley prime 23 ≡ −1 (mod 24), the antipode unit = the σ involution of the
   perspective key. Both generators are Paley tournament gauges (THM-480/481), so the whole
   Type II weight-enumerator ring is generated at the two scales of the involution modulus.
3. **The two code families are the two order-8 groups.** Doubling/d⁺ is 𝔽₂^m-graded
   (Walsh/Sylvester, additive, elementary-abelian); border/eQR is ℤ/q-graded (QR,
   multiplicative, cyclic). They coincide at ê₈ = RM(1,3) = eQR(8) = d₈ (unique Type II [8,4])
   and split upward (THM-482: doubling → d₁₆⁺). The cyclic-vs-elementary-abelian fork of the
   two order-8 groups is the additive/multiplicative temperature axis (S720/729) at the
   self-dual-code level — and the puzzle's cyclic 8-cycle vs the involution group's 𝔽₂³ are
   the two poles made concrete.
4. **Squaring is the shared operation, with two faces**: crystalline mod 24 (x² ≡ 1, maximal
   self-duality, the doubly-even shadow) and chaotic on digits (the unique 8-cycle attractor).

## What is honest coincidence (flagged, not buried)

The cycle LENGTH 8 = φ(24) is NOT structural. A base-by-base computation shows the digit-square
map's nontrivial cycle lengths run 1, {1,2}, {1,3}, {1,4}, {1,2,3}, {1,2,3,10}, … and a unique
8-cycle occurs ONLY at bases 6 and 10 (base 6 = lcm(2,3), the repo's recurring modulus — itself
a wink, not a theorem). So the puzzle is a memorable avatar of the cyclic order-8 pole, not a
derivation of 24. Saying so is the point: the discipline is to keep the genuine bridge
(involution modulus → Gleason orders → tournament gauges) separate from the evocative rhyme.

## The real progress (OPEN-Q-061 / HYP-2415): a famous open problem, reframed

Tracking exactly WHERE the eQR gauge ladder is extremal forced contact with a famous problem.
THM-481's codes are extremal Type II at q = 7, 23, 31, 47 (lengths 8, 24, 32, 48; minimum
distances 4, 8, 8, 12 = the Gleason bound 4⌊n/24⌋+4, all verified) and FIRST FALL SHORT at
**q = 71**: eQR(72) has d = 12 < the extremal 16. The extremal Type II **[72,36,16]** is one of
the most famous open problems in coding theory (existence still unknown, 2026). Since order
72 ≡ 8 (mod 16), the tournament gauge of EVERY skew-Hadamard matrix of order 72 is a Type II
[72,36] code — so a sufficient route to the famous code is: find an order-72 skew-Hadamard whose
gauge minimum distance is 16. Paley (the most symmetric one) gives only 12; thousands of other
order-72 skew-Hadamard matrices exist (Đoković–Kotsireas). The tournament-gauge line thus lands
on a real, named open problem with a concrete computational handle (t-0120), and a sharp
conjecture: the very symmetry that makes a tournament's gauge code computable caps its minimum
distance, so the extremal code — if it exists — is asymmetric.

## Honest scope

Parts 1-2 of THM-484 are classical arithmetic (the bridge/framing is the contribution). The
[72,36,16] reframing is a SUFFICIENT route, not an equivalence (an extremal code need not be a
gauge). All distances ≤ 48 verified by exact enumeration; d=12 at 72 cited (QR tables). The
cannonball/η²⁴/Leech "deep why 24" is cited context, not claimed.
