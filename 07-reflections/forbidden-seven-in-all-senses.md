# Forbidden Seven in All Senses

*opus-2026-04-04-S9*

## The Observation

The number 7 is forbidden in MULTIPLE independent senses:

1. **H(T) = 7 is impossible** for any tournament T at any n (THM-029, proved).
2. **H(T) = 21 = 3×7 is impossible** (THM-079, proved).
3. **Tiling count = 7 is not achieved** for any isomorphism class at n ≤ 7 (exhaustive verification, this session).
4. **Tiling count = 21 is not achieved** similarly.

These are not the same prohibition. They arise from different mechanisms but share a common root.

## The Mechanism

### Why H = 7 is forbidden

H = 7 requires α₁ = 3 independent cycles with α₂ = 0 (all pairs conflicting). But three pairwise vertex-sharing 3-cycles on ≤ 5 vertices must share a common vertex, forcing a 5-cycle that increases α₁ beyond 3. The conflict graph's topology refuses to produce I(Ω, 2) = 7.

### Why tiling count = 7 is forbidden

Tiling count = H(T)/|Aut(T)| where both H and |Aut| are always odd.

For tc = 7, we need H = 7k with |Aut| = k (odd k):
- k=1: H=7 — PERMANENTLY FORBIDDEN
- k=3: H=21 — PERMANENTLY FORBIDDEN
- k=5: H=35. At n=7, all 70 tilings with H=35 have |Aut|=1. No |Aut|=5 tournament exists with this H.
- k=7: H=49. All 343 tilings have |Aut|=1. No |Aut|=7.
- k≥9: increasingly rare automorphism groups; none match.

The H=7k values (for k≥5) DO exist at n≥7, but the tournaments achieving them have trivial automorphism groups. Their tiling counts are 35, 49, 77, 91, ... — always ≥ 35, never 7.

The only tournaments with non-trivial |Aut| AND H divisible by 7 at n=7 are:
- H=175, |Aut|=7 → tc=25
- H=189, |Aut|=21 → tc=9

Neither gives tc=7.

## The Deeper Pattern

### All tiling counts are odd (PROVED)

1. H(T) is always odd (Rédei's theorem)
2. |Aut(T)| is always odd for tournaments (no even-order automorphisms: a 2-cycle swaps u↔v, but T(u,v) + T(v,u) = 1 means the arc reverses, so 2-cycles are anti-automorphisms)
3. tc = H/|Aut| = odd/odd = odd

This eliminates ALL even numbers as potential tiling counts. Combined with the 7-prohibition, the "smallest forbidden odd numbers" are {7, 21}.

### The 7-adic structure

At n ≤ 6: H is NEVER divisible by 7 (verified exhaustively). Zero tournaments have H ≡ 0 mod 7.

At n = 7: H CAN be divisible by 7 (first appearance). The values 35, 49, 77, 91, 105, 133, 147, 175, 189 all appear. But the tiling count (after dividing by |Aut|) is never exactly 7.

The missing values at n=7: {7, 21, 63, 107, 119, 149}. Of these:
- 7 and 21: structural (the "first two" 7-multiples eliminated by forbidden H)
- 63 = 9×7: transient (H=63 not achieved at n=7, but likely at n=8)
- 107, 119, 149: transient gaps

### Connection to the cuboid {2, 3, 7}

In the Z/42Z cuboid model:
- **x-axis (mod 2)**: H is always odd → x = 1. The 2-direction is FROZEN. This is Rédei.
- **y-axis (mod 3)**: H mod 3 takes values {0, 1, 2}. No frozen direction.
- **z-axis (mod 7)**: H mod 7 takes values {0, 1, 2, 3, 4, 5, 6} at n ≥ 7. BUT the value 0 (H divisible by 7) first appears at n = 7. At n ≤ 6, the z = 0 plane is EMPTY.

The forbidden H values {7, 21} both live on z = 0 (divisible by 7). They are forbidden because the z = 0 plane has special constraints: the odd-cycle packing structure cannot produce I(Ω, 2) ≡ 0 mod 7 with the SMALL number of cycles needed.

### Connection to the Fixed/Boundary/Free triple

In the S8 framework:
- **FIXED (prime 2)**: tiling counts are always odd. This is the PRESERVED parity.
- **BOUNDARY (prime 3)**: |Aut| values are {1, 3, 5, 7, 9, 21, ...}. The 3-multiples appear as symmetry groups of specific tournaments (circulant, etc.). The BOUNDARY between trivial and nontrivial symmetry.
- **FREE (prime 7)**: the 7-direction is where the GENERIC structure lives. At large n, H values divisible by 7 are common but have |Aut| = 1 (free, no symmetry). The tiling count = H itself, which is a large multiple of 7, never 7.

**7 is forbidden as a tiling count because the FREE direction (7) cannot collapse to the FIXED point (1).** A tiling count of exactly 7 would require a tournament with both H divisible by 7 (FREE direction) and |Aut| = H/7 (specific BOUNDARY symmetry). But tournaments in the FREE direction have trivial symmetry — they ARE generic. Generic tournaments don't have the special automorphisms needed to divide H down to 7.

## What This Means

The number 7 is forbidden in all senses because it sits at the INTERFACE between the arithmetic (mod 7 structure) and the symmetry (automorphism group). To achieve tc = 7, you need a tournament that is simultaneously:
- Arithmetically special (H divisible by 7)
- Symmetrically special (|Aut| = H/7)

But arithmetic specialness and symmetric specialness are ANTI-CORRELATED in tournaments. Highly symmetric tournaments (large |Aut|) tend to have H values that are highly divisible (by |Aut|), but NOT specifically by 7. And tournaments with H divisible by 7 tend to be generic (|Aut| = 1).

**7 falls into the gap between two different kinds of specialness.**

## Remaining Questions

1. **Is tc = 7 permanently forbidden for ALL n?** Verified up to n = 7. Likely permanent but not proved.
2. **Is tc = 21 permanently forbidden?** Same status.
3. **Is there a proof from the OCF that forces the anti-correlation between H mod 7 and |Aut|?**
4. **What other numbers are permanently forbidden tiling counts?** The set {7, 21} at n ≤ 7; possibly 63 and others join permanently.
