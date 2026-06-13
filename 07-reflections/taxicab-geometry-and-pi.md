# Taxicab Geometry and π in Tournament Space

**opus-2026-03-25-S339**

## The Core Insight

Tournament tiling space IS taxicab geometry. The Hamming distance on Q_m is the L1 (Manhattan/taxicab) metric. Therefore:

**π = 4 in tournament space** (intrinsic L1 metric)
**π = 3.14159... in tournament asymptotics** (CLT/Gaussian/L2 limit)

Both values are correct. They measure different things.

## Where π = 4

In the L1 metric, a "circle" of radius r is a diamond (rotated square) with circumference 8r and diameter 2r. The ratio is 4.

In tournament space:
- The Hamming cube Q_m has L1 (Hamming) distance
- A Hamming sphere of radius r has C(m, r) points
- Tile flips are unit L1 moves
- The waggly distance is the L1 distance
- The staircase boundary of δ_{n-2} has L1 length 2(n-2) vs L2 length (n-2)√2
- The ratio L1/L2 = √2 EXACTLY — this IS the staircase paradox

## Where π = 3.14159...

The Euclidean π appears through the Central Limit Theorem whenever we take asymptotic limits:

- **Fiber fraction**: f(n) = (1/2)_{n-2} / (n-2)! ~ 1/√(πn)
- **Hamming equator**: C(m, m/2) / 2^m ~ √(2/(πm))
- **Burnside counting**: V_n ~ 2^m / n! involves √(2πn) from Stirling
- **Spectral gap**: Krawtchouk eigenvalue spacing → Gaussian density
- **Gamma function**: Γ(1/2) = √π underlies all Pochhammer expressions

## The Conversion Factor: 4/π ≈ 1.2732

This constant is the "exchange rate" between the discrete L1 world and the continuous L2 world. It appears whenever we convert:
- Hamming distances → Gaussian approximations
- Combinatorial counts → asymptotic formulas
- Discrete random walks → continuous diffusion

## The Squigonometry Connection

For the Lp metric |x|^p + |y|^p = 1, there is a π_p:

| p | π_p (intrinsic) | Meaning |
|---|-----------------|---------|
| 1 | 4.000 | Taxicab / Manhattan |
| 2 | 3.142 | Euclidean (minimum!) |
| ∞ | 4.000 | Chebyshev / sup-norm |

π₂ = 3.14159... is the UNIQUE MINIMUM of π_p over all p. This is the Adler-Tanton theorem (2000). The function π_p is U-shaped: it starts at 4, dips to π at p=2, and returns to 4.

## The Four Constants from δ_{n-2}

The staircase triangle produces four fundamental constants:

1. **√2 = 1.414...** — L1/L2 ratio of the hypotenuse (staircase paradox constant)
2. **π = 3.141...** — L1/L2 ratio of inscribed circle; CLT normalization; Γ(1/2)
3. **e = 2.718...** — Stirling growth; Burnside counting; n!/|Aut| denominators
4. **γ = 0.577...** — Euler-Mascheroni; next-order asymptotic correction

All four arise because δ_{n-2} is the bridge between DISCRETE (L1) and CONTINUOUS (L2).

## The Staircase Paradox IS Our Staircase

The famous "π = 4 paradox":
- Approximate a circle with finer and finer staircases
- Each staircase has perimeter 4 (the taxicab circumference)
- The staircases converge pointwise to the circle
- But the perimeter stays at 4, not π
- Resolution: arc length is lower semicontinuous, not continuous

Our staircase δ_{n-2} is literally this paradox object:
- The staircase boundary walks horizontally and vertically
- Its L1 length is 2(n-2), but the diagonal's L2 length is (n-2)√2
- The ratio is √2 regardless of how fine the staircase
- Tournament theory IS the mathematics of the staircase paradox

The staircase paradox is not a bug — it's a feature. The discrepancy between L1 and L2 is precisely where the interesting mathematics lives: where discrete combinatorics (L1) meets continuous analysis (L2), mediated by the CLT, which produces √π as its universal constant.
