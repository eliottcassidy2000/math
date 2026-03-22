# Why π — The Circle Inside the Tournament

**opus-2026-03-22-S197**

## The Question

The fiber fraction f(n) = C(2k,k)/4^k ≈ 1/√(πk) measures the probability that a random arc flip preserves the tournament's score sequence. Why does π — the ratio of a circle's circumference to its diameter — appear here?

## The Answer

π appears in tournaments for the same reason it appears in circles: **both involve sums of squares in a flat space**.

### The chain

1. A tournament is a set of **independent** binary choices (arcs)
2. A score is a **sum** of these independent choices
3. An arc flip changes one score by ±1 — a **random walk step**
4. The probability of returning to the same score is C(2k,k)/4^k
5. By the **CLT**, the walk converges to a Gaussian N(0,σ²)
6. The Gaussian normalizes with ∫e^{-x²}dx = √π
7. This integral is computed via **polar coordinates** on x²+y²=r²
8. x²+y²=r² is the equation of a **circle**
9. The angular integral ∫dθ = 2π **IS** the circle

### The common root

Both circles and tournaments reduce to the quadratic form x²+y²:

- **Circle**: x²+y²=1 defines the shape. Circumference = 2π by measuring along it.
- **Gaussian**: e^{-(x²+y²)} has circular level curves. The angular integral gives 2π.
- **Tournament**: Var(ΣXᵢ) = ΣVar(Xᵢ) = Σ1 = k. Additive variances ARE the Pythagorean theorem. CLT gives Gaussian, which gives π.

### The deepest point

**π is the signature of independence.** When random variables are independent:
- Variances add (probabilistic Pythagorean theorem)
- The CLT applies (sums → Gaussian)
- The Gaussian normalization involves √π (circular symmetry of e^{-r²})

If arc flips were correlated, variances would not add, the CLT would not apply, and π would not appear. The fiber fraction would be something else entirely.

**π is not about circles.** π is about the quadratic form x²+y². Circles happen to use this form for distances. Gaussians happen to use it for exponents. Tournaments happen to use it through variance additivity. All three are the same computation viewed from different angles.

## The four faces of π in tournaments

1. **Geometry**: f(n) = (1/π)∫₀^π sin^{2k}(t)dt — sine traces the unit circle
2. **Probability**: f(n) ~ 1/√(πk) — CLT Gaussian normalization
3. **Analysis**: ∫e^{-x²}dx = √π — polar coordinates on the 2D Gaussian
4. **Information**: structural bits ~ n^{3/2}/(2√π) — circular symmetry limits score capture

All four reduce to the angular integral ∫₀^{2π}dθ = 2π around the circle x²+y²=r².

## What this means

Tournament space computes π through its fiber structure. The Wallis product

π/2 = (2/1)·(2/3)·(4/3)·(4/5)·(6/5)·(6/7)·...

is exactly the reciprocal of the Pochhammer ladder:

f(n+1)/f(n) = (2n-3)/(2n-2)

Each time we add a vertex to the tournament, the fiber fraction thins by one step of the Wallis product. The infinite product converges to π because it IS the angular integral, decomposed into sine moment ratios.

This is not a metaphor. It is a theorem.
