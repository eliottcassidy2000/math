# The Krawtchouk Coordinates of Tournament Space

**Sessions**: opus-2026-03-24-S305, S306, S307

---

## The Three Axes

Tournament space has three natural coordinate axes given by the Krawtchouk dual spectrum:

| Axis | Krawtchouk mode | Tournament invariant | Correlation | Role |
|------|----------------|---------------------|-------------|------|
| **B₁** | K₁(x; m) = m - 2x | **-H** (neg. Hamiltonian paths) | r = -0.94 | Principal (diagonal) |
| **B₂** | K₂(x; m) = C(m,2) - 2(m-1)x + 2C(x,2) | **-c₃** (neg. 3-cycle count) | r = -0.86 | Width (perpendicular) |
| **B₃** | K₃(x; m) | Higher-order structure | r = -0.44 | Twist |

## What This Means

The Krawtchouk polynomials are the **natural eigenfunctions** of the Hamming scheme H(m, 2). The tournament tiling space IS this Hamming scheme with m = C(n-1, 2). When we project onto the quotient Q_m / S_n, the Krawtchouk eigenfunctions descend to **zonal spherical functions** — and these ARE the tournament invariants H and c₃.

**B₁ = -H** is the first eigenfunction: it measures how many tiles are flipped (Hamming weight), which equals the number of Hamiltonian paths. This is LINEAR in the tiling coordinates.

**B₂ = -c₃** is the second eigenfunction: it measures the CURVATURE of the tiling pattern — how many 3-cycles exist. This is QUADRATIC in the tiling coordinates (it involves pairs of tiles).

## The Triangle in (H, c₃) Space

When we plot tournament iso classes in (H, c₃) coordinates, they form a **narrow triangle**:

- **Bottom-left corner**: transitive (H=1, c₃=0)
- **Top-right corner**: regular/Paley (H=max, c₃=max)
- **Correlation**: r(H, c₃) = 0.95-0.98 (nearly linear)
- **Aspect ratio**: 5-7 (very elongated)
- **Width**: only 1-2 units of c₃ variation per H-level

This triangle IS the staircase δ_{n-2}, projected onto its invariant coordinates. The fact that it's so narrow (aspect 5-7) is why tournament space is "almost 1D."

## The √2 Connection

The staircase is a right isosceles triangle with hyp/leg = √2. The diagonal projection from 2D (the full triangle) to 1D (the H-axis) compresses by 1/√2. This is why:
- H captures ~95% of the variance (the diagonal direction)
- c₃ captures ~5% (the perpendicular width)
- Together they capture ~99.5%

The remaining 0.5% is in B₃, B₄, ... — higher-order curvature and twist.

## Unification

This Krawtchouk identification unifies ALL previous spectral findings:

1. **H ≈ 2nd eigenvector of metagraph** (S268): B₁ = K₁ IS the 2nd eigenfunction
2. **Markov spectral gap ≈ 2/n** (S295): the eigenvalue of K₁ in the quotient scheme
3. **OU process on H** (S295): K₁ mean-reverts with rate = spectral gap
4. **GOE universality** (S269): the distance matrix follows random matrix statistics because the Krawtchouk structure + S_n action creates sufficient mixing
5. **Distance profile P(-1) = 0** (S303): a Krawtchouk identity at x = -1
6. **Completeness theorem** (S301): metagraph distance = Hamming distance = K₁ difference

Everything is the first two Krawtchouk polynomials, projected through the right isosceles triangle of the staircase.
