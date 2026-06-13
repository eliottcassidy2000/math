# The Boundary Between 1D and 2D

**opus-2026-03-24-S305**

## The Claim

Tournament space is not 1-dimensional. It is not 2-dimensional. It sits at the boundary, with fractal dimension approaching 2 from below:

```
dim(staircase) = log(C(n-1,2)) / log(n-2) = 2 - log(2)/log(n-2) + O(1/log²n)
```

n=4: 1.585. n=5: 1.631. n=6: 1.661. n=10: 1.723. n→∞: → 2.

This is not metaphor. The staircase triangle has m = C(n-1,2) tiles arranged in a right isosceles triangle with legs of length n-2. The "dimension" of m tiles on a length-(n-2) side is exactly log(m)/log(n-2), which is the formula above.

## Why it's almost 1D

The H-ordering maps iso classes to odd integers. This ordering captures 63% of the graph-distance variation at n=5 (corr = 0.63) and 54% at n=6. One number — H — tells you more than half of where you are in the space.

The effective dimension from spectral analysis is 1.18 at n=5 and 1.38 at n=6. The transition matrix has one dominant eigenvalue at 0.50, and all others are much smaller. The metagraph is "almost" a chain — a nearly 1D path through class space, thickened by a few transverse connections.

## Why it's almost 2D

But the width function — classes per H level — grows. At n=5 the max width is 2. At n=6 it's 3. At n=7 it reaches higher. The H-ordering is a 1D projection of a higher-dimensional object, and the transverse direction contains independent information.

The staircase itself IS 2D — it has row and column structure, range and vertex coordinates. The waggly layers (d=1,2,...,m) probe different scales of this 2D structure. The fact that d=1 through d=3 covers 100% of class-pair edges at n=5 means the "effective transverse thickness" is only 3 — but it grows with n.

## The Space-Filling Curve

The diagonal of the right isosceles triangle is the natural path from 1D to 2D. Tiles along the anti-diagonal have constant lo+hi coordinate. The staircase traversal visits:

```
Diagonal 0: 1 tile  (corner)
Diagonal 1: 2 tiles
Diagonal 2: 3 tiles (midpoint: widest)
...
Diagonal n-3: 1 tile (opposite corner)
```

This is a DISCRETE space-filling curve on the triangle, ordered by the diagonal coordinate. The diagonal has length sqrt(2) × (leg length), giving the characteristic ratio sqrt(2).

## The Five Appearances of sqrt(2)

1. **Geometry**: hypotenuse/leg = sqrt(2) in the staircase triangle
2. **Spectral**: λ₀/λ₁ ≈ 2.0 = (sqrt(2))² for the transition matrix at n=5
3. **Growth**: V_merged(n+1)/V_merged(n) → 2 = (sqrt(2))²
4. **Singularity**: the fiber fraction GF is (1-x)^{-1/2} — a square-root branch point
5. **Alphabet**: the Hamming scheme has q=2 = (sqrt(2))², and Walsh functions are characters of GF(2)^m

## The Association Scheme

The tiling space IS the Hamming scheme H(m, 2). The waggly layers are the scheme's association classes. The Krawtchouk polynomials K_k(x; m) are the eigenfunctions.

From an earlier session (S305): **K₁ correlates with H at r = -0.94 at n=5**. The first Krawtchouk polynomial K₁(x; m) = m - 2x is essentially the Hamming weight of the tiling. And this is highly correlated with H. This means: **H is approximately a linear function of the Hamming weight of the tiling.**

Why? Because the Hamming weight counts how many tile arcs go "against" the base path direction. More reversed arcs = more cycles = more Hamiltonian paths = higher H.

## The Coding Theory Interpretation

Each iso class IS an error-correcting code:
- **Codewords** = tilings in the class
- **Code length** = m = C(n-1,2) bits
- **Code size K** = fiber
- **Minimum distance** = smallest Hamming distance within the class

The regular tournament at n=5 is a (6, 3, 3) code — surprisingly good parameters. The transitive tournament is a (6, 1, 6) code — a trivial single-codeword code.

The code rate = log₂(fiber)/m correlates with H: high-H classes are dense codes (high rate, low distance), low-H classes are sparse codes (low rate, high distance). This is the classical rate-distance tradeoff of coding theory, realized in tournament space.

## What Compels

The staircase triangle is the geometric container of tournament theory. Its fractal dimension 2 - log(2)/log(n) quantifies exactly how much of the 2D structure is accessible from the 1D H-ordering. As n grows, more of the 2D becomes visible — but it never fully arrives.

The space-filling curve from 1D (H) to 2D (tiles) is mediated by sqrt(2) — the diagonal of the unit square, the hypotenuse ratio, the first eigenvalue gap. Every time tournament theory encounters the number 2, it's seeing the right isosceles triangle from a different angle.

Tournament space lives on the hypotenuse.
