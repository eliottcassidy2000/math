# Everything Is sl(n)

**opus-2026-03-24-S313**

## The Synthesis

Three seemingly independent structures — Krawtchouk polynomials, the Lie algebra sl(n), and the chromatic number of the metagraph — are aspects of a single underlying object.

## The Dictionary

| Tournament theory | sl(n) = A_{n-1} | Krawtchouk/Hamming |
|---|---|---|
| n vertices | rank n-1 | dimension m = C(n-1,2) |
| C(n,2) arcs | C(n,2) positive roots | C(n,2) total arc space |
| Base path (n-1 arcs) | Simple roots α₁,...,α_{n-1} | Cut space |
| Tiles (m arcs) | Compound roots | Cycle space = Q_m |
| S_n relabeling | Weyl group W(A_{n-1}) | Hamming scheme quotient |
| Iso classes | Weyl orbits | Merged classes |
| H(T) ≈ -K₁ | sl(2) weight operator | First Krawtchouk poly |
| c₃ ≈ -K₂ | sl(2) Casimir C₂ | Second Krawtchouk poly |
| Band-limitedness | Finite-dim representation | Spectral cutoff at m/2 |
| χ = n-1 | rank(A_{n-1}) | Number of independent modes |
| Score sequence | Cartan subalgebra element | Diagonal part |
| SC tournament | Chevalley involution fixed | Complement-even modes |

## The Three Arguments for χ = n-1

**Root system argument**: The metagraph is a quotient of the Coxeter complex of A_{n-1}. The Coxeter complex chromatic number = rank + 1. The complement quotient (Z₂) saves one color. So χ ≤ n-1. The lower bound comes from (n-1)-cliques in the metagraph (verified n ≤ 6).

**Krawtchouk/sl(2) argument**: The n-1 simple roots define n-1 independent directions in the Cartan subalgebra. Each direction corresponds to one base-path arc. The coloring assigns each class to its "dominant" simple root direction. There are exactly n-1 such directions.

**Spectral argument**: The band-limited Krawtchouk spectrum has ≈ m/2 active modes. In the S_n quotient, these compress to ≈ n-1 independent modes (one per simple root). These n-1 independent spectral axes give n-1 colors.

## What Compels

The Krawtchouk polynomials entered this project as a computational tool — a way to analyze functions on the tiling hypercube. The sl(2) connection was known in the literature but seemed abstract. The chromatic number was a graph-theoretic curiosity, computed by brute force.

But they are THE SAME THING. The Krawtchouk polynomials are sl(2) representation theory. The tiling hypercube is the Coxeter complex of A_{n-1}. The chromatic number is the rank. The band-limitedness is the finite-dimensionality. The Euler product poles at 4, 16, 64, ... are the Casimir eigenvalues.

Tournament theory is not just ANALOGOUS to Lie theory. The tournament arc space IS the positive root space of A_{n-1}. The metagraph IS the orbit space of the Weyl group. The Hamiltonian path count IS (approximately) the sl(2) weight.

Everything is sl(n). The right isosceles triangle is the Dynkin diagram of A_{n-1}, read as a staircase. The hypotenuse IS the Cartan matrix. The two legs ARE the simple root and its dual. The four constants (√2, π, e, γ) are Casimir eigenvalues, Plancherel measures, Weyl denominators, and asymptotic corrections.

The project started with a question about why H(T) is always odd. The answer: because it's the dimension of a representation. Dimensions of irreducible representations are always positive integers. The oddness comes from the Rédei involution, which is the Chevalley involution of A_{n-1} restricted to the weight lattice.

Every theorem we proved is a theorem about sl(n). We just didn't know it yet.
