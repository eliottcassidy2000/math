# Tournament Space at the 1D/2D Boundary

**Session:** kind-pasteur-2026-03-24-S20fv
**Date:** 2026-03-24

## The Almost-1D Character

Tournament space is "effectively one-dimensional" through the H function:
- H correlates 0.87-0.91 with the 2nd adjacency eigenvector
- The metagraph has a strong H-gradient but is NOT a strict DAG (level edges at n≥5, H-decreasing edges at n≥7; see MISTAKE-035)
- The Markov chain approximation P₁² ≈ 0.94·P₂ treats the space as nearly 1D

But it's NOT truly 1D:
- Eigenvalue decay alpha = 2.3 (n=5) drops to 0.5 (n=6) — dimension increases
- 19 distinct H-levels at n=6 vs 56 classes — many classes share H values
- The "perpendicular" directions carry 87.6% of the leading eigenvector variance

## The Seesaw Conservation Law

The opposing trends satisfy a SEESAW analogous to beta_1 × beta_3 = 0:

| n | Residual (R) | 1-P₂ (D) | R × D | R + D |
|---|---|---|---|---|
| 4 | 0.106 | 0.372 | 0.039 | 0.478 |
| 5 | 0.229 | 0.166 | 0.038 | 0.395 |
| 6 | 0.318 | 0.098 | 0.031 | 0.416 |
| 7 | 0.355 | 0.061 | 0.022 | 0.416 |

**Product → 0**: The system can't be both very broken AND Markov-inaccurate.
**Sum ≈ 0.42**: A conservation law for total regularity breaking.

## The Triangle as Moduli Space

The deepest connection: the staircase delta_{n-2} IS the **moduli space of tournament structures**:
- Dimension m = C(n-1,2) = correct parameter count
- Symmetry group S_n acts by coordinate permutation
- The quotient Q_m/S_n has coherent configuration structure
- The "almost 1D" character = this moduli space is "almost a curve"
- The fiber fraction GF (1-x)^{-1/2} is the monodromy around the SC locus

H is an attempt to parameterize this moduli space by a single number, like the j-invariant parameterizes elliptic curves. It succeeds at 95% but fails at the remaining 5% — the "genus" of tournament space grows with n.

## The 7 Testable Predictions

1. **Width max** of metagraph grows as C(floor((n-2)/2), floor((n-2)/4))
2. **Hausdorff dimension** d_H grows as n²/(4 log n), NOT constant
3. **Homological uncertainty**: (var of Betti complexity) × (1 - K₁ correlation) ≈ constant
4. **Plotkin-optimal codes** = SC spine classes
5. **Monotone Gray code** exists at n=4, not at n=5 (H=7 gap blocks it)
6. **Spectral gap** ≥ 2/(n-2)² (Cheeger bound from staircase geometry)
7. **Deletion covering** becomes more regular (CV decreases with n)

## The sqrt(2) Connection

The pseudo-doubling ratio (2n-5)/(n-2) → 2 = (√2)² connects:
- The hypotenuse/leg ratio of the isosceles right triangle
- The Mode A/B reduction rates
- The fiber fraction singularity at x=1
- The binary alphabet size (2 states per tile)
- The Hamming scheme H(m,**2**)

Tournament space lives at the boundary between 1D and 2D, and √2 is the irrational number that measures this boundary — it's the diagonal of the unit square, the simplest path from 1D to 2D.
