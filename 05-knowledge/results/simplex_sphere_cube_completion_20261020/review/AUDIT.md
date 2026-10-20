# Independent scoped audit: geometric character completion

Reviewed 2026-10-20 UTC. **PASS: no blocking mathematical defect found in the scoped claims.** This verdict applies to the frozen 3,708-check version, including §5.3 and the two parent-approved dimensional/polarity wording corrections below. It is a human-readable proof audit plus bounded independent exact controls, not formal verification or a new priority claim.

## Reviewed identity

Main packet SHA-256:

    d3e800332a5bad27ba3918bf17d626e989f9356ccef0e9400101773c328b0b15

Primary verifier SHA-256:

    75fc4e2c0c6d12ce215d4a6833becfb94670668958d3fa017dfa5ba337f1672e

Both primary receipt files have SHA-256:

    ed849ab3ff9f6c27b1ccfaa832b2455e2a4db6792f262ec3d788605860daeae9

All four supplied manifest entries match. The two supplied receipts agree byte-for-byte, and their counters sum to 3,708. The full supplied verifier was inspected but not rerun. An earlier initial read preceded the author's addition of §5.3; the old 3,400-check version is not the version certified here. `REVIEWED_INPUTS.sha256` records the final source and inherited dependencies.

## Mathematical findings

1. **Fisher normalization and dimension: correct.** For the stated metric, the map is exactly `p -> 2 sqrt(p)`, with radius 2 and pullback `sum(u_i v_i/p_i)`. The positive great-circle arc yields `2 arccos(sum sqrt(p_i q_i))`. Its completed diameter is π, attained exactly at disjoint supports. A three-atom law gives a two-dimensional spherical triangle on S²; a two-atom law gives an arc on S¹. The interior metric and its boundary completion are distinguished correctly.

2. **Minimal-observer theorem: correct.** The recorded constant ensures every kernel vector has coordinate sum zero, so both sufficiently small uniform perturbations remain probability laws. This proves the necessary row-span condition without a boundary-support loophole. The arbitrary-expectation minimum is the quotient rank gap. Restriction to individual Fourier coordinates instead requires the union of missing character supports; these are different optimization problems. Independent rank controls verify the four-root pair recorder has rank 7 and the star adds exactly one dimension.

3. **Maximally indistinguishable laws: correct.** Every nontrivial character of the even-character group is balanced. The parity laws have identical surviving character moments and disjoint supports. Along their affine family, the Fisher speed and the arcsine distance formula agree with the chosen radius. Endpoint consumer separation is `2 |fhat(A)|`, giving the stated minimax lower bound and the four-root star's cost values 1 and 3/2.

4. **Complete moment matrix and triangle: correct.** The even subsets form a group of order `N=2^(m-1)` under symmetric difference. Its character table diagonalizes `M` with eigenvalues `N lambda_x`. Thus PSD is exactly atom positivity after normalization, and determines the orbit law, not an arbitrary nonsymmetrized law. For three roots the full Klein-four lift uses the same three scalar pair moments. The −1/2 equicorrelation passes ordinary 3×3 PSD but has a −1/2 full-matrix eigenvalue; the −1/3 q17 triangle passes the complete local test.

5. **One-missing-row hostile: correct, including the scalar distinction.** With `c=-1/(N-2)`, every proper principal test has nonnegative equicorrelation eigenvalues, while the full matrix has eigenvalue `-1/(N-2)`. The masses are exactly `-1/[N(N-2)]` and `(N-1)/[N(N-2)]`. At four roots the seven-row product table already contains the quartic scalar; the obstruction is the missing row, not an absent scalar moment. The stated pair projection really is realizable with quartic +1/2. For the one-negative-atom family the all-ones eigenvalue of an `(N-1)`-row test is `[1-N(N-2) epsilon]/(N-1)`, proving the sharp stated threshold (with epsilon nonnegative as intended).

6. **Added equicorrelation parity boundary: correct.** PSD gives `c>=-1/(m-1)`. Since `S=sum X_i` has the parity of m, `E S²=m+m(m-1)c` is at least 0 for even m and 1 for odd m. Uniform balanced words attain the resulting bounds; symmetric all-equal mixing fills the interval to 1. The odd-m gap is exactly `1/[m(m-1)]`. Independent sign enumeration checks m=2 through 7.

7. **q17 quantitative separation: correct as inherited-bound consequences.** The stronger addendum gives `Cut<=1730`; the base alone gives only 1732. The packet cites both correctly. Hence realizable edge correlations sum to at least −859 versus the witness's −867. Hölder/Cauchy–Schwarz yields lower bounds 8, 8/51, and 8/2601 for the edge-coordinate ℓ¹, ℓ², and ℓ∞ distances. Neither equality nor exact maximum cut is claimed. A fresh graph reconstruction confirms 289 vertices, degree 18, 2,601 edges, 1,734 triangles, and two triangles per edge. The inherited compressed proof certificate's hash matches, but its tree was not replayed in this bounded review. The q17 conclusions remain dependent on that separately reviewed finite certificate and the inherited PSD witness.

8. **Walsh stochastic composition defect: correct.** Squaring the matrix-product entries produces precisely the displayed cross terms. Orthogonality makes each squared-entry matrix doubly stochastic, so both defect margins vanish. The normalized symmetric Walsh involution gives coherent identity versus classical uniform mixing; point-mass total-variation discrepancy is `1-1/N`. Retaining signed amplitudes or the full covariance matrix restores the intended coherent composition. No claim that stochastic composition itself is invalid is made.

9. **Geometry and transfer types: sound.** Coordinate complement is an involutive order anti-isomorphism of `[0,1]^d`; it does not identify its face lattice or Euclidean polar with itself. Positive norm sections, probability amplitude charts, polar bodies, and ellipsoid relaxations are kept separate. The stated arithmetic and planar-Jacobian caveats correctly forbid transfer without a map preserving the actual integer/complex-algebraic predicates. No Collatz, Jacobian, or LRC conclusion follows from these results.

## Independent controls and source checks

`independent_controls.py` imports neither the source verifier nor repository code. It passed **77 exact controls** in ordinary Python and with `-O`; the resulting JSON files agree byte-for-byte. Its checks include actual sign-tuple enumeration, rational Gaussian elimination, character Gram reconstruction, every proper equicorrelation size in m=3 through 6, the missing-row sharp threshold, parity-boundary attainment, a fresh q17 adjacency/triangle reconstruction, and rational orthogonal composition examples.

Reproduce from the staging root:

    python -B session11-geometric-review/independent_controls.py
    python -B -O session11-geometric-review/independent_controls.py

The script intentionally reads the sibling staging packet and the inherited q17 certificate under `math/`; adjust those roots if publishing into a different layout. It prints a receipt and writes no files itself. All review artifacts are outside the repository; no Git mutation or repository write was performed.

The [Manifolds.jl probability-simplex documentation](https://juliamanifolds.github.io/Manifolds.jl/stable/manifolds/probabilitysimplex/) was checked for the stated radius-2 convention, metric, and distance. [Laurent, Lemma 8.15 and the ±1 note following it](https://homepages.cwi.nl/~monique/files/moment-ima-update-new.pdf) support the full finite Boolean moment criterion. [Ball, Lecture 3, Theorem 3.1](https://math.uchicago.edu/~shmuel/AAT-readings/Combinatorial%20Geometry,%20Concentration,%20Real%20Algebraic%20Geometry/ball.pdf) supports the John contact-point certificate. These classical citations are contextual; the scoped completion and defect proofs were checked directly.

## Optional editorial precision, not release blockers

- Section 5's heading calls the normalized simplex a matrix cone; strictly it is a normalized affine slice of the full character-positive cone. The theorem itself correctly imposes `mu_empty=1`.
- In §2.2, explicitly repeating `a_i>0` and symmetric `Q>0` would make the inherited polarity formulas standalone. The inherited October 9 source states the necessary convex-body and positive-definiteness assumptions.

Two parent-approved wording corrections were applied to the staging manuscript before final hashing: §2.2 now explicitly says polarity exchanges vertices with facets, and §2.5 restricts the infinite extreme-point count to dimension at least two while describing the one-dimensional interval exception. No verifier or mathematical-result changes were made. All final hashes and both independent receipts were regenerated afterward. The two optional points above were left unchanged.
