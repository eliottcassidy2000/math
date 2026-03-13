# THM-163: H Walsh-Hadamard Fourier Structure

**Status:** VERIFIED (n=3 through n=6, exhaustive)
**Session:** kind-pasteur-2026-03-13-S61

## Statement

Let T be a tournament on n vertices with skew-symmetric encoding sigma in {-1,+1}^m (m = C(n,2)). Define the Walsh-Hadamard expansion:

H(sigma) = sum_{S subset [m]} hat(H)(S) * prod_{i in S} sigma_i

Then:

### 1. Even-Degree Vanishing

hat(H)(S) = 0 for all |S| odd.

**Proof:** H(T) = H(T^op) (path reversal bijection). T^op corresponds to sigma -> -sigma, which multiplies hat(H)(S) by (-1)^{|S|}. Invariance forces odd-degree terms to zero.

### 2. Degree-2 Coefficient Formula

For edges e_1 = (a,b) and e_2 = (c,d) sharing exactly one vertex v:

hat(H)({e_1, e_2}) = (-1)^{pos_1 + pos_2 + 1} * (n-2)!/2^{n-2}

where pos_i = 0 if v is the smaller endpoint of edge i, 1 if larger.

For edges sharing NO vertex: hat(H)({e_1, e_2}) = 0.

**Verified exhaustively at n = 3, 4, 5, 6.**

### 3. Vertex Balance Function Z_v

Define Z_v = sum of signed products of sigma over edge pairs sharing vertex v:

Z_v = sum_{a<v<b} sigma_{av} * sigma_{vb} - sum_{a<b<v} sigma_{av} * sigma_{bv} - sum_{v<a<b} sigma_{va} * sigma_{vb}

Then Z_v depends ONLY on the score s_v = out-degree of v:

**Z_v = -2(s_v - (n-1)/2)^2 + (n-1)/2**

This is a downward parabola centered at the median score (n-1)/2, with maximum value (n-1)/2 achieved at the regular score.

**Verified exhaustively at n = 3, 4, 5, 6, 7.**

### 4. H_2 as Score Regularity

The degree-2 part of H is:

H_2 = (n-2)!/2^{n-2} * sum_v Z_v
    = (n-2)!/2^{n-2} * n * ((n-1)/2 - 2*Var(score))

where Var(score) = (1/n) * sum_v (s_v - (n-1)/2)^2.

**Key consequences:**
- H_2 is maximized when Var(score) = 0 (regular tournaments)
- H_2 is minimized at the transitive tournament (max variance)
- corr(score_variance, H) = -0.973 at n=5

### 5. Maximum Fourier Degree

The maximum non-zero even degree is exactly 2 * floor((n-1)/2):
- n=3,4: max degree 2 (purely quadratic in sigma)
- n=5,6: max degree 4
- n=7,8: max degree 6 (predicted)

**Verified exhaustively at n = 3, 4, 5, 6.**

### 6. Degree-4 Structure

All degree-4 coefficients have support on 5-vertex subsets with graph degree type (2,2,2,1,1). At n=5: all 60 non-zero coefficients have magnitude exactly 1/8.

At n=6: two types appear:
- 360 coefficients with |coeff| = 1/8 (5-vertex support)
- 90 coefficients with |coeff| = 1/4 (6-vertex support, type (2,2,1,1,1,1))

## Energy Distribution

| n | Degree 0 | Degree 2 | Degree 4 | Degree 6+ |
|---|----------|----------|----------|-----------|
| 3 | 75.00%   | 25.00%   | 0%       | 0%        |
| 4 | 75.00%   | 25.00%   | 0%       | 0%        |
| 5 | 75.95%   | 22.78%   | 1.27%    | 0%        |
| 6 | 77.59%   | 20.69%   | 1.72%    | 0%        |

The Fourier energy is DOMINATED by the constant and degree-2 terms. At n=3,4, H is exactly a degree-2 polynomial. The degree-4 correction is < 2% of the total energy.

## Connection to Regular Tournament Maximization

For regular tournaments:
- Z_v = (n-1)/2 for all v (maximum value)
- H_2 = (n-2)!/2^{n-2} * n * (n-1)/2 = n!/2^{n-1} = E[H]
- Therefore H_0 + H_2 = 2 * E[H]
- At n=5: degree-4 terms vanish exactly, giving H = 2 * E[H] = 15

## Connection to Vitali Set Analogy

The Fourier structure reveals the hidden higher-dimensional structure:
- Tournament space is {-1,+1}^m, a Boolean hypercube of dimension m = n(n-1)/2
- H is an EVEN function (invariant under global sign flip = tournament reversal)
- The effective dimension is much lower: H depends primarily on score variance (1 DOF!)
- The degree-4 terms encode cycle-level interactions beyond score structure

## Verification

- vitali_tournament_structure.py, H_fourier_structure.py, H_fourier_formula.py, H_fourier_deep.py
- All results in 05-knowledge/results/
