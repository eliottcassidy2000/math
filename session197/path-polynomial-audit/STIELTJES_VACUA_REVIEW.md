# Independent review: Stieltjes fractions for the two physical vacua

**PASS. No correction required.** Review date2026-11-05. Frozen source files in `session197-odd-cycle-path-polynomials/`:

- `STIELTJES_VACUA.md`: `dbe7c0e68e77e5e14e389c57221578191a22381b086cb4c5864c7b2d1abc6fe9`
- `verify_stieltjes_vacua.py`: `a1700d8806eccb0307fe65431690307c17a057064e4e392da7f9afa1580cc098`
- Normal/optimized receipt: `030f493471bf02a865e172effb824e81e506c543c933ad791352ac8a952ec216`

## All-order proof review

The Jacobi parameters agree with the previously audited physical wrapper operator. The odd pivot estimate is a direct valid induction. For the even sequence, the nine displayed boundary pivots are correct, the final one exceeds76=f_8, and the exact residual identity is `(3n+5)(n-9)`. Its nonnegativity for n>=9 and positive denominator n^2-5 establish all remaining pivots. Thus every finite symmetric Jacobi corner is genuinely positive definite; the argument does not confuse positive off-diagonal entries with positive spectrum.

The stated Stieltjes normalization is correct: alpha_(2n-1)*alpha_(2n)=b_n and alpha_(2n)+alpha_(2n+1)=a_n, with alpha_1=a_0. All alpha values are positive rational numbers. Formal continued-fraction contraction is legitimate coefficient by coefficient, without an analytic convergence or infinite-operator assertion.

For any requested finite Hankel minor, a sufficiently large positive-definite Jacobi corner reproduces its moments. The irreducible symmetric tridiagonal matrix has simple positive eigenvalues and strictly positive spectral weights at e0. Cauchy-Binet then expresses every minor as a positive sum of products of two generalized Vandermonde determinants. The proof of their positive sign is valid: induction using Rolle bounds positive zeros of a nonzero sparse polynomial, and the coalescing-node Wronskian has the positive exponent-Vandermonde sign. This establishes strict total positivity of every minor with increasing row and column index lists, not merely the leading minors.

The shifted-leading determinant product also passes. The N-corner reproduces all required moments through2N-1, since an excursion reaching the omitted levelN and returning needs at least2N steps. The Krylov matrix is triangular with squared determinant product b_k^(N-k), and its Gram matrix for the J inner product supplies det(J_N)=product ell_k.

All assertions are correctly restricted to the two vacuum starting graphs K1 and Tr2. The Tr4/Tr5 negative-Hankel controls remain valid and exclude an arbitrary-input extension. No independent physical operation is attributed to an individual rational Stieltjes weight.

## Independent checks

The independent script in `stieltjes-audit/` imports no producer code. It verifies201 exact LDL pivots in each parity; computes25 moments independently from the tridiagonal walk rule and by truncated rational-series inversion of the Stieltjes fraction; checks nine shifted determinant products; and checks2,940 nonprincipal minors per parity. It independently counts the first five actual wrapper graphs in each parity by walk inclusion–exclusion, plus both negative starting-state examples. The induction polynomial identity is checked coefficientwise.

Copied-source normal and optimized replays reproduce the frozen receipt exactly. Original files remain unchanged. The finite tests corroborate the all-order finite-matrix argument rather than replace it.
