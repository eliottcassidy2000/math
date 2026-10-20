# Independent audit of the commutator-clock mixing proof

2026-10-20. **PASS for the stated finite-state estimates.** Reviewed `COMMUTATOR_CLOCK_MIXING.md` and `verify_clock_mixing.py`. The numerical eigenvalue checks are controls only; the proof below establishes the general bounds.

## Algebra and constants

The convention `U_g f(x)=f(g^(-1)x)` is homomorphic. Thus the affine commutator `EOE^(-1)O^(-1):x↦x-1/4` gives precisely its permutation operator. Since `D` is coprime to six, the translation step is a unit and its mean-zero squared displacement has lower bound `4 sin²(pi/D)`.

For `W=UVU^(-1)V^(-1)`, telescoping the four factors and using unitary prefixes gives

`||(W-I)f|| <= 2(||(U-I)f||+||(V-I)f||)`.

Squaring and applying Cauchy–Schwarz yields the stated factor eight. Consequently the sum of the two generator residual squares is at least `(1/2)sin²(pi/D)||f||²` on mean-zero functions.

The symmetric operator `S=(U+U*+V+V*)/4` has Dirichlet form one quarter of that residual sum, hence gap at least `sin²(pi/D)/8`. Laziness divides the gap by two and puts the full spectrum in `[0,1]`. This proves the displayed lazy gap `sin²(pi/D)/16` and the total-variation estimate with prefactor `sqrt(D-1)/2`. The prefactor is correct for the uniform-normalized L2 norm of a point-mass density minus one.

For a forward eigenvector of `P=(U+V)/2`, center both generators at `zeta=lambda/|lambda|`. Their commutator is unchanged. The residual sum is exactly `4(1-|lambda|)||f||²`, so the same comparison yields `|lambda|<=1-sin²(pi/D)/8`. The `lambda=0` case is immediate.

The fixed space consists of constants: equality in the unitary average forces invariance under both generators and their transitive translation. For **algebraic** simplicity, add the elementary fact that `||P^t||<=1` excludes a nontrivial Jordan block at eigenvalue one. This is a minor proof-explicitness improvement, not a counterexample.

## Scope checks

- The explicit finite-time mixing estimate applies only to the self-adjoint lazy symmetric chain. A forward eigenvalue bound alone does not establish the same operator-norm bound for a nonnormal operator.
- The symmetric chain includes inverses and holding. The forward random-choice chain averages over compositions. Neither is the fixed-length, fixed-number-of-odd-steps law required for the carry problem.
- The coefficient of `z^p` in `(U_E+zU_O)^q` is the correct retained fixed-composition object; noncommutative order must remain.
- The comparison between the sufficient time `D² log D` and the word length `q` must stay conditional on the stated large-clock regime. Entropy alone does not control `|2^q-3^p|`; no new uniform denominator lower bound is supplied here.
- A commutator with several orbits leaves additional invariant functions; the one-cycle hypothesis or an explicitly restricted invariant subspace is load-bearing.

## Independent exact controls

`check_clock_energy.py` checks the homomorphic permutation convention and the phase-centered energy inequality using Gaussian integers, not floating point. It covers every modulus 5–80 coprime to six, twenty deterministic mean-zero vectors per modulus, and phases `1,i,-1,-i`: **3,120 explicit checks, PASS**. Normal and `python3 -O` outputs are byte-identical.

The supplied clock verifier correctly labels its eigenvalue comparisons as numerical, with stated tolerances, and its permutation checks as exact. Its explicit exception gates are not disabled by optimization. No numerical bound is substituted for the proof.

## Reviewed follow-on: 3-smooth Fourier rectangle and failure of a uniform lazy gap

The parent's proposed upper-bound construction is also **mathematically valid**, with the following explicit carrier checks.

Take `L>=1`, `6^L<=sqrt(D)`, and normalized additive Fourier modes at integer frequencies `a_ij=2^i3^j`, `1<=i,j<=L`, weighted by the product of the normalized first Dirichlet sine vectors. For `D>4`, all source frequencies and both required shifted supports are positive, distinct, and below `2·6^L<=2sqrt(D)<D`. Unique factorization therefore rules out hidden modular collisions.

In the declared representation, `U_E` moves `(i,j)` to `(i+1,j)`. `U_O` moves it to `(i+1,j-1)` and multiplies by `exp(-2πia_ij/(3D))`. The latter expression is justified by **`j>=1`**, which makes the integer frequency divisible by three; it is not a casual replacement of a modular inverse by real division.

With `theta=pi/(L+1)`, the one-coordinate sine overlap is `cos(theta)`, while the two-coordinate diagonal shift overlap is `cos²(theta)`. Hence the exact squared residuals are `2(1-cos theta)` for `E` and `2sin²(theta)` for the unphased `O` shift. The actual phase adds norm error at most `delta=2π6^L/(3D)`. The lazy Dirichlet form is thus bounded by

`[2(1-cos theta)+(sqrt(2)sin theta+delta)^2]/8`

`<= 5pi²/[8(L+1)^2] + pi²6^(2L)/(9D²)`.

The coarse constant follows directly from `2(1-cos theta)<=theta²`, `sin theta<=theta`, and `(a+b)^2<=2a²+2b²`. Choosing `L=floor(log D/(2log 6))` for sufficiently large `D` coprime to six makes this tend to zero as `O((log D)^-2)`. Thus the lazy symmetric spectral gap genuinely has no uniform positive lower bound over all these moduli. This is compatible with the previously proved lower bound of order `D^-2`.

The argument uses a Rayleigh quotient for a self-adjoint operator. It does **not** promote approximate vectors to forward eigenvalues for the nonnormal directed operator, and it does not address fixed-composition carry cancellation.
