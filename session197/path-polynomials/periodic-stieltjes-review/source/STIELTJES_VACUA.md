# Stieltjes fractions and strict Hankel total positivity for two physical vacua

Session197. **PROVED strengthening restricted to two physical wrapper families.** This is not a universal positivity statement for arbitrary starting tournaments. The negative-Hankel witnesses in `PROOF.md` remain valid.

## 1. Physical moments and the inherited Jacobi matrices

Let W(T)=C3[T,K1,K1] and put

    mu_n^odd = H(W^n(K1)),
    mu_n^even = H(W^n(Tr2)).

The exact central-factorial calculation in `PROOF.md` gives, for either parity, mu_n=(J^n)_(0,0), where the symmetric tridiagonal matrix has diagonal a_n and neighboring entries sqrt(b_n), with

    odd:  a_n=n^2+6n+3,  b_n=2n^2(2n+1);
    even: a_n=n^2+6n+5,  b_n=2n(n+1)(2n+1), n>=1.

The earlier packet proved the formal Jacobi fraction and positive unshifted leading Hankel determinants. Positive Jacobi entries alone do not imply support on the positive half-line. We prove the additional matrix positivity here rather than assuming it.

## 2. Every finite Jacobi corner is positive definite

Define the rational LDL pivots

    ell_0=a_0,
    ell_n=a_n-b_n/ell_(n-1), n>=1.                          (1)

Positive pivots imply positive definiteness of every finite principal corner by exact LDL factorization.

For the odd matrix, induction gives

    ell_n >= (n+1)^2 > 0.                                  (2)

The base is ell_0=3>=1. If ell_(n-1)>=n^2, then

    ell_n >= n^2+6n+3-2(2n+1)=(n+1)^2.

For the even matrix the first nine pivots are exactly

    5, 48/5, 59/4, 1216/59, 4185/152, 10052/279,
    16762/359, 503016/8381, 1597341/20959.

All are positive, and the final value exceeds76. Let f_n=n^2+2n-4. For n>=9, f_(n-1)=n^2-5>0 and

    (a_n-f_n)f_(n-1)-b_n
      =(4n+9)(n^2-5)-2n(n+1)(2n+1)
      =(3n+5)(n-9) >= 0.

Consequently ell_(n-1)>=f_(n-1) implies ell_n>=f_n. Starting at n=8 proves positive pivots for all remaining indices. This finite boundary check plus the displayed polynomial identity is an all-order proof.

## 3. Positive Stieltjes fractions, with the normalization stated

Put

    alpha_1=ell_0,
    alpha_(2n)=b_n/ell_(n-1),
    alpha_(2n+1)=ell_n, n>=1.

Every alpha is a positive rational number. The ordinary generating series of either physical sequence has the formal Stieltjes fraction

    sum_(n>=0) mu_n z^n
      = 1/(1-alpha_1 z/(1-alpha_2 z/(1-alpha_3 z/(...)))).   (3)

Indeed contracting pairs of Dyck steps gives the Jacobi parameters

    a_0=alpha_1,
    a_n=alpha_(2n)+alpha_(2n+1),
    b_n=alpha_(2n-1)alpha_(2n),

which are exactly (1). Equivalently these identities verify the formal two-level continued-fraction contraction. Each coefficient depends on finitely many levels, so no analytic convergence is asserted or needed.

For example the odd alpha sequence starts3,2,8,5,14,9,21,96/7,205/7. The coefficients are counting weights after a rational change of representation; they are not independent physical graph operations.

## 4. Every Hankel minor is strictly positive

For either sequence, every finite minor

    det[mu_(r_i+c_j)]_(i,j=1..k)

is strictly positive when0<=r_1<...<r_k and0<=c_1<...<c_k. In particular every shifted Hankel matrix[mu_(i+j+s)]_(i,j=0..N-1), s>=0, is positive definite. This is strict Hankel total positivity for these two sequences only.

**Finite-dimensional proof.** Choose a finite Jacobi corner J_L large enough to reproduce all moments appearing in the minor and with L>=k. A walk reaching the omitted level must use at least L steps, so choosing L greater than every required exponent certainly suffices. Section2 makes J_L positive definite. Its eigenvalues theta_1<...<theta_L are positive and distinct, and every normalized eigenvector has a nonzero first coordinate: if that coordinate vanished, the tridiagonal recurrence with nonzero neighboring entries would force the entire vector to vanish. The same recurrence implies simple eigenvalues. Thus the reproduced moments have a finite representation

    mu_t=sum_(l=1..L) w_l theta_l^t, w_l>0.

Cauchy-Binet expresses the desired minor as a sum over k-element node subsets of

    (product w_l) det[theta_l^(r_i)] det[theta_l^(c_j)].

Both generalized Vandermonde determinants are positive on positive increasing nodes. Here is an elementary sign justification. A nonzero linear combination of k distinct nonnegative integer powers has at most k-1 positive roots: divide out the lowest power, differentiate, and induct using Rolle's theorem. Hence the determinant never vanishes on the connected region of ordered positive nodes. As the nodes coalesce, its first nonzero term has the sign of the Wronskian of those powers, which equals a positive power of the common node times the positive Vandermonde product in the increasing exponents. Its sign is therefore positive everywhere. Every Cauchy-Binet summand is positive, proving the claim.

For the first shifted leading determinants there is also an exact product:

    det[mu_(i+j+1)]_(i,j=0..N-1)
      = (product_(k=1..N-1) b_k^(N-k))
        (product_(k=0..N-1) ell_k).                         (4)

The vectors e0,J e0,...,J^(N-1)e0 form a triangular matrix whose squared determinant is the first product. Their Gram matrix for the J inner product is the shifted Hankel matrix; det(J_N) is the second product. The N-corner reproduces moments through2N-1, so this calculation is exact.

## 5. Essential scope and checks

Starting from Tr4 instead gives the actual wrapper counts1,17,241, whose2-by-2 Hankel determinant is-48. Starting from Tr5 gives1,33,645 and determinant-444. Their positive central-factorial initial states do not turn them into the vacuum state e0 for the same symmetric realization. Nothing above extends total positivity to those inputs, arbitrary rooted contexts, or arbitrary positive algebraic states.

`verify_stieltjes_vacua.py` checks the exact boundary pivots, both induction identities, positive pivots through100, Stieltjes Dyck moments against the physical operator, shifted determinant products and nonprincipal minors, plus literal small wrapper graphs and both negative starting-state controls. Run normally and with `-O`; explicit exception gates remain active. The all-order claims rest on Sections2-4, not on a finite determinant census.
