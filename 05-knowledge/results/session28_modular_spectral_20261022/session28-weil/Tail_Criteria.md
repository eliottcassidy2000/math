# Quantitative finite positivity versus global Weil positivity

These are standard operator/form sufficient criteria, not new RH claims. The missing work is establishing their hypotheses for the full desired Weil form. All inner products use one fixed Hilbert-space norm; changing coordinates requires tracking its Gram matrix.

## 1. Three different approximation operations

An archimedean integration cutoff T at fixed frequency dimension N changes entries of one matrix. Increasing N admits new test directions. Increasing the prime/spatial cutoff c changes the allowed test-function support. A sign certificate for the first limit alone does not settle either of the other two.

For fixed (c,N), assume exact matrices satisfy 0 <= Q_(N,infinity)-Q_(N,T) <= b_(N,T) I. Then

lambda_min(Q_(N,T)) <= lambda_min(Q_(N,infinity)) <= lambda_min(Q_(N,T))+b_(N,T).

This follows by evaluating Rayleigh quotients and taking infima. Certified finite-cutoff positivity survives the integration tail. A negative value below -b is genuinely negative. Anything between -b and zero is inconclusive. For an interval/numerical matrix with norm error epsilon, replace the lower endpoint by lambda_min(approximation)-epsilon and the upper endpoint by lambda_min(approximation)+epsilon+b.

By contrast, if Q_N are exact compressions of ONE form to nested subspaces in ONE norm, their minimum Rayleigh quotients are nonincreasing in N. A new direction can introduce a negative value even when every earlier tested compression was positive. Example: diag(1,...,1,-1) after any prescribed finite initial segment. Restrictions to different test dictionaries need not be nested, so nesting itself must be proved.

## 2. Joint cutoff scheduling, a sufficient quantitative rule

The recent finite Guinand-Weil source proves a fixed-(c,N) positive archimedean tail with the bound below (Theorem3.2 and Corollary3.3). This is a scoped cited result, not an independently rerun interval certificate:

b_(N,T) <= [2(2N+1)rho/pi²] [log(T)/(T-rho N) + log(T/(T-rho N))/(rho N)],

where rho=2pi/log(c), N>=1 and T>max(rho N,7). See https://arxiv.org/html/2607.02828v1 .

Our direct simplification: if T>=2rho N, then

b_(N,T) <= [4(2N+1)rho/pi²] (log(T)+1)/T.

Proof: 1/(T-rho N)<=2/T and -log(1-u)<=u/(1-u)<=2u for u=rho N/T<=1/2. Thus at fixed c the explicit sufficient schedule T_N/(N log T_N) -> infinity forces this certified error budget to zero. For instance T_N=N(log N)^2 works eventually. This schedule is sufficient for the displayed trace-based bound, NOT a necessary cutoff growth condition for actual operator-norm convergence. It does not furnish positive minimum eigenvalues or a coercivity constant. The latter can shrink faster than the error budget.

## 3. A full-space cross-term sufficient criterion

Let a self-adjoint block form on H=E direct-sum E-perp be

Q = [[A,B],[B*,C]],

with A>=a I, C>=gamma I, a,gamma>0, and bounded B with ||B||<=beta. Then

Q >= kappa I, where
kappa = (a+gamma-sqrt((a-gamma)^2+4beta²))/2.

Proof: Q[(x,y)] >= a||x||² -2beta||x||||y|| + gamma||y||², and the smallest eigenvalue of the real2x2 comparison matrix is kappa. In particular a gamma>beta² gives uniform positivity. Equality only gives nonnegativity. When C is positive invertible, the exact condition is A-B C^-1 B*>=0, by completing the square. Hence omitted directions require BOTH their own positive lower bound and payment for cross terms. A finite principal-block certificate by itself supplies neither.

For unbounded forms, require a common dense form domain compatible with the decomposition, a closed semibounded extension, and the indicated cross-form bound. Do not treat arbitrary formal matrices as a self-adjoint operator by declaration.

## 4. A countable, checkable sufficient target

Suppose q is a Hermitian form on finitely supported sequences with entries q_ij and positive diagonal d_i=q_ii. Let K_ii=0 and K_ij=q_ij/sqrt(d_i d_j) for i!=j. If positive weights w_i and eta<1 satisfy

sum_(j!=i) |K_ij| w_j <= eta w_i for EVERY i,

then the weighted Schur test gives ||K||_(ell²)<=eta and

q[x] >= (1-eta) sum_i d_i |x_i|².

Proof of the bound can avoid any operator assumption: apply
2|u_i||u_j| <= |u_i|² w_j/w_i + |u_j|² w_i/w_j,
sum against |K_ij|, and use Hermitian symmetry. This holds first for finite support. The diagonal form is closed; the relative perturbation with norm less than1 gives a closed extension and preserves the bound. If inf_i d_i>0 this is coercive in ell². Otherwise it is positive in the diagonal energy norm, not necessarily in the ambient norm. Establishing the infinite row inequalities for an actual Weil dictionary would be substantive new work, not something supplied by a finite experiment.

## 5. Signed channels must be justified, not discarded

A decomposition into even/odd, irreducible, score/cycle, or local prime/archimedean pieces is valid only with the actual action and observables. If a form commutes with an involution, even and odd blocks have no cross term, but positivity of one block says nothing about the other: diag(1,-1) is already a counterexample. If the form does not commute, dropping mixed terms can destroy positivity: [[1,2],[2,1]] has positive diagonal blocks but a negative eigenvalue. Representation multiplicity spaces can carry nontrivial Hermitian matrices even when irreducible types are known. Symmetry therefore reduces a problem, not automatically solves it.

The tournament Hessian provides a concrete instructive signed decomposition: scores have eigenvalue c(2-n), cycles2c. Replacing it by its absolute value produces a positive operator but loses the original second variation. Similarly, replacing a signed arithmetic explicit-formula channel by its absolute value cannot retain an exact zero-side identity without a new correction theorem.

## 6. What all-dimensional positivity WOULD prove

If a closed semibounded form q has a form core equal to the union of nested finite spaces E_N, and q>=0 on EVERY E_N, continuity in the form norm implies q>=0 on its full domain. A common q>=kappa||x||² bound gives coercivity. Neither conclusion follows from checking merely a finite number of N. For the Weil criterion one must additionally prove that the chosen core or dictionary covers the required test class and that the normalization is the intended explicit-formula form. There is no claim here that these hypotheses have been established.
