# Integral Fourier coordinates: admissibility, information loss and bad-prime collapse

**Status:** elementary proofs below; independent audit requested. Classical Walsh and Smith ingredients are not claimed new. This is a common integral observer for twist profiles, boundary tables and signed-orbit lattices, not an identification of their target predicates.

## 1. The exact integral lattice

Let G=F2^r, N=2^r and H_(a,x)=(-1)^(a dot x). For an integer histogram n put A=Hn. Orthogonality gives H²=NI, so n=HA/N over Q.

**Theorem F1.** The Smith factors of H are 2^j with multiplicity binom(r,j), j=0,...,r. Thus Z^N/HZ^N is the direct sum of (Z/2^j)^(binom(r,j)), and |det H|=2^(r 2^(r-1)) for r>=1.

**Proof.** H1=[[1,1],[1,-1]] is integrally equivalent to diag(1,2) by elementary row and column operations. H_r is its r-fold tensor power. Tensor powers of the unimodular equivalence matrices remain unimodular. The resulting diagonal entries are products of r choices from1 and2; sorting by divisibility proves the assertion. QED.

For64 coordinates the index is2^192, not1. The transform is invertible over every odd-characteristic field but has rank1 modulo2. Over Z/2^k the full Smith exponents, not only the field rank, control its defect. The tournament Fourier matrix has the same shape, but its action and interpretation differ.

## 2. Positivity and integrality are separate exact conditions

**Theorem F2.** An integer vector A is the Fourier transform of a nonnegative integer histogram iff every coordinate of HA is a nonnegative multiple of N. Define the real symmetric matrix C_A by (C_A)_(a,b)=A_(a+b). Its Walsh eigenvectors have eigenvalues (HA)_x=N n_x. The criterion is equivalently

    A belongs to H Z^N, and C_A is positive semidefinite.

This follows immediately from H²=NI and direct diagonalization. It characterizes abstract histograms, not which histograms a chosen algebraic arrangement, unit-distance exterior or tournament can realize.

For r=2, n=(1,1,1,-1) has transform A=(2,2,2,-2). It satisfies A0=2, |Aa|<=A0 and all Fourier integrality congruences, but C_A has a negative eigenvalue. Scalar bounds and correct lattice membership do not imply effectivity.

## 3. Arbitrary pinned continuation requires the full table

For nonnegative tables f,g on G, shared-state gluing is pointwise multiplication; total mass is sum_x f_x g_x. If all indicator contexts g=1_(x=x0) are allowed, identical continuation responses are equivalent to equality of every table entry. No proper linear quotient preserves all those named queries. This concerns the declared context bank, not geometric realizability of every context.

For A=Hf, B=Hg, direct orthogonality gives

    H(fg)_c = (1/N) sum_a A_a B_(c+a).

The denominator is exact on valid integral Fourier tables, because the left side transforms an integer table. It cannot be discarded after reducing arbitrary Fourier coordinates modulo2. A mean, selected low modes or an energy norm need not be closed under this operation.

For a Kummer torsor, n_x counts base points with square-class profile x. Twist point counts are N n_x and character sums are Hn. Their mean A0 is the base count. The full profile and its marked origin recover a specified twist.

## 4. Full two-point correlations still forget phase

Translation n'_x=n_(x+h) changes Aa to (-1)^(a dot h)Aa. It preserves A0, every Aa² and the full translation autocorrelation, but may change named counts. Covariance is stronger than a mean, weaker than a marked response table.

Exact algebraic controls: on U_lambda=P1 minus{0,1,lambda,infinity}, the Legendre double cover y²=x(x-1)(x-lambda) has two twist counts (q-3)+S_lambda and(q-3)-S_lambda, where S_lambda is the quadratic character sum of the displayed cubic.

- OverF7, lambda2 gives S=0 and counts(4,4); lambda3 gives S=-4 and counts(0,8). The same four-point branch matroid and base count do not determine variance:0 versus16.
- OverF5, lambda2 gives S=2 and counts(4,0); lambda3 gives S=-2 and counts(0,4). Mean and variance agree but marked twists differ. This is an origin/rational-twist effect, not a claim that their unmarked twist families differ.

Thus additive branch coefficients, signs and chosen twist labels cannot be replaced by incidence, cardinality or unsigned energies alone.

## 5. Proposed interpretation

A faithful integral observer retains three pieces: its integral lattice and action; the operation/context bank; and the positive or geometric model cone. The same matrix may describe twist histograms, coloring boundary counts or tournament Fourier functions, but the theorems above transport only the declared linear/composition data. They do not transport colorability, rational points, torification or tournament moves without another proved map.

Run check_fourier_effectivity.py normally and under Python optimization. It checks Smith factors throughr6, all81 ternary nonnegative tables onF2² and all6561 pair products, phase translations and the four explicit Legendre controls. Finite tests supplement the all-r proofs.
