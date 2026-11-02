# Arithmetic second proof of the tensor compensation lemma

Status: independent complete proof; exact controls in
`verify_eisenstein_bernstein.py`. This note is separate from the shorter
inductive proof and from all frozen cutoff packets.

## Theorem

Let P(x_1,...,x_r) be a multiaffine polynomial with integer coefficients.
Assume its coefficient of x_1...x_r is odd and all other coefficients are
even. Then, throughout |x_i|<=1,

    sum_(S subset [r]) (partial_S P(x))^2 >= product_i(1+x_i^2).

The empty derivative means P itself. For r=0 this says an odd integer's
square is at least1. Equivalently, if b_empty is odd, all other b_U are
even, and

    f_S=sum_(J subset S) b_(S\J) product_(j in J) x_j,

then sum_S f_S^2>=product_i(1+x_i^2), since f_S is the mixed derivative
complementary to S of P=f_[r].

## Proof

This independently found proof is longer, but confirms the parity mechanism
without using the preceding induction.

Let P be an integer multiaffine polynomial whose top coefficient is odd and
all other coefficients are even. Its mixed-derivative jet norm is

    J(P;x)=sum_(S subset [r]) (partial_S P(x))^2.

Changing signs of coordinates preserves the coefficient parity and the
jet norm, so it suffices to use0<=x_i<=1. For the coefficient vector c in
the monomial basis,

    J(P;x)=c^T [tensor_i Q(x_i)] c,
    Q(x)=[[1,x],[x,1+x^2]].

The degree-two Bernstein coefficient matrices on[0,1] are

    Q0=I,
    Q1=[[1,1/2],[1/2,1]],
    Q2=[[1,1],[1,2]]=M^T M,  M=[[1,1],[0,1]].

The corresponding scalar coefficients of1+x^2 are1,1,2. Fix a tensor
Bernstein coefficient with k coordinates of type2. Applying the integral
matrix M in each such coordinate sends the parity vector(0,1) to(1,1).
The identity matrices in type0 coordinates separate coordinate slices.
Retain the slice with top index1 in all type0 coordinates. There are2^k
retained slices indexed by the type2 coordinates; each retained slice has
an integer coefficient vector with only its top coefficient odd. All other
slices contribute nonnegative quadratic forms and may be discarded.

It remains to show c^T(Q1 tensor ... tensor Q1)c>=1 for an odd-top/even-rest
integer coefficient vector. Put zeta=exp(i*pi/3), so

    zeta^2-zeta+1=0,        conjugate(zeta)=1-zeta.

The quadratic form is the average of

    |P(zeta_1,...,zeta_m)|^2,

where each zeta_i is independently zeta or its conjugate. Each evaluation
is an Eisenstein integer. Modulo2, the polynomial z^2-z+1 becomes
z^2+z+1, which is irreducible over F2; thus the quotient is F4. Each zeta_i
has nonzero residue, and the odd-top/even-rest coefficient condition makes
the residue of P equal the nonzero product of these residues. Therefore
P(zeta_1,...,zeta_m) is nonzero. Its squared absolute value is the positive
integer a^2+ab+b^2 for integers a,b, hence is at least1.

Each of the2^k retained slices contributes at least1. Thus every tensor
Bernstein coefficient of J(P;x)-product_i(1+x_i^2) is nonnegative. This
proves(1), since f_S is the mixed derivative complementary to S of P=f_[r].

This second proof uses a specific mod2-to-F4 nonvanishing argument. It does
NOT establish an arbitrary-modulus replacement, nor any Riemann-hypothesis
conclusion.


## Why this is a genuine compensation statement

The assertion is not the false pointwise inequality P(x)^2>=product x_i^2.
For example P(x)=3x-2 vanishes at x=2/3. Its complete derivative jet instead
contains the coefficient square9, which pays for that cancellation. The
arithmetic proof does this simultaneously for an entire hypercube, with no
assumption that a binary complementary Pfaffian has magnitude1.

For matching uncertainty in a tournament, each suitably grouped hypercube
of principal Pfaffians has exactly the polynomial parity pattern in this
theorem. The coefficient parity follows from counting perfect matchings
modulo2 after deleting disjoint prescribed edges. Summing the resulting
bounds by vertex parity gives the absolute matching margin and the nonlinear
retained-baseline gap described in `INDEPENDENT_TENSOR_MATCHING_REVIEW.md`.
The present note independently proves only the tensor lemma; it does not
replace the determinant identities or the separate Pfaffian grouping audit.

## Exact controls and limitations

`verify_eisenstein_bernstein.py` checks the three Bernstein matrices, the
integral Q2 Gram factor, Eisenstein multiplication/conjugation/norm identities,
all nonzero residues of F4, the phase-average Gram identity in small dimensions,
and mixed tensor Bernstein coefficients for seeded integer odd-coset vectors.
It also rejects P=x-1 as a control when the even-proper-coefficient hypothesis
is removed. These finite tests audit algebra and implementation; the proof
above, not the finite sample, supplies arbitrary dimension.

The arithmetic input is specifically the nonvanishing of the relevant coset
in F4=F2[z]/(z^2+z+1). It provides no automatic analogue for every modulus,
and has no implication about zeros of the Riemann zeta function.
