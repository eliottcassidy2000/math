# Independent structural review: tensor compensation and nonlinear matching gap

Status: independently reconstructed; PASS on the mathematical argument below.
This note does not modify the frozen n>=748 or n>=526 packets. Exact matrix
controls are being checked separately by the tensor-verification agent.

The conclusion is the all-order matching-uncertainty refinement

    u-|v| >= product_i(1-q_i)                     if n>2r,
    u-|v| >= product_i(1-q_i)-product_i q_i        if n=2r,

and consequently

    Delta(A) >= E Delta(T)+1-product_i(1-q_i)
                  -1_{n=2r} product_i q_i.

Here r is the number of disjoint uncertain pairs, all other edges are
binary, q_i=p_i(1-p_i), u=det(I+A), and v=(-1)^n det A. Empty products
are1. The zero-vertex remainder has margin0 and is handled separately.
This is an absolute margin bound, NOT either false conditional product
factor documented in `LOCAL_FACTOR_HOSTILES.md`.

## 1. The elementary tensor lemma checks

For integer coefficients b_U indexed by subsets of [r], suppose b_empty
is odd and every b_U with U nonempty is even. Define

    f_S(x)=sum_(J subset S) b_(S\J) product_(j in J) x_j.

Then, for |x_i|<=1,

    sum_(S subset [r]) f_S(x)^2 >= product_i(1+x_i^2).    (1)

The induction is valid in every dimension. On splitting off one coordinate,
the vector of all f_S has the form (a,b+x a). The coefficient array for a
has the required odd/even parity; every coefficient of b is even. Thus a,
a+b, and a-b each satisfy the preceding-dimensional hypothesis. With

    L=product_previous(1+x_i^2),
    A=||a||^2, C=||b||^2, D=<a,b>,

we have A>=L and A+C+/-2D>=L, so |2D|<=A+C-L. Therefore

    A+||b+x a||^2
      >=(1-|x|+x^2)A+(1-|x|)C+|x|L
      >=(1+x^2)L.

Every coefficient multiplying A or C is nonnegative for |x|<=1. No
coordinatewise unit-Pfaffian assertion is used. The r=0 base is the square
of an odd integer, at least1. In particular the previous hostile3x-2 is
compensated by its coefficient square: 3^2+(3x-2)^2>=1+x^2.

## 2. Independent second proof via Bernstein coefficients and F4

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

## 3. Pfaffian hypercube grouping is legitimate

Label each uncertain matching pair consecutively. Fix which other pairs
contribute exactly one endpoint, choose those endpoints, and fix a subset
of the binary isolated vertices. In the observer-containing sum, include
the observer as a fixed vertex. Call this binary base B0. Only choices
with even |B0| contribute Pfaffians; adding full pairs never changes parity.

For every remaining pair, allow either neither or both endpoints. For a
set S of such full pairs, let f_S be the corresponding principal Pfaffian.
Because the pair endpoints are adjacent in the induced order,

    partial_(x_i) Pf = Pf_delete_pair

with positive sign. Thus these f_S have exactly the convolution form of
Section1. The constant b_empty is the binary Pfaffian of B0, hence an odd
integer, including the empty Pfaffian1.

For nonempty U, the constant b_U is the Pfaffian on B0 plus the full pairs U,
with those internal entries set to zero. Modulo2 every other off-diagonal
entry is1. Inclusion-exclusion over the |U| disjoint forbidden edges gives

    b_U mod2 = sum_(J subset U) 1 =0.

Every residual complete even graph has an odd number of perfect matchings,
including the empty residual graph. Hence all required b_U are even, and
(1) applies to the entire hypercube of principal Pfaffians. This groups
coefficient energy with the terms it controls, rather than bounding an
individual Pfaffian square falsely.

## 4. Normalization and parity signs check

Write x_i=2p_i-1, so q_i=(1-x_i^2)/4, and t=n-2r. After applying(1) to each
hypercube group, its none/full choice has aggregate weight1+x_i^2; its
single-endpoint choice has weight2 and odd vertex parity. Thus the parity
polynomial for all groups is

    product_i[(1+x_i^2)+2z] (1+z)^t,

with z tracking parity. The accepted augmented-skew identities express
u+v and u-v as2^(1-n) times the respective even and odd parity sums.

If t>=1, the z=-1 evaluation vanishes. Both margins are at least

    2^(-2r) product_i(3+x_i^2)=product_i(1-q_i).

If t=0, the two lower bounds are

    product_i(1-q_i) +/- (-1)^r product_i q_i.

Their minimum is the stated absolute margin bound. At n=2,r=1 this is
1-2q, exactly correct. At n=3,r=1 it is1-q, also sharp.

## 5. The retained-baseline nonlinear gap follows without new sign assumptions

Round matching edges in any fixed order, retaining the terminal expectation
B=E Delta(T). The accepted one-edge inequality gives correction q_i times
the minimum margin of the pair-deleted remainder on every conditioning branch.
All cross edges of that remainder remain binary, so the new absolute margin
bound applies regardless of the orientations already rounded.

For n>2r, every remainder has a binary isolated vertex; the total correction is

    sum_i q_i product_(j>i)(1-q_j)=1-product_i(1-q_i).

For n=2r, only the first pair-deleted remainder lacks a binary isolate.
Its extra subtraction is q_1 product_(j>1)q_j=product_i q_i. Later steps
have the two vertices of a previously rounded pair as binary isolates.
This proves the claimed nonlinear bound with B retained once.

The proof needs neither a guessed determinant sign nor a fixed-sign surplus:
the absolute-margin theorem is strong enough that the already accepted
minimum-at-each-edge interpolation telescopes correctly. Thus the suggested
signed surplus remains potentially useful elsewhere, but is not a missing
assumption here.

## Scope

The theorem is all-order on matching uncertainty with arbitrary binary cross
orientations. It is not an all-order sharp stability theorem for arbitrary
uncertainty graphs. Any improved universal cutoff obtained from this new
input requires its own proof and endpoint audit; the frozen n>=526 theorem
is unchanged by this review.
