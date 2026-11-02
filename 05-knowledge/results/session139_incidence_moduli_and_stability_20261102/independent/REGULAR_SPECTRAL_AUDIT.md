# Independent regular-spectrum and exact slice audit

Status: VERIFIED / PROVED in the exact scopes below. The polynomials and
Hamiltonian counts were checked independently using custom standard-library
Fraction polynomial arithmetic, without importing the producer's verifier.

## General regular weighted tournament geometry

Let A be a real weighted tournament matrix of order n>=2 with every row sum
d=(n-1)/2. Put K=A-A^T. Then K1=0 and

    A|_(1-perp)=-I/2+K/2.

The real skew restriction is normal with imaginary spectrum i theta. Thus
all non-Perron eigenvalues of A have real part -1/2, and the multiset is
stable under lambda -> -1-lambda. The Perron eigenvalue d is separate.
Writing P(z)=det(I-zA) and Q(z)=P(z)/(1-dz), direct factor pairing yields

    Q(z)=(1+z)^(n-1) Q(-z/(1+z)).

This is a polynomial identity after clearing denominators, including z=-1
by polynomial continuation. Each reduced determinant root z=1/lambda obeys
|z+1|=|1+lambda|/|lambda|=1. The unreduced Perron pole is not on this circle.
This is automatic finite-dimensional skew symmetry, not a Riemann-hypothesis
argument and not a sufficient Hamiltonian-count invariant.

## Explicit regular cospectral hostile

Let K have directed circulation t on the triangle0->3->4->0 and all other
entries0. Set A=(J-I+K)/2 and conjugate by Q=I-ww^T/3 with
w=(1,1,-2,0,0). Both A and QAQ^T are valid regular weighted tournaments
for every |t|<=1. Their common determinant factorization is

    (1-2z)(1+z/2)^2[(1+z/2)^2+3t^2 z^2/4].

The independent verifier establishes coefficientwise

    H(A)=15/2+9t^2/4,
    H(QAQ^T)-H(A)=-t^4/27,
    M23(QAQ^T)-M23(A)=-t^4/54.

At t=1/4 the Hamiltonian values are489/64 and52811/6912, a difference
-1/6912. V=157/64, tau=163/64, u=57/256, v=-19/128 agree. All-one rooted
resolvents agree structurally because Q1=1. Thus the exact critical line,
pole circle, functional equation, full spectrum, and all-one scalar rooted
response together still do not recover H. Arbitrary individual vertex-root
responses are not preserved: (B^2)11-(A^2)11=5t^2/18. This is an explicit
firewall against calling the whole vertex-rooted system invariant.
For every modulus coprime to6,
the rational difference remains a nonzero unit after valid denominator
reduction; no reduction claim is made at moduli where denominators fail.

## Exact3+2 regular rank-two slice

Let K_i4=x_i, K_i5=-x_i fori=1,2,3, K45=t, with
sum x_i=t, sum x_i^2=t^2, |t|<=1. The skew matrix is v wedge w, where

    v=(x1,x2,x3,-t/2,-t/2),    w=(0,0,0,1,-1).

They are orthogonal with squared norms3t^2/2 and2, so rank is at most2,
with nonzero skew eigenvalues +/-i sqrt(3)t. Row sums vanish, and
|x_i|<=|t|<=1 ensures all adjacency weights are valid. For t=0 the entire
matrix is zero and must be handled before normalizing x/t.

For t!=0, put s_i=x_i/t. The rational parameter

    s=(-u,u(u+1),u+1)/(u^2+u+1)

covers the circle sum s_i=sum s_i^2=1 except the limit point(0,1,0).
Its denominator is strictly positive. Independent polynomial enumeration
cleared denominators and verified on this entire rational parameter

    M23=5/8+t x1 x2 x3/8,
    H=15/2+9t^2/4+t x1 x2 x3/4.

Because these are polynomial identities after clearing a positive
polynomial denominator, they hold for every real parameter. The omitted
point follows directly or by continuity, proving the full slice identity.

Let r=s1s2s3. The s_i are the three real roots of z^3-z^2-r; its
discriminant is -r(4+27r). Therefore -4/27<=r<=0. Conversely each r in
this interval gives three real roots with the required elementary sums,
so every intermediate value is attained. The exact H range is

    [15/2+9t^2/4-t^4/27, 15/2+9t^2/4].

For t!=0, the lower endpoint consists of permutations of
(-1/3,2/3,2/3) in s coordinates; the upper consists of permutations of
(1,0,0). At t=0 the endpoints coincide. This classifies this3+2 slice only,
not the entire regular orthogonal orbit or every weighted tournament.

## Receipts and hostile controls

`audit_regular_spectrum.py` checks the family, interval validity, regularity,
complete determinant factorization, reduced functional equation, and exact
H/M23 formulas. `audit_regular_circle.py` independently checks the full
rational parameter identity. Both normal and -O outputs agree byte-for-byte.
The false H coefficient-1/28 is rejected in both modes by explicit exceptions.
