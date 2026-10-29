# Supersingular Frobenius: coarse invariants versus named two-torsion

**CITED arithmetic interface; proved finite algebra; finite-exact certificates.**
This note does not prove Ogg's all-prime classification anew. It builds a
reproducible check of the arithmetic condition, and identifies exactly why
level parameters cannot replace the invariant in that condition.

For an odd prime p, the Legendre elliptic curve

    E_lambda: y²=x(x-1)(x-lambda)

is nonsingular for lambda not0,1. Its supersingular parameters are the roots
of the Hasse polynomial

    H_p(lambda)=sum_(i=0)^((p-1)/2) binom((p-1)/2,i)^2 lambda^i.

This is the classical Hasse-invariant criterion; see [Chai's Geometry and
Numbers lectures](https://www2.math.upenn.edu/~chai/lectures/geometry_and_numbers_v6_slides.pdf).
For the Ogg property and its modular interpretation, see [Aricheta,
Supersingular Elliptic Curves and Moonshine](https://arxiv.org/abs/1809.07421).
The fifteen absolute supersingular primes concern **all supersingular
j-invariants in characteristic p**, not supersingular reduction of one
chosen rational elliptic curve.

## 1. A finite algebra certificate

For p>=5 work in A_p=F_p[lambda]/H_p. The polynomial H_p is separable and
has neither0 nor1 as a root. Consequently

    J=256(1-lambda+lambda²)^3 / [lambda²(1-lambda)²]

is a well-defined element of A_p. All supersingular j-invariants lie in F_p
if and only if

    J^p-J=0 in A_p.                                      (1)

Indeed A_p is a product of finite separable fields after decomposition, so
vanishing is equivalent to vanishing at every geometric Hasse root.
Every elliptic curve over an algebraic closure of characteristic not2 has a
Legendre presentation, so this covers all supersingular j-invariants.

The minimal polynomial S_p(X) of J in A_p has as its distinct roots exactly
those j-invariants. Thus a second exact calculation is

    rational_j_count=degree gcd(S_p(X),X^p-X).

The checker constructs S_p by incremental linear algebra on1,J,J²,... in
A_p and verifies the relation directly. It does not assume a factorization
of H_p or approximate its roots. The difference between degree S_p and the
rational count is twice the number of quadratic Frobenius pairs, using the
classical fact that supersingular j-invariants lie in F_(p²).

The generators lambda->1-lambda and lambda->1/lambda permute the level-two
markings and preserve J. The checker verifies these substitutions directly.
This is the exact observation being quotiented. A rational coarse invariant
need not have a rational individual marking above it.

## 2. Two decisive small controls

At p=5,

    H_5=lambda²-lambda+1.

It has no F5 root, but its J element is zero. Thus all supersingular
j-invariants are rational, while no supersingular Legendre parameter is
rational. Replacing(1) by lambda^p-lambda=0 would reject an absolute
supersingular prime. This is a failure of marking descent, not supersingularity.

At p=37 the checker obtains

    S_37(X)=X³+23X²+5X+11
           =(X-8)(X²-6X-6) modulo37.

The quadratic factor is irreducible. Hence there is one rational j and one
nonrational conjugate pair. This provides an explicit certificate for the
first prime outside Ogg's list. The constant in that quadratic is **minus6**;
its roots are3 plus/minus sqrt15 over F_(37²).

## 3. A genuine, but carefully typed, occurrence of seven and twenty-one

At the last Ogg prime71, the supersingular j polynomial is

    X^7+48X^6+9X^5+17X^4+21X^3+65X.

All seven j-invariants lie in F71. Exactly21 supersingular Legendre
parameters lie in F71, distributed over just four of those j-invariants:

    j=17:6 parameters; j=24:3; j=40:6; j=48:6.

Here24 is1728 modulo71, explaining its smaller level-marking orbit.
The other three rational j-invariants have no rational Legendre parameter.
These are actual arithmetic counts in a single quotient map: seven coarse
objects, twenty-one rational marked presentations, and nonconstant fibres.
They do not prove an equivalence with forbidden tournament Hamilton-path
counts7 and21, and the21 parameters are not the21 unordered pairs of seven
objects. Their provenance is explicitly different.

The adjacent Richelot note uses genus-two two-torsion and15 Lagrangian
kernels. The present test concerns elliptic j-invariants and their level-two
markings. Both require Frobenius-compatible source data, but neither a shared
number nor a bare bijection identifies their arithmetic objects.

## 4. Scope and reproduction

`check_hasse_quotient.py` tests every prime5<=p<=199. It recovers exactly
5,7,11,13,17,19,23,29,31,41,47,59,71 in that finite scope; the additional
absolute primes2,3 require their own small-characteristic descriptions.
Completeness beyond199 is the cited Ogg theorem, not an inference from this
scan. The output records each exact polynomial, rational invariant count,
quadratic-pair count and rational level-parameter count. All arithmetic is
modular integer arithmetic using the standard library.

The script compares two quotient-Frobenius calculations and verifies S3
marking substitutions, separability, the classical supersingular degree
formula, and the p5/p37/p71 hostiles. These finite tests do not infer
rank over Q, count rational Richelot isogenies without Frobenius guards, or
settle any general distribution or rank conjecture.
