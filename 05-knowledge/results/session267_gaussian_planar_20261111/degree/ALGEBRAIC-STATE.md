# How many Gaussian scalar tests can an all-row repair need?

**PROVED over C.** This is an obstruction to universal polynomial
zero-testing, together with an exact Gaussian encoding that retains the
necessary growing state. It does not assert that actual Keller sources
fill the ambient row space.

## 1. The algebraic zero-test lower bound

Let V be an N-dimensional complex vector space, and let U be a linear
subspace of codimension q. Suppose polynomial observables
`sigma_1,...,sigma_r` have the exact test property

`v in U  iff  sigma_1(v)=...=sigma_r(v)=0`.

Then **r≥q**.

Proof. Identify C[V] with C[z_1,...,z_N]. The common zero set is U, so
the radical of `(sigma_1,...,sigma_r)` is the prime ideal I(U), of
height q. Krull's height theorem bounds the height of every minimal
prime of an r-generated ideal by r. Therefore q≤r. Equivalently, r
equations cannot cut an affine linear space by codimension greater
than r over an algebraically closed field. QED.

The bound is attained by q independent linear coordinates on V/U.
The same proof applies when the test value is sigma(v)=sigma(0), by
subtracting constants. Thus polynomial nonlinear observations do not
evade this particular lower bound.

The hypothesis is important: arbitrary nonalgebraic encodings are not
covered. A Hermitian norm `sum |d_j|²` tests q complex debts with one
nonnegative real number, but it uses conjugation and is not a complex
polynomial observable. In contrast, `sum d_j²` has nonzero complex
zeros already at `(d_1,d_2)=(1,i)`.

## 2. Application to the inherited planar row flag

For row m≥5, set `d_m=floor((m-1)/2)+2` and let

`V_m=C[x]_(≤m)`,
`U_m=L_m(C[x]_(≤d_m))`,
`L_m(theta)=-x theta+(x²+6)theta'/(2m)`.

Session254 proves the source-depth bound, and the nonzero leading
coefficients of L_m on monomials of degree≤d_m show that its image
has dimension d_m+1. Consequently

`codim_(V_m) U_m=m-d_m=floor(m/2)-1`.

Any universal polynomial scalar-test bank recognizing this image
therefore requires at least floor(m/2)-1 equations. The same-row
Hamiltonian source packet pays exactly one quotient direction, leaving
a repaired image of codimension floor(m/2)-2. Recognizing that enlarged
image still needs at least floor(m/2)-2 polynomial equations.

At row16 these lower bounds are seven before source payment and six
after it. Their divergence with m precludes a fixed finite bank of
polynomial scalar moments that decides unrestricted row repair at
every level. This is stronger than a linear rank count but narrower
than a no-compression theorem: source restrictions, rational generating
functions, auxiliary structure, and a growing marked state can all help.

It also distinguishes two bounds in this session. The 2m+1 lower bound
in [the rotation note](../planar/PROOF.md) concerns retaining EVERY
future multiplication observable on all binary forms. The bound here
concerns recognizing ONE specified repair subspace. Different tasks
have different information requirements.

## 3. An optimal marked Gaussian ledger

Let Z be a normalized circular complex Gaussian and W=conjugate(Z),
so `E[Z^a W^b]=a!` if a=b and is zero otherwise. A q-coordinate debt
vector `(d_0,...,d_(q-1))` has the encoding

`D(Z)=sum_(j=0)^(q-1) d_j Z^j`.

Then

`E[W^j D(Z)]/j! = d_j`,
`E[exp(sW) D(Z)] = D(s)`.

The second identity is formal: after expectation only finitely many
terms survive. Thus a single *marked generating polynomial* retains
the entire ledger exactly, while it still has q independent coefficients.
It is not a one-scalar compression. Multiplication of two debt encodings
gives the ordinary convolution of their coefficient arrays, compatible
with a recursive composition tree when that multiplication is the
correct native operation.

The simplest hostile is a nonzero D with no constant term. Every
unmarked moment E[D(Z)^m] vanishes. Hence even ALL power moments of
this encoding fail to detect any of its nonconstant debts. GMC(2)
correctly says it lies in the one-sided nullcone; it does not say
the encoded debt vector is zero. The marked observers are exactly
the information that an invalid bridge would discard.

## 4. A useful conjectural target

The promising target is not bounded universal scalar memory. It is an
explicit *actual-source generating carrier* whose coefficients are the
whole debt ledger, whose operations are proven to preserve realizability,
and whose all-order regularity can be decided by finitely many algebraic
conditions on its numerator and denominator. The rational-carrier
extension in the planar lane pursues precisely this alternative.

The height bound is general commutative algebra; its application and
the Gaussian ledger are proved here. These results impose requirements
on a GMC-to-JC bridge, not an obstruction to JC itself.
