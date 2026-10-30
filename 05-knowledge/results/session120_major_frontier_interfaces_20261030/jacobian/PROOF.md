# A complete affine-fibre Hamiltonian image classification

**Status: PROVED elementary, characteristic zero; FINITE-EXACT regression checks.**
This is a restricted planar polynomial calculation, not a proof of JC(2),
not a reduction of arbitrary Keller maps to this form, and not a novelty claim.
No purported higher-dimensional Jacobian counterexample is used.

**Inheritance correction after targeted canon search:** the nonzero-A
trichotomy and rational-primitive lemma below are reconstructed from existing
THM-3326 (`linear-in-z-unit-response-trichotomy-and-jet-torsion`); the general
one-root valuation without a smoothness assumption is already in THM-3348
(`linear-z-generic-puncture-response-and-one-root-valuation`). Their notation
P=f(x)+g(x)z differs only by swapping variables (and bracket sign). Thus
Sections 1–5 are an attributed unified reconstruction, NOT a new Jacobian
case or a new response theorem. Sections 7–8 record the zero-coefficient
boundary, kernel jump, and a precise failure of regular parameter lifting.
No literature-priority claim is made for those elementary extensions either.

## 1. Exact statement

Let K be an algebraically closed field of characteristic zero, A,B in K[y],
A nonzero, and P=A(y)x+B(y). Put D_P(Q)={P,Q}=P_x Q_y-P_y Q_x.
Identify K[T] with K[P], and define

    I_P = {R(T) in K[T] : D_P(Q)=R(P) for some Q in K[x,y]}.

This is an ideal because D_P(P)=0. The following list is complete.

1. **A=a is constant, a nonzero.** I_P=K[T], and all companions are

       Q=(y/a) R(P)+H(P),       H in K[T].

2. **A=a(y-b)^m, m>=2.** Write z=y-b and beta=B(b). If B is nonconstant,
   set s=ord_b(B-beta) and e=ceil((m-1)/s). If B is constant, set e=1.
   Then

       I_P=(T-beta)^e K[T],
       Q= -R(P)/(a(m-1)z^(m-1)) + H(P),    H in K[T].       (1)

   The displayed expression is polynomial exactly for R in that ideal.

3. **Every other nonconstant A.** I_P=0. This includes a linear A and
   every A having at least two distinct roots, even if every root is multiple.

In all cases ker D_P=K[P]. Thus the formulas also give all solutions, rather
than only an example for each admissible right-hand side.

## 2. Proof: separate integration from polynomial descent

The substitutions T=A(y)x+B(y), x=(T-B(y))/A(y) identify
K(y)[x] with K(y)[T]. The chain rule gives

    D_P = A(y) (partial/partial y at fixed T).                 (2)

Write a polynomial Q in these rational fibre coordinates as
Q=sum_j q_j(y)T^j, with each q_j rational. Its derivative is zero iff all
q_j are constant, proving the kernel assertion.

If R is nonzero and D_P(Q)=R(P), choose any nonzero coefficient r_j of R.
Equation (2) gives q_j'=r_j/A. Hence 1/A must have a rational primitive h.
Conversely, once h exists, all solutions in K(y)[T] are

    Q=h(y)R(T)+H(T),       H in K[T].                         (3)

This coefficientwise argument is important: H cannot hide rational y-poles.
It remains to classify rational primitives and determine when (3) descends
back to K[x,y].

## 3. Reciprocal-polynomial rational primitives: one root or none

**Lemma.** For nonzero A in K[y], 1/A has a rational primitive iff A is
constant or A=a(y-b)^m with m>=2.

The stated cases have primitives y/a and
-1/(a(m-1)(y-b)^(m-1)). Conversely, a rational derivative cannot have a simple
pole: differentiating a Laurent series creates no exponent -1. Thus a
primitive is impossible if A has a simple root. Suppose every root has
multiplicity m_i>=2, let r be the number of distinct roots, and N=sum m_i.
A rational primitive h would have poles exactly at those roots, of orders
m_i-1, and no other finite poles. Its total finite pole order is N-r.
At infinity h'=1/A has order N in the coordinate 1/y. Integrating its
leading Laurent term in characteristic zero shows that h has a finite
limit h_infinity and h-h_infinity has a zero of order exactly N-1 at
infinity. A nonzero rational function with total pole order N-r cannot
have a zero of order greater than N-r. Therefore N-1<=N-r, so r<=1.
This proves the lemma, including the impossibility of cancellation among
several multiple-root poles. The degree-one case already has a simple pole.

## 4. The exact descent exponent

In the only nonconstant surviving case, A=a z^m, the primitive's denominator
is z^(m-1). Therefore (3) is polynomial iff

    ord_z R(a z^m x+B(y)) >= m-1.                             (4)

Reduction modulo z^(m-1) removes the a z^m x term. Thus (4) is equivalent
to ord_z R(B(y))>=m-1. If B is nonconstant and R has root multiplicity k
at beta=B(b), this order is exactly sk (and is zero for k=0). Hence
k>=ceil((m-1)/s). If B is constant, any positive k suffices, since then
P-beta=a z^m x; k=0 fails. This gives e=1 and proves the entire list.

The infinity argument and the valuation argument are distinct gates:
first rational integration, then descent to the polynomial ring. Merely
checking that every root of A is multiple does not pass the first gate.

## 5. Sharp consequences and examples

**Keller corollary.** If {P,Q}=c in K* and P is affine in x, then A=a is
constant and Q=(c/a)y+H(P). Conversely these are all Keller pairs in this
class. They are automorphisms: given (U,V)=(P,Q), recover

    y=(a/c)(V-H(U)),       x=(U-B(y))/a.

No degree bound on B or Q is required. Finding coordinates which make an
arbitrary Keller component affine in one variable is NOT supplied.

**Same A, changing the vanishing order of B changes the exact image ideal.**
For A=y^5 and B=y^s, the ideal is (T^ceil(4/s)) for every s>=1.
For B=0 it is (T), not the zero ideal. Thus the defect is not determined
by the support of A or its root multiplicity alone.

**A non-Keller self-scaling Hamiltonian.** For P=y^2 x+y, take
Q=-xy-1. Then {P,Q}=P. Nonzero Hamiltonian image on K[P] does not imply
that a constant belongs to the image.

**Repeated roots are not enough.** For P=y^2(y-1)^2 x+B(y), for every B,
there is no Q with {P,Q}=R(P) for any nonzero R. All roots of A are multiple,
yet the rational-primitive obstruction kills the entire ideal.

**A sharp threshold family.** For P=y^m x+y, m>=2, the least positive
power of P in the image is P^(m-1), and its explicit preimage is
-P^(m-1)/((m-1)y^(m-1)). Every smaller nonnegative power is excluded.

## 6. Connection to the inherited frontier and limitations

Session 113's `jacobian/README.md` gives the omitted-character-bracket
identity and a hostile example to Reynolds projection. Session 114's
`jacobian/README.md` classifies one-coordinate mixed-sign torus symmetry.
The inherited THM-3326/3348 sector used here assumes no torus symmetry: A and B may have arbitrary polynomial
coefficients, but the displayed affine-in-x condition is essential.
This adds a different exact image criterion (rational integration plus
polynomial descent), not a way around either prior hypothesis.

The current frontier in `00-navigation/CURRENT-FRONTIER.md` lists unresolved
planar source and chart-entry problems. None is claimed closed here.
The mixed-weight local discussion in THM-4727 is not a dependency, and
no higher-dimensional claims in that file are endorsed or used.

`check_image.py` checks exact Poisson identities, complete bounded image
intersections by independent coefficient-matrix ranks, sharp missing powers,
and residue controls. The proof above, not finite checks, establishes the
all-parameter classification. Run under Python and Python -O; no assert gates.

## 7. The zero-coefficient boundary: a critical-value radical

The hypothesis A nonzero is meaningful: the kernel jumps at A=0. For a
nonconstant polynomial P=B(y), put

    C_B(T)= product (T-v),

where v ranges over the DISTINCT critical values B(a) at roots a of B'.
An empty product is 1. Then the complete classification is

    I_B=C_B(T)K[T],       ker D_B=K[y],
    Q= -x R(B(y))/B'(y) + H(y),       H in K[y].             (5)

Indeed D_B=-B'(y) partial_x. A polynomial solution exists iff B' divides
R(B). If a root a of B' has multiplicity m, characteristic zero gives
ord_a(B-B(a))=m+1. Thus divisibility at a is equivalent to R(B(a))=0:
any positive multiplicity of that target root suffices. Repeated critical
values are counted once, not once per source critical point. The displayed
formula and kernel follow by integration in x. This proves (5).

Equivalently C_B is the monic squarefree part of the resultant
Res_y(B(y)-T,B'(y)), up to the irrelevant nonzero scalar in the resultant.
This is a computational formula, not an additional hypothesis.

Examples: B=y^n, n>=2, gives C_B=T, regardless of the derivative's
multiplicity n-1. For B=(y^2-1)^2 the critical points -1,0,1 have only the
two critical values 0 and 1, so C_B=T(T-1), not T^2(T-1). The companion
for the generator is -x B(B-1)/B'. A nonzero constant determinant is
possible precisely when B is linear, as expected.

Combining Sections 1 and 7 classifies the image intersection for EVERY
nonconstant P with degree_x(P)<=1. It does not classify maps of higher
x-degree. The constant P case has zero Hamiltonian derivation and is
excluded from the identification K[T]=K[P].

## 8. Fibrewise solvability does not give a regular parameter companion

This is a concrete obstruction to a missing implication, not a JC result.
Treat t as a parameter, and use the Poisson bracket in x,y only. Put

    P_t=t y^5 x+y^3.

For every lambda in K there is a polynomial Q_lambda in K[x,y] such that
{P_lambda,Q_lambda}=P_lambda^2. For lambda nonzero take

    Q_lambda= -lambda*y^6*x^2/4 - y^4*x/2 - y^2/(4*lambda),

and for lambda=0 take Q_0=-x*y^4/3. These identities follow by direct
polynomial differentiation. Nevertheless there is NO polynomial
Q in K[t,x,y] satisfying {P_t,Q}=P_t^2.

**Proof of nonlifting.** Work over K(t). The complete companion formula in
Section 1 forces every companion to have the form

    Q= -t*y^6*x^2/4 - y^4*x/2 - y^2/(4*t) + H_t(P_t),
    H_t(T) in K(t)[T].

Expand the finitely many coefficients of H_t as Laurent series at t=0.
If its smallest t-order is -r<0, its leading coefficient is a nonzero
polynomial h(T). Substitution T=P_t gives leading coefficient h(y^3),
which is nonzero because y^3 is transcendental over K. Thus substitution
cannot increase the t-order of H_t. Polynomiality of Q consequently
requires r<=1; at order t^-1 it requires

    h(y^3)=y^2/4.

This is impossible: every monomial in K[y^3] has degree divisible by 3.
If H_t has no pole, the displayed y^2/(4t) pole remains instead. This proves
nonexistence. The argument even excludes companions regular in t at zero
with rational parameter coefficients, not just polynomial coefficients.

**A second specialization warning: image ideals need not be nested.**
Let P_t=t y^5 x+(y^2-1)^2. At every nonzero t, Section 1 gives

    I_t=((T-1)^2),

since B(0)=1 and ord_0(B-1)=2. At t=0, Section 7 gives

    I_0=(T(T-1)).

Neither ideal contains the other: (T-1)^2 does not vanish at T=0, while
T(T-1) has only a simple zero at T=1. This is a statement about these
marked image ideals, not a blanket assertion about semicontinuity of
finite-dimensional ranks. A specialization-based argument must retain
regularity of the companion and the correct kernel, not just solvability
on each fibre or the generic annihilator.
