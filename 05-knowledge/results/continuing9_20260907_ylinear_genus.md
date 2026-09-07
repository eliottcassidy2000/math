# Exact genus for every polynomial multiplier linear in y with constant slope

**Status: PROVED with cited curve inputs; INDEPENDENTLY AUDITED.**
Let K be any characteristic-zero field, let gamma be a nonzero element of K,
and let A be any polynomial in K[p]. In the actual cusp source, put

`D=p³-y²`, `I=p²D(A(p)+gamma y)`.

The generic invariant field K(p,y)/K(I) is geometrically integral. Its
smooth proper curve has genus 6 if A has degree at most one (including
A=0), and genus `2 deg(A)+3` otherwise. Every nonconstant polynomial
Hamiltonian f(I) has a nonrational specialization at every nonzero scalar
time: its two coordinate images cannot both be rational. General source
carriers, variable coefficients gamma(p), and compositions of different
flows remain OPEN. Planar JC remains OPEN.

## 1. Inheritance and the new operation

The closest proved mechanism is the completed universal carrier and its
explicit infinite-order displacement in
`05-knowledge/results/planar_jc_long_20260906_nonrational.md`.
The immediate input is the independently audited affine, univariate and
monomial genus theorem in
[continuing9_20260907_flow_carrier_genus.md](continuing9_20260907_flow_carrier_genus.md).
The source modules retain **THM-3989 /
`01-canon/theorems/THM-3989-cusp-log-laurent-conductor-and-nondividing-depth-reduction.md`**
and **THM-4308 /
`01-canon/theorems/THM-4308-source-normal-bracket-hasse-truncation-through-row-eight.md`**.

The canonical hostile is the non-LND rational flow of x²t. The corrected
near miss is replacing the primitive invariant I by f(I), whose generic
fibre can split. The least-used sidecar is the finite set of critical
values of a polynomial, which replaces a growing parameter-dependent
critical-point resultant without discarding the original invariant.

The live concepts are: actual carrier entry; primitive invariant; finite
critical values; ramification at zero; ramification at infinity; and
the nonzero first flow displacement. Increasing deg A changes the
finite branch count, while preserving the two boundary ramification
types and the field over which the flow must act. This is the reason
the affine calculation extends to every degree.

The connection maps an actual completed source flow to a selfmap of its
generic invariant curve. It preserves the invariant and the iterate;
it loses the affine presentation and formal completion. Geometric
integrality, genus, and the fixed-input comparison retain the information
needed to infer failure of rationality. The cheapest decisive test is
the full discriminant degree together with both boundary charts.

## 2. Finite critical values in characteristic zero

**Lemma.** A nonconstant polynomial J in K[p,y] has only finitely many
critical values over an algebraic closure of K.

Work over the algebraic closure. The algebraic set J_p=J_y=0 has finitely
many irreducible components and no two-dimensional component, since in
characteristic zero that would make J constant. A zero-dimensional
component is a point. On a curve component, the differential of the
restriction of J to its function field is zero. If that restriction
were transcendental, the function field would be finite separable over
K(J), and the derivation d/dJ would extend to it. The chain rule would
give both dJ=0 and dJ=1, a contradiction. Thus J is constant on each
component. There are only finitely many resulting values.

The lemma excludes critical points on the fibre J=c with c transcendental
over K. It does not claim that the critical points themselves are finite.
This distinction is essential when coefficient specialization creates
critical curves.

## 3. The complete finite branch divisor

Divide I by gamma and replace A by A/gamma, absorbing the scalar in c.
The normalized generic equation is

`F=(y²-p³)(y+A(p))+c/p²=0`.                             (1)

At p=0, every root has valuation -2/3: when the valuation is nonnegative,
the pole c/p² cannot cancel; when it is negative, y³ is the unique lowest
term except when it balances c/p². Since rational functions in p have
integral p-valuation, the cubic has no rational root even after extending
the constant field to an algebraic closure of K(c). It is therefore
irreducible and defines a geometrically integral degree-three cover.

The complete discriminant numerator is the polynomial identity

`N=p⁴ disc_y(F)`

`=4p⁷(p³-A²)²+4Ac p²(9p³-A²)-27c²`.                  (2)

It satisfies N(0)=-27c². Let m=deg A for nonzero A. If m<=1 or A=0,
its degree in p is 13 with leading coefficient 4. If m>=2, its degree
is 4m+7 with leading coefficient four times the fourth power of the
leading coefficient of A. No leading terms cancel: 2m cannot equal 3.

Every root of N is simple. A triple root of the cubic would imply
`y=-A/3` and `A²+3p³=0`. The latter is never the zero polynomial, because
p³ is not a polynomial square, so p and then y would be algebraic over K.
The value of I at such a point cannot be the transcendental c.

At a double but not triple root, translate y to express the specialized
cubic as `y²(y+q)`, q nonzero. Under a first-order perturbation epsilon V,
the discriminant changes by `-4q³ V(0) epsilon`. Thus a repeated
discriminant zero would force F_p=F_y=0 at F=0. Since
`F=-(I-c)/p²` and p is nonzero there, this would be a critical point of
I with value c, excluded by the lemma. This proves simplicity uniformly
for every polynomial A, not only generic coefficients.

The finite nonzero branch divisor therefore has n simple index-two
points, where n=13 for m<=1 and n=4m+7 for m>=2. Elsewhere away from zero
and infinity, (1) is a monic cubic with invertible discriminant, so there
are no missing finite branches.

## 4. The two boundary types do not change with degree

At zero substitute `p=r³,y=r^-2Y` and multiply (1) by r⁶. Its initial
equation is Y³+c=0, with three simple roots permuted transitively by the
cubic deck action. They descend to one index-three point, contributing
two to the ramification divisor.

For m<=1 the inherited infinity chart is

`(Y²-1)(Y+r³A(r^-2))+c r^13=0`,

obtained from `p=r^-2,y=r^-3Y`. Its initial roots are 1,-1,0; the first
two descend to an index-two point and the third is unramified.

For m>=2 multiply that chart by r^(2m-3). If a_m is the leading
coefficient of A, the equation becomes

`(Y²-1)(r^(2m-3)Y+r^(2m)A(r^-2))+c r^(2m+10)=0`.       (3)

Its initial polynomial is `a_m(Y²-1)`, with two simple roots. The deck
involution `(r,Y)->(-r,-Y)` exchanges their formal branches while fixing
the original p,y. Thus they give exactly one index-two point at infinity.

The third point can be seen without a fractional-power base change:
put `p=q^-1,y=q^-mV` in (1) and multiply by q^(3m). The initial polynomial
is `V²(V+a_m)`. Its root -a_m is simple, so it gives an unramified point
with y asymptotic to -A(p). The two charts account for degrees two and
one of the degree-three cover; no further point lies above infinity.
The ramification contribution at infinity is always one.

The total ramification is n+3. Characteristic-zero Riemann--Hurwitz gives

`2g-2=-6+(n+3)`, hence `g=(n-1)/2`.

This is genus 6 for m<=1 and genus 2m+3 for m>=2, proving the stated
classification. In particular increasing the degree of A creates more
finite branch points; it does not create a genus-zero or genus-one escape.

## 5. Actual flow consequence and stopping boundary

Every nonconstant f(I) lies in the inherited universal source carrier
`K+p²D K[p,y]`. The completed exponential and its scalar-time group law
therefore exist with the actual depth and source-preservation properties.
The invariant retained for the curve argument is I itself.

In logarithmic source coordinates `p=s²+tau,y=sp`, one has

`I=tau W(s)+O(tau²)`,

`W(s)=s⁸(A(s²)+gamma s³)!=0`.

The even polynomial A(s²) cannot cancel the nonzero odd term gamma s³.
If the first nonconstant term of f is f_k z^k, the original Hamiltonian
bracket gives

`exp(lambda delta)(p)-p`

`=2lambda k f_k s W(s)^k tau^k+O(tau^(k+1))`.

Every positive iterate has nonzero displacement. A rational specialization
would induce a nonconstant selfmap of the smooth proper generic curve
over K(I). [Stacks, Theorem 53.2.6](https://stacks.math.columbia.edu/tag/0BY1)
The genus is at least six, so Riemann--Hurwitz makes the map an automorphism.
[Stacks, Section 53.12](https://stacks.math.columbia.edu/tag/0C1B)
The geometric automorphism group is finite, contradicting the infinite
order just proved. [Stacks, Section 109.7](https://stacks.math.columbia.edu/tag/0DST)

The actual source comparison and the non-LND rational hostile remain
exactly those in the audited primary carrier-genus theorem. No conclusion
is drawn about arbitrary coefficients gamma(p), higher powers of y in
the multiplier, compositions of different flows, or arbitrary Keller maps.
Those require their own branch divisors or a different retained invariant.

## 6. Exact reproduction

The standalone source reconstructs the full discriminant and triple-root
eliminant with an unevaluated A, checks the double-root derivative formula,
and keeps c transcendental in nine exact control polynomials of degrees
zero through six, including A=0 and repeated factors. For each it checks
the entire finite branch degree and squarefreeness, both boundary types,
and the original logarithmic leading response at outer orders one and two.
The all-degree proof is the preceding argument, not an extrapolation of
these nine controls.

    python3 04-computation/continuing9_20260907_ylinear_genus.py
    python3 -O 04-computation/continuing9_20260907_ylinear_genus.py

Both modes pass 97 always-active exact gates with identical raw LF output.
The final session manifest pins source, transcript, and proof bytes.
