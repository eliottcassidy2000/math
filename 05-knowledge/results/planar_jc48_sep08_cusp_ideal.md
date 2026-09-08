# Poles close the elliptic boundary of the cusp-ideal flow

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [independent analytic and source audit](planar_jc48_sep08_cusp_ideal_audit.md)
accepts the complete proof and byte-identical normal/optimized replays.
JC(2) remains OPEN. This is a statement about specified scalar-time formal
Hamiltonian operations, not arbitrary Keller maps or compositions with
unrelated invariants. No theorem ID or external priority is claimed.

## 1. Statement and inheritance

Let K have characteristic zero. Put

    Delta=p^3-y^2,
    p=t(1+x^2t),   y=xtp,
    {F,G}=F_x G_t-F_t G_x.

For every nonconstant polynomial

    S=c0+I,       0!=I in Delta K[p,y],                 (1)

its literal scalar-time Hamiltonian exponential on K[x][[t]] exists.
The theorem is: **at every nonzero scalar time lambda in K,
the two coordinate images cannot both be rational**. Equivalently,

    not (exp(lambda delta_S)(p) in K(p,y)
         and exp(lambda delta_S)(y) in K(p,y)),         (2)

with rationality compared to the literal source flow as in Section 5.
This excludes rational source automorphisms at such times. It does not
assert that either individual image is transcendental.

The geometric input is slightly different from the previous carrier:
every geometric generic component of I has genus at least one. If any
component has genus one, the induced rational vector field on each such
component has a genuine pole. The latter information pays the elliptic
case that genus alone leaves open.

The closest mechanism is the actual exceptional-cover argument of
[the universal genus theorem](planar_jc48_sep07_universal_genus.md), with
[its independent audit](planar_jc48_sep07_universal_genus_audit.md).
Its source-invariant consumer extends the already proved
[incoming nonrational-time theorem](planar_jc_long_20260906_nonrational.md).
The earlier universal result requires p Delta for its genus bound and
p^2 Delta for its completed source-preserving flow. We now retain the
missing affine divisor p=0 instead of imposing its factor in I.

The named hostile is S=Delta: its generic proper curve has genus one, so
it refutes a genus-at-least-two assertion for all of (1). It does not
refute (2). The corrected near miss is to invoke finite automorphisms of
an unmarked elliptic curve, or to assume a genus-one selfmap has degree
one. Translations and duplication prevent those inferences. The least-used
sidecar is the pole divisor of the Hamiltonian vector field at the
otherwise discarded affine points p=0.

The five live concepts are the exceptional cyclic cover, geometric
components, the boundary vector field, finite relative constants, and
literal scalar-time transport. The connection sends a putative rational
time to a selfmap of a generic invariant curve. It preserves both the
invariant and commutation with its vector field. Discarding the affine
model loses p=0; retaining the vector-field pole divisor restores the
information needed at genus one. The cheapest tests are Delta, powers of
Delta, and an actual elliptic duplication map.

The proof is first over C. Section 7 gives descent over arbitrary K,
including a field containing an independent boundary-time parameter.

## 2. The exceptional cover now has genus at least one

The least nonzero weight-(2,3) part of I has the form

    C p^A y^B product_(i=1)^r(p^3-lambda_i y^2)^e_i,     (3)

where A,B>=0, r>=1, e_i>=1, C!=0, and the lambda_i are distinct nonzero
constants, one equal to 1. The latter condition follows from Delta|I.
The removed p hypothesis means that A may now vanish. Set

    E=sum e_i,  M=2A+3B+6E,
    n0=A+B+2E,  ninf=A+B+3E,
    d=gcd(A,B,e_1,...,e_r).

For completeness, the inherited map is the ordinary blowup chart
p=v^2w,y=v^3w. It gives exactly

    I=v^M(Phi(w)+v Psi(v,w)),
    Phi=C w^n0 product_i(w-lambda_i)^e_i.

Let Sigma be the sphere with disjoint small open disks removed about
0,infinity,lambda_1,...,lambda_r. On a uniform neighborhood of Sigma, Phi is
a unit. The binomial coordinate

    vtilde=v(1+v Psi/Phi)^(1/M)

is single-valued and fibrewise invertible near v=0: choose the root near
1, with its argument uniformly in the unit disk. The derivative is
uniformly close to 1, so a uniform disk can be chosen on which the map is
injective and its image contains another uniform disk. Therefore, for
all sufficiently small c!=0, all M roots of

    vtilde^M Phi(w)=c,       w in Sigma                (4)

form an actual compact bordered surface in I=c. Blowdown is injective
there, since v,w are nonzero and v=y/p,w=p^3/y^2. This is an embedding of
a surface with boundary; the compact curve obtained by capping its
boundary is used only to compute its genus.

The cyclic increments are -n0,-e_i,+ninf. They sum to zero and generate
d connected components. The equality

    gcd(M,n0,e_i)=gcd(A,B,e_i)=d

follows by using M-2n0=B+2E and then n0-B-2E=A. The compactification of
each component has genus g given by

    2d(g-1)=rM-gcd(M,n0)-gcd(M,ninf)-sum_i gcd(M,e_i).  (5)

This is Riemann--Hurwitz for the cyclic cover. Divide all A,B,e_i by d;
the component genus is unchanged and the new gcd is 1.

If B>0, the elementary bounds

    gcd(M,n0)<=B+2E,  gcd(M,ninf)<=B,
    sum_i gcd(M,e_i)<=E

make the right side of (5) at least

    rM-2B-3E >= 2A+B+3E >= 4.

If B=0, then ninf=M/2 and gcd(M,n0)<=2E. The right side is at least

    (2r-1)A+6(r-1)E.

This is positive whenever A>0 or r>=2, and is even, so is at least 2.
The only remaining type is A=B=0,r=1. Its primitive data are

    e_1=1, M=6, n0=2, ninf=3,

and (5) equals zero: g=1. Since one lambda is 1, this exceptional initial
form is exactly C Delta^e before primitive normalization.

Thus every central component has genus >=1, and has genus >=2 unless
the initial form is a pure power of Delta. This last condition is
necessary, not sufficient, for the whole generic curve to have genus one.

The local/global transfer is unchanged by A=0. Resolve the polynomial
pencil on a proper surface and take its finite Stein base. Away from
finitely many base values the proper component curves are smooth and
have constant genus. Choose small c avoiding those values; an embedded
genus-g bordered surface forces its containing proper component to have
genus at least g, by its nondegenerate handle intersection form. The
polynomial I-c is irreducible over C(c) before algebraic extension:
it is primitive and irreducible, being linear in c in C[p,y,c]. Its
geometric components are conjugate. Consequently all their smooth proper
genera agree, and the bound applies to every geometric generic component.
This retains outer powers and does not assume geometric integrality.

In particular p|I gives genus >=2. Delta=c is the actual sharp genus-one
boundary, since y^2=p^3-c is a squarefree cubic for transcendental c.
Each geometric component of Delta^e=c has the same elliptic genus.

## 3. The missing affine divisor supplies a pole on every component

Suppose p does not divide I. The polynomial

    f(y)=I(0,y)=-y^2 R(0,y),       I=Delta R,

is nonzero, nonconstant, and vanishes at y=0. Over the algebraic closure
of C(c), every root of f(y)=c is nonzero and simple. Indeed c is
transcendental, so it is neither zero nor a critical value of the fixed
polynomial f; equivalently gcd(f-c,f')=gcd(f-c,y)=1.

At a point (p,y)=(0,y0) of the generic fibre, I_y=f'(y0)!=0. The fibre
is smooth there and p is a uniformizer on it. The *actual* Poisson
identity and vector-field coordinate are

    {p,y}=-Delta/p,
    delta_I(p)=-Delta I_y/p.                           (6)

The numerator of (6) at that point is

    y0^2 f'(y0) != 0.

Hence the rational vector field has a simple pole in the local parameter
p. The factor 1/p is not a removable coordinate artifact: p is itself a
uniformizer on the smooth curve. For I=Delta, the leading coefficient is
-2y0^3, providing the simplest concrete check.

At least one geometric generic component contains such a point. All
components are conjugate over C(c), and both p and the derivation are
defined over that field. Conjugation transports the point and its pole
to every other component. Thus *each* geometric generic component has
a nonempty pole divisor whenever p does not divide I.

Combining Sections 2 and 3, every geometric generic component falls in
one of the sufficient cases

    genus >=2; or genus 1 with a nonempty vector-field pole divisor. (7)

No elliptic classification of all possible R is needed.

## 4. A pole-preserving elliptic map has finite order

Here is the curve lemma used in (7). Let C be a smooth proper connected
curve over an algebraically closed characteristic-zero field. Let V be
a nonzero meromorphic vector field, and let h:C->C be a nonconstant map
satisfying

    dh(V)=h^*V.                                       (8)

If g(C)>=2, Riemann--Hurwitz gives deg(h)=1 and the automorphism group
is finite. If g(C)=1, Riemann--Hurwitz first gives that h is étale;
it does **not** yet bound its degree. The cited formula and its tame
ramification hypotheses are in
[Stacks, Section 53.12](https://stacks.math.columbia.edu/tag/0C1B).

Since h is étale, dh is an everywhere invertible map of tangent lines.
Taking orders in (8) therefore gives

    div(V)=h^* div(V).

If P is the effective pole divisor, this implies P=h^*P. When P!=0,
its degree is positive, so deg(P)=deg(h)deg(P) forces deg(h)=1.
Now h belongs to the automorphism group preserving the nonempty finite
set Supp(P), a finite group. To see finiteness without suppressing the
elliptic translations, take the kernel of its permutation action on
Supp(P). It fixes a chosen point P0. The point stabilizer is a closed
subgroup of the finite-type automorphism group, with tangent space
H^0(C,T_C(-P0))=0: on a genus-one curve T_C has degree zero, and twisting
by -P0 gives negative degree. Thus the stabilizer is zero-dimensional
and has finitely many geometric points; the finite permutation quotient
does also. These are the automorphism-scheme and tangent-space facts of
[Stacks, Section 109.7](https://stacks.math.columbia.edu/tag/0DST), applied
to the additional point-stabilizer condition. This proves the lemma.

The pole hypothesis is load-bearing. An elliptic curve with a holomorphic
translation-invariant vector field admits translations of infinite
order commuting with that vector field. Duplication is a second warning:
on y^2=p^3-1 its p-coordinate is

    (p^4+8p)/(4(p^3-1)),

and it has degree four. It does not preserve the cusp-flow pole support:
p=-2 maps to p=0, although the original vector field is regular at the
two points with p=-2. The finite verifier checks this actual curve map.

## 5. Literal convergence, log comparison, and infinite order

For any I in (1), its pullback has t-order at least three, because

    Delta=t^3(1+x^2t)^2.

Thus I_t has order >=2 and I_x has order >=3. The derivation

    delta_I=I_t partial_x-I_x partial_t

raises t-adic order by at least two on K[x][[t]]. Its exponential converges
coefficientwise for every scalar lambda, defines a continuous ring
automorphism, obeys the additive scalar-time group law, and has inverse
at -lambda. This proves existence in the literal source completion even
for generators outside the universal source carrier.

For the rationality argument use a different, explicitly compared chart:

    p=s^2+tau,  y=sp,  Delta=tau p^2,
    {F,G}=tau(F_s G_tau-F_tau G_s).

Write I=Delta^e T, with e>=1 and Delta not dividing T. Then

    I=tau^e W(s)+O(tau^(e+1)),
    W=s^(4e) T(s^2,s^3) != 0.                         (9)

Nonvanishing follows from the cusp restriction kernel being (Delta).
The derivation raises tau-order by at least e on K(s)((tau)), for every
integer initial order. It therefore defines scalar-time automorphisms
there, preserves K[s][[tau]], fixes I, and commutes with delta_I.
The first displacement is

    exp(lambda delta_I)(p)-p
      =2lambda e s W(s) tau^e+O(tau^(e+1)).             (10)

Every second or later iterate starts at order >=2e>=e+1. In particular
(10) is nonzero for lambda!=0 and remains nonzero after replacing lambda
by any positive integer multiple. Every nonzero scalar time has infinite
order on the embedded rational functions unless rationality itself fails.

The fixed-input comparison is exact. The map

    sigma:K[[s,tau]] -> K[x][[t]],  s->xt, tau->t

is injective: s^i tau^j becomes x^i t^(i+j), with distinct exponent pairs.
There are finitely many pairs at any fixed t order. The chain rule gives

    (sigma F)_x=t sigma(F_s),
    (sigma F)_t=x sigma(F_s)+sigma(F_tau),

hence sigma intertwines the log bracket with the literal bracket. It
therefore intertwines their convergent exponentials on p,y. If a literal
image were a rational function of x,t, substitute x=s/tau,t=tau and clear
powers of tau to write that rational function as g(s,tau)/h(s,tau), with
polynomial g,h and h!=0. Multiplying by h and using sigma's injectivity
forces the corresponding log image to be that rational function.
This is only a comparison for the fixed inputs, not an identification of
completions or a claim of injectivity for the whole old Delta-adic carrier.

Finally K(x,t)=K(p,y), with

    x=yp/Delta,  t=Delta/p^2.

Since the formal maps are injective and extend to fractions, simultaneous
rationality of x,t is equivalent to simultaneous rationality of p,y.

## 6. Relative constants and the rational-time contradiction

Work over C and suppose both images in (2) are rational. The log
exponential then restricts to an injective field endomorphism

    sigma:L=C(p,y) -> L,

which fixes E=C(I) and commutes with delta_I. Rationality of its inverse
has not been assumed. Let E0 be the algebraic closure of E in L. This is
a finite extension, and sigma restricts to an E-automorphism of E0:
its image is contained in E0 and injectivity plus finite dimension gives
surjectivity. Some positive power sigma^n fixes E0 pointwise. Also delta_I
kills E0, since it kills E and algebraic minimal polynomials are separable.

The extension L/E0 is a regular one-variable function field. After an
algebraic closure of E0 it is the function field of a smooth proper
geometrically connected curve. The endomorphism sigma^n induces a
nonconstant selfmap h of that curve. This is the standard contravariant
[curve/function-field correspondence, Stacks Theorem 53.2.6](https://stacks.math.columbia.edu/tag/0BY1).
Commutation with delta_I becomes exactly (8), as can be checked by
applying both sides to any rational function on the curve.

Section 2 gives genus >=1. If its genus is one, p cannot divide I, so
Section 3 supplies the required nonempty pole divisor. The curve lemma
therefore makes h finite-order in both cases of (7). A further positive
power of sigma is the identity on L. The scalar-time group law would
then make a positive integer multiple of lambda act trivially on p,
contradicting (10). This proves (2) over C.

The finite E0 step is necessary for composite invariants. For example
I=Delta^e has e elliptic geometric generic components, not an integral
geometric generic fibre. The proof preserves the component field rather
than wrongly applying a curve theorem over C(I) without separating it.

## 7. Arbitrary characteristic-zero fields and the fixed-H consequence

The topological genus proof over C extends to every characteristic-zero
field K by descending the finitely many coefficients of I to a field
K0 finitely generated over Q and embedding K0 into C. Smooth proper
geometric component genera are unchanged by algebraically closed field
extension. For a putative rational time, add to K0 the scalar lambda and
the finitely many coefficients of the two proposed rational images.
After clearing denominators, the formal-flow identities are coefficientwise
identities over K0 obtained by field operations and factorial denominators.
They survive the embedding into C; lambda and the polynomial denominators
remain nonzero. This contradicts the complex case. One does not embed the
entire possibly much larger K into C.

An independent boundary-time parameter zeta may belong to K, or K may
be extended by it. Then a nonzero j(zeta) is a scalar time in the field.
It is unrelated to the source coordinate s and is not a source-dependent
time. Specializations at which j vanishes are identity times, excluded
by hypothesis; the proof does not use a scalar group law for variable
source-dependent times.

Here is the exact connection to the fixed-source problem. The inherited
[fixed-H carrier theorem](planar_jc48_sep06_hamiltonian.md), Section 3,
for G_H=-u/2+H(p,y), u=y^2/Delta, says

    {G_H,S} in K[p,y]
    iff S=c0+Delta R and p divides J_(p,y)(H,S).        (11)

Every nonconstant infinitesimal preserver in (11) is covered by (2),
without a condition on H. This does not claim that its full flow preserves
the affine source form, all depth constraints, or the old completed carrier.
The literal K[x][[t]] flow exists by Section 5; those further preservation
properties are separate obligations.

There is a useful narrower genus strengthening. If H_p(0,0)!=0 and (11)
holds, the exceptional initial form C Delta^e is impossible. For such an
initial form, higher weighted terms give

    S_y(0,y)=2e C(-1)^e y^(2e-1)+O(y^(2e)),
    S_p(0,y)=O(y^(2e)).

Indeed a higher-weight monomial p y^b has 2+3b>6e and hence b>=2e;
a higher pure y term has exponent >=2e+1. Thus the coefficient of
y^(2e-1) in J(H,S)|_(p=0) is

    2e C(-1)^e H_p(0,0) != 0,

contradicting (11). Such fixed-H preservers have every generic component
of genus >=2. The actual reconstructed source has H_p(0,0)=-3.
The condition is essential for this stronger *genus* conclusion:
H=0,S=Delta satisfies (11) and has genus one, while still failing all
nonzero rational scalar times by the main result.

The theorem already includes every polynomial outer composition f(S)
with nonconstant f, since f(S)-f(c0) remains a nonzero Delta multiple.
It makes no assertion about rational expressions outside K+Delta K[p,y],
arbitrary formal infinite tails, or compositions with different invariants.
The rational Hamiltonian flow of x^2t remains an outside-carrier hostile.
Constant S and time zero give the identity and are the stated exceptions.

## 8. Reproducible finite controls and audit boundary

The [standalone source](../../04-computation/planar_jc48_sep08_cusp_ideal.py)
imports no inherited mathematical implementation. Its explicit universe is
975 multiplicity tuples: A,B=0..4, one to three distinct binomial factors,
and each multiplicity 1..3. It checks the true component gcd, primitive
normalization, and the exact genus-one exceptional type. A separate 54
literal sheet-permutation covers reconstruct orbit sizes and ramification.

Further controls check actual blowup substitutions, bracket signs, source
and log finite flows for two generators through order 8, scalar group
law, inverse, commutation, the leading displacement, eight actual generic
p-zero pole families including composite Delta powers, the fixed-H
coefficient obstruction, a fixed-H preserver outside the p ideal, and an
elliptic duplication map that fails to preserve the pole support.
The non-LND rational-flow hostile is retained with its actual symplectic
identity. None of these finite checks supplies the local/global surface
transfer, the curve-map lemma, or arbitrary-field descent by itself.

    python3 -B 04-computation/planar_jc48_sep08_cusp_ideal.py
    python3 -B -O 04-computation/planar_jc48_sep08_cusp_ideal.py

The matching [.out](planar_jc48_sep08_cusp_ideal.out) records 6,956 always-active
gates. Normal and optimized output are byte-identical; the frozen output is the same
513 bytes. The semantic digest is
`a191e56d20f8f29a1d7ad4597fd5b8032b496039ccbaf30b27e5dca405426a20`.

- Source SHA256: `a7f9ef420b8782026c5df71fa361a981126283e14e5f67440e1ed6ba2ff5c1e1`.
- Output SHA256: `08a1fb64c65b3a4c1330f2d0b3df082ada0220e99ab33f6d97f73d1c8b7d69fc`.

The independent audit accepts this result. The frozen source and output
retain their original RESERVED text as a record of their production status;
proof status is governed by this note and the linked analytic audit.
