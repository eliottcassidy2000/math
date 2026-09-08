# Independent audit of the full cusp-ideal nonrational-time theorem

**Status: INDEPENDENT ANALYTIC + SOURCE + EXACT REPLAY AUDIT PASS.**
I read the complete [primary proof](planar_jc48_sep08_cusp_ideal.md) and
standalone implementation, independently rederived its new curve-map and
fixed-input arguments, and replayed normal and optimized modes. The
all-parameter theorem is accepted for every nonconstant S in
K+(p^3-y^2)K[p,y], in characteristic zero. Each nonzero scalar time fails
simultaneous rationality of the two source images. This includes the
entire infinitesimal fixed-H carrier; it does not integrate its additional
source-preservation or depth constraints, exclude mixed-invariant
compositions, or prove JC(2).

## 1. The local surface and all geometric components

The inherited universal-genus proof was already independently audited.
I checked the exact change here: allowing A=0 in the weighted initial
factorization leaves the actual embedded-cover construction intact.
There is still at least one binomial factor because Delta divides I,
M is positive, and the compact punctured sphere avoids every zero of Phi.
The near-one analytic binomial root makes the local equation exactly
vtilde^M Phi(w)=c on a uniform neighborhood. All M roots are small there;
there is no loss of sheets from the actual fibre. The blowdown inverse
v=y/p,w=p^3/y^2 proves embedding off the excluded axes.

The common divisor is gcd(A,B,e_i), including A=B=0. Dividing all
multiplicities by it divides the Riemann--Hurwitz numerator by the same
integer. For B>0 the primitive numerator is at least 2A+B+3E>=4.
For B=0 it is at least (2r-1)A+6(r-1)E, positive unless A=0,r=1.
The numerator is even, so positivity forces genus at least two. The
exceptional primitive cover has M=6 and exponents2,1,3, giving genus one.
Since a mandatory binomial is Delta, its initial form before primitive
normalization is exactly C Delta^e. This is necessary for whole-fibre
genus one, not sufficient: further global handles may remain.

The proof of the genus lower bound is about an embedded bordered surface,
whose handle intersection pairing persists in its ambient compact curve.
It does not embed a capped compact curve in an affine fibre. After a
proper pencil resolution and finite Stein factor, one can choose the
small value outside the finitely many exceptional values. The smooth
proper component genus is constant on the connected Stein base's smooth
open. Equivalently I-c is irreducible over C(c) before extending constants,
and its geometric components are conjugate. The bound therefore applies
to every geometric generic component, including the powers Delta^e.

## 2. Every possible elliptic component has an intrinsic pole

If p divides I, the already audited stronger genus bound applies.
Otherwise f(y)=I(0,y) is a nonzero polynomial divisible by y^2, so is
nonconstant. For transcendental c, f-c has only simple nonzero roots.
At each point (0,y0), the condition I_y(0,y0)!=0 makes p a local
parameter on the smooth curve I=c. This is the indispensable type
check on the claimed pole: a denominator in an arbitrary coordinate
would not by itself establish a vector-field singularity.

Direct chain-rule differentiation of the actual source coordinates gives
{p,y}=-Delta/p. Consequently delta_I(p)=-Delta I_y/p has a nonzero
coefficient y0^2 f'(y0) at p-order minus one. For I=Delta it is
-2y0^3. Thus this is a genuine simple pole of the vector field in a
uniformizing coordinate. At least one component contains such a point;
conjugation over C(c) sends its pole to every component, since both
the rational p and the derivation are defined over that field. A pole
cannot vanish on another component through an untracked primitive-invariant
choice. No claim that these are the vector field's only poles is used.

## 3. The elliptic curve-map argument is valid in the required direction

For a nonconstant selfmap h of a proper genus-one curve in characteristic
zero, Riemann--Hurwitz says h is etale. It does not give degree one.
Commutation with V is the equality dh(V)=h^*V of meromorphic sections
of h^*T_C. Since dh is everywhere invertible, taking orders gives
div(V)=h^*div(V). Etaleness also ensures that pulling back this divisor
retains the positive and negative parts separately. Thus the nonempty
effective polar divisor P satisfies P=h^*P. Its degree forces deg h=1.

The subgroup of automorphisms preserving the finite set Supp(P) is
finite. Its kernel in the permutation group fixes a chosen point P0.
That stabilizer is a closed subgroup of the finite-type automorphism
scheme and has tangent space H^0(T_C(-P0))=0, since its line bundle
has degree -1. A zero-dimensional finite-type group scheme has only
finitely many geometric points; the finite permutation quotient does too.
This pays the finiteness step without falsely invoking finiteness of the
whole unmarked elliptic automorphism group. For genus at least two the
usual Riemann--Hurwitz and finite-automorphism argument already applies.

I re-accessed the primary inputs on September 8: [Stacks 53.12](https://stacks.math.columbia.edu/tag/0C1B)
for the separable curve formula, [Stacks 109.7](https://stacks.math.columbia.edu/tag/0DST)
for finite-type automorphism schemes and their derivation tangent spaces,
and [Stacks 53.2.6](https://stacks.math.columbia.edu/tag/0BY1) for the
contravariant curve/function-field correspondence. The new pole-divisor
and point-stabilizer deductions above supply the additional hypotheses;
they are not quoted as an existing theorem of those pages.

The two elliptic hostiles have distinct roles. Translations show that
commutation with a holomorphic vector field need not have finite order.
Duplication shows genus one alone does not force degree one. The source
checks the actual duplication map and its coprime degree-four rational
p-coordinate. Its extra inverse image p=-2 of the pole at p=0 is a
regular point of the cusp-flow vector field. This prevents duplication
from satisfying the required pole-preservation condition.

## 4. Actual scalar time, rational comparison, and relative constants

For the literal pullback of I, t-order is at least three. The two
derivation coefficients I_t and -I_x have orders at least two and three,
respectively. Acting on any term of K[x][[t]] therefore raises its order
by at least two. At each t-order the exponential is a finite sum of
polynomials in x. It is a continuous automorphism with additive time
law and inverse at negative time for every scalar. This argument
requires neither the universal carrier nor a local-nilpotence assertion.

In the log chart, I=tau^e W(s)+higher order, W!=0, and the derivation
raises tau-order by at least e, even on K(s)((tau)). On p=s^2+tau its
first displacement is 2 lambda e s W tau^e; every later iterate has
order at least2e>=e+1. No nonzero positive integer multiple of lambda
can yield the identity in characteristic zero.

The map K[[s,tau]] to K[x][[t]] taking s to xt and tau to t is defined
on the whole displayed domain and is injective: its monomial map is
(i,j) to (i,i+j), and only finitely many pairs contribute at each t-order.
The chain rule gives the exact bracket intertwining. Both log images of
p,y have polynomial-in-s coefficients at each tau-order, so lie in this
domain. If a literal image is rational in x,t, substitution x=s/tau and
clearing a finite power of tau gives g/h with g,h polynomial and h!=0.
Multiplying before applying injectivity is legitimate even when h is
not a unit in the power-series ring. It proves the required fixed-input
rationality comparison without identifying completions.

The resulting injective field endomorphism fixes E=C(I) and commutes
with delta_I. The relative algebraic closure E0 of E in L=C(p,y) is
finite. The endomorphism acts as an E-automorphism on E0, so a positive
power fixes it pointwise. Separability also forces delta_I to kill E0.
After algebraic closure of E0, regularity of L/E0 yields a geometrically
integral curve and the induced nonconstant map. Commutation becomes
exactly the differential condition audited above. Its finite order
contradicts the explicit leading scalar displacement. The rationality
of the inverse endomorphism was never assumed.

Finally the birational inverse x=yp/Delta,t=Delta/p^2 transfers
simultaneous rationality between the two coordinate systems. Arbitrary
characteristic-zero fields are handled by the finite field of definition
of I, lambda and both alleged rational images. Embedding that field
into C preserves the cleared formal identities and nonzero coefficients.
The proof does not assume that the entire original field embeds into C.

## 5. Fixed-source scope and replay

The fixed-H infinitesimal condition is precisely S=c0+Delta R and
p|J(H,S) in the inherited carrier theorem. All such nonconstant
Hamiltonians are covered. When H_p(0,0)!=0, a pure initial C Delta^e
would leave the coefficient 2e C(-1)^e H_p(0,0) in J(H,S)|p=0.
Higher-weight terms of S_p have y-order at least2e, as do the remaining
terms of S_y, so they cannot cancel it. Thus every generic component
then has genus at least two. H=0,S=Delta checks the exact genus boundary
when that derivative condition is omitted. Full fixed-H source-form
preservation at later orders is not proved or assumed.

I read every gate in the standalone source. Its 975 multiplicity tuples
and 54 literal permutation covers check separate mechanisms; the latter
reconstruct actual signed puncture orbits and ramification without the
gcd count. Eight generic p-zero families include repeated Delta factors,
and literal/log jet engines use their separate brackets. The finite
coefficient controls support the analytic proof rather than extrapolating
it. Both independent replays are byte-identical to the frozen report.

    python3 -B 04-computation/planar_jc48_sep08_cusp_ideal.py
    python3 -B -O 04-computation/planar_jc48_sep08_cusp_ideal.py

Each mode passes6,956 always-active gates. Frozen SHA256:

    source a7f9ef420b8782026c5df71fa361a981126283e14e5f67440e1ed6ba2ff5c1e1
    output 08a1fb64c65b3a4c1330f2d0b3df082ada0220e99ab33f6d97f73d1c8b7d69fc

No producer source or output bytes were changed by this audit. Their
original RESERVED text records the status when the controls were frozen;
the accepted analytic scope is the primary proof plus this audit.
