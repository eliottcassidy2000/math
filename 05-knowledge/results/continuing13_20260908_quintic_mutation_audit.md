# Independent referee: the univariate supplier and its quintic mutation

**Status: PASS — independent analytic audit and FINITE-EXACT controls.**
The two separately pinned primary packets are accepted without mathematical
repair. The supplier has the claimed exact whole-ansatz conditions and
fixed pointed Weyl module. Its polynomial mutation is a global submersion
with the stated full fibre, collision, unit-response and pair-degree laws.
The explicit quintic has non-isotrivial eight-punctured original-source
fibres; their projective completions have genus zero.

## 1. Targets, inherited types, and independent route

Mutation target `continuing13_20260908_quintic_mutation`:

    report a91d2aca021335eeb182d637fc61050486ce557b021a64fcdfeb7b55cf348c4a
    source 3fa9e69e59453836cc12c79ad3b0695db95ffe7abc28b9ae6d404bfb7eddc7bc
    output 185d702971314b232a076d1ae92ff64d3461927ad1c29e793a0fdc3da610aa34
    cert   a3ccce3c5aae0bf601b65bea5f7e9152d26348a45889df86c640ccb174f22b55

Supplier target `continuing13_20260908_univariate_fixed_unit`:

    report 57f9dc1edc0325afcce632611a0fa11b90ac7c52c49bea208f682a557290dee0
    source 8f0752834accd88b91fde2348f627cd6e95f97f59c3dead20e55181c042d60ef
    output 02084442a557d2cf7db9f9a61d45f266a506dfed597903f4664e18f0fa25c5a1
    cert   3e96a8c597d6328455d31556e385240678b85de65bfc5f923a57649b12f1ba01

The actual component-principal-part supplier and canonical connection in
`planar_jc48_sep06_torsion.md`, Sections2–3, remain the mathematical
dependency. They concern the original ring C[x,t], require a polynomial
unit gradient and rational constants exactly C(first function), and keep
all labelled fibre components before quotienting by the common diagonal.
Every use below pays those conditions. Neither O(W2) nor a localized
source ring replaces the response ring.

The companion imports no producer. It independently checks literal chart
identities, source jets, rational and polynomial brackets, whole fibre
inverses, actual factor divisibility, parameter collisions in exact
quadratic quotient fields, and formal Laurent-module operations. The
all-degree and all-parameter conclusions rest on the proofs below.

## 2. Independent audit of the complete univariate supplier

Put u=x-h, z=u-2h+u³t and A=u f(z), with f nonconstant. On the added
chart x=1/r,t=-r²-r⁴b, one has

    z=rM, M=-h²(3-hr)-b(1-hr)³, u=(1-hr)/r.

The unique possible negative boundary term of A is f(0)/r. Therefore
globality is exactly f(0)=0. If f=sum a_j z^j with a_0=0, its complete
extension is (1-hr) sum_(j>=1) a_j r^(j-1)M^j, whose boundary
tangential derivative is -f'(0).

In the original source u=0, the derivative A_x is f(-2h). On u!=0,
the coordinate Jacobian J(u,z)=u³ is a unit, and the derivatives of A
in (u,z) are f(z),u f'(z). Their simultaneous vanishing is exactly a
repeated root of f, realized by actual source points. Consequently, in
the global ansatz, submersion is equivalent to f squarefree and
a0=f(-2h)!=0. These conditions force h!=0 and f'(0)!=0; the boundary
test then pays every remaining point. There is no hidden generic-point
qualification in this iff.

Let Q'=f² and Q(0)=0. The original-source Jacobian is

    J(A,Q)=u³ f Q'=A³,
    J(A,Q/A³)=1.

For every c!=0 the full source fibre has ring C[z,f^-1], with u=c/f
and t=(z-u+2h)/u³. These are entire fibre inverses, because u and f
are already units when their product is c. Thus every nonzero fibre is
irreducible. The full field is C(A)(z), with nonzero derivation on z,
so its rational constants are exactly C(A).

The zero fibre has precisely q+1 components: E_u and E_rho for every
root rho of f, where q=deg f. Each equation z-rho is primitive linear
in t, because at u=0 its constant is -2h-rho!=0. The components are
irreducible, reduced and pairwise disjoint. The full source torsion
therefore has q arms, all at target value zero.

For the actual primitive Q/A³, complete principal parts are

    E_u: Q(-2h)/A³+a0/A²,
    E_rho: Q(rho)/A³.

At a simple f-root, Q-Q(rho) has order at least three. At E_u, impose
the actual source jet z=-2h+u+u³t. Then

    A=a0u+f'(-2h)u²+O(u³),
    Q(z)-Q(-2h)=a0²u+a0 f'(-2h)u²+O(u³).

Subtracting a0*A leaves order at least three, proving the entire local
formula, including absence of a simple pole. The component E_0 has
Q(0)=0 and is regular; it must remain in the component quotient.

For q>=2, at least one nonzero f-root has Q(rho)!=0. Otherwise every
root would be a triple zero of Q, since Q'=f². Squarefreeness then gives
f³|Q, contradicting 3q>2q+1=deg Q. Thus the two coefficient vectors
C3,C2 in theta=C3/g³+C2/g² are independent: C2 is supported only on
E_u and is nonzero, while C3 has a nonzero root-component coordinate.
The exact scalar order is three, and the unit generates precisely two
full principal-part arms for every admitted q>=2.

In the free C[partial]-basis C3/g,C2/g, the columns g theta and g²theta
are (-partial,1) and (1,0), of determinant -1. Direct differentiation
gives the relation R=1+partial*g+(partial²/2)g². Modulo Dg³, every
operator has unique PBW remainder a0(partial)+a1(partial)g+a2(partial)g².
Subtracting a0R eliminates its constant term. Independence of the two
pivot columns then forces both remaining coefficients to vanish if it
annihilates theta. Hence the exact all-operator ideal is

    Ann_D(theta)=Dg³+D[1+partial*g+(partial²/2)g²].

The pointed module is the same throughout this supplier, while the
ambient arm count q and rational pair degree deg Q=2q+1 vary. For
q=1, the degree inequality becomes equality and the component space has
one arm; that boundary is not part of the two-arm assertion. The
supplier's root-critical-value collision is also paid independently in
Section8 below and does not invalidate the coefficient independence.

## 3. Mutation and all-point global submersion

Now assume q>=2 and retain the supplier's admissible f,h. Let lambda
satisfy lambda*(lambda+a0)!=0, and set

    T=lambda*A+Q(z), G=1/(2A²).

Since J(A,Q)=A³ and d(1/(2A²))=-A^-3 dA, the exact original-source
Jacobian is J(T,G)=1. In particular T is submersive on A!=0, where G
is regular. At an f-root, Q'=0 and dT=lambda*dA!=0. At u=0, the
actual coordinate relation, not independent (u,z) coordinates, gives

    T_x=a0*(lambda+a0), T_t=0.

This pays every source point. Lambda=0 makes every f-root component
critical; lambda=-a0 makes the entire E_u line critical. The rational
Jacobian identity remains an identity of fields on those excluded
parameters, so it does not override their source poles and critical loci.

Both A and z are global. Since f has a simple zero at zero, Q has
order three there. Thus Q(rM) and its first boundary derivatives vanish
at r=0. The complete added-chart T has boundary tangential derivative
-lambda*f'(0)!=0. Every added point is paid, so T is a global
submersion on fixed W2.

The degree of Q is 2q+1 with leading coefficient lead(f)²/(2q+1).
The term lambda*A has source t-degree q, while Q(z) has degree 2q+1
and leading coefficient lead(f)²*u^(3(2q+1))/(2q+1). They cannot
cancel. The mutation has exactly the claimed original source degree.

## 4. Whole rational field, every source fibre and every extra component

The identities

    u=(T-Q(z))/(lambda*f(z)),
    t=(z-u+2h)/u³

recover the entire function field C(x,t)=C(T)(z). The derivation on z
is lambda*u³*f(z), which is nonzero. A nonzero multiple of ordinary
z differentiation on C(T)(z) has exactly C(T) as its rational constants.

Let z0=-2h and let S be the distinct values among Q(z0) and all Q(rho)
for f-roots rho. A=0 is the union of exactly q+1 actual source
components, and these values are their T-images. For every c, the
locus A!=0 in T=c has the exact coordinate ring

    C[z, f(z)^-1, (c-Q(z))^-1].

All inversions are legitimate in that locus; conversely its displayed
inverse formulas recover every original source point. The ring is a
nonzero localization of C[z], so is a nonempty irreducible affine curve.
For c outside S it is the entire fibre, because A=0 has no point there.

For c in S its closure is one irreducible component R_c. Every other
fibre component must lie in A=0 and is exactly one of the q+1 source
components assigned to c. No extra component is hidden by localization.
Source submersion makes every fibre reduced and its distinct components
disjoint: an intersection of two polynomial fibre factors would force
both derivatives of T to vanish. In particular R_c cannot acquire an
intersection with an A=0 component on closure; it is precisely the
regular irreducible component represented by the localized ring. G is
regular on it. This proves the complete special-fibre description even
when several root values coincide with each other or with Q(z0).

Thus if e_c is the number of A=0 components assigned to c, the fibre
has e_c+1 components and contributes e_c full source torsion arms. The
sum over all special targets is q+1. No global fibre component or
boundary point has been substituted for an original affine component.

## 5. Non-isotrivial punctures, with constant completed genus

For c outside S, the polynomial c-Q has degree 2q+1 and simple roots:
a multiple root would satisfy Q'=f²=0 and put c among the root values
already excluded. Those moving roots are disjoint from the q f-roots.
The exact fibre ring therefore describes P1 minus q fixed f-roots,
the fixed point at infinity, and 2q+1 moving roots. The puncture count
is 3q+2, and the smooth projective completion always has genus zero.

Here is a proof of non-isotriviality that does not label the punctures.
Fix any one reference puncture set P, of size r=3q+2. Any isomorphism
of the punctured curves extends to an automorphism of their smooth
projective completions, hence to a Mobius transformation. Choose three
fixed points from the f-roots and infinity; q>=2 guarantees three.
Their images in P have at most r(r-1)(r-2) choices and determine the
transformation. Thus one isomorphism class can contain at most that many
different puncture sets from this family. Distinct parameters c give
distinct sets, because after removing the fixed punctures their remaining
sets are the disjoint root sets of Q=c. Each isomorphism class therefore
contains only finitely many parameters, proving non-isotriviality on the
ordinary locus.

This concerns original affine fibres. It does not claim varying
projective genus or automatically make the same puncture count for W2
fibres. The q=1 boundary illustrates the role of the fixed-point sidecar:
for f=z, the moving roots of z³/3=c differ only by scaling, and all
ordinary punctured curves are isomorphic. The three-fixed-point argument
correctly does not apply there.

## 6. Complete scalar parts and the collision law

At a root component E_rho, Q-Q(rho)=O(A³), because Q'=f² has a
double zero and A is a local uniformizer. Thus, for g_c=T-c,

    g_c=lambda*A+O(A³).

At E_u the actual source expansion from Section2 gives
Q-Q(z0)=a0*A+O(A³), so g_c=(lambda+a0)A+O(A³).
For either nonzero linear coefficient kappa, reparametrization
g=kappa*A+O(A³) gives

    1/(2A²)=kappa²/(2g²)+a regular remainder.

There is no simple pole: equivalently g²/(2A²) has constant
kappa²/2 and no linear term in either local uniformizer. The complete
principal parts are therefore lambda²/(2g_c²) on every root component,
(lambda+a0)²/(2g_c²) on E_u, and zero on R_c.

The regular R_c coordinate must remain before taking the common-diagonal
quotient. At each support c the coefficient vector has at least one
nonzero entry, its R_c entry is zero, and there is only one coefficient
vector, at order two. Thus it is nonzero in the e_c-dimensional
component quotient and generates exactly one principal-part arm under
the actual canonical Weyl action. It does not generate e_c arms when
e_c>1. Multiplication by g_c and differentiation produce all negative
powers in this one direction. Polynomial CRT at the full squared
support ideals separates the unit's pieces at distinct target values.

It follows exactly that total ambient arms=q+1, unit-generated arms=|S|,
and the unit generates all torsion iff the q+1 component values are
distinct. The scalar annihilator is

    ( product_(c in S)(T-c)² ).

Its degree is 2|S|; it should not be abbreviated to a global scalar order
two without specifying each primary support. The upper bound is also an
actual polynomial repair: if P(T)=product_(c in S)(T-c), every reduced
irreducible factor of A divides P(T), so P(T)²/(2A²) is polynomial in
the original source. The exact local nonzero double poles and the
regular R_c components show that no exponent can be lowered, using the
paid constants C(T) to account for every possible rational correction.

The divisibility statement is typed in C[u,t]. It is not generally
true in the independent-coordinate ring C[u,z]. The referee first divides
P(lambda*u*f+Q) by f(z) in C[u,z], then checks that the quotient at
u=0,z=z0 is P(Q(z0))/f(z0)=0. Only after the actual relation
z=z0+u+u³t is imposed does the remaining factor u divide. This is a
useful hostile against losing the E_u component in a localized chart.

After identifying the formal target T, the pointed unit-generated module
depends only on S: each primary piece is one standard principal-part
arm with distinguished element g_c^-2, after rescaling its nonzero
coefficient vector. This preserves all Weyl operations on the unit but
does not recover e_c or the complementary torsion directions.

## 7. Intrinsic rational-pair degree and the actual quintic

Put k=C(T). The element A=(T-Q(z))/lambda is transcendental over k,
and the polynomial map z to Q(z) has degree 2q+1. The irreducibility
of Q(Z)-Y over k(Y), by the prime polynomial ideal and Gauss's lemma,
gives [k(z):k(A)]=2q+1. Since [k(A):k(A²)]=2 and G=1/(2A²),

    [C(x,t):C(T,G)]=2(2q+1).

The degree-two step is not lost because the polynomial defining the
outer map is a square; that simply factors the field degree through
the intermediate field k(A). Every rational mate G' with nonzero
constant Jacobian is alpha*G+R(T), alpha!=0, by ker D_T=C(T).
It generates the same embedded pair field, so the degree is intrinsic
to this first function among all rational constant-Jacobian mates.

For f=z(z-1), h=1 and lambda=1, Q=z⁵/5-z⁴/2+z³/3, a0=6, and

    S={0,1/30,-256/15}.

All three are distinct. The actual source t-degree is five, with leading
coefficient u^15/5; the rational pair degree is ten; ordinary original
fibres have eight punctures. Complete ambient and unit-generated torsion
both have three arms, and the exact scalar annihilator is
[T(T-1/30)(T+256/15)]². The nonzero unit rules out a polynomial mate.
No claim that degree five is minimal for non-isotrivial affine fibres
is inferred from this construction.

## 8. Independent collision and boundary controls

For the explicit quintic f, choose z0 satisfying 6z0²-15z0+10=0,
h=-z0/2 and lambda=1. The polynomial Q is zero modulo this quadratic,
while gcd(f(1+f),6z²-15z+10)=1. Both complex choices therefore retain
all source and boundary submersion conditions and make Q(z0)=Q(0)=0.
Q(1)=1/30 remains different. Ambient torsion stays at three arms while
the unit generates only two. This is an admissible collision, not an
excluded singular parameter.

A separate exact root/root collision uses f=z(z-1)(z-a) with
7a²-7a+2=0 and h=lambda=1. The referee checks Q(1)=Q(0)=0,
Q(a)!=0, Q(-2)!=0, Q(a)!=Q(-2), and f(-2)*(1+f(-2))!=0 in this
exact quotient field. It has four ambient arms and three unit-generated
arms after mutation. Before mutation, it also independently verifies the
supplier's claim that a repeated critical value does not destroy its
two coefficient directions: Q(a) remains nonzero.

The declared polynomial bank q=2,...,6 uses f=product_(j=0)^(q-1)(z-j),
h=lambda=1. Every value is retained, rather than replacing the complete
tuple by its maximum pole. Additional controls pay lambda=0,
lambda=-a0, the q=1 isotrivial boundary, the complete actual Eu jet,
all root remainders, and exact squared-support CRT projectors.

## 9. Frozen verification and stopping scope

The combined verifier is a single independent engine for these two
separately pinned targets. Its always-active normal and optimized runs
are compared as raw-LF stdout and certificate bytes. It imports neither
producer and writes its certificate beside itself before installation,
or into 05-knowledge/results after placement in 04-computation.

Both modes pass **341 always-active exact gates**, with identical raw-LF
stdout and certificate bytes. Frozen SHA256:

    source 1a3046ead7e21b5edf1b6722eb55dfc3fa290aa928bdecb6fb0dc12777bd230e
    output dd8ba9799db9fb3618523e1cc413d867c8544e2c9ffe43885a6f831bbfa0119f
    cert   c3a57cdd3e499074ca3e410b9ee3ae38a720a10fe08a4f3f0f47ab29f9e6e820

    python -B 04-computation/continuing13_20260908_quintic_mutation_audit.py
    python -B -O 04-computation/continuing13_20260908_quintic_mutation_audit.py

The full proofs of all-q scope, entire fibres, non-isotriviality and
arbitrary operator equality are analytic. The finite controls support
those proofs rather than replacing them. The old higher-degree A+H
exploration was not frozen as another theorem or used as a dependency.
All response statements remain in C[x,t], with rational primitives
allowed to have poles. Arbitrary cubic/quintic classification, a global
regular mate, and general Jacobian-conjecture claims remain outside scope.
