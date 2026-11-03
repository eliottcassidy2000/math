# Nonlinear orbit scaling already occurs under composition

Status: **PROVED**, with exact symbolic controls; same-agent audit.
This corrects the unrestricted linear-scaling language of
[THM-4727, mixed-weight torus reduction](../../01-canon/theorems/THM-4727-mixed-weight-torus-reduction-of-keller-maps-ramification-exponent-twisted-planar-jacobian-equation-and-the-dimension-profile.md)
and settles the existence clause of HYP-9343(c)(2). It does not classify
all Keller counterexamples or settle JC(2).

**Priority/integration:** concurrent
[Session144](../session144_jacobian_eth_mrb_20261103/jacobian/equivariant_composition_tower.md)
published this same composition tower and degree law before our checkpoint,
with stronger disjoint-divisor results. The calculations here are independent
confirmation. Its subsequent gauge correction, incorporated below, settles
the raw degree2..6 question too. No new primitive species is claimed.

## 1. The actual map and the retained quotient

Reorder the published counterexample's outputs as K=(R,Q,P), where

    R=2x-3x^2y-x^3z,
    Q=y+3x(1+xy)^2z+3xy^2(4+3xy),
    P=(1+xy)^3z+y^2(1+xy)(4+3xy).

Then det DK=2. It is equivariant for the SAME source and target weights
(1,-1,-2), so iteration preserves these coordinates. Put v=xy,t=x^2z and

    l=2-3v-t,
    A=v+3(1+v)^2t+3v^2(4+3v),
    B=(1+v)^3t+v^2(1+v)(4+3v).

On x!=0 the map is (xl,A/x,B/x^2). These expressions are polynomial
everywhere: A belongs to (v,t), B to (v^2,t). The quotient is

    Fbar(v,t)=(l A,l^2 B),   det D Fbar=2l^2.            (1)

The three rational points

    q=(0,0,-1/4),  (1,-3/2,13/2),  (-1,3/2,13/2)

all map to q, which is fixed. Thus EVERY K^k is noninjective, with these
same visible collisions, and has constant Jacobian2^k. Padding with the
identity preserves this witness in every dimension>=3.

The starting map and its provenance are independently displayed in
[Meng--Yang, equations(1),(2)](https://arxiv.org/html/2607.22198v2).
Our checker recomputes the determinant and collision over Q, rather than
assuming the conjecture's status from a news report.

## 2. Composition is a multiplicative cocycle, not constant scaling

For two polynomial maps of the same weights, with data (l_F,A_F,B_F),
(l_G,A_G,B_G), direct substitution on x!=0 yields

    l_(G o F) = l_F * l_G(Fbar),
    A_(G o F) = A_G(Fbar)/l_F,
    B_(G o F) = B_G(Fbar)/l_F^2,
    (G o F)bar = Gbar o Fbar.                           (2)

The apparent divisions cancel in the invariant polynomials: substitution
of (l_F A_F,l_F^2 B_F) into A_G in(v,t) or B_G in(v^2,t) supplies the
required powers. The same weighted ideals ensure polynomiality at x=0.
The identity extends across that plane because the original maps are
polynomial. No assertion that a rational expression is automatically
integral is needed.

For K squared specifically,

    l_2=l(2-3lA-l^2B).                                (3)

Its total degree in v,t is7, with30 monomials. The corresponding A_2,B_2
are polynomials of total degrees13,16, satisfy the weighted ideals, and
obey the exact twisted equation

    l_2 {A_2,B_2}+A_2 {l_2,B_2}+2B_2 {A_2,l_2}=4.       (4)

Thus the older question whether ANY solutions have deg l>=2 has an
explicit positive answer. No coordinate change or finite numerical search
is involved. Requiring an inequivalent composition-indecomposable map would
be a different question; (3) does not answer it. Raw nonlinear degree alone
is already settled by the coordinate shear described in Section5 below.

## 3. All-order degree formula: keep three different clocks

Let (v_j,t_j)=Fbar^j(v,t). Equation(2) gives

    l_k=product_(j=0)^(k-1) (2-3v_j-t_j).                (5)

Writing s=v(3v+t), the top homogeneous parts of v_1,t_1 are -3s^2,s^3,
of degrees4,6. If deg t_j>deg v_j, the unique highest weighted terms in
the next substitution are -3v_j^2t_j^2 and v_j^3t_j^3. Their coefficients
are nonzero. Induction proves

    deg v_j=4*5^(j-1),  deg t_j=6*5^(j-1)  (j>=1).

Consequently

    delta_k:=deg l_k=(3*5^(k-1)-1)/2,
    delta_1=1,  delta_(k+1)=5delta_k+2.                 (6)

The initial degrees are1,7,37,187,937. This is ordinary polynomial degree
on the quotient, not the geometric degree of K^k. The latter is3^k by
the field-extension tower law and the inherited generic degree3 of K;
the determinant is2^k. The triple (2^k,3^k,delta_k) keeps three different
invariants that a single degree label would conflate.

The symbolic proof of the top-degree law is supplemented by exact
specialization v=t=u modulo1009 through k=5. This yields matching lower
bounds on degree; it is not being used to infer a rational identity from
one specialization. Equation(4) is checked over Q in full.

## 4. What a complete classification would still require

[THM-1330, Keller monoid necessary-structure atlas](../../01-canon/theorems/THM-1330-keller-monoid-exact-picture-inverse-jelonek-cusp-rule.md)
already distinguishes necessary covering data from a classification. Its
composition law is the mechanism missing from the later linear-scaling
claim. The inherited degree spectrum in dimensions>=3 is{1,3,4,5,...};
degree1 consists of automorphisms and degree2 is excluded by the cited
Galois theorem. This is a value spectrum, not a list of equivalence classes.

A full classification must still resolve composition atoms, equivalence
within a degree, realizability of non-properness/monodromy data, and possible
maps outside the displayed equivariant families. Even an atlas covering
all weight signatures does not enumerate all polynomial maps in a signature.
[Gao's primary construction](https://arxiv.org/html/2608.00222v1) produces
families and explicitly defers some full fibre stratifications.
[Shaska's discussion](https://arxiv.org/html/2607.20210v2#S11) retains
specific orbit and weight-type questions. These sources support neither
an unrestricted complete classification nor a proof of JC(2).

A second direction correction is necessary: Meng--Yang establish
HC(4)=>JC(2). They do not establish the converse. A JC(2) counterexample
would refute HC(4) by doubling, but an HC(4) counterexample does not follow
backward to JC(2) from this implication. The affected legacy descriptions
are repaired together with the composition error.

## 5. Incoming stronger classification boundaries

The independently audited
[equivariant gauge hostile](../session144_jacobian_eth_mrb_20261103/jacobian_followup/equivariant_gauge_hostile.md)
is the elementary automorphism T_k=(x,y+x^(2k-1)z^k,z), k>=2. Its quotient
is(v+t^k,t), so K composed with T_k has scale2-3v-t-3t^k of degree k,
with unchanged determinant2 and generic degree3. It is equivalent to K.
Thus even primitive maps exhibit every nonlinear displayed degree; asking
for the smallest such degree was the wrong next question. A meaningful
classification must fix an equivalence relation and use invariants such
as divisor components and boundary fibre/inertia data.

Session144 proves the factors l composed with Fbar^j have pairwise disjoint
smooth zero divisors, and the kth iterate has at least k components in its
scale-zero divisor. Unlike polynomial degree, this component count survives
equivariant source/target automorphisms. It distinguishes the composite
tower from any affine-scale map in that category, but not automatically
under arbitrary coordinate changes that alter the chosen torus action.

The concurrent
[weighted-family classification, W1--W3](../jacobian_eth_prime_wait_20261102/WEIGHTED-CLASSIFICATION.md)
goes further in a different direction. For its specified admissible R,
off the target coordinate plane C=0, put T(w)=R(w)-BCw+cAC^2 and D=gcd(T,T').
The exact finite-fibre count is n-deg D-deg rad D: only simple roots
reconstruct finite points. On C=0 there are separate two/three-sheet rules;
the off-plane formula must not be substituted there. Repeated nonzero
roots supply separate inertia cycles; a root of multiplicity mu at zero
instead supplies a cycle of length mu-1. This marked-zero offset is the
essential missing coordinate in a root-multiplicity-only classification.

That proof constructs at least p(n-3)+p(n-2)-1 pairwise inequivalent
determinant-one primitive maps of generic degree n>=4, where p is the
partition function. The first bounds are2,4,7,11,17,25. They all have
global monodromy S_n and are separated by the unordered boundary-component
ledger (finite sheets, inertia type). These are inherited proved family
results, not a count of all Keller maps and not new results of this packet.
