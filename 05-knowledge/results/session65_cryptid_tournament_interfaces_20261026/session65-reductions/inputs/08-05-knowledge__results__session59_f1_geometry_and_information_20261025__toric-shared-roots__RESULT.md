# The shared-root surface is a twisted toric surface: what counts conceal

**Status:** PROVED elementary results, with FINITE-EXACT checks in `check.py`.
**Novelty:** new deductions for this repository, not historical novelty claims.
The cross-polytope toric variety and blueprint constructions are classical.
No integer-square existence theorem, universal F1 no-go, or RH theorem follows.

## 0. Inheritance and six-concept board

Closest inherited mechanism: Session58's exact shared-root base
`hourglass-fiber-product/fiber_product.md`, together with its finite locally
free hourglass maps and its C2-coupled point count. Canonical hostile example:
equal marginal data do not determine the glued count (q=5 versus q=7 have
opposite covariance signs). Corrected near miss: a toric support carrier and
q=1 evaluation do not automatically descend the additive square equations.
Least-used sidecar: the *Galois action on boundary strata and singular points*,
not just the total point polynomial.

| Concept | Object / operation | Lost coordinate / decisive probe |
|---|---|---|
| AP roots | r pairs with common center root | Change variables over Q(i), not Q |
| Toric carrier | cross-polytope semigroup | Preserve Galois swap of opposite generators |
| Point counts | open torus + boundary | Check the strata separately before summing |
| F1 conventions | rank-zero, monoid Hom, blue presentations | They give different numbers |
| Root gluing | H+ tensor over B with H− | Retain shared roots and all additive relations |
| Morphism obstruction | branch divisor inside dense torus | Test one ramified point, not polynomial counts |

Anchor: explicit geometric model of B. Niche: hidden Frobenius cancellation.
Wildcard: higher-dimensional AP-root bases, which reveal why r=2 is special.

## 1. Cross-polytope model in every dimension

Let r>=1 and let k have characteristic different from two. Define

    B_r = Proj k[c,a_1,b_1,...,a_r,b_r]
                 /(a_j^2+b_j^2-2c^2 : 1<=j<=r).

The inherited shared-root surface is B_2, with
`(a_1,b_1,a_2,b_2)=(a,b,d,e)`.

### Theorem 1 (explicit split toric model)

If i in k satisfies i^2=-1, the invertible linear coordinate change

    u_j=(a_j+i b_j)/(1+i),   v_j=(a_j-i b_j)/(1-i)

identifies B_r with

    X_r = Proj k[c,u_1,v_1,...,u_r,v_r]/(u_j v_j-c^2).

This is the normal projective toric variety of the cross-polytope
`conv(±e_1,...,±e_r)`. Its dense split torus is c!=0. Its boundary strata are
indexed by a nonempty subset S of {1,...,r} and one sign for each element of
S; the stratum is G_m^(|S|-1). Thus there are

    2^s binom(r,s)

boundary strata of dimension s-1, and one r-dimensional open torus.

**Proof.** The denominators are units since `(1+i)(1-i)=2`. Their product
gives `u_j v_j=(a_j^2+b_j^2)/2=c^2`; the inverse is obtained by adding and
subtracting `(1+i)u_j` and `(1-i)v_j`.

Map the homogeneous coordinate ring to the Laurent polynomial ring by

    c -> t, u_j -> t x_j, v_j -> t x_j^(-1).

Reducing each occurrence of u_j v_j to c^2 gives the unique normal monomial

    c^(n-sum |m_j|) product_j u_j^(max(m_j,0)) v_j^(max(-m_j,0)),

indexed by `n>=sum |m_j|`. Distinct reduced monomials have distinct Laurent
images, so the displayed binomials generate the entire kernel. The image is
the semigroup algebra of

    S_r={(n,m) in Z x Z^r : n>=sum |m_j|}.

This semigroup is saturated: if a positive multiple of an integral lattice
point satisfies the inequality, that point satisfies it. Its group is
Z^(r+1), since differences of its degree-one generators give `(0,±e_j)`.
The standard saturated-semigroup criterion gives normality; alternatively
the usual lattice-cone proof of that criterion applies directly here.
The degree-one slice is exactly the lattice points of the cross-polytope.

On c!=0 set c=1; then v_j=u_j^(-1), giving G_m^r. At c=0 each equation says
u_j v_j=0. A nonzero projective point therefore selects a nonempty set of
indices and at most one coordinate from each opposite pair. Fixing exactly
which selected coordinates are nonzero gives a projective torus of dimension
|S|-1. These are the claimed toric strata. QED.

The toric model is defined over Z by the binomial presentation. Its comparison
with the original square-root presentation is over Z[1/2,i], **not over Z**.
For the original Q-form, conjugation exchanges u_j and v_j simultaneously.
This explicit lattice involution is part of the descent data.

## 2. Exact finite-field counts and the exceptional two-pair cancellation

Put

    D_r(T)=sum_(s=1)^r 2^s binom(r,s)(T-1)^(s-1)
          =((2T-1)^r-1)/(T-1),
    P_r^+(T)=(T-1)^r+D_r(T),
    P_r^-(T)=(T+1)^r.

The quotient defining D_r is a polynomial; the summation is also valid at
T=1. Superscripts here mean split/nonsplit, not the two hourglasses.

### Theorem 2 (all odd prime powers)

For an odd prime power q and epsilon=chi_q(-1),

    #B_r(F_q)=(q-epsilon)^r + (1+epsilon)D_r(q)/2.

Equivalently it is P_r^+(q) if -1 is square and P_r^-(q) otherwise. More
precisely the open c!=0 has `(q-epsilon)^r` points. The boundary has D_r(q)
points in the split case and none in the nonsplit case.

**Proof.** Normalize c=1 on the open. The number of pairs a,b satisfying
`a^2+b^2=2` is q-epsilon. One elementary proof uses the standard identity
`sum_a chi(a^2-d)=-1` for d!=0: count solutions of `y^2=x^2-d` by
`(x-y)(x+y)=d`, which gives q-1 solutions and hence that character sum.
Consequently

    sum_a (1+chi(2-a^2))=q-epsilon.

The r pairs are independent after this normalization. On c=0, if -1 is
nonsquare, a_j^2+b_j^2=0 forces both coordinates zero for every j, leaving no
projective point. In the split case apply the signed-support strata of
Theorem 1. QED.

An independent cone expression, useful for checking the normalization, is

    #B_r(F_q)=((q-1)(q-epsilon)^r
               +(q+epsilon(q-1))^r-1)/(q-1).

Here the last subtraction removes the cone vertex. It is not an extra
projective point.

### Corollary 3 (exact cancellation boundary)

`P_r^+=P_r^-` if and only if r=1 or r=2. Specifically

    #B_1(F_q)=q+1,             #B_2(F_q)=(q+1)^2

for every odd q. For r>=3 there is no one polynomial P(T) whose value P(q)
equals #B_r(F_q) for all odd prime powers q. There is also no such extension
count polynomial over a fixed nonsplit finite field.

**Proof.** Direct expansion proves r=1,2. The coefficient of T^(r-1) in
P_r^+ is `2^r-r`, whereas in P_r^- it is r. Their difference `2^r-2r` is
strictly positive for r>=3. If a common polynomial existed, evaluation at
q=5^m for all m would force it to equal P_r^+, whereas q=3^(2m+1) would
force it to equal P_r^-. For a fixed nonsplit q the even and odd extensions
give the same contradiction on two infinite sets. QED.

First failure:

    P_3^+=T^3+5T^2-T+1,
    P_3^-=T^3+3T^2+3T+1,
    P_3^+-P_3^-=2T(T-2).

Thus r=2 is an exceptional cancellation, not evidence for a general
polynomial-count principle. In characteristic two the original equations
become `(a_j+b_j)^2=0`; the reduced scheme is P^r and its field-point count
is `1+q+...+q^r`. In particular B_2 has 7 points at q=2, not 9. Even its
odd-field polynomial does **not** count the original integral model at all
primes. Nonreduced structure is deliberately not detected by field points.

## 3. The four nodes recover the lost splitting sidecar

### Theorem 4 (surface, forms, and zeta cancellation)

Over an algebraically closed field of odd characteristic, X_2 has exactly
four singularities, each locally `wz=c^2`, and exactly four boundary lines.
The boundary lines form a cycle. In the original Q-form B_2, conjugation
swaps the singularities in two pairs and the boundary lines in two pairs.
In particular B_2 has no rational singular point and is not isomorphic over
Q to a split toric surface under **any** split toric structure.

For a fixed odd F_q, B_2 and the split X_2 have the same counts over every
finite extension and the same zeta function

    Z(B_2/F_q,t)=1/((1-t)(1-qt)^2(1-q^2 t)).

If -1 is nonsquare in F_q they are nevertheless not isomorphic over F_q:
the singular locus of B_2 has no F_q-points and that of X_2 has four.

**Proof.** In the chart u_1=1, eliminate v_1=c^2. The other equation is
u_2 v_2=c^2, whose gradient vanishes exactly at its origin. It is a quadratic
ordinary double point (an A1 surface singularity). The other three noncentral
coordinate charts give the other three nodes. The open c!=0 is smooth and
these charts cover X_2, so there are no more singularities. At c=0 choose
one zero in each opposite pair; this gives four projective lines, meeting
as the stated cycle at the four coordinate points.

Conjugation swaps u_j,v_j, hence its node permutation is `(u_1 v_1)(u_2 v_2)`.
It reverses both signs of a boundary line's two-sign index, giving two line
orbits. The original coordinates also show this directly: a node has c=0,
one pair identically zero, and the other pair in ratio ±i.

For the split-toric obstruction, the finite geometric singular locus of a
toric surface is stable under its dense torus. A connected torus acts
trivially on a finite set, so each singular point is torus-fixed. In a split
toric variety every zero-dimensional torus orbit is Spec k, hence rational.
This contradicts the rational singular locus of B_2. This argument concerns
a global split torus action, **not arbitrary geometric torifications**.

The count and zeta assertions follow from Corollary 3 and
`Z=exp(sum_n #B_2(F_(q^n))t^n/n)`. To expose the cancellation, let epsilon be
chi_q(-1), so the extension value is epsilon^n. The open and boundary counts
are respectively

    U_n=q^(2n)-2(epsilon q)^n+1,
    D_n=2q^n+2(epsilon q)^n.

Therefore

    Z(U,t)=(1-epsilon q t)^2/((1-q^2t)(1-t)),
    Z(D,t)=1/((1-qt)^2(1-epsilon q t)^2).

The epsilon-dependent factors cancel exactly in their product. If epsilon=-1,
the singular subscheme has zeta `(1-t^2)^(-2)`, versus `(1-t)^(-4)` for X_2.
This separately marked finite scheme recovers information that the total
zeta function has lost. QED.

The same total zeta is also that of P^1 x P^1, but X_2 is singular and
P^1 x P^1 is smooth. Thus even splitness and singularity data are not
recoverable from the displayed total point polynomial.

## 4. A precise obstruction to extending the toric bridge to hourglass maps

Use the inherited equations over an algebraically closed field of
characteristic zero:

    B=B_2: a^2+b^2=2c^2=d^2+e^2,
    H+ -> B: h^2=F=a^2+d^2-c^2,
              h'^2=G=3c^2-a^2-d^2.

The second hourglass uses `J=c^2+a^2-d^2`, `K=c^2-a^2+d^2` instead. These are
the actual finite locally free rank-four maps, not a newly chosen cover.

### Theorem 5 (interior-branch obstruction)

With the fixed dense torus `T={c!=0}` on B supplied by Theorem 1, neither
H+ -> B nor H− -> B can be a finite toric morphism from a normal toric
source extending that target torus structure. The same applies to M -> B.

**Proof.** At

    (c,a,b,d,e)=(1,0,sqrt(2),1,1)

one has F=0, G=2, J=0, K=2. This is inside T: after the coordinate change
every u_j,v_j is nonzero because their product is c^2=1. On the smooth open
T, F has a simple zero at this point. For example take tangent variation
`delta a=delta b=0, delta d=1, delta e=-1`; it respects both AP equations
and `delta F=2`. The square-root equation h^2=F is therefore ramified along
the divisor F=0 there. The relative differential at h=0 is already nonzero;
the other new square root is a unit, so it does not remove that ramification.
The same calculation for J has `delta J=-2`.

A finite dominant toric morphism between normal toric varieties has finite
lattice index; the preimage of the target dense torus is the source dense
torus. In characteristic zero its restriction is an isogeny of tori, hence
etale: in Laurent coordinates the Jacobian is the nonzero lattice-index
determinant times units. It cannot ramify over an interior point of T. This
contradicts the exhibited branch point for either hourglass and for their
combined cover. QED.

**Exact failure boundary.** The conclusion is relative to this specified
toric base structure and finite toric morphisms. It does not rule out every
torification of H±, blue models, alternative additive geometries, or all
forms of F1 descent. The surviving positive statement is the genuine toric
base after adjoining i, together with the original exact finite flat
root-gluing diagram. Branch loci are a necessary morphism sidecar.

## 5. A weak but exact blue-presentation repair

Here a blueprint presentation consists of a pointed multiplicative monoid,
with declared equalities between formal finite sums. Ringification replaces
those declarations by ordinary polynomial relations over Z. This is the
framework of Lorscheid's *The geometry of blueprints, Part I*, especially
Sections 1.6 and 1.14, Proposition 1.12:
https://arxiv.org/pdf/1103.1745v2 . We use the presentation universal property,
not an unqualified identification between every affine blue scheme and its
coordinate blueprint.

Let C_B be the free pointed monoid on c,a,b,d,e with preaddition generated by

    a^2+b^2 ≡ c^2+c^2,
    d^2+e^2 ≡ c^2+c^2.

Form C_+ by adding h,h' and

    h^2+c^2 ≡ a^2+d^2,
    h'^2+a^2+d^2 ≡ c^2+c^2+c^2.

Form C_- by adding s,t and

    s^2+d^2 ≡ c^2+a^2,
    t^2+a^2 ≡ c^2+d^2.

Let C_M contain all nine generators and all six relations.

### Proposition 6 (blue algebra gluing and literal F1 assignments)

There is a pushout in the category of presented blueprints

    C_M = C_+ tensor_(C_B) C_-.

Ringification gives exactly the inherited affine cone diagram, including
`R_M=R_+ tensor_(R_B) R_-`. Each of the four blueprints C_B,C_+,C_-,C_M has
exactly two morphisms to the initial pointed-monoid blueprint F1={0,1}: the
all-zero variable assignment and the all-one variable assignment. Its
nonzero homogeneous-coordinate assignment count is therefore one.

**Proof.** A morphism from C_M to any blueprint is an assignment of all nine
generators satisfying all six relations. This is precisely a pair of
morphisms from C_+ and C_- agreeing on the five shared generators. That is
the pushout universal property. Ringification preserves that property, or
one can read off the ordinary quotient polynomial rings explicitly.

In F1 with no additional preaddition, equality of sums of ones means equality
of their nonnegative integer multiplicities; in particular `1+1` is NOT 1.
Write A,B,C,D,E in {0,1} for generator images; squares have the same images.
The first two relations give `A+B=2C=D+E`. If C=0 all five images are zero;
if C=1 all five are one. The additional equations then force all new images
to be the same value. Both resulting assignments exist. QED.

This is a genuine, explicitly specified **blue algebra presentation**, but
its existence is a weak statement: integer polynomial equations can always
be recorded by putting positive and negative terms on opposite sides.
We do not claim flat descent, preferred F1 rational-point conventions,
Frobenius lifts, or an RH input. Projective blue gluing can be developed with
the corresponding graded theory; it is not needed for this proposition.
The usual projective scheme gluing remains exactly the already audited
Session58 shared-chart proof.

For comparison, the split monoidal model of X_r has one monoid point per
toric orbit, hence `1+sum_s 2^s binom(r,s)=3^r` literal Hom-points from the
specified pointed F1. Its rank-zero points number 2r, which is P_r^+(1).
At r=2 these are 9 and 4; the nonzero assignments of the original additive
blue presentation number 1. These are **different models and conventions**,
not inconsistent answers to one specified functor. The blue presentation
base-changes to the original integral B_r; the monoidal model base-changes
to X_r; their comparison here needs 2 inverted and i adjoined.

## 6. Connection contracts and strongest survivors

1. **AP roots -> toric X_r.** Map: explicit linear change over k(i).
   Preserves the entire scheme after that base extension. Loses the k-form
   if conjugation is discarded. Sidecar: simultaneous opposite-coordinate
   swap. Decisive test: number of rational singular points when r=2.
2. **Toric strata -> finite-field polynomial.** Map: sum of stratum counts.
   Preserves total counts. Destroys splitting and cancellation location.
   Sidecar: open/boundary counts or their separate zetas. Decisive test:
   r=3 or the marked four-node subscheme.
3. **B toric model -> H± cover.** Attempted map: retain the actual
   coordinate-forgetting map as a toric morphism. Fails at interior branch.
   Strongest survivor: a finite flat cover of a toric base, not a finite
   toric map. Sidecar: the four additive radicand divisors.
4. **Root relations -> blue presentations.** Map: formal additive
   relations and pushout. Preserves labels, equations, and exact cone
   gluing under ringification. Does not preserve a guessed q=1 point
   interpretation. Sidecar: category, base extension, and point functor.

Board revisit: the toric coordinate change solves the carrier question but
creates a Galois sidecar; the count theorem shows that sidecar can disappear
from all total extension counts; the branch obstruction shows even a
correct carrier does not make its maps toric; the blue pushout repairs
equation bookkeeping without manufacturing the missing arithmetic or
cohomological inputs.

## 7. Verification and source scope

Run `python3 check.py` in this directory. The script checks the formulas by
independent pair enumeration over listed prime fields and quadratic
extensions, checks direct projective enumeration for small fields, checks
the toric coordinate transformation and branch tangent symbolically, and
checks every {0,1} assignment of the presented blue diagrams. It saves exact
data in `results.json`. Finite checks support rather than replace the proofs.

Primary framework references (definitions/classical constructions, not
novelty sources for this packet's model-specific deductions):

- Anton Deitmar, *F1-schemes and toric varieties*, Section 4 and Theorem 4.1,
  https://arxiv.org/pdf/math/0608179v10 . We do not infer the original
  Q-form has a split monoidal model from its geometric toricity.
- Javier Lopez Pena and Oliver Lorscheid, *Torified varieties and their
  geometries over F1*, https://arxiv.org/abs/0903.2173 . The distinction
  between a torus decomposition and compatible morphisms is retained.
- Oliver Lorscheid, *The geometry of blueprints. Part I*, Proposition 1.12,
  https://arxiv.org/pdf/1103.1745v2 . Only the stated algebraic presentation
  and pushout machinery is imported.

No broad geometric claim is inferred from a bare count. The original
integer-square hourglass/full-square problems and RH remain untouched.
