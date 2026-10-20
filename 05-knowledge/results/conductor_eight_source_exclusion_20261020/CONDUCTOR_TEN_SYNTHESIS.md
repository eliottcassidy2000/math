# The primitive monic Keller source-conductor bound is at least ten

20 October 2026. **PROVED; exact-dependency synthesis independently audited PASS.**

## Exact theorem

Put R=C[x], K=C(x). Let P,Q belong to R[y], with P monic of positive
y-degree, nonzero constant Jacobian, and

    K(P,Q)=K(y).

Set A=K[P,Q], B=K[y], and delta=dim_K(B/A). If (P,Q) is not a
polynomial automorphism, then

    delta >= 5,
    deg_y c >= 10,                                 (T)

where cB is the normalization conductor of A in B.
The same conclusion holds if P has nonzero constant leading coefficient
instead of being monic, by constant rescaling.

This is a lower bound on the generic source-conductor degree, not on the
function-field degree of the planar map. It does not assert constant
special-fiber delta, does not classify conductor ten, and does not resolve
the planar Jacobian conjecture.

All proof dependencies below have their own independent PASS audits.
This note combines them without changing their hypotheses. Exact source
hashes are recorded in `SYNTHESIS_DEPENDENCY_HASHES.json`.

## 1. Algebraic setup and the zero-defect case

Since P has positive degree, y is integral over K[P]. Hence B is a finite
K[P]-module. The primitive fraction-field hypothesis and the inclusion
A subset B give finite codimension delta and a nonzero conductor.
The inherited trace-dual/norm-square argument proves

    deg_y c = 2 delta.

It also proves that the degree semigroup Gamma(A) is numerical, has
delta gaps, and is symmetric with largest gap 2 delta-1 when delta>0.
These facts are read from the normalization variable y; they do not
identify deg_y c with the degree of C(x,y) over C(P,Q).

If delta=0, then A=B, so y is a finite polynomial expression in P,Q
with coefficients in C(x). Clearing its finitely many denominators yields

    D(x)y = sum r_ij(x) P^i Q^j

with D nonzero and r_ij in C[x]. Choose a in C with D(a)!=0.
Then y belongs to C[P(a,y),Q(a,y)], so the restriction of (P,Q)
to the affine line x=a is injective.
Gwozdziewicz's one-line injectivity theorem makes (P,Q) a polynomial
automorphism. Therefore the noninvertible case has delta>0.

## 2. Every positive genus at most four is excluded

The inherited low-delta classification and the independently audited
genus-four classification give the complete list:

- delta=1: <2,3>;
- delta=2: <2,5>;
- delta=3: <2,7> or <3,4>;
- delta=4: <2,9>, <3,5>, or <4,5,6>.

Each exclusion is supplied by a precise audited theorem:

1. Every <2,m> with odd m>=3 is impossible for an actual monic Keller
   source. This is the all-quadratic-semigroup theorem, including composite
   m. In particular it excludes <2,3>, <2,5>, <2,7>, and <2,9>.
   Its proof retains nonprincipal alternate embeddings such as (6,9):
   an even terminal coordinate determines a polynomial center, and the
   Jacobian makes the centered critical line injective.
2. The constant-leading conductor-six theorem excludes <3,4>.
   Its direct endpoint handles both possible terminal pairs (3,4) and
   (4,6), without requiring a general relative ambient lifting theorem.
3. For <3,5>, the only nonprincipal terminal pair is (3,5).
   Constant-leading reductions take place over C and preserve a nonzero
   constant Jacobian. The direct cubic/quintic addendum excludes that
   entire terminal degree pair.
4. For <4,5,6>, the only nonprincipal terminal pair is (4,6).
   The exact degree-separated module calculation forces its coefficient
   normal form. The actual Jacobian then makes the quartic's linear
   coefficient a unit and gives the contradiction in the classification
   packet, Section 5.

Thus no positive delta<=4 remains. Combined with Section 1, a
noninvertible source has delta>=5, and deg_y c=2delta>=10.
This proves (T).

## 3. Why the obvious shortcuts are not the proof

- All three genus-four semigroups genuinely occur for immersed polynomial
  plane curves. The exact classification is not an abstract obstruction
  to plane realizability.
- The alternate terminal (6,9) pair is genuinely realized:

      f=t^6+4t^2, g=t^9+6t^5+6t,
      g^2=f^3+8f+4t^2.

  Its algebra has semigroup <2,9> and four ordinary nodes. It is excluded
  from a Keller source by the all-m line argument, not by pretending it
  is a quadratic coordinate pair.
- A polynomial family with constant cubic/quintic first integral need
  not be stationary. The explicit family a=u^2,b=u^3/3 has H=0 but
  J=-(35/27)u^6u'. The final constant-Jacobian identity is essential.
- Monicity of generic chosen generators is never assumed. The actual
  Keller pair is reduced by constant target operations, and its own
  coefficients certify each needed polynomial source shear.
- The assertions concern the generic source algebra. No all-fiber
  conductor bound or arbitrary nonmonic relative-order theorem follows.

## 4. Exact dependencies and reproducibility

The theorem uses these audited mathematical results:

1. [Low-delta plane source theorem](../low_delta_source_conductor_exclusion_20261020/LOW_DELTA_QUADRATIC_AND_RELATIVE_GLOBALIZATION.md),
   Sections 1--4: conductor degree 2delta, symmetric degree semigroup,
   and classification through genus three.
2. [Constant-leading conductor-six exclusion](../constant_leading_conductor_six_exclusion_20261020/CONSTANT_LEADING_GLOBALIZATION.md):
   top-coefficient target reduction and exclusion of generic <3,4>.
3. [All quadratic-semigroup exclusion](../quadratic_semigroup_keller_obstruction_20261020/QUADRATIC_SEMIGROUP_KELLER_OBSTRUCTION.md):
   all odd m>=3, including the genuine terminal-(6,9) possibility;
   its source audit pins Gwozdziewicz's one-line theorem.
4. [Genus-four classification and quartic/sextic exclusion](CONDUCTOR_EIGHT_CLASSIFICATION.md),
   Sections 1--5, and [its independent audit](review/AUDIT.md).
5. [Direct cubic/quintic exclusion](CUBIC_QUINTIC_ADDENDUM.md),
   and [its separate independent audit](review/CUBIC_QUINTIC_AUDIT.md).

The source for the one-line theorem is J. Gwozdziewicz, *Injectivity on one
line*, Theorem 1.1, original primary PDF
https://arxiv.org/pdf/alg-geom/9305008v1 .
The characteristic bounds and primitive-entry rule are in
Sathaye--Stenerson, *On Plane Polynomial Curves*, Section 3.1,
https://www.ms.uky.edu/~sohum/ma561/affgeom/papnew.pdf .
The individual audits record exact use and applicability.

The classification checker has 57 exact controls; its independent checker
has 42. The cubic/quintic checker has 31 exact controls; its independent
checker has 32. Both ordinary and optimized Python outputs agree for all
four checkers. The sextic/nonic realization has 11 further exact controls.
The all-m quadratic packet and conductor-six packet retain their own
independent receipts. These computations verify identities and bounded
numerical universes; the universal exclusions are the written proofs.

Reproduce the new packet:

    python3 verify_frontier.py
    python3 -O verify_frontier.py
    python3 verify_cubic_quintic.py
    python3 -O verify_cubic_quintic.py
    python3 review/independent_checks.py
    python3 -O review/independent_checks.py
    python3 review/independent_cubic_quintic.py
    python3 -O review/independent_cubic_quintic.py
    python3 sextic_nonic/verify_realization.py

The bound ten is a synthesis of these exact results, not a new finite
coefficient search and not a claim of historical priority.
