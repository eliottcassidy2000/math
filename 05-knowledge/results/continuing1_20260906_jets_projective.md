# Complete projective higher jets: attained precision and an intrinsic p31 packet

**Status: PROVED relative to the named inherited inputs; FINITE-EXACT;
INDEPENDENTLY AUDITED.** See the [independent root audit](continuing1_20260906_jets_audit.md).
This note extends the declared domain of an existing inverse-precision theorem
and gives a bracket formulation of an existing intermediate Smith packet.
It does not re-claim the affine theorem, the new two-jet projective law, or the
p31 factor formulas as new discoveries. No full p31 partition or general
higher-jet metric law is asserted.

## 1. Inheritance and the connection being made

The immediate incoming supplier is
`05-knowledge/results/creative_20260906_smith_bridge.md`, proved and independently
audited. It gives primitive projective two-jet precision, complete-jet covariance,
and a residue-class splitting argument. Its explicit open question asks for
an intrinsic bracket sidecar transporting higher intermediate ideals.

The affine arbitrary-multiplicity precision theorem is already **THM-4443**,
`01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md`,
proved in `05-knowledge/results/hermite-higher-jet-unit-boundary-overnight-hexagon-sep05.md`.
The normalized shell-jet mechanism was also recovered from
`05-knowledge/results/overnight6_20260906_jets_sidecar.md`. The present largest-
factor formula retains those reciprocal jets and moves them to a homogeneous
integral coefficient lattice, where no common integral affine chart is assumed.

The proved p31 suppliers are
`05-knowledge/results/overnight13_20260906_jets_p31_intermediate.md` and
`05-knowledge/results/overnight14_20260906_jets_p31_adjacent_even.md`, with their
independent audits. They settle three full determinantal ideals and two factors
of the standalone48-dimensional three-node16-jet observer. Their indices are
kept separate from the larger projective observer constructed below.

The six-concept board is: primitive homogeneous coefficient lattices; complete
Hasse banks; Gauss content; reciprocal jets; projective residue splitting; and
the p31 relative-position residue. The canonical hostile is THM-4443's isometric
three-jet pair with unequal worst precision. The corrected near miss is assuming
that a chart or representative change is integral without retaining its units.
The least-used sidecar is the content of a homogeneous cardinal factor, followed
by a transverse direction in a bracket cross ratio.

The source-to-target map preserves the entire Smith module only under integral
unit changes. Its reciprocal coefficient vector retains higher-jet cancellation.
Taking its largest denominator forgets intermediate factors; retaining the p31
bracket residue recovers the specified intermediate ideal decision. Individual
selected minors are not asserted to transform as scalar invariants.

## 2. Exact projective precision for arbitrary complete multiplicities

Let v_i=(a_i,b_i) be primitive integer vectors with pairwise distinct rational
directions, and let m_i≥1. Set M=Σm_i and D=M−1. Choose integer w_i with
det(v_i,w_i)=1. On the rank-M lattice of homogeneous binary forms of degreeD,
observe all coefficients

    [T^r] F(v_i+T w_i),                 0≤r<m_i.                 (1)

These are complete local Hasse jets. Write

    Δ_ij=det(v_i,v_j),
    Q_i(Z)=∏_(j≠i) det(Z,v_j)^(m_j),
    q_i(T)=Q_i(v_i+T w_i),
    h_(i,r)=[T^r] q_i(T)^(-1),          0≤r<m_i.                 (2)

**Theorem A.** The observer is nonsingular. Its largest integer Smith factor is

    s_max=lcm_(i,r<m_i) denominator(h_(i,r)),                      (3)

where every rational denominator is reduced and positive. At any primep,
the exact extra coefficient-recovery precision is

    L_p=max_(i,r<m_i) {−v_p(h_(i,r))},          v_p(0)=∞.           (4)

The maximum is nonnegative because h_(i,0) is the reciprocal of a nonzero
integer. For N≥1, all observations modulo p^(N+L_p) determine all coefficients
modulo p^N; when L_p>0 one fewer digit fails for actual integral coefficient
vectors. In fact **a perturbation in a single value observation suffices for
sharpness, with every derivative observation unchanged**.

### Homogeneous cardinal proof and attainment

Put s_i(Z)=det(Z,w_i), t_i(Z)=det(v_i,Z). These are integral unimodular
coordinates, with s_i(v_i+Tw_i)=1 and t_i(v_i+Tw_i)=T. The cardinal form for
jetk at nodei is

    Φ_(i,k)(Z)=Q_i(Z) Σ_(r=0)^(m_i−1−k)
                      h_(i,r) s_i(Z)^(m_i−1−k−r) t_i(Z)^(k+r). (5)

It is homogeneous of degreeM−1. At nodei, multiplying the truncated reciprocal
series gives T^k modulo T^(m_i). At every other nodej, the factor
det(Z,v_j)^(m_j) kills its entire observed jet. Thus (5) gives every column
of the inverse observer matrix and proves nonsingularity.

Every linear form det(Z,v_j) is primitive. Products of primitive integer
polynomials are primitive: reducing modulo any prime gives a product of nonzero
polynomials in an integral domain. Consequently Q_i is primitive. Multiplication
by Q_i preserves the least rational coefficient denominator of any homogeneous
factor; primewise, its minimum coefficient valuation is zero and the minimum
valuation of a product is the sum of the two minima. The unimodular change
from(X,Y) to(s_i,t_i) also preserves that denominator, since substitution and
its inverse have integer matrices on each homogeneous degree.

The exact denominator of Φ_(i,k) is therefore the lcm of the denominators
of h_(i,0),…,h_(i,m_i−1−k). The value column k=0 attains the entire local
list. This supplies projective denominator attainment without a monic common
affine coordinate. Taking the lcm over all columns proves (3): the least
denominator of the inverse of a nonsingular integer matrix is its largest
Smith factor, by integral unimodular Smith changes.

For the stronger sharpness statement, choosei whose value-column denominator
d_i has v_p(d_i)=L_p>0. The integer form d_i Φ_(i,0) has at least one coefficient
not divisible byp; otherwise its p denominator would not be minimal. The form
p^(N−1)d_i Φ_(i,0) is integral, has a coefficient of valuationN−1, and its only
nonzero observation is the value at i, of valuationN+L_p−1. It is invisible at
the proposed one-less data precision but not modulo p^N in the source.

For one direction, Q_i=1 and all denominators equal one. Thus arbitrary
singleton multiplicity has no loss, including a direction at infinity.

### An integral recurrence and a finite local precision bound

Write q_i(T)=q₀+q₁T+⋯ and define

    N₀=1,     N_r=−Σ_(j=1)^r q_j q₀^(j−1) N_(r−j).              (6)

Then h_(i,r)=N_r/q₀^(r+1). The exact local value-column denominator is

    d_i=|q₀|^(m_i) /
          gcd_{0≤r<m_i}(q₀^(m_i−1−r) N_r).                    (7)

The gcd is nonnegative and nonzero because its r=0 operand is q₀^(m_i−1).
It divides q₀^(m_i), so (7) is an integer. This is an exact coefficient-content
formula, not a replacement by valuations of separate summands.

If h_ij=v_p(Δ_ij), A_i=Σ_(j≠i)m_j h_ij, and f_i=max_(j≠i)h_ij, then

    L_p≤max_i[A_i+(m_i−1)f_i].                                (8)

For a singleton use f_i=0. To prove (8), expand each factor
(Δ_ij+T det(w_i,v_j))^(−m_j); its normalized linear coefficient has valuation
at least−h_ij. A degree-r term has valuation at least−A_i−rf_i, and sums
cannot lower that bound. The verifier uses (8), not a guessed largest factor,
to choose sufficient modulus precision for complete local Smith peeling.

## 3. Covariance, weighted determinant, and full residue splitting

The observer retains complete higher banks under the lawful operations already
identified in the incoming projective note. Replacing w_i by w_i+k_i v_i gives

    F(v_i+T(w_i+k_i v_i))
       =(1+k_iT)^D F(v_i+[T/(1+k_iT)]w_i).                     (9)

Through orderm_i−1 this is triangular with diagonal one. For g∈GL₂(Z_p), use
v_i'=g v_i and w_i'=(det g)^(-1)g w_i. Form substitution is unimodular, and
the order-r target row gets a unit factor(det g)^(-r). Replacing a local
primitive representative by u_i v_i and its tangent by u_i^(-1)w_i multiplies
that row by u_i^(D−2r), again a unit. OverZ use g∈GL₂(Z) and representative
units±1. These preserve the full Smith module, not just the largest factor.

The full weighted homogeneous determinant is

    |det E|=∏_(i<j)|Δ_ij|^(m_i m_j).                           (10)

For proof on a rational chart b_i≠0, write x_i=a_i/b_i and P(X)=F(X,1).
The jet transformation to affine Hasse derivatives is triangular with diagonal
entries (−1)^r b_i^(D−2r). Its block determinant has absolute factor
b_i^(m_i(M−m_i)). The ordinary confluent determinant has factors
(x_i−x_j)^(m_i m_j); substituting Δ_ij/(b_i b_j) cancels every power of b_i.
A determinant-one integer shear can make all second coordinates nonzero;
only finitely many shear parameters are forbidden. This proves the identity
over the rationals, including infinity, without claiming that this chart is
integral at a chosen prime.

Partition the directions by reduction in P¹(F_p). For a class C, put M_C=Σ_(i∈C)m_i
and multiply its degreeM_C−1 source by

    R_C(Z)=∏_(j∉C)det(Z,v_j)^(m_j).                            (11)

This gives a source map B from the direct sum of the class lattices to the
degreeM−1 full lattice. Its observations vanish outsideC. WithinC, multiplication
by R_C acts on the full jet at nodei by a triangular matrix with unit diagonal
R_C(v_i), since every cross-class bracket is a unit. Hence EB is a direct sum
of the class observers up to target units. By (10), the full determinant
valuation is the sum of the class determinant valuations. Taking determinants
therefore gives v_p(det B)=0. Thus B is unimodular and the **entire p-Smith
module splits by projective residue class for arbitrary multiplicities**.

This fixed-degree gluing argument remains valid when all p+1 projective residue
classes occur. A single integral affine chart then cannot contain all directions:
its excluded residue class is occupied, and GL₂(F_p) only permutes the classes.
In particular a class containing one direction of multiplicitym contributesm
unit Smith factors, not a nontrivial precision block.

## 4. The p31 intermediate packet in intrinsic bracket coordinates

Take three local primitive directions v₀,v₁,v₂, each with multiplicity16, and
assume all three bracket valuations at31 equal e≥1. Choose any primitive
transverse z with Δ₀z a31-adic unit. The other two Δ_iz are then units too.
Set

    U=31^(−e) Δ₀₂ Δ₁z,     V=31^(−e) Δ₀₁ Δ₂z,
    ā=U/V mod31.                                             (12)

Both U,V are units. Plücker's identity gives
U−V=31^(−e)Δ₁₂Δ₀z, also a unit, so ā avoids0 and1.

The residue is independent of the transverse choice. Indeed write another
transverse vector z'=αv₀+βz, with α∈Z₃₁ and β a unit. Both U and V change
to β times themselves minus the same quantity
α31^(−e)Δ₀₁Δ₀₂, which is divisible by31^e. Consequently their quotient is
unchanged modulo31^e, in particular modulo31. This does **not** assert equality
of the exact rational cross ratios. Unit changes of representatives and GL₂(Z₃₁)
cancel directly in the bracket ratio. Permuting the triple gives the six usual
transforms a,1/a,1−a,1/(1−a),a/(a−1),(a−1)/a on its residue.

Define κ=1 when ā belongs to{3,11,15,17,21,29}, the cross-ratio orbit of3,
and κ=0 otherwise. Thus κ is an intrinsic, unordered relative-position sidecar
for this equilateral projective triple. Pairwise bracket valuations do not
determine it.

For a literal binary packet, recover the inherited primitive polynomial

    q₁₄(a)=1/19380 · Σ_(j=0)^14 binom(14+j,j)binom(28−j,14−j)
                                     a^(14−j)(a−1)^j,
    Q₁₄(U,V)=V^14 q₁₄(U/V).                                 (13)

The proved all-lift packet theorem says
min(v₃₁(Q₁₄(U,V)),v₃₁(Q₁₄(V,U)))=κ. Its divided companion caps the common
cancellation; a single selected packet can cancel deeper and is not substituted
for this pair. Formula(13) writes the intermediate decision directly in brackets.

To transport the full ideals lawfully, use coordinates
T(Z)=det(v₀,Z), S(Z)=det(Z,z). Their coordinate matrix has unit determinant.
After unit representative normalization, the three affine nodes are
0, x₁, x₂ with x_i=Δ₀i/Δ_iz and v₃₁(x₁)=e. Write x₁=31^e u with u a
unit, and scale the affine coordinate by u⁻¹. This is still a unit change;
it gives exactly0,31^e,31^e a with a=x₂/x₁. No division of the observer
coordinate by31^e is performed. Complete-jet covariance now identifies its
whole local Smith module with the inherited affine observer.

Therefore, on the **standalone48-dimensional cluster module**, letting E_r
denote the sum of the first r exponents after its16 zero factors,

    E₂₈=588e+2,     E₂₉=631e+1+κ,     E₃₀=675e+1,              (14)
    λ₂₉=43e−1+κ,               λ₃₀=44e−κ.                    (15)

These are the already proved affine formulas, now transported by an intrinsic
projective packet. They concern full determinantal ideals and ordered factors;
they do not assert covariance of the individual minimum minors used in their
original proof or supply the other p31 factors.

If other projective residue classes contain one direction apiece, with total
multiplicityR, they add exactlyR unit factors by Section3. Thus (14) concerns
the full ideals D_(R+44),D_(R+45),D_(R+46), and the factor pair (15) occupies
positions R+45 and R+46 in the full observer. Arbitrary nontrivial additional
clusters can interleave positive factors; this index statement is not made
for them, although the direct-sum module statement still holds.

## 5. Exact controls, failures, and reproduction

The principal new no-common-chart control has directions

    (1,0), (1,31^e), (1,31^e a), and (r,1), 0≤r<31,

with multiplicities16,16,16 on the first three and one on each remaining
direction. All32 projective residue classes occur. The source has dimension79;
R=31, so there are47 initial zero factors. At e=1, a=3 versus4 gives

| a | v₃₁(D₇₅) | v₃₁(D₇₆) | v₃₁(D₇₇) | factors76,77 |
|---|---:|---:|---:|---|
| 3 | 590 | 633 | 676 | 43,43 |
| 4 | 590 | 632 | 676 | 42,44 |

Literal full79×79 matrices have kernel orders31^762 and31^761 modulo31^43.
This is a FINITE-EXACT full-observer consequence; the all-depth theorem is
(14)–(15). Two additional full matrices use a=879 and−20 at e=2. The script
also checks transverse changes and homogeneous packet valuations, including
the inherited deep individual-cancellation hostile879.

For heterogeneous complete-residue examples, take multiplicities3 and4 on
(1,0),(1,p^e), and multiplicity one on each(r,1),0≤r<p, at p=2,3,5 and
e=1,2. Every full partition agrees with the two-direction partition plusp
unit factors. The general splitting follows from the proof, not these cases.

The higher-jet metric firewall is retained: the affine three-jet triples
(0,4,8) and(0,4,12), transported to directions(1,−x), have dyadic losses18
and19. The inherited labeled metrics agree. Thus Theorem A does not extend
the incoming metric-only two-jet terminal law to higher multiplicity.

Primitivity and completeness are independently tested. For degree-four forms,
the primitive directions(1,0),(0,1), multiplicities3,2, give all unit factors.
Replacing only(1,0) by(2,0) while retaining tangent(0,1) gives dyadic factors
0,0,2,3,4. This violates the declared primitive/unimodular-frame hypothesis.
Conversely, a local unit coordinate change with its tangent correctly adjusted
by a unit denominator preserves the full p-Smith module. Finally, keeping four
observations on cubics but replacing one first-order datum by a second-order
datum creates a singular matrix in the two-axis example; equal cardinality
does not replace the complete-bank condition.

[Standalone source](../../04-computation/continuing1_20260906_jets_projective.py)
and [output](continuing1_20260906_jets_projective.out):

```
python -B 04-computation/continuing1_20260906_jets_projective.py
python -B -O 04-computation/continuing1_20260906_jets_projective.py
```

Both runs pass **1,645 always-active gates** with identical LF output. The
source uses exact standard-library integer/rational operations plus SymPy for
independent integer Smith forms and determinants. It imports no prior producer.
The finite universe is14 small global configurations,42 prime comparisons,
42 complete integer Smith symmetry checks, nine integral value-only sharpness
hostiles, six heterogeneous complete-residue configurations, and four full
79-dimensional p31 matrices. Literal homogeneous jet rows are constructed
independently of the cardinal formulas; modular unit-pivot peeling is checked
against full integer Smith forms on the small configurations.

Source SHA256:
`aba09c132b7d7f5b96bdf239869362ed173986444ad1e68389dad823bf2ec6e3`.
Output and optimized-output SHA256:
`53532c90b6204058f8ff5adea492af94c88933e822ff90ad7b70a835ae82acd2`.
No Git, navigation, scarce identifier, or inherited artifact was modified.
