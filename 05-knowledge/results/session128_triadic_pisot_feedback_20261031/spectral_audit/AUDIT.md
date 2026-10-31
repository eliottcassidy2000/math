# Session128 independent audit: equivalence, Pisot graph, and feedback maps

**Verdict: SOUND WITH THE CORRECTIONS BELOW APPLIED.**
The author packets are `../triadic/` and `../pisot/`; the bounded final
extension is `../resolvent_norm/`. This audit writes no repository files and
makes no claim of a new pure-spectrum, tiling, or Amanda-loop theorem.

## 1. The triad is a genuine fixed-degree equivalence

For a fixed irreducible `f_k`, the categories have all their morphisms:
finite free K-modules, rational block-companion intertwiners, and rational
shift-commuting maps of bilateral k-nacci streams. Multiplication-matrix
realization and the full stream `s_n=e_0^T C^n v` are mutually invertible
on the named skeletal objects and Hom spaces. The inverse reads each
commuting block's first column and reconstructs coefficients from k samples.
The explicit observability inverse and unimodularity are correct.

Irreducibility at k=3,4,5 is independently certified here over F_3,F_2,F_5,
respectively, by exhausting all possible monic factors of degree at most
floor(k/2). The all-k Rouche argument in the spectral packet also checks out.

The paper's reversed-stream state and the substitution-length coordinates
are different. With the author notation, `U B=C^-1` and `U^-1=B C`.
For `v=(-103,56,0)`, they are respectively `(56,0,-103)` and
`(-47,56,0)`. The draft's erroneous last coordinate -103 in the latter was
self-corrected by its author, independently confirmed, and covered by a
literal exact control. Negative coordinates are not positive-word counts.

Norm scalars, one stream sample, and different degrees are not equivalent
objects. The norm-determinant identity for K-linear matrices is correct;
its static feedback consequence solves `x=u+gamma D x`. It is not a theorem
that arbitrary additive dynamical feedback is conjugate to the original
evolution, nor a stability theorem. This scope correction was applied.

## 2. Prefix edges, base fields, and exact contraction

Use edges `j -> i` for occurrences of i in sigma(j), with column j the
source. For sigma=(12,13,1), natural lengths are
`(1,beta-1,beta^-1)` and their contracting images are
`(1,alpha-1,alpha^-1)`. Its five edge digits give exactly

    R1=alpha(R1 union R2 union R3),
    R2=alpha R1+1, R3=alpha R2+1.

The real beta field and nonreal alpha field are isomorphic as abstract
fields but are distinct embedded subfields. In the real plane with basis
`(1,alpha)`, coefficients `p=beta^-1,s=1-beta` belong to the real beta field,
and `alpha^2=s alpha-p`. This is not identifying beta with alpha.

The base norm is complex modulus, not automatically the standard Euclidean
norm on the nonorthogonal coordinates. Its Gram matrix is

    Q=[[1,s/2],[s/2,p]],   A^T Q A=p Q.

For `S_h=[[1,h],[0,1]]`, the transported Gram matrix
`Q_h=S_h^-T Q S_h^-1` satisfies `B_h^T Q_h B_h=p Q_h` exactly.
This correction was added to the spectral packet. It also avoids using
`G_h` both for the feedback gain and the positive-definite metric.

The module `Lambda=Z+Z alpha+Z alpha^2` is rank three, not a planar lattice.
Its real intersection is Z. Since `S_h alpha-alpha=h` and
`S_h alpha^2-alpha^2=hs`, its inclusion in itself forces h and hs integral,
hence h=0. Rational coefficients give the stronger corresponding assertion
for `Q(alpha)`. Nonzero shears require the transported arithmetic module.

## 3. Affine branch conjugacy: a necessary coincidence correction

For common invertible S and component translations b_i, the equations

    B S=S A,
    t'_e=S t_e+b_i-B b_j, e:j -> i,

are exactly necessary and sufficient for ambient labeled branch conjugacy.
They transport path codes and the compact component tuple. However, with
different b_i, they preserve address equality only within a fixed target
component, or on the disjoint labeled union. They need not preserve a
spatial coincidence between different components.

**Minimal hostile:** two components, each with its self-loop `x -> x/2`,
have R1=R2={0}. Taking b1=0,b2=1 gives transformed components {0},{1}.
Both branch-conjugacy equations hold; the spatial cross-component equality
does not survive. The author has incorporated this repair and counterexample.
For the advertised shear family all b_i=0, one common ambient map acts on
every component, so the full spatial coincidence relation does survive.

## 4. Twisted substitution and Fourier gauge directions

For sigma'=(21,13,1), only digits 1->1 and 1->2 change, to alpha-1 and0.
The obstruction to a common-linear, edge-preserving affine conjugacy is
stronger than needed: its digit equations alone are inconsistent, even
before imposing that S commute with A.

Indeed unchanged edges 2->1 and3->1 force `b1=A b2=A b3`. Subtracting
the equations for 2->3 (digit1 to1) and1->2 (digit1 to0) yields
`1=A(A-I)b2`. The self-edge1->1 then demands
`alpha-1=A(I-A)b2=-1`, impossible since alpha!=0.
This excludes precisely the specified ambient, edge-preserving common-linear
affine conjugacy, not arbitrary relabelings or measurable identifications.

With positive Fourier phase and column-source convention,

    H'_n(y)=D_b(y) H_n(S^T y) D_b((B^T)^n y)^-1

has the correct signs, transpose directions, and endpoint labels. The
induction on path digits is

    t'_path=S t_path+b_target-B^n b_source.

The independent checker tests all indexed paths of lengths0 through6, five
gains (including10), and distinct nonzero b_i: 1,135 endpoint identities.
It does not discard path multiplicities. This verifies the actual telescoping
claim rather than only its b=0 special case. Pure dynamical spectrum does
not follow from this covariance identity.

## 5. The added resolvent-norm join

Let gamma=a/b be reduced with b>0 and

    lambda=b beta/(b-a beta),
    C_gamma=(I-gamma C)^-1 C.

Then `Q(lambda)=Q(beta)`, since `beta=lambda/(1+gamma lambda)`.
The nonzero denominator is automatic for irrational beta and rational gamma.
The norm and coprimality give

    N(lambda)=(-1)^(k+1)b^k/F_k(b,a),
    gcd(F_k(b,a),b)=1,
    lambda algebraic integer iff F_k(b,a)=+/-1.

Necessity uses integer norm, not the false general converse that an integer
norm implies integrality. Here coprimality forces the denominator to be a
unit. Conversely norm +/-1 makes b-a beta a unit already in Z[beta], so
lambda belongs to that order. The monic minimal polynomial is
`F_k(bY,b+aY)/F_k(b,a)`; the preserved degree proves it is minimal.

**CITED completeness, independently reopened:** Knapp, Lemos and Neumann,
[*Integral values of generating functions of recurrence sequences*](https://math.loyola.edu/~mknapp/papers/General%20Quadratic%20-%202-10-23.pdf),
Lemma2.4 p7 and Example4.1 pp10-11, gives exactly the relevant cubic
norm-one solutions. Its coordinates are reversed from F_k(b,a). Applying
the sign symmetry and b>0 leaves exactly gamma=0,1/2,56/103. This is a
cited global classification, not an extrapolation of bounded enumeration.

The corresponding cubic roots are

    beta, 2 beta^4, -103 beta^17.

The first two are positive Pisot numbers. The third is negative; its negative
is Pisot, so it must not be called a positive Pisot number. Both nontrivial
internal contractions follow from beta>9/5:
`2 beta^-2<1` and `103 beta^(-17/2)<1`. Their exact monic polynomials are

    Y^3-22Y^2-20Y-8,
    Y^3+3249959Y^2+2875039Y+1092727.

The field is retained but the generated order need not be. The conjugate
difference formula

    lambda_i-lambda_j =
      b^2(beta_i-beta_j)/((b-a beta_i)(b-a beta_j))

gives `disc(lambda)/disc(beta)=b^(2k(k-1))/F_k(b,a)^(2(k-1))`.
In the integral cases, `[Z[beta]:Z[lambda]]=b^(k(k-1))`, respectively
64 and103^6=1194052296529 for the two nontrivial cubic gains.
These formulas are independently checked in exact cubic arithmetic.
They do not assert equality with the maximal order.

The resolvent changes the spectrum and is not the earlier isospectral
shear absorption. It may define a new contracting internal graph system,
but does not identify that system with the original substitution attractor.

An additional exact join in the author's final checker is also sound:
the constant power-basis coefficient of lambda is
`U_k(gamma)=gamma^(k-1)/(1-gamma-...-gamma^k)=b a^(k-1)/F_k(b,a)`.
To avoid a convergence assumption, compare coefficients in
`(1-gamma beta)lambda=beta`. They give `c0=gamma c_(k-1)` and, for i>=1,
`c_i=gamma c_(i-1)+gamma c_(k-1)+delta_(i,1)`, hence the displayed rational
identity. Coprimality with both a and b makes integrality of this one scalar
equivalent to integrality of the entire operator **inside this specific
rational-gain family**. It is not a generic equivalence from one sample.

**Final-byte corollary review:** the fixed-S extension is sound: replacing
integer norm by a Z[S^-1]-integral norm forces F_k(b,a) to be an S-unit;
the adjugate proves the converse. Under this gate the operator is
S-unimodular exactly when b is also an S-unit. The fixed-S finiteness
statement remains explicitly CITED Thue--Mahler input, not a computation.
The rational operator-category change of generator is likewise sound:
`T -> T(I-gamma T)^-1` and `D -> D(I+gamma D)^-1` have identity actions on
intertwiner arrows, and inverse object maps. Irreducibility excludes the
linear poles. This changes the time generator, not its similarity class.
The claimed integral limitation is witnessed by the proper lattice Z[lambda],
which is lambda-stable but not alpha-stable when the computed order index
exceeds one. These two final paragraphs require no change to the checker.

## 6. Nonunit finite-place addendum

The separate `../nonunit_feedback/` addendum is sound without further
correction. Repeating each sigma^n(j) exactly q times gives incidence qM^n,
expansion lambda=q beta^n, and internal contraction q alpha^n. Its ordered
prefix translations are the originals plus r alpha^n ell_j(alpha), so word
order is specified, not inferred from the matrix. The cases are (q,n)=(2,4)
and(103,17); the latter realizes the sign-flipped positive operator, not
the raw negative resolvent.

The planar determinant is delta=q^2 beta^-n. If the pieces in each target
equation are essentially disjoint and m denotes finite component areas,
then m=delta N m. A strictly positive left Perron vector gives
`ell m=delta lambda ell m=q^3 ell m`. Thus all areas vanish, or some pair
of edge pieces has positive-measure overlap. No area value is thereby
computed. The distinction between within-equation and cross-component
intersection is correctly maintained.

Since beta is a unit in the maximal order as well, the finite places
dividing q beta^n are exactly those dividing q. Their normalized Haar
moduli multiply to q^-3, so the complete internal scaling is
`delta q^-3=lambda^-1`. This fixes the necessary measure balance but
does not prove tiling or purity. The representation-space normalization
was directly rechecked against [Minervino--Thuswaldner, Section3 and
Lemma3.1](https://arxiv.org/pdf/1402.2002), independently of the author audit.

## 7. Verification and remaining boundary

`verify_independent.py` uses only the standard library, no author imports,
exact rational arithmetic in the real beta field, and active exceptions.
`normal.json` and `optimized.json` are byte-identical. Its finite controls
also include a nonintegral hostile gain1/3, all cited cubic Thue pairs and
their sign reversals, the Mobius inverse, and independent minimal polynomial
and power-basis determinant computations. General proofs, not finite counts,
establish the all-real-h and all-rational-gamma statements.

The precise residual claims not established are: symbolic conjugacy,
pure dynamical spectrum, Rauzy tiling/disjointness, a dynamical equivalence
with Amanda's threshold loop, arbitrary-feedback absorption, and any
equivalence between number fields of degrees3,4,5.
