# One extra failure cannot repair the mixed stress: the 29-adic bound rises to 445

Date: 2026-10-19. Status: **PROVED**, with **FINITE-EXACT** controls and independent audit PASS.

## Outcome and boundary

For the norm-one graph over F_(29^2), every section to its unramified
precision-two graph fails on at least **445 of the 12,615 edges**.

More specifically, the three direction families

    u=(1,0), v=(0,1), w=(3,9),   Q(x,y)=x^2+xy+y^2

force at least **89 failures**, improving the preceding proved bound of 88.
There is no section realizing 88 failures on these three families. Neither
89 nor 445 is asserted to be sharp.

The mechanism is sparse-error stability. With just one failure beyond the
87 disjoint-cycle lower bound, the extra error is a two-point dipole. The
mixed stress confines all nonzero mixed second differences to an eight-point
cube. Its line sums force a very specific form. That form can be removed by
an explicit local infinitesimal motion, leaving the impossible straight-seam
configuration from the preceding checkpoint.

No new chromatic bound or Euclidean realization is claimed. This continues
only the prescribed-section edge-defect problem.

## 1. Inherited algebra, stated explicitly

Put p=29 and V=F_p^2. For d in {u,v,w}, write

    D_d f(x)=f(x+d)-f(x),
    ell_d(h)=(2d_x+d_y)h_x+(d_x+2d_y)h_y.

For fixed residue representatives, the first carry c and the error vector of
any section satisfy

    z_d=c_d+ell_d D_d h.

Thus the admissible errors form the affine coset c+im(R), where
(Rh)_d=ell_d D_d h. The error sum around every straight p-cycle is -2.
There are p such cycles per direction family.

The previous mixed-stress checkpoint proved the identity

    -3 D_v D_w z_u -9 D_u D_w z_v +D_u D_v z_w=0.     (1)

Both the rigidity term and the carry satisfy it. In brief, if
X_(x+d)-X_x=d-p k_d(x), then

    c_d=constant-ell_d k_d,   D_e k_d=D_d k_e,
    ell_w=3 ell_u+9 ell_v.

Commuting differences therefore proves (1), including its carry term.

The balanced-digit implication will be used repeatedly:

    A,B,C in {-2,-1,0,1,2},  A-3B-9C=0 mod 29
        implies A=B=C=0.                              (2)

Indeed, the left side has absolute value at most 26 and hence is zero as an
integer. Reducing modulo 3 and then dividing by 3 gives the assertion.

The previous checkpoint also proves that a binary one-error-per-cycle
support whose three mixed derivatives vanish over the integers is a triple
of straight seams. All eight possible orientations have horizontal period
obstructions

    16,14,4,2,27,25,15,13 modulo 29,                   (3)

independent of seam offsets. In particular, if J_d is binary with one point
on every d-cycle and

    D_vD_w J_u=D_uD_w J_v=D_uD_v J_w=0,

then -2J does not belong to c+im(R).

For reference, the orientations are encoded by eta_u,eta_v,eta_w in {0,1}:
J_u is v-parallel or w-parallel, J_v is u-parallel or w-parallel, and J_w is
u-parallel or v-parallel, respectively. Their period is

    16-18 eta_u+17 eta_v-2 eta_w modulo 29.

This obstruction concerns the full carry coset, not just (1).

## 2. Normal form for an 88-failure candidate

Suppose a section had exactly 88 failures on these three families. Of the
87 straight cycles, exactly one has two failed edges; each other cycle has
one. Let e be the exceptional direction. Let P,Q be the two failed edge
origins on its exceptional cycle, and choose their order so that

    Q=P+k e,   1<=k<=14.

Choose a binary field I_d containing one selected failure on every d-cycle,
using P on the exceptional cycle. The single failures have error -2. If the
exceptional errors at P,Q are alpha,beta, then alpha+beta=-2 and beta!=0.
Consequently, for some nonzero t,

    z_d=-2 I_d+t 1_(d=e) K,
    K=delta_P-delta_Q.                               (4)

Set

    S=sum_(j=1)^k delta_(P+j e).

Our translation convention gives D_e S=K. Define the alternating cube

    C=D_u D_v D_w S.

It has eight distinct support points and values +1 and -1. It is the double
difference of the dipole K in the two directions other than e.

For each d, let H_d be the mixed derivative of I_d in the other two
directions. For example, H_u=D_vD_w I_u. Substituting (4) into (1), with

    (s_u,s_v,s_w)=(-3,-9,1)

used here only as stress coefficients, gives

    -2 sum_d s_d H_d+t s_e C=0 modulo p.              (5)

Outside the cube support, (2) forces all
three integer H_d to be zero. It remains to understand them on eight points.

## 3. Eight-point cube lemma

**Lemma.** In the preceding situation there are integers
lambda_u,lambda_v,lambda_w in {-2,-1,0,1,2} such that

    H_d=lambda_d C over the integers,
    sum_d s_d lambda_d=t s_e/2 modulo p.              (6)

### Geometry behind the lemma

The determinant magnitudes of the three direction pairs are 1,3,9 modulo
29. They are nonzero and distinct up to sign. In particular, for any one
direction, the projections of the other two onto the quotient by that
direction are neither equal nor opposite.

The cube is the set of subset sums of three axes: the exceptional axis is
k e, and the other two are the corresponding unit directions. Translation
and signs do not matter. All eight points are distinct. Indeed, a collision
would make a nonzero multiple of e equal to a sum of the other two directions
with coefficients in {-1,0,1}; independence excludes a zero coefficient,
and the unequal-up-to-sign quotient projections exclude the remaining case.

For any pair of directions, form the bipartite incidence graph whose vertices
are their parallel lines meeting the cube, and whose edges are the cube's
points. The two faces in that pair of directions form two four-cycles. They
share at most one incidence vertex:

- If the pair includes e, the two faces have disjoint line vertices on one
  side. This follows because the third direction's coordinate on the other
  unit direction is neither 0 nor +/-1. Their two-element line sets on the
  other side share at most one vertex, since the remaining coordinate and k
  are nonzero and p is odd.
- If the pair excludes e, both coordinates of the face displacement k e are
  nonzero. Sharing a line vertex on both sides would require both coordinates
  to be +/-1, contrary to their ratio being different from +/-1.

A function on the cube with zero sum on each of those parallel lines is
therefore a linear combination of the alternating functions on the two
faces. This can be checked by eliminating the degree-two vertices of the
two four-cycles, even when they share a vertex. There are exactly two free
face coefficients. This argument holds over the integers and over F_29.

### Coupling the three face pairs

Every H_d has zero sum on every line in each of its two derivative directions,
by telescoping. The preceding incidence argument describes it by two face
coefficients, indexed by the cube coordinate in direction d.

At a cube vertex indexed by (i,j,l) in {0,1}^3, divide (5) by the alternating
cube sign. The result is a sum of three functions, each depending on just one
of i,j,l, equal to a constant. Each function must be constant. As every s_d
is nonzero, this says H_d=lambda_d C modulo p, and (5) gives the second part
of (6).

Every H_d value lies in [-2,2], since I_d is binary; every nonzero C value is
+1 or -1. Taking lambda_d from one cube vertex therefore places it in
[-2,2]. The equality modulo 29 is then an equality of integers at every
vertex: the possible difference has absolute value at most four. This proves
the lemma.

The proof is uniform in all three choices of e and all nonzero k. The
short-separation condition k<=14 is needed in the next step, not here.

## 4. A short-column binary-normalization lemma

**Lemma.** Let p>=5. Suppose an integer array J on a p by p grid satisfies

    J(s,t)=F(s)+G(t).

Suppose I is binary, has total weight p, and I-J is supported in at most two
columns, in at most p-2 positions of each. If sum(I-J)=0, then J is binary
and has weight p.

**Proof.** The weight statement is immediate. If F is constant, each column
has an unmodified position, so its constant J value is 0 or 1.

Otherwise an unmodified column shows that, after changing the additive
constant, F is a nonconstant binary function. Every unmodified column must
then have the same values F: shifting a nonconstant binary function by a
nonzero constant cannot keep it binary. If m is the number of ones of F,
these at least p-2 columns contribute (p-2)m ones to I. Since I has only p
ones and p>=5, m=1.

An exceptional column has at most p-2 modified positions. Among the p-1
zero positions of F, at least one is unmodified, so its G value is 0 or 1.
All unmodified columns have G=0. But sum J=p and sum F=1 imply sum G=0.
Thus both exceptional G values are zero. Hence J=F is binary. QED.

The strict support condition matters. At p-1 changes per exceptional column,
a nonbinary additive array can be concealed by a binary array of weight p.
The verifier includes an explicit counterexample at that boundary.

## 5. Remove the dipole by a local motion

By the cube lemma, define

    J_d=I_d-lambda_d D_d S.                           (7)

Its mixed derivative in the other two directions vanishes over the integers,
so it has the additive form in the binary-normalization lemma.

For d!=e, D_d S is the signed difference of two parallel e-segments, each of
length k<=14. In coordinates of the two derivative directions, these occupy
two columns with at most 14 modified positions each. For d=e, D_e S=K has
just two points, occupying two distinct columns in either derivative grid.
Moreover sum D_d S=0 and sum I_d=p. The lemma therefore applies to every
J_d, making all three binary with weight p.

The correction D_d S also sums to zero on each d-cycle. Thus J_d still has
one selected point per d-cycle, not merely the right total weight. Each J_d
is consequently a straight seam of one of the orientations listed in (3).

Now put

    L_d=-2 lambda_d+t 1_(d=e).

Equation (6) implies

    -3L_u-9L_v+L_w=0 modulo p.

Choose a constant vector H with ell_u(H)=L_u and ell_v(H)=L_v; it exists
because the projection matrix [[2,1],[1,2]] has nonzero determinant 3.
Then ell_w(H)=L_w too. Equations (4) and (7) become

    z_d=-2J_d+L_d D_d S,
    z=-2J+R(H S).                                    (8)

This is an explicit local motion, not merely a vanishing-stress assertion.
Since z belongs to c+im(R), (8) would put -2J in that same affine coset.
Equivalently, subtracting H S from the section correction would construct a
section with just 87 failed edges, already excluded by the previous theorem. But
all straight-seam orientations fail the nonzero period (3). Contradiction.

**Theorem.** An 88-failure section on {u,v,w} does not exist. Combined with
the previously proved 88 lower bound, every section has at least 89 failures
on these three direction families.

## 6. The full graph requires at least 445 failures

Rotate the triple by all 15 norm-one direction pairs modulo sign. Norm-one
rotations lift modulo p^2, so the theorem applies to each rotated triple.
Every direction family occurs exactly three times. Therefore, if M is the
full section defect,

    3M >= 15*89,
    M >= 445.

No disjointness of the rotated triples is needed. This remains a bound on
first-carry section defect in the unramified quadratic norm graph.

## 7. Exact controls and a useful hostile

Run

    python3 verify_one_extra.py

It uses only the standard library. Its controls include:

- All 42 short dipole geometries: three exceptional directions and fourteen
  separations. Every cube has eight points; every incidence pair is two
  four-cycles sharing at most one vertex; the three pair kernels each have
  dimension two; the coupled eight-point system has rank four and exactly
  the two homogeneous degrees of freedom predicted by (6).
- Eighty-four explicit 88-error arrays, with both choices of exceptional seam
  orientation for every geometry. They pass every straight-cycle check and
  every compact mixed-stress check. An explicit local motion carries each
  back to a straight seam, and its global period is directly checked nonzero.
- The boundary hostile to binary normalization, using p-1 changes per
  exceptional column.

The 88-error arrays are **not sections**. Their role is important: merely
finding an 88-error vector satisfying the compact stress, or merely adding
local checks around the exceptional edges, would not establish attainment.
The preserved nonzero carry period is the final obstruction.

The cube argument is proved above; the rank checks are independent finite
controls rather than a substitute for the incidence proof. No solver
infeasibility result is used.

This is a separate continuation of the reviewed [balanced-stress checkpoint](../balanced_digit_section_defect_20261019/README.md)
that raised 435 to 440. The earlier packet is left unchanged. The ingredients
inherited from it are explicitly identified in section 1; the new ingredients
are the dipole normal form, cube rigidity, short-column binary normalization,
and the explicit removal by a local motion. No broad literature-priority
claim is made.

## Independent reproduction

Read the [independent proof audit](review/REVIEW.md). Run the submitted checker with assertions enabled: `python3 -B verify_one_extra.py`. Do not pass `-O` to this submitted checker. Compare stdout with `verify_one_extra.out` and its generated JSON with `one-extra-exact-results.json`.

Run `python3 -B review/independent_controls.py`, then `python3 -B -O review/independent_controls.py`. Compare stdout and generated JSON with the corresponding normal/optimized review receipts. These 328,679 independent checks use explicit failures that remain enabled under optimization. Both checkers use only Python 3 and its standard library. `SOURCE_PROVENANCE.json` preserves reviewed and published hashes; no mathematical proof or verifier logic was changed.
