# Two extra errors: a three-failure cycle is impossible, and disjoint cubes can interact

Date: 2026-10-19. Status: **PROVED, independently audited, with a finite-exact p=29 lemma**.

## Outcome

For the three norm-one direction families

    u=(1,0), v=(0,1), w=(3,9),   Q(x,y)=x^2+xy+y^2 over F_29,

a section to the unramified precision-two norm graph cannot have exactly
89 failed edges concentrated as one failure on every straight cycle except
for **one cycle with three failures**.

Consequently, any section attaining 89 failures on this triple must have
**two distinct cycles with two failures each**, and one failure on every
other cycle. We also prove an explicit small-rank correction criterion that
excludes many of these remaining configurations.

There is a real limitation: two eight-point cube supports can be disjoint
while their coupled incidence kernel has dimension six, rather than the
four dimensions supplied by the two original dipoles. An explicit hostile
below exhibits and completely describes its two additional connecting modes.

The established numerical lower bounds remain **89 on the triple and 445
on the full norm graph**. This packet proves neither 90 nor 450, constructs
no attaining section, and makes no chromatic-number claim.

## 1. Setup and inherited period obstruction

Put p=29 and V=F_p^2. For d in {u,v,w}, define

    D_d f(x)=f(x+d)-f(x),
    ell_d(h)=(2d_x+d_y)h_x+(d_x+2d_y)h_y,
    (Rh)_d=ell_d D_d h.

The section error vector has the form z=c+Rh, where c is the first carry.
The error sum around every straight p-cycle in every direction is -2.
There are p cycles per direction family, hence 87 cycles on the triple.
Each cycle must contain at least one nonzero error.

We use the following previously proved facts from the mixed-stress and
one-extra-error checkpoints, restated here to specify the dependency:

1. The error of every section obeys the compact mixed stress

       sum_d sigma_d D_a D_b z_d=0,
       sigma=(-3,-9,1),   {a,b}={u,v,w}\{d}.           (1)

   This includes the carry: if X_(x+d)-X_x=d-p k_d(x), then
   c_d=constant-ell_d k_d, D_a k_d=D_d k_a, and
   ell_w=3ell_u+9ell_v. Commuting differences proves (1).

2. For A,B,C in [-2,2] integral,

       A-3B-9C=0 modulo 29  implies  A=B=C=0.        (2)

   The integer expression has absolute value at most 26, so it vanishes
   integrally. Reduction modulo 3, twice, then proves the assertion.

3. A binary J_d with one point per d-cycle and D_aD_b J_d=0 integrally
   is a straight seam parallel to a or b. The eight possible seam triples
   have nonzero horizontal carry periods

       16,14,4,2,27,25,15,13 modulo 29.              (3)

   More explicitly their period is
   16-18 eta_u+17 eta_v-2 eta_w, with eta recording the two seam orientations
   in each family. Seam offsets do not affect this period. Therefore -2J
   does not lie in c+im(R).

4. Every section on the triple has at least 89 failures; rotation and
   multiplicity-three counting give at least 445 on the full graph.

This packet continues the prescribed-section defect problem in the
anisotropic unramified quadratic norm plane. Fact (3) is a global carry
obstruction: the compact stress alone is not sufficient for a section.

## 2. Incidence spaces and a complete twelve-point certificate

For a finite support T in V, let W_d(T) be the F_29 vector space of functions
supported on T whose sums vanish on every line parallel to each of the two
directions other than d. Define

    L_T : W_u(T) + W_v(T) + W_w(T) -> F_29^T,
    L_T(H_u,H_v,H_w)=sum_d sigma_d H_d.              (4)

Here the domain is an external direct sum. The dimensions below concern
this explicitly specified linear map, not a graph-coloring search.

Fix e in {u,v,w} and three distinct points P,Q,R on one e-cycle. Set

    K_1=delta_P-delta_Q,   K_2=delta_P-delta_R,
    C_j=D_aD_b K_j,        {a,b}={u,v,w}\{e},
    T=supp(C_1) union supp(C_2).

**Finite-exact twelve-point lemma.** T has 12 points. Each C_j has eight
points with values +/-1, their supports intersect in four points, and they
agree on that intersection. Moreover

    dim ker L_T=4,
    ker L_T = { H_d=lambda_d1 C_1+lambda_d2 C_2 :
                sum_d sigma_d lambda_dj=0 for j=1,2 }.               (5)

### Universe and certificate

Translation lets P=0. Write Q=q e and R=r e. Interchanging Q and R changes
no assertion, so the complete universe is

    e in {u,v,w},   1<=q<r<=28.

There are exactly 3*binomial(28,2)=1,134 geometries. No further symmetry
quotient or unexamined generic-position assumption is used.

For each geometry the accompanying standard-library verifier:

- forms the integer sparse arrays C_1,C_2 directly using D_d;
- checks the support and sign assertions above;
- constructs the line-sum matrix for each W_d, using cross(d,x) as its
  line label, and computes its nullspace by exact elimination in F_29;
- computes the rank of L_T on those nullspace bases;
- separately forms all line and coupling constraints on the 36 entries
  (H_u,H_v,H_w), and checks that this direct matrix has rank 32;
- verifies the four independent model vectors
  (C_j,0,3C_j) and (0,C_j,9C_j), j=1,2.

Thus (5) follows by inclusion and equal dimension. Every pivot inverse is
computed modulo 29; there are no floating-point tolerances, external
packages, optimization solvers, or infeasibility inferences.

The dimension histogram is:

| e | (dim W_u,dim W_v,dim W_w) | number | rank L_T |
|---|---|---:|---:|
| u | (3,4,4) | 372 | 7 |
| u | (3,4,5) | 3 | 8 |
| u | (3,5,4) | 3 | 8 |
| v | (4,3,4) | 372 | 7 |
| v | (4,3,5) | 3 | 8 |
| v | (5,3,4) | 3 | 8 |
| w | (4,4,3) | 372 | 7 |
| w | (4,5,3) | 3 | 8 |
| w | (5,4,3) | 3 | 8 |

In particular, the exceptional extra line-incidence modes in the 18
nongeneric geometries do not survive the coupled stress. This lemma is a
finite-exact result at the stated prime and directions; no uniform-prime
version is claimed.

## 3. Two binary-normalization lemmas

All arrays and sums in this section are over the integers.

### Column version

**Lemma.** Let I be binary on a p by p grid with total weight p, and let
J(s,t)=F(s)+G(t) be integral. Suppose sum(I-J)=0. If I-J is supported in at
most r columns, each containing at most p-2 modified positions, and p>2r,
then J is binary with weight p.

**Proof.** The total weight is immediate. If F is constant, every column
retains an unmodified position, so its constant J value is 0 or 1.
Otherwise an unmodified column permits normalization with F nonconstant
binary. Every unmodified column then has G=0: a nonzero integral shift of
a nonconstant binary column cannot be binary. If F has m ones, the at least
p-r unmodified columns imply (p-r)m<=p. Since p>2r and m>=1, m=1.
Every exceptional column retains an unmodified zero of F, because F has
p-1 zeros and at most p-2 positions are changed. Its G value is therefore
0 or 1. But sum J=p and sum F=1 imply sum G=0. All G vanish. QED.

The p-2 condition cannot be replaced by p-1 on this information alone.
An explicit weight-p binary I masking a nonbinary additive J at exactly
p-1 changes in each of two columns is replayed by the verifier.

### Crossing rows and columns

**Lemma.** Suppose I,J satisfy the same weight, additivity, and zero-total-
correction conditions. Let I-J be supported in the union of r designated
rows and s designated columns. Suppose

    p>2r, p>2s, (p-r)(p-s)>p;

and that each designated row has at most p-s-2 modified positions outside
the designated columns, while each designated column has at most p-r-2
modified positions outside the designated rows. Then J is binary.

**Proof.** On the clean (p-r) by (p-s) rectangle, J is binary and additive.
Therefore at least one of the two additive factors is constant there.
The rectangle cannot be all ones, since its area exceeds the total weight.

If both factors are constant, the rectangle is all zero. Normalize both
factors to zero on its row and column index sets. Each designated row and
column retains a clean unmodified point, so all remaining factor values
are 0 or 1. Since sum J=p, the sum of all factor values is one. Their
additive array is therefore binary.

Otherwise, suppose the row factor F is nonconstant on the clean rectangle;
the other case is symmetric. Normalize F there to values 0,1 and G=0 on
clean columns. If there are m clean rows with F=1, then (p-s)m<=p; p>2s
forces m=1. Each designated row has an unmodified clean-column position,
so its F value is 0 or 1. Among the p-r-1 clean zero rows, each designated
column has an unmodified position by the assumed p-r-2 bound. Its G value
is 0 or 1. Again sum J=p says sum F+sum G=1, forcing all G=0 and no
additional ones of F. QED.

At p=29, two designated rows and two designated columns with at most 14
changes on each strip satisfy these hypotheses: the thresholds are 25,
and the clean rectangle has area 27^2>29.

## 4. Excluding the one-three-failure-cycle case

Suppose a section has 89 failures and one e-cycle has three of them. Every
other cycle has one. Choose one exceptional origin P and call the other two
Q,R. Let I_d select one failed origin on every d-cycle, using P on the
exceptional cycle. The straight-cycle sums give

    z_d=-2 I_d+1_(d=e)(t_1 K_1+t_2 K_2),             (6)

where t_1,t_2 are nonzero in F_29 and K_j are the dipoles from Section 2.
Let H_d=D_aD_b I_d for the other two directions. Each H_d is integral with
values in [-2,2]. Substituting (6) in (1) gives

    sum_d sigma_d H_d = (sigma_e/2)(t_1 C_1+t_2 C_2).                (7)

Outside T, (2) forces H_u=H_v=H_w=0 integrally. Every H_d has zero sums on
its two derivative line families, by telescoping. A particular solution of
(7) is H_e=(t_1C_1+t_2C_2)/2 and the other components zero. The twelve-point
lemma consequently implies

    H_d=lambda_d1 C_1+lambda_d2 C_2 modulo 29,
    sum_d sigma_d lambda_dj=sigma_e t_j/2.            (8)

The Q-face of C_1 is disjoint from C_2, and the R-face of C_2 is disjoint
from C_1. Reading H there chooses each lambda_dj as an integer in [-2,2].
On T the right side of (8) then has absolute value at most four; its
difference from H_d has absolute value at most six, less than 29. Thus the
first equality in (8) holds over the integers everywhere.

Now choose the baseline P more specifically: take it immediately after a
largest gap in the cyclically ordered three failed origins. Their three
positive gaps sum to 29, so a largest gap is at least 10. The remaining arc
has length at most 19. Write Q=P+k_1 e, R=P+k_2 e within that arc, and set

    S_j=sum_(l=1)^(k_j) delta_(P+l e),   1<=k_j<=19.

Then D_e S_j=K_j and D_uD_vD_w S_j=C_j. Define

    J_d=I_d-sum_(j=1)^2 lambda_dj D_d S_j.           (9)

Its mixed derivative in the other two directions vanishes integrally, so
in the grid generated by those directions it has the form F(s)+G(t).
Its total weight is p, since each derivative has total zero.

For d!=e, the two S_j are nested initial segments on one e-line. The union
of the correction supports lies in two e-parallel columns with at most 19
positions each. For d=e, it consists of P,Q,R, which lie in three distinct
columns of either derivative grid, with one position each. The column
lemma applies in every case and makes J_d binary. Also D_d S_j sums to
zero on each d-cycle, so J_d still has exactly one point per d-cycle.

For j=1,2, put

    L_dj=-2 lambda_dj+t_j 1_(d=e).

The second part of (8) says sum_d sigma_d L_dj=0. The matrix of ell_u and
ell_v is [[2,1],[1,2]], with determinant 3 nonzero in F_29. There is a
constant vector h_j satisfying ell_d(h_j)=L_dj for all three d, because
ell_w=3ell_u+9ell_v. Equations (6) and (9) now give the explicit identity

    z=-2J+R(h_1 S_1+h_2 S_2).                       (10)

As z lies in c+im(R), so would -2J. But J is a binary one-per-cycle triple
with vanishing mixed derivatives, contradicted by the global period (3).

**Theorem.** A section with 89 failures on this triple cannot have one
three-failure cycle. Since every cycle must contain a failure, the only
remaining 89-failure distribution is two distinct double-failure cycles.

## 5. A reusable two-dipole correction criterion

The remaining distribution has exceptional cycles in directions e_1,e_2,
possibly equal. Choose the order of the two failed origins on each cycle so
that Q_j=P_j+k_j e_j with 1<=k_j<=14. Set

    K_j=delta_(P_j)-delta_(Q_j),
    S_j=sum_(l=1)^(k_j) delta_(P_j+l e_j),
    C_j=D_uD_vD_w S_j.

Choosing I at P_j on each exceptional cycle gives

    z_d=-2I_d+sum_j t_j 1_(d=e_j)K_j.                (11)

**Criterion.** Assume the two eight-point supports supp(C_1),supp(C_2) are
disjoint, and their union T satisfies dim ker L_T=4. Assume also that, in
each component's derivative grid, the union of the supports D_d S_1 and
D_d S_2 satisfies either binary-normalization lemma of Section 3. Then no
section with the error pattern (11) exists.

**Proof.** Each C_j lies in every W_d(T). The four model vectors
(C_j,0,3C_j), (0,C_j,9C_j) are independent, because the supports are disjoint.
They therefore span ker L_T by its assumed dimension. Apply (2) off T, and
then use the inhomogeneous version of (5), to obtain

    H_d=sum_j lambda_dj C_j,
    sum_d sigma_d lambda_dj=sigma_(e_j)t_j/2.        (12)

Disjointness lets each coefficient be read at an exclusive cube point, so
lambda_dj lies in [-2,2] and (12) lifts to integer equality. Normalize
J_d=I_d-sum_j lambda_dj D_d S_j using the specified support hypotheses.
Its d-cycle sums remain one. Taking
L_dj=-2lambda_dj+t_j 1_(d=e_j) constructs constant h_j and identity (10)
just as before. The same period obstruction finishes the proof. QED.

Useful sufficient cases for the normalization part are:

- If e_1!=e_2, each component's correction is contained in at most two rows
  and two columns, each with at most 14 modified positions outside their
  intersections. If d equals an exceptional direction, that dipole is only
  two points and can be assigned to two rows or columns. Thus the crossing
  lemma applies automatically at p=29.
- If e_1=e_2=e, then for d!=e there are at most four e-parallel columns.
  The column lemma applies whenever each contains at most 27 positions in
  the union of the correction supports. Failure of this sufficient support
  condition requires two length-14 segments sharing a column with disjoint
  position sets, covering 28 positions. For d=e there are at most four
  individual points, so the normalization condition is automatic.

Disjoint derivative-line incidence is a simple sufficient way to obtain
the rank condition: if no line parallel to u, v, or w meets both cube
supports, the line systems and the pointwise coupling split by cube. Each
single cube contributes homogeneous dimension two, hence the union has
homogeneous dimension four. The verifier directly checks this single-cube
statement for all 42 short geometries (three directions, fourteen lengths). This is stronger than necessary. The direct
small-matrix test in (4) permits configurations with shared lines as well.

For a concrete nonvacuous example, take (e_1,P_1,Q_1)=(u,0,u) and
(e_2,P_2,Q_2)=(v,(7,5),(7,6)). The supports are disjoint, the pair kernel
dimensions are (4,4,4), and L_T has rank eight, hence nullity four. The
criterion excludes this configuration. Translation of the whole pattern
has the same conclusion.

## 6. A disjoint-support hostile, with the missing modes explicit

Let

    C_j=D_uD_vD_w delta_(j v),   j=0,1,2,
    B={0,-u,-w,-u-w}.

Then

    supp(C_0)=B union (B-v),
    supp(C_2)=(B+2v) union (B+v),
    supp(C_1)=(B+v) union B.

The four layers are disjoint, so the first and third cubes have disjoint
eight-point supports. Nevertheless the intermediate cube C_1 is supported
in their sixteen-point union. It is independent of C_0,C_2: the outer
layers B-v and B+2v force both coefficients of any purported combination
to vanish.

For T=supp(C_0) union supp(C_2), exact elimination gives

    (dim W_u,dim W_v,dim W_w)=(6,4,6),
    rank L_T=10,   dim ker L_T=6.                   (13)

The two additional vectors are explicitly

    (C_1,0,3C_1),   (0,C_1,9C_1).

Together with the four endpoint-cube vectors they form a verified basis
of ker L_T. Thus (13) is fully explained by connecting local-motion modes;
it is not merely an unexpected numerical rank.

The endpoint cubes themselves arise from short dipoles on two different
u-cycles: use S_0=delta_0 and S_2=delta_(2v), with failed origins (-u,0)
and (2v-u,2v). The example therefore lies in the relevant two-double-cycle
geometric universe.

This hostile refutes the claim that **disjoint cube supports alone imply
the four-dimensional correction form**. It is a linear-incidence hostile,
not a binary error vector and not a section. For example, one displayed
mode has H_w values +/-3, outside the binary second-difference range. The
additional modes are themselves local-motion modes, so their existence
does not prove that correction or exclusion is impossible; it identifies
what a stronger multi-dipole normalization argument must handle.

## 7. Replay, boundary, and next obstruction

Run

    python3 verify_two_extra.py

The verifier writes `two-extra-exact-results.json` alongside itself. The
JSON gives every one of the 1,134 support geometries and both rank checks,
all largest-gap path lengths and normalization support bounds, the explicit
six-dimensional hostile arrays, and the passing separated example.
`verify_two_extra.out` records the successful run.

The research artifacts are this note, `verify_two_extra.py`, its stdout
receipt, and `two-extra-exact-results.json`. Exploratory probes are not proof dependencies.

The unresolved 89-failure candidates now have two distinct double-failure
cycles. The rank-four criterion excludes the separated configurations that
meet its normalization hypotheses. Remaining possibilities include
interacting cube geometries with extra connecting modes, overlapping cube
supports, and the narrow same-direction normalization crowding case. No
claim is made that any of these remaining possibilities yields a section.

The immediate reusable problem is to control binary normalization after
including all connecting local-motion potentials, rather than assuming
that the two original dipole potentials generate the whole correction.

## Independent audit and reproduction

This continues the [445 section-defect bound](../section_defect_stability_20261019/README.md). Read the [independent proof review](review/REVIEW.md).

Run the submitted `verify_two_extra.py` with assertions enabled; do not use `-O` as a verification mode. For the independent implementation, use Python 3 and its standard library from this package directory:

    python3 -B review/independent_controls.py
    python3 -B -O review/independent_controls.py

Compare stdout with the corresponding normal/optimized `.out` receipts and the generated `review/independent-results.json` with the corresponding JSON receipts. The independent checker uses explicit exceptions and remains effective under optimization. `SOURCE_PROVENANCE.json` records original and published hashes; only status, navigation, reproduction instructions and provenance paths changed. The theorem and mathematical controls are unchanged.
