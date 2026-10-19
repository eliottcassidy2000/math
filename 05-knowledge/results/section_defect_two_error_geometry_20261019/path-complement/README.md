# Path complements remove the same-direction normalization exception

Date: 2026-10-19. Status: **PROVED**, independently audited, with finite-exact support controls.

This is an addendum to the [two-extra-error packet](../README.md). Its
conclusion strengthens the two-dipole criterion there:

> At p=29, two disjoint mixed-cube supports with coupled homogeneous
> incidence dimension four cannot support an 89-failure section. No
> additional sparsity hypothesis is needed.

The numerical lower bounds remain 89 on the mixed triple and 445 on the
full norm graph. The dimension-six disjoint-support hostile and the general
two-double-cycle case remain unresolved by this criterion.

## 1. Notation and the only former exception

Use the directions

    u=(1,0), v=(0,1), w=(3,9) in F_29^2,
    D_d f(x)=f(x+d)-f(x).

For a dipole with endpoints P_j,Q_j=P_j+k_j e_j, choose their order so
1<=k_j<=14, and let

    S_j=sum_(l=1)^(k_j) delta_(P_j+l e_j),
    K_j=D_(e_j) S_j=delta_(P_j)-delta_(Q_j),
    C_j=D_uD_vD_w S_j.

The preceding packet proves that if the two cube supports are disjoint
and their coupled incidence kernel has dimension four, then all candidate
binary mixed derivatives have the form

    H_d=lambda_d1 C_1+lambda_d2 C_2,
    lambda_dj in {-2,-1,0,1,2},

with the projection relations needed to remove the dipoles by a local
motion. The remaining step is to show that the additive integer arrays

    J_d=I_d-lambda_d1 D_d S_1-lambda_d2 D_d S_2         (1)

are binary. They have total weight 29, and all correction terms sum to zero
on every d-cycle.

The r-column normalization lemma from that packet applies whenever the
correction is supported in at most r columns, with at most p-2 positions
in each, and p>2r. For distinct exceptional directions, the crossing-row/
column lemma already always applies to the short paths above.

For equal exceptional directions e_1=e_2=e, the two exceptional e-cycles
are distinct. Writing their affine lines as L_1,L_2, the correction in a
component d!=e lies on

    L_1, L_1-d, L_2, L_2-d.                          (2)

There are at most four columns. The only potential failure of the p-2=27
condition is a shared column containing two disjoint 14-point segments,
whose union has 28 positions. We now remove that exception.

## 2. At most one component can have a shared column

Let bars denote images in the one-dimensional quotient V/<e>. Since
L_1!=L_2 and d is not parallel to e, a coincidence among the four columns
in (2) occurs exactly when

    L_2-L_1 = +bar(d) or -bar(d).                    (3)

For a fixed d, there can be at most one shared column: two coincidences
would imply 2 bar(d)=0, impossible in characteristic 29.

The quotient coordinates of the two possible choices of d!=e, measured
by cross(e,d), are respectively

    e=u:  1,9;
    e=v:  28,26;
    e=w:  20,3.

In each case they are nonzero and unequal up to sign. Thus (3) can hold
for at most one of the two nonexceptional components. The other component
has four unshared correction columns.

## 3. Complement repair

If no column contains 28 positions, leave both paths unchanged. Otherwise,
let d be the unique component with a crowded column. Both paths must have
length 14 and their two segments on that shared line must be disjoint.

Replace the first path by the signed complement

    S'_1 = S_1 - 1_(L_1).                           (4)

This is minus the indicator of the 15 positions of L_1 outside S_1.
Since a whole e-line is invariant under translation by e,

    D_e 1_(L_1)=0,
    D_e S'_1=D_e S_1,
    D_uD_vD_w S'_1=D_uD_vD_w S_1                    (5)

over the integers. Thus the endpoint dipole and the mixed cube are exactly
unchanged, as are the incidence spaces and all coefficients lambda_dj.

In the crowded component, (4) replaces the first 14-point segment on the
shared line by its 15-point complement. Because the second 14-point
segment was disjoint from the first, it is contained in this complement.
The union support on that column now has 15 positions. Every other column
of this component is unshared, and has length at most 15.

The other nonexceptional component had no shared columns by Section 2.
Its four columns remain unshared after (4); their lengths are 14 or 15.
For the exceptional component d=e, the correction is unchanged by (5)
and has at most four individual points. In either derivative grid these
occupy at most four columns with at most four positions each.

Consequently, after this optional single complement operation, every
component meets the column-normalization hypotheses: at most four columns,
at most 27 positions per column, and 29>2*4. When a repair was necessary,
all nonexceptional components actually have at most 15 positions per column.

This argument is independent of the cube supports being disjoint and of
any incidence-rank assumption. It is a support-normalization result for
any two same-direction dipoles on distinct cycles.

## 4. Strengthened exclusion criterion

For equal exceptional directions, use S'_1 from (4), when necessary, in
(1). Its mixed cube is the same C_1, so the mixed derivative of J_d still
vanishes exactly. All total sums and d-cycle sums of the correction remain
zero, because it is still a d-derivative. The preceding section therefore
makes J_d binary with one point per d-cycle.

The original local-motion calculation now applies with the revised
potential. If h_j are the constant vectors determined by the projection
relations, it reads

    z=-2J+R(h_1 S'_1+h_2 S_2).

It would place the straight-seam error -2J in the carry coset, contradicting
the established nonzero period of every seam orientation. For distinct
exceptional directions, retain the original short paths and the crossing
normalization lemma.

**Corollary.** In the p=29 mixed triple, any two-double-cycle pattern whose
two mixed cubes have disjoint supports and whose coupled homogeneous
incidence kernel has dimension four is impossible for a section.

The qualification about kernel dimension is essential: the previous
packet's explicit disjoint-cube example has dimension six. The corollary
does not exclude that geometry by itself and does not establish a 90-error
or 450-error lower bound.

## 5. Exact controls and complete reduced universe

Run

    python3 verify_path_complement.py

Only potential failures of the column bound need enumerating. If either
short path has length at most 13, the union in any shared column has size
at most 27. If there is no shared column, every column has at most 14.
Hence the only remaining universe has k_1=k_2=14 and line difference
+/-d for one component d!=e.

Translate P_1 to zero. Every geometry in that universe has the unique form

    e in {u,v,w},  d in {u,v,w}\{e},
    sign in {-1,+1},  h in F_29,
    P_2=sign*d+h*e.

There are 3*2*2*29=348 cases. Unequal quotient projections ensure that the
chosen shared component and sign do not omit or duplicate a line class.
The verifier checks all 348 cases using integer sparse arrays. Exactly
24 have a 28-position crowded column, corresponding to h=14 or h=15 in
each of the twelve direction/sign configurations. In each of those cases
it applies (4), checks the two equalities (5), and checks every component's
column bounds and every d-cycle correction sum. It also checks the 324
cases requiring no repair.

The JSON records all 348 before/after support certificates. The script
uses only the Python standard library. These controls confirm the support
argument; the proof does not infer anything from solver infeasibility.

## Independent audit and reproduction

Read the [independent proof review](review/REVIEW.md). Use Python 3 and its standard library. Run the submitted `verify_path_complement.py` with assertions enabled; `-O` is not a verification mode for that checker. From this package directory, run the independent implementation in both effective modes:

    python3 -B review/independent_controls.py
    python3 -B -O review/independent_controls.py

Compare stdout and the generated `review/independent-results.json` with the corresponding normal/optimized receipts. Besides the 348 reduced cases, the independent implementation checks all 477,456 ordered short-path geometries up to translation, with exactly 24 complement repairs. `SOURCE_PROVENANCE.json` retains reviewed and published hashes. No mathematical proof or test logic changed.
