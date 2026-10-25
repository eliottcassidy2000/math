# Sixteen Fano sign quotients as an affine four-bit carrier

2026-10-25. **PROVED**, conditional only on the explicitly cited established
root-surface theorem for the inherited surface; the elementary counts and
local quotient calculations below are proved here. **FINITE-EXACT** controls
are supplementary. Geometry is over an algebraically closed field of
characteristic zero. No literature-priority claim or classification up to
isomorphism of surfaces is made independently here.

**Incoming recovery, commit c24c123d86 (2026-10-25).**
[Session60's seven-line theorem](../../session60_six_seven_geometry_20261025/seven-lines/RESULT.md)
and its [independent geometric audit](../../session60_six_seven_geometry_20261025/audits/campedelli/AUDIT.md)
arrived while this packet's independent derivation was running. They already
prove the 16/12/2 count, the singularity profiles, the 4+12 geometric orbits,
and the stronger Campedelli/properly-elliptic/Enriques classification.
Those parts below are independent corroboration and explanation, not a new
frontier. The additional results here are the explicit affine four-bit
W-action in Section 5 and the no-6+4+6 subgroup-orbit obstruction.

## 1. Inheritance and exact object

Use the seven branch forms, in cyclic hexagon order,

    L0=E, L1=E+x, L2=E+x+z, L3=E+z,
    L4=E-x, L5=E-x-z, L6=E-z.

Their root surface X in P6 is the inverse image of their parameter plane
under coordinate squaring. The established
[Session59 root-surface theorem T1](../../session59_f1_geometry_and_information_20261025/toric-shared-roots/TORIFICATION_OBSTRUCTION.md)
proves normality, 24 A1 nodes, canonical sheaf O(1), K_X^2=16,
p_g=7, and q=0. Its three triple intersections have label triples

    T1={0,1,4}, T2={0,2,5}, T3={0,3,6}.                 (1)

No other triple or higher crossing occurs. Each triple has eight nodes
above it. The projective sign group is G=F2^7/<1111111>, of order 64.

Closest mechanism: the independent-root sign cover and its local A1 model.
Canonical hostile: a sign subgroup can be free on every coordinate line
and double crossing but acquire stabilizers at a triple crossing. Corrected
near miss: K_X^2/8 is the canonical square of the singular quotient, not
automatically of its resolution. Least-used sidecar: the labelled incidence
matroid of actual branch lines, retained alongside the binary code.

## 2. Why precisely thirty, and why the crossing histogram is 16,12,2

A Fano plane F on the seven fixed labels has seven three-element lines.
For each such line A take the sign mask of its complement A^c, of weight
four. Together with zero these form a three-dimensional binary simplex
code C_F. Its image H_F in G is a subgroup of order eight. Equivalently,
its seven nonidentity classes have unique weight-three representatives,
precisely the Fano lines. The even-weight representatives give a canonical
lift of H_F back to C_F: complementing a mask changes parity since 7 is odd.

Every labelling of the seven nonzero vectors of F2^3 gives a Fano plane;
two labellings agree precisely up to GL(3,2), of order 7*6*4=168. Thus
there are 7!/168=30 distinct labelled planes and distinct H_F. This is the
specified simplex-code family, not all order-eight subgroups of G.

For any fixed point 0, the three Fano lines through it partition the other
six labels into a perfect matching M. Conversely each matching has exactly
two Fano extensions. To see the latter, write its pairs {a_i,b_i}, i=1,2,3.
The other four lines must choose one point from each pair. Encode a choice
by a bit triple. Every two-coordinate projection must contain all four
possibilities exactly once, forcing the four triples to have either even
or odd parity. Both choices work. There are 15 matchings and hence 30 planes.

Let M0={{1,4},{2,5},{3,6}} be the matching in (1). The number b(F) of
forbidden triple crossings belonging to F is |M intersect M0|. Two shared
matching edges force the third; b=2 is impossible. The matching census is:

* b=3: one matching, giving two Fano planes;
* b=1: choose the shared edge in three ways, then one of two other
  matchings of the remaining four labels: six matchings, twelve planes;
* b=0: the remaining eight matchings, sixteen planes.

Equivalently each forbidden triple belongs to six planes, while every pair
of forbidden triples belongs to the same two planes. Inclusion-exclusion
gives 30-(3*6-3*2+2)=16. Thus the three tests are strongly dependent; treating
them as independent loses the reason for the sixteen.

## 3. All thirty quotient singularity profiles

Put Z_F=X/H_F, and let S_F be its minimal resolution of singularities.
This means minimal among resolutions of Z_F, not a minimal model in the
birational classification of smooth surfaces.

At a point of X, a projective sign fixes it exactly when the signs on all
nonzero root coordinates agree. Choose that common sign to be positive;
the remaining sign support is contained in the vanishing coordinates.
Those coordinates have size at most three. H_F has no nonzero representative
of weight one or two, so its action is free except possibly at the nodes.
At a triple A it has stabilizer of order two iff A is a Fano line, and
otherwise trivial stabilizer. In the order-two case its generator changes
all three vanishing root signs simultaneously.

The local A1 equation is w^2=alpha*u^2+beta*v^2. Over the present field it
is the quadratic cone, equivalently A2/{+1,-1}. Simultaneous negation of
its three linear coordinates lifts to scalar multiplication by i on A2.
Its quotient is therefore the cyclic singularity 1/4(1,1). Consequently:

| b | Labelled subgroups | A1 points on Z_F | 1/4(1,1) points on Z_F | K_(S_F)^2 | p_g(S_F),q(S_F) |
|---:|---:|---:|---:|---:|---:|
| 0 | 16 | 3 | 0 | 2 | 0,0 |
| 1 | 12 | 2 | 2 | 0 | 0,0 |
| 3 | 2 | 0 | 6 | -4 | 0,0 |

Indeed eight nodes modulo H_F give one orbit if the action is free and
two if the stabilizer has order two. This proves all singularity counts;
no new singularity arises on the smooth free locus.

The finite map X->Z_F is unramified in codimension one. Hence
K_X=pi^*K_Z and K_Z^2=16/8=2. An A1 resolution is crepant. The minimal
resolution of 1/4(1,1) has one exceptional rational curve D of square -4:
one may view the singularity as the degree-four Veronese cone over P1,
whose resolution is the zero section of O_P1(-4). Write
K_S=f^*K_Z+aD. Adjunction gives K_S.D=-2-D^2=2, so a=-1/2.
Each such point therefore changes K^2 by a^2 D^2=-1. There are 2b of
them, proving K_S^2=2-2b.

For clarity, p_g=q=0 does not assume free action on the resolution.
Use the canonical even-weight lift C_F. Its four quadric equations are
invariant, and the determinant of each sign matrix is +1. Adjunction's
residue linearization on H^0(K_X)=<x0,...,x6> is therefore exactly the
coordinate sign representation. Every coordinate takes value -1 on some
weight-four codeword, so no coordinate character is trivial and this
seven-dimensional representation has no invariants. The seven characters
are in fact all the nontrivial characters of H_F. Thus H^2(O_X)^H=0 by
Serre duality. Finite quotient and exact averaging in characteristic zero
give H^i(O_Z)=H^i(O_X)^H. In particular H^1(O_Z)=0 and H^2(O_Z)=0.
The A1 and 1/4(1,1) singularities are rational: their local resolutions have
single exceptional P1 with conormal O(2) or O(4), whose successive powers
have H^1=0; formal functions gives R^1 f_*O_S=0. Hence the same statements
hold on S.

For b=0, the lifted H_F action on the smooth resolution of X is free:
there is no stabilizer at a node and the exceptional curves are permuted
freely. K_S is nef and big, being the pullback of ample K_Z, so S has no
(-1)-curves and is a numerical Campedelli surface. Its singular canonical
model Z has exactly three A1 points. The other fourteen quotients have
the displayed noncanonical singularities. The incoming Session60 proof
settles their minimal models, as reviewed in Section 6: the twelve b=1
resolutions are already minimal properly elliptic surfaces; each b=3
resolution contracts four disjoint (-1)-curves to an Enriques surface.
These conclusions require the retained central-curve/pencil data;
K_S^2=0 alone does not make S a K3, and K_S^2=-4 alone does not determine
rationality.

## 4. Which symmetry really splits sixteen?

The full permutation group preserving the crossing star (1) is
W=C2 wr S3, of order 48: permute the three opposite pairs and swap their
endpoints. It acts transitively on the eight matchings disjoint from M0.
For one such matching M, the union M0 union M is a six-cycle. The
stabilizer of both coloured matchings has order six. It exchanges the
two parity extensions: for the representative M={{1,2},{3,4},{5,6}},
the permutation (1 2)(3 6)(4 5) preserves M and M0 and exchanges those
extensions, as is seen by its action on the even/odd transversal triples.
Transitivity on matchings transports this check. Thus W is transitive on all sixteen good
planes, and every stabilizer has order three.

The actual hexagon symmetries of the branch forms supply a natural
dihedral subgroup D6 of order twelve. In cyclic labels 1,...,6 its
generators are rotation i->i+1 and reflection i->7-i (fixing 0).
They are induced by linear transformations of x,z. For example rotation
is x->x+z,z->-x. Reflection x->-z,z->-x reverses the displayed hexagon.

A good matching consists either of three adjacent hexagon edges, or of
one adjacent edge and two distance-two edges. There are two matchings of
the first kind and six of the second. The former give four planes, the
latter twelve. D6 is transitive within each collection: it is transitive
on the two/six matchings, and the reflection fixing a representative
matching exchanges its two parity extensions. Representatives and the
exchange permutations are recorded by the checker. Therefore the exact
geometric orbit decomposition is **4+12**.

There is a stronger obstruction to reading 6+4+6 as a symmetry-orbit
partition here. For any subgroup K of W, a point stabilizer has size 1 or
3, since it is contained in a W-stabilizer of order three. Every K-orbit
therefore has size |K| or |K|/3. Orbit sizes four and six cannot both occur:
they would require |K| in {4,12} and in {6,18}, respectively. No subgroup
preserving these three crossings has orbit partition 6+4+6 on the sixteen.

This does not prohibit an externally chosen partition into those sizes.
It rules out the claimed intrinsic mechanism from the specified symmetry.
The recovered norm-graph 6,4,6 consists instead of dimensions of three
line-incidence spaces on sixteen support points; its maps and equations
must be transported before any identification is available.

## 5. A lossless four-bit encoding, with the action retained

The sixteen are not merely equipotent with F2^4. There is an explicit
affine action model for every permutation preserving the crossing star.
Index the opposite pairs P_i={i,i+3}, i=1,2,3. Write
W=N semidirect S3, with N=F2^3 swapping the members of each P_i and S3
permuting the pair indices while preserving the selected first members.
For v in N let s_v be its sign-swap permutation, and let

    rho=(1 2 3)(4 5 6),    tau=(1 2)(4 5).

Choose the explicitly labelled Fano plane

    F0={015,026,034,123,146,245,356}.

It avoids (1), is fixed by rho, and its W-stabilizer is exactly <rho>:
the previous transitivity/stabilizer argument gives order three. Define

    Psi(v,epsilon)=s_v tau^epsilon F0,
        v in F2^3, epsilon in F2.                         (2)

These sixteen planes are distinct and exhaust the good planes, because
they are representatives for W/<rho>. More strongly, for g=(n,sigma) in W,

    g : (v,epsilon) -> (n + sigma(v), epsilon + parity(sigma)).  (3)

Here sigma(v) permutes three coordinates, and all additions are in F2.
To prove (3), conjugate s_v past sigma, then note that <rho>=A3 is the
kernel of parity in S3. Thus sigma*tau^epsilon and
tau^(epsilon+parity(sigma)) have the same right coset modulo <rho>.
Every group element acts affinely on four bits. The map (2) reconstructs
the entire Fano line system and simplex sign subgroup; it loses no code
or port information inside the sixteen-good-code universe.

This gives a specific meaning to “three bits plus one extra parity”:
the first three record independent exchanges within opposite branch
pairs; the fourth records the orientation of the permutation of those
three pairs. It does not assert a preferred origin or a natural vector
space on the sixteen before a labelled seed and pair ordering are chosen.
All 48*16 affine transport identities are checked independently against
literal permutations of seven-point Fano line sets.

The actual branch-geometry D6 action is obtained by restricting (3) to
the twelve realized linear symmetries. Its orbit sizes remain 4+12.
The larger W is an incidence symmetry of the three crossing tests;
most of its elements are not automorphisms of the complex branch
arrangement. Keeping this distinction is essential.

## 6. Review of the incoming minimal-model proof

The fresh Session60 argument adds the coordinate previously missing from
the numerical table. Blow up the three collinear triple centers in P2,
writing H for the pulled-back line and E_j for the exceptional curves.
Normalize the rank-three quotient cover over this base. The strict seven
lines have total class 7H-3*sum(E_j). An exceptional E_j is branched
exactly when the three inertia labels have nonzero sum. At every new
branch crossing the two labels are distinct nonzero vectors, hence
independent, so the resulting cover g:S->B is smooth. Hurwitz gives the
exact linear equivalence

    2K_S ~ g^*(H - sum_(bad j) E_j).                       (4)

This reproduces the resolution K-squares above, but retains substantially
more information. If b=1, H-E_bad is the free pencil of lines through
that center. A generic line has four branch points labelled by the
complement of a Fano line. These four labels span F2^3, so its degree-eight
cover is connected, and Riemann--Hurwitz gives genus one. Equation (4)
makes K nef and supplies h^0(2mK)=m+1. This proves minimality and Kodaira
dimension one, rather than inferring them from K^2=0.

If b=3, H-sum(E_j) is the strict central line B0, with B0^2=-2.
Its residual cover splits into four disjoint rational components R_i,
each mapping isomorphically, and g^*B0=2*sum(R_i). Projection formula
gives R_i^2=-1. Contract these four curves; the canonical blowup formula
and (4) give 2K_T~0 on the resulting smooth surface T. Its unchanged
p_g=q=0 force K_T nontrivial, so T is Enriques. The six exceptional
(-4)-curves meet the four central components according to K4; after
contraction they become (-2)-curves. These multiplicities and the exact
linear equivalence check the mechanism independently of the original
singular-discrepancy calculation.

**Review conclusion: PASS for the incoming scoped trichotomy over C.**
The finite checker here does not reconstruct that geometric proof by
enumeration. Its independent corroboration concerns the codes, local
stabilizers and canonical characters. The displayed geometric review
checks the relevant extra deduction before inheriting its status.

## 7. Connection ledger and next target

The useful bridge is **kernel versus local stabilizer**, not an equality
between the numbers six, seven, and sixteen:

    branch incidence triples + sign code -> stabilizer intersections
        -> quotient singularities -> correction to K^2.

Preserved: labelled roots, code subgroup, crossing incidence and exact
local action. Lost when keeping only the number 16: which crossings can
be regenerated, subgroup orbit data, singularity type and discrepancies.
Needed sidecar: the actual incidence triples and their images in G/H.
Cheapest hostile: change one ordinary crossing into a code-dependent
triple; two 1/4 points replace one A1 point, and resolution K^2 drops two.

This is thematically compatible with movable-cut repairs: a quotient
certificate must check its interaction with the actual local forbidden
supports. It is **not a map to tournament blackness**, nor a connectivity,
Collatz or Riemann-hypothesis argument. The geometric minimal-model question
has already been answered by the incoming trichotomy. A precise combinatorial generalization
is to classify crossing-avoidance codes for other seven-line incidence
hypergraphs; the present three-star count need not persist.

## 8. Sources and reproduction

Classical context: Mendes Lopes--Pardini--Reid,
[Campedelli surfaces with fundamental group of order 8](https://arxiv.org/abs/0805.0006),
identifies degree-eight etale covers with intersections of four quadrics.
Alexeev--Pardini,
[Explicit compactifications, v4, Sections 2.2 and 3](https://arxiv.org/html/0901.4431v4),
supplies the standard abelian-cover and labelled-seven-line conventions.
The arithmetic hourglass specialization, thirty-code count and displayed
local quotient/discrepancy computations are derived above rather than
claimed as consequences of those citations alone.

Run `python check.py` and `python -O check.py` from this folder. The checker
uses two independent finite descriptions: all 7! label permutations of
the binary plane, and direct exact-cover enumeration of Steiner triples.
It checks each code, crossing stabilizer, quotient-node orbit, canonical
character and group orbit, and includes controls distinguishing the bad
cases. It proves no surface-classification statement by enumeration.
