# Exact physical-state hulls and structural graft images

**Status: PROVED, independently audited with FINITE-EXACT controls.** All-order mathematical claims
are stated separately from finite controls. Chan–Pak all-odd completeness
remains OPEN. No literature-priority claim is made.

## 0. Inheritance, conventions, and scope

A tournament is finite and nonempty unless explicitly deleted to the empty
tournament. Directed paths are ordered vertex words, and H(empty)=1.
For a rooted tournament (A,a), use the Session225 state

    x(A,a)=(D,E,U,d),

where D=H(A), E counts Hamilton paths ending at a, d=H(A−a), and U counts
ordered spanning pairs (P,R), with P possibly empty and R nonempty and
starting at a. There is one empty path in an allowed empty bin.

The inherited graft G_B(A,a) retains A and nonempty rooted B, directs
B→a and A−a→B, and marks B's root r. Its four-state matrix K_B is proved in
`session225_graft_recovery_20261108/PROOF.md`, Sections 1–3. Session226
`session226_structural_continuation_20261108/graft/PROOF.md`, Sections 1–3,
proves det K_B=seΔ, Δ=qQ−Hb, for |B|≥2, and provides its exact rational
inverse. Rational state inversion is not a graph inverse.

The inequalities E≤d and U−d≥D−E are inherited endpoint-defect facts, not
new claims: see `adaptive_metagraph_reformulation_20261103/ADAPTIVE.md`,
Section 3, `adaptive_state_induction_20261103/BLOCK_SWITCH.md`, B2, and
Session226 `root_descent/PROOF.md`, Section 1. The positive coordinates in
those sources do not claim exactness of the closed cone of all physical
states. The present refinement is exactness, its extreme-ray witnesses,
the exact affine integer lattice hull, and their transport to graft images.
The cycle-with-transitive-block count H=1+2^m is also inherited; for example
`nonzero_preparation_20261103/algebra/ONE_EXIT_AND_PERSISTENCE.md`, Section 5.
Our use derives its deletion-ratio signs rather than claiming a new family.
The independent H37 graph-image hostile in Section 6 is not the older
H37 arbitrary-root-switching hostile in ADAPTIVE.md, Section 4.

Closest mechanism: four-state invertible grafts. Canonical hostile:
state equality can lose graph modules. Corrected near miss: an integral
point in the exact cone/lattice hull is still not a physical state.
Least-used sidecars: a cyclic transitive block as an extremal-ray sequence,
and the named homogeneous gadget module. Concept board: physical cone,
affine lattice, actual graph module, deletion-ratio signs, scalar coverage.

## 1. Exact closed cone of all rooted tournament states

Let S be the set of all physical four-state vectors, over all finite orders,
and let cone(S) mean finite nonnegative real linear combinations. Then

    closure(cone(S)) = C,
    C = { (D,E,U,d): E≥0, d−E≥0, D−d≥0,
                       U−D−d+E≥0 }.                    (1)

This is a simplicial cone. Its four extreme rays, in coordinate order
(D,E,U,d), are

    r1=(1,1,1,1), r2=(1,0,2,1),
    r3=(1,0,1,0), r4=(0,0,1,0).                      (2)

**Necessity.** E≥0 is immediate. Deleting a terminal root injects those E
paths into the d deletion paths. Every Hamilton path of A−a has an insertion
position for a; choose a fixed deterministic such position. Deleting a from
the resulting path recovers the original, so d≤D. Finally split every
Hamilton path not ending at a immediately before a. This injects the D−E
paths into the U−d covers whose rooted path is nontrivial. Thus all four
inequalities hold, also for a singleton source.

For completeness, insertion exists by the usual binary-transition argument:
insert first if a beats the first vertex, last if the last vertex beats a;
otherwise an in-neighbor followed by an out-neighbor supplies an interior
position. This also proves Hamilton-path existence inductively.

**Sufficiency for the closed hull.** Every x∈C has the unique expansion

    x=E r1 +(d−E)r2 +(D−d)r3 +(U−D−d+E)r4.          (3)

The singleton realizes r1. The source-rooted directed edge realizes r2.
For m≥1, let B_m consist of the transitive m-block X and singletons u,v,
with X→u→v→X. Mark u. Section 3 proves

    x(B_m,u)=(2^m+1,1,2^m+1,1).

Dividing by 2^m+1 tends to r3. The source-rooted transitive tournament on
n≥2 vertices has state (1,0,2^{n−1},1); division by 2^{n−1} tends to r4.
Thus all rays belong to the closed hull, proving (1). Linear independence
in (3) proves the extreme-ray assertion.

**Consequence.** Every homogeneous linear inequality valid on all physical
states is a nonnegative linear combination of the four displayed facets.
This is exactly the dual-cone statement from the invertible coordinates
in (3). It concerns homogeneous linear inequalities over all orders, not
nonlinear, fixed-order, affine, or discrete constraints. The closure is
essential: every nonzero finite nonnegative combination of physical states
has d>0, whereas the nonzero rays r3 and r4 have d=0.

## 2. Exact affine integer lattice hull, and a strict realizability gap

The smallest affine integer lattice containing S is

    L={ (D,E,U,d)∈Z^4 : D≡d≡1 (mod 2), U≡E (mod 2)}. (4)

Its difference lattice has index 8 in Z^4.

**Necessity.** Rédei parity gives D and d odd, with H(empty)=1. The inherited
singleton-graft formula gives H(J(A,a))=D−E+U, also odd. Subtracting D proves
U≡E (mod 2). Rédei parity is a classical cited dependency, catalogued at
`01-canon/theorems/THM-001-redei.md`.

**Generation.** The following physical states suffice:

    s0=(1,1,1,1)       rooted singleton;
    s1=(1,0,2,1)       source-rooted edge;
    s2=(1,0,4,1)       source-rooted transitive triple;
    s3=(3,1,3,1)       rooted directed triangle;
    s4=(5,1,9,3)       B_2 rooted at the first X vertex.

Differences and integer combinations give

    s1−s0=(0,−1,1,0),
    s2−s1=(0,0,2,0),
    (s3−s0)−(s2−s1)=(2,0,0,0),
    (s4−s0)−2(2,0,0,0)−4(0,0,2,0)=(0,0,0,2).

These generate precisely the difference lattice in (4), of determinant
magnitude 8. Hence no additional affine subgroup congruence holds on every
physical state. This is a statement about the affine lattice hull; nonlinear
residue restrictions are not ruled out.

The vector z=(1,0,6,1) lies in C∩L but is not physical. Indeed Hamilton-path
count is monotone under induced vertex addition by the insertion injection.
A tournament with a directed triangle therefore has at least three Hamilton
paths. A tournament without a directed triangle is transitive (a shortest
cycle would otherwise have length three). Thus H=1 forces transitivity.
In a transitive rooted tournament every root-starting path is an arbitrary
subset of its out-neighbors in their forced order; the complementary path
is uniquely ordered. Therefore U=2^{outdegree(a)}, never 6.
This proves a genuine integer realizability gap even after both exact hulls.

## 3. Full-rank gadgets of both signs at every order at least four

For B_m=X→u→v→X above, m≥1, let H_m=H(B_m), Q_m be its ordered two-bin cover
count with empty bins allowed. Then

    H_m=2^m+1,    Q_m=2·3^m+2^{m+2}−2.               (5)

**Hamilton paths.** A path in which u,v are consecutive in that order is
obtained by choosing which X vertices lie before u; all remaining X lie
after v. Each chosen subsequence is in its transitive order, giving 2^m
paths, including the empty extreme choices. The only path with v before u
is v, all X in order, u. If u precedes v they must be consecutive: u has
only the out-neighbor v. These cases prove H_m.

**Two bins.** If u,v lie in the same bin, choose the bin (two choices) and
its X subset, of size k. Its Hamilton count is 2^k+1 for k≥1 and 1 for k=0;
the other bin is transitive. The total is 2(3^m+2^m−1). If u,v lie in distinct
bins, every assignment of X gives two transitive paths, contributing 2^{m+1}.
This proves Q_m.

For root u, s=e=1, q=1 and b=2^{m+1}, since deletion leaves a transitive
(m+1)-tournament. Also U=H_m: the rooted bin is either {u}, contributing 1,
or begins u,v and takes any subset of X, contributing 2^m. The same endpoint
counts s=e=1 hold for root v. Consequently, for either singleton root,

    Δ=2(3^m+2^m−1−4^m).                             (6)

This is zero for m=1 and negative for every m≥2. One elementary proof is
to set f_m=4^m−3^m−2^m+1: f_1=0 and
f_(m+1)=4f_m+3^m+2·2^m−3>0 for m≥1.

For a root at the i-th X vertex, indexed i=0,...,m−1, with m≥2, deletion
leaves B_(m−1). Hence

    q=2^{m−1}+1, b=2·3^{m−1}+2^{m+1}−2,
    Δ=(2^m+4)3^{m−1}+3·2^m>0.                       (7)

The endpoint counts are s=2^{m−1−i}, e=2^i. To start at the root in a path
with u before v, its initial X run can take an arbitrary subset of the
later X vertices, and its suffix after v contains every other X vertex;
this gives s. The reverse argument gives e.

All B_m are strong: the three nonempty blocks form a cyclic directed quotient.
The inherited determinant theorem therefore implies invertibility for every
mark when m≥2. Singleton-root determinant is (6); an X-root determinant is
2^{m−1} times (7). Thus each order n≥4 has strong full-rank gadgets with
both signs of Δ. The order-three boundary is singular. No infinite strong
singular family and no scalar-coverage conclusion is claimed.

## 4. Exact image cone and affine lattice, not a physical inverse

Fix any invertible gadget B. Its physical output states are S_B=K_B S.
Because K_B is a linear homeomorphism,

    closure(cone(S_B))=K_B C,
    affine_integer_lattice(S_B)=K_B L.              (8)

Therefore the four facet tests on K_B^(-1)y are the complete homogeneous
linear tests for the graft's physical-state image; inverse membership in L
is the exact affine-lattice-hull test. Their conjunction is necessary, not
sufficient, for y to lie in S_B. The explicit z from Section 2 gives K_B z
in both image hulls but outside S_B by injectivity. Its output vector is not
claimed to be a physical tournament state.

The difference lattice in (8) has index 8|det K_B| in Z^4. More generally,
for a finite prescribed chain P=K_(B_k)...K_(B_1) of invertible grafts, the
same conclusions hold with P, and the index is 8 product_j |det K_(B_j)|.
Every physical source has a literal graph realizing this chain; this is why
these are equalities of hulls, rather than merely formal matrix containments.

For the inherited B4=[8,1,3,6], write y=(A,B,C,F) to avoid primed notation.
Its exact inverse is

    D=C+5B−(A+13F)/2,
    E=C+4B−(5A+11F)/4,
    U=F−B,
    d=(A+3F−8B)/8.                                 (9)

Thus K_B4 Z^4 consists exactly of integer y with A+3F≡0 (mod 8).
Sufficiency follows by substituting this congruence into every denominator
in (9), and necessity follows from the last line. Its lattice index is 8,
consistent with det K_B4=−8. The affine physical-lattice hull has index 64:
use (9) together with D,d odd and U≡E modulo 2. The image-cone test is the
four inequalities in (1) applied to (9).

### Exact scalar affine congruence and an infinite factor obstruction

For any fixed nonempty gadget, even a singular one, the affine integer
lattice generated by its actual scalar outputs is exactly

    H+g Z,    g=gcd(2H,Q).                           (10)

Indeed the first output coordinate at s0 is H. Apply the first matrix row
to the difference generators (2,0,0,0), (0,−1,1,0), (0,0,2,0), (0,0,0,2)
from Section 2. The respective increments are

    2H, Q, 2(Q−H), 2(2H−Q),

whose gcd is g. Integer affine hull commutes with projection, proving (10).
This is the strongest affine scalar congruence shared by all actual source
states, not a statement that all numbers in the progression occur.

For every odd m≥3, B_m is a full-rank gadget at every root, yet 3 divides
both H_m and Q_m. Consequently every graft using it has Hamilton count
divisible by 3. Full rank and unbounded gadget orders therefore do not by
themselves remove even this elementary scalar arithmetic obstruction.
For B4, g=gcd(10,32)=2, so this affine congruence imposes only oddness;
the full physical scalar image is still not identified.

## 5. A complete graph-level image certificate

A nonempty subset M of a tournament T is a module if each vertex outside M
uniformly beats all of M or loses to all of M. Fix a nonempty rooted gadget
(B,r), and let (T,t) be a rooted tournament. Then T is isomorphic to an
actual G_B(A,a), preserving the marked root, if and only if there is a subset
M satisfying all of the following:

1. t∈M and the rooted induced tournament (T[M],t) is isomorphic to (B,r).
2. M is a module of T.
3. Exactly one vertex a outside M loses to M; every other outside vertex
   beats M. In particular M is proper and its complement is nonempty.

**Proof.** The actual gadget copy has precisely these properties by the
cross-edge definition. Conversely take A=T−M, rooted at the unique a.
The induced graphs and all cross-arcs are then exactly the graft definition.
This is an if-and-only-if graph certificate, valid even for singleton or
nonstrong gadgets. It reconstructs the source once M is supplied; uniqueness
of M is not asserted. Fixed-gadget recognition can simply inspect all
root-containing subsets of size |B| and test these three finite conditions.

For a prescribed gadget chain, recursively apply this certificate in reverse
order, marking the uniquely specified old root each time and retaining each
selected module. A complete successful module sequence is equivalent to an
actual graph derivation from its remaining nonempty rooted source. If gadget
matrices are invertible the recovered source states also agree with the
successive rational inverses, but those numerical inverses cannot replace
the module sequence.

## 6. A physical numerical inverse does not recognize a graph image

Let B4=[8,1,3,6], root 0, as inherited. Consider the rooted six-tournaments

    T_good=[18,32,3,7,14,29], root 0,
    T_bad =[18,40,3,5,14,29], root 0.

Both are strong and have the identical physical state

    (D,E,U,d)=(37,7,43,9).

T_good has the module M={0,2,3,4}, with unique external loss vertex a=1.
Its induced rooted module is isomorphic to B4 via gadget vertices
(0,1,2,3)→(0,2,3,4). The remaining rooted source is the edge 1→5, rooted
at 1. Thus the common numerical inverse is (1,0,2,1), genuinely physical.
Nevertheless T_bad has no four-vertex module containing 0, and hence no
rooted B4 graph-image certificate. All ten four-subsets containing 0 fail
the module condition, before the rooted-isomorphism test is even needed.

This is FINITE-EXACT evidence for the explicit pair; the general separation
statement follows immediately from the pair. It says that even equality to
a physical graft-output state does not imply graph-level membership. It
makes no claim that a scalar count fails to occur by another construction.

## 7. Verification and boundary

The portable standard-library checker uses explicit exceptions, remains
active under python -O, and separates finite controls from the all-order
proofs. `inherited_verify.py` is an unchanged copy of Session225's literal
subset DP, permutation counters, graft constructor, and matrix utilities.
The present checker independently compares profiles with permutation counts.

The unfiltered profile universe is every labelled tournament of orders 1
through 5 and every root: 5,405 rooted presentations. Further controls are
all 5,625 pairs of labelled source/gadget tournaments through order 4 with
root 0; 75 literal two-B4 chains; all 75 rooted cyclic family examples for
m=1,...,10, including 18 independent permutation-statistics controls through
m=4; 11 transitive-ray controls; 6,561 output vectors in [-4,4]^4 for the
exact B4 lattice test; and all 75 small gadget scalar-lattice projections.

Graph recognition is independently tested on every labelled rooted-at-0
tournament of orders 5 and 6, totaling 33,792 graphs. The recognizer tests
modules; its comparison bank is separately generated by every appropriate
source and every root-preserving relabelling of the constructed graft.
The graph-image counts are 24 and 240, respectively. The explicit collision
pair is additionally counted by direct permutations and its module lists
are checked. A finite check does not establish any all-order quantifier.

Reproduce in this directory:

    python3 -I -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B -O verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.optimized.json

The two checked-in receipts are identical. Independent audit and its separate
checker are in the sibling image-audit packet supplied by the coordinator.
Neither cone/lattice exactness, graph-image recognition, nor the full-rank
families prove an unrestricted scalar coverage range or the all-odd
Chan–Pak conjecture. The latter remains OPEN.
