# Parallel-coherent permutohedra and a Grassmannian lift of Hamilton paths

Status: PROVED elementary statements below; independently audited;
FINITE-EXACT controls are separately specified. This does not prove that the
only missing odd Hamilton counts are 7 and 21. No literature-priority claim
is made.

## 0. Inheritance, target, and connection contract

The target is a natural graded carrier for the tournament Hamilton-path
count spectrum, rather than a renaming of odd integers. The closest inherited
mechanisms are the transitive-tournament permutohedron (THM-207), the classical
one-reversed-arc formula H(NT_n)=1+2^(n-2) (Moon, as attributed by THM-4748),
and the already proved gaps 7 and 21 (THM-1370). Source paths and immutable
revision are in SOURCES.json. The scalar one-arc formula and the gaps are NOT
new results of this packet. The forward-adjacency distribution is also already
present in THM-062; our face-incidence interpretation uses that observable.

The hostile is the strong four-vertex tournament whose five local minima
have feedback energies 1,2,2,2,3, so local minima are not globally optimal
rankings. The corrected near miss is ordinary convex optimization: a generic
linear functional on a convex polytope has only one sink. The underused sidecar
is the feedback energy of each Hamilton path, rather than its count alone.

Concept board: parallel edge classes / adjacent-swap sinks / all-pairs
feedback energy / face factorization / binary inversions / finite-field
row reduction. Anchor: an exact geometric carrier. Niche: which relaxations
lose the forbidden counts. Wildcard: lifting a binary path family to subspaces.

The full carrier map is a bijection, not a compression: labelled tournaments
map to parallel-coherent orientations of a fixed permutohedron, preserving
order and Hamilton paths as sinks. Its class definition is geometric, but
it retains every original arc, so this equivalence by itself is not a solution
of the spectrum problem. The additional nontrivial bridge is an explicit
energy-refined bijection for NT_n: its nontrivial Hamilton paths index binary
words, and field decorations on their inversions parametrize ALL subspaces
of a finite-dimensional vector space. Forgetting the decorations leaves the
Schubert-cell index, not the subspace itself. The cheapest checks are a
four-vertex local/global hostile and row-reduction over F_2 and F_3.

## 1. An exact geometric carrier

Let P_(n-1) be the standard permutohedron, with a vertex x(pi) for every
permutation pi of [n]. Coordinate x_v is the position of v in pi. Its graph
joins permutations differing by one adjacent interchange. An edge swapping
labels a,b is parallel to e_a-e_b. These are exactly the binomial(n,2)
parallel edge classes.

Call an orientation parallel-coherent if every edge in one such class is
oriented in the same vector direction. This definition uses the embedded
polytope and does not mention tournaments.

**Theorem 1.** Labelled tournaments on [n] biject parallel-coherent
orientations of P_(n-1). Every such orientation is acyclic. Its sinks biject
the tournament's directed Hamilton paths. The correspondence commutes with
relabeling and takes converse to reversal of all oriented edges.

**Proof.** Given T, when a->b direct every adjacent interchange toward the
ordering with a before b. Its vector direction is e_b-e_a, independent of
where the interchange occurs. Conversely one coherent direction on each
parallel class uniquely specifies the orientation of each pair a,b.

For an ordering pi define its all-pairs feedback energy

    E_T(pi) = #{i<j : pi_j -> pi_i in T}.

Interchanging adjacent a,b changes the relative order of no other pair.
Our directed interchange therefore decreases E_T by exactly one. A directed
cycle is impossible. A vertex is a sink exactly when every consecutive pair
follows a tournament arc, which is exactly the Hamilton-path condition.
Relabeling and converse have the stated effects. QED.

All positive edge weights can replace the unit costs in E_T without changing
any directed interchange or sink. A decreasing interchange then changes the
energy by the positive weight of its swapped pair. This shows that the
landscape's sinks depend on signs, not the magnitudes of pairwise penalties.
The polynomial grading below uses UNIT weights.

Under coordinate-permutation equivalence, the carrier also identifies
unlabelled tournaments with equivalence classes of these orientations.
It is not the native fixed-order tournament metagraph: the vertices here
are rankings of one fixed tournament, and its edges swap adjacent ranks.

## 2. Faces, exact product law, and incidence polynomial

A nonempty face is an ordered partition B1|...|Bk of [n]. Its vertices list
B1 first, then B2, and so on, in arbitrary internal order. Its dimension is
n-k and its graph is the Cartesian product of the permutohedron graphs of
the blocks. The whole polytope is included as the one-block face.

**Theorem 2.** For the orientation induced by T,

    sinks(B1|...|Bk) = product_j H(T[Bj]).                         (1)

Let f_T(pi) denote the number of forward ADJACENCIES of pi, not the number
of forward pairs. Then

    sum_nonempty_faces F z^dim(F) sinks(F)
      = sum_permutations pi (1+z)^f_T(pi).                       (2)

In particular,

    sum_F (-1)^dim(F) sinks(F) = H(T).                          (3)

**Proof.** A sink of a Cartesian product orientation is precisely a tuple
of sinks. Concatenating one Hamilton path from each induced block proves
(1). For (2), count sink/face incidences at a fixed vertex pi. A face
containing pi is specified by choosing which adjacent positions stay in
the same block. A choice of d positions gives dimension d, and pi is a
sink of that face exactly when every chosen adjacency is forward. Its
contribution is (1+z)^f_T(pi). Sum over pi. At z=-1 only orderings with no
forward adjacency survive; reversing them bijects Hamilton paths. QED.

Equation (2) is a geometric refinement of the already known tournament
forward-step polynomial. It distinguishes this local adjacent statistic
from E_T, which sums over all pairs. For even n, separating out the
whole-polytope term in (3) gives

    (H(T)-1)/2 = (1/4) sum_proper_nonempty_faces F
                         (-1)^dim(F) (sinks(F)-1).              (4)

Here the ordinary alternating face sum is 2 because the proper boundary
has even dimension n-2. Formula (4) is signed and does not manufacture a
positive set of size (H-1)/2 or explain its forbidden values by itself.

## 3. Two sharp boundaries: convex linearity and arbitrary acyclicity

**Theorem 3.** A parallel-coherent orientation is induced by decreasing a
generic linear functional on the standard embedded permutohedron if and
only if its tournament is transitive. Such an orientation has one sink.

**Proof.** Write the functional as sum_v c_v x_v. Genericity on all edges
means the c_v are pairwise distinct. A directed edge e_b-e_a decreases it
exactly when c_a>c_b. Therefore the tournament direction is a->b exactly
when c_a>c_b, a transitive ranking. Conversely distinct coefficients in
any transitive order induce its coherent orientation. A generic linear
functional on a polytope has one minimizing vertex and every other vertex
has a descending incident edge, so it has one sink. QED.

**Theorem 4.** If parallel coherence is dropped, the sink-count spectrum of
ALL acyclic orientations of the same permutohedron graph is exactly

    {1,2,...,n!/2} for n>=2, and {1} for n=1.                    (5)

**Proof.** An acyclic orientation of a finite nonempty graph has a sink;
its sinks are independent. Adjacent interchange changes permutation parity,
and interchanging the first two positions gives a perfect matching.
Consequently every independent set has at most n!/2 vertices.

For any 1<=k<=n!/2 choose k even permutations S. They are independent.
Orient the graph by decreasing the lexicographic potential

    (distance to S, permutation in a fixed total order).

This is acyclic. Every vertex outside S has a neighbor at smaller distance,
so is not a sink. Every vertex in S has all neighbors at distance one, hence
is a sink. This realizes exactly S and proves (5). QED.

Thus ordinary acyclic orientations already realize 7 sinks at order 4
and 21 sinks at order 5; these are the smallest possible orders for those
sink counts by (5). The verifier supplies deterministic orientations and
identifies the parallel classes in which coherence fails. The forbidden
tournament counts do not arise from the permutohedron graph or acyclicity
alone. Nor can a convex-linear model accommodate counts greater than one.
The exact carrier is strictly between those two overly coarse models.

This sharp boundary localizes a necessary structural constraint; it is not
an independent proof of the inherited global exclusions 7 and 21.

## 4. A Gaussian-binomial energy refinement of the one-arc family

Let NT_n, n>=3, have vertices 0,1,...,n-1 and all arcs i->j for i<j except
that n-1->0 is reversed. Set m=n-2, with interior vertices 1,...,m.

The path P0=(0,1,...,n-1) has E(P0)=1. For each subset A of the interior,
put B=[m] minus A and define

    P_A = (A in increasing order, n-1, 0, B in increasing order).

**Theorem 5.** The Hamilton paths are exactly P0 and the 2^m paths P_A.
Moreover, with inv(A,B)=#{(a,b):a in A,b in B,a>b},

    E(P_A)=m+inv(A,B).                                         (6)

Writing [m choose r]_q for the Gaussian-binomial polynomial,

    Z_n(q) := sum_Hamilton_paths P q^E(P)
      = q + q^m sum_(r=0)^m [m choose r]_q.                     (7)

**Proof.** Every descending adjacent pair in a Hamilton path must be the
only backward arc n-1->0. A path with no descent is P0. A path using that
descent has exactly two increasing runs and is uniquely P_A. All choices
of A work because no remaining increasing adjacency is the reversed pair.

For P_A the wrong-way pairs consist of: 0 versus A (|A| pairs), B versus
n-1 (|B| pairs), and exactly the cross pairs a in A,b in B with a>b.
All other pairs agree with the tournament. This gives (6).

For fixed r=|A|, the polynomial sum_A q^inv(A,B) is [m choose r]_q:
split according as m is outside or inside A to obtain the recurrence

    G(m,r)=G(m-1,r)+q^(m-r)G(m-1,r-1),

with G(0,0)=1 and zero out of range. These are precisely the Gaussian-
binomial recurrence and initial values. Summing (6) proves (7). QED.

At q=1 we recover the CLASSICAL H(NT_n)=1+2^m, attributed by THM-4748 to
Moon. The path/Boolean-set bijection is the familiar increasing-run proof;
the retained all-pairs energy is what enables the next, richer map.

**Corollary 5.1 (local search hostile).** For n>=4, P0 is the unique globally
minimum feedback ranking and has energy 1. Every other Hamilton path is a
strict adjacent-swap local minimum of energy at least n-2. The largest such
local energy is exactly

    m+floor(m^2/4) = floor(n^2/4)-1.                             (8)

**Proof.** A global minimizer must be a local sink; (6) excludes the other
sinks when m>=2. For |A|=r one has inv(A,B)<=r(m-r), with equality when A
is the set of the largest r interior labels. Optimize over r. QED.

Thus adjacent-swap local optimality can have a quadratic approximation ratio
even when the global optimum is one. With arbitrary positive weights the
ratio can be unbounded already at n=4: choose P_A with A={2}, B={1}; it
violates 0->2, 1->2 and 1->3, so give one such arc weight M and all others
weight 1. Its energy is at least M, it remains a local sink, while P0 still
has energy 1. This uses the weighted observation after Theorem 1.

The smallest unit-weight local/global mismatch is order 4. The lex-pair
code 4 (bit 1 means i->j for i<j) is the converse of NT_4 and has sink
energies {1:1,2:3,3:1}. At orders at most 3, a tournament is transitive or a
cyclic triangle possibly with fewer vertices; transitive tournaments have
one sink and the triangle's three sinks all have energy 1.

## 5. Explicit finite-field bijection: paths to subspace cells

Now q denotes a prime power and F_q a field of that size. Decorate P_A by
one freely chosen field element lambda_(a,b) for every inversion pair
(a,b) with a in A, b in B, b<a. For each a in A form a row in F_q^m:

    v_a = e_a + sum_(b in B, b<a) lambda_(a,b) e_b.              (9)

Map the decorated path to the subspace spanned by these rows.

**Theorem 6.** This is a bijection from decorated non-P0 Hamilton paths of
NT_(m+2) to ALL linear subspaces of F_q^m. Paths with |A|=r map to r-dimensional
subspaces. Consequently

    (Z_(m+2)(q)-q)/q^m
      = # {linear subspaces of F_q^m}
      = sum_(r=0)^m [m choose r]_q.                            (10)

**Proof.** In (9), pivot columns A form the identity matrix, and a row with
pivot a has only zero entries to the right of a. This is the unique reverse-
column reduced row-echelon form: reverse the columns to obtain ordinary
RREF, and order the rows by their reversed pivot positions. Every subspace
has exactly one such form. Its rightmost pivot set supplies A, the
complement supplies B, and its free entries supply the decorations.
Conversely all freely chosen entries in precisely those positions give
such a form. This proves both injectivity and surjectivity and preserves
rank. Equation (10) follows from (6). QED.

In standard language these pivot classes are the Grassmannian Schubert
cells. The argument has supplied their coordinates explicitly, so no
abstract identification or merely equal-cardinality argument is needed.
Over F_q the cell indexed by A has q^inv(A,B) points. Over a general field
the same matrices give an affine cell of dimension inv(A,B).

For example m=3 gives 2^3=8 non-P0 Hamilton paths, one for every subset, plus P0 for a total of 9. Their decorated
objects over F_2 are the 16 subspaces of F_2^3, in dimensions 0,1,2,3 with
counts 1,7,7,1. The distinction between eight cell indices and sixteen
points is essential. It is not an H=16 tournament construction.

The map uses the known transitive order, exceptional backward arc, and
field-coordinate decorations. It is not invariant under forgetting those
choices, and it does not identify arbitrary tournament Hamilton paths
with arbitrary subspaces. The odd-to-natural index on this family is
(H-1)/2=2^(n-3), so this family alone covers powers of two, not all
nonnegative integers other than 3 and 10.

## 6. Verification and boundaries

Run `python3 verify.py > replay.json` and compare with RESULT.json; repeat
with `python3 -O verify.py`. The checker has explicit exceptions and no
assertion-dependent logic. It imports no repository data or third-party
package.

- All 1,099 labelled tournaments through order 5: subset-DP Hamilton counts
  independently equal permutation local-sink counts; every adjacent-swap
  energy difference is checked; all face-incidence coefficients and the
  alternating identity are checked. The complete order supports are
  {1}, {1}, {1,3}, {1,3,5}, {1,3,5,9,11,13,15}.
- NT_n through n=12: every constructed path and energy, every fixed-rank
  Gaussian coefficient, and the sharp maximum (8) are checked. Through
  n=8 a separate exhaustive permutation scan proves finite completeness.
- The field map is fully enumerated for q=2,m<=5 and q=3,m<=4. Standard
  RREF checks pairwise distinct images and rank counts. For m<=3, an
  independent closure construction adjoining every field vector generates
  all subspaces and confirms surjectivity without using the pivot recipe.
- Deterministic ordinary acyclic orientations with 7 and 21 sinks are
  checked at orders 4 and 5, respectively, and their coherence violations
  are explicitly reported.

All-order claims follow from the proofs, not extrapolation of those tests.
The infinite spectrum-completeness conjecture, namely surjectivity of
H -> (H-1)/2 onto the nonnegative integers other than 3 and 10, and a
universal constructive realization grammar remain OPEN. This packet gives a faithful geometric carrier, sharp no-go tests
for two tempting relaxations, and an explicit nontrivial finite-field lift
of one classical subfamily, rather than a claimed global classification.
