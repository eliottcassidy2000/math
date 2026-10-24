# Converse-merged refinement-error spectrum at order ten

24 October 2026. **PROVED:** the stable-cell obstruction, error-refined SCC
identities, and the order-ten extraction from the inherited exhaustive
classification. **FINITE-EXACT:** full arc-isomorphism, automorphism, matrix-orbit,
and bounded SCC-word certificates. No claim of a full metagraph wiring law,
full order-ten tournament census, or new anti-involution group theorem is made.

## 1. Object, inherited inputs, and result

Begin directed color refinement on a finite tournament T with one color. At
every round retain the old color and the outneighbor count into every current
color. The eventual parts are the **actual stable cells**. If O ranges over
Aut(T)-vertex orbits and C(O) is its stable cell, define

    e(T) = sum_O ( |C(O)| - |O| ).

This is the pointed-refinement excess. It is not a Hamiltonian-path count,
blue-line count, black-line count, or the defect of a different root matrix.
An ordinary class identifies only directed isomorphisms. A converse-merged node
identifies an ordinary class with its reverse T^op.

The Session51 order-ten packet supplies the exhaustive strong positive-error
classification, with ordinary histogram

    4 t^14 + 32 t^18 + 84 t^40.

The inherited Session50 order-nine and regular-nine packets, together with SCC
additivity, supply the complete ordinary positive-error order-ten histogram

    P_10(t) = 12 t^14 + 46 t^18 + 4 t^36 + 84 t^40 + 14 t^72.

These statements are recorded in the pinned source theorems and receipts. No
unrestricted order-ten atlas or numerical total count a_10 is an input.

**Theorem.** Exactly 24 of these 160 ordinary classes are self-converse. All
24 are strong, rigid, and have stable cells 5+5 and e=40. There are four in
the crossing-degree-one sector and twenty in crossing degree two. Hence

    P_SC,10(t) = 24 t^40,
    P_merged,10(t) = (P_10(t) + P_SC,10(t))/2
                   = 6 t^14 + 23 t^18 + 2 t^36 + 54 t^40 + 7 t^72.

The positive-error slice therefore contains 92 converse-merged nodes. Its
node-weighted error sum is 3234 and its raw second moment is 133908. The
ordinary and self-converse error sums are respectively 5508 and 960, so
3234=(5508+960)/2. The factor is coefficientwise Burnside, not naive halving.

A self-converse node need not be a pure-blue point of symmetry. Neither these
counts nor e(T) determines any blue/black incidence or metagraph wiring.

## 2. Converse preserves the actual stable partition and the error

For a vertex v and a current cell D,

    out_(T^op)(v,D) = |D| - 1_{v in D} - out_T(v,D).

Because refinement retains the current cell, two vertices in that cell have
equal outcount vectors in T exactly when they do in T^op. Induction from the
one-cell partition shows that the partitions agree as sets of vertex subsets
at every round, although arbitrary names for their colors can differ.
Aut(T)=Aut(T^op) as permutation groups. Consequently e(T)=e(T^op), so the
error is a well-defined function on converse-merged nodes.

Every stable cell C induces a regular tournament: equitability makes the
internal outdegree constant, and counting its internal arcs gives
|C|(|C|-1)/2=|C|d. Thus every stable cell has odd size.

## 3. Stable-cell anti-involution obstruction

**Lemma.** For a self-converse tournament of order n, some anti-automorphism q
satisfies q^2=id. Every such q has exactly n mod 2 fixed vertices and exactly
n mod 2 fixed actual stable cells. All other cells occur in equal-sized pairs.
For n=0 or 1, q may be the identity; “involution” here means order dividing two.
For n>=2 an anti-involution has order exactly two.

**Proof.** A tournament automorphism cannot have even order: an appropriate
power would be a nonidentity involution, which exchanges some pair and hence
reverses its arc. Given an anti-automorphism a, its square is an automorphism
of odd order m. Then q=a^m is anti because m is odd and q^2=(a^2)^m=id. This
odd-power mechanism is explicitly inherited from THM-647,
`THM-647-anti-redei-theorem.md`, and is not claimed as new here.

An anti-automorphism cannot fix two vertices, since it would reverse their
arc while leaving both endpoints fixed. For q^2=id, its number of fixed
vertices has the parity of n, so it equals n mod 2. Naturality of refinement
and Section 2 show that q permutes the actual stable cells. A fixed cell is
an odd-sized set acted on by an involution, so contains a fixed vertex. A
fixed vertex conversely lies in a fixed cell. There is at most one such
vertex, hence exactly n mod 2 such cells, each with exactly one fixed vertex.
The remaining cell orbits have length two and equal cell sizes. ∎

In particular, at even order, **every stable-cell size must occur an even
number of times**. At odd order exactly one size has odd multiplicity, namely
the size of the unique fixed cell. These are necessary, not sufficient,
conditions: the packet contains an explicit non-self-converse ten-tournament
with actual stable cells 5+5.

For the strong order-ten classification, 7+1+1+1 and 3+3+3+1 both violate the
even-order condition. Only 5+5 can contribute self-converse positive error.

**Divisibility corollary.** Every tournament has even e. Every even-order
self-converse tournament has e divisible by four. Indeed Aut(T) has odd order
(the absence of an involution plus Cauchy's theorem), so all vertex orbit
sizes are odd. If a stable cell C contains r automorphism orbits, then r is
odd and its error contribution is |C|(r-1), an even number. An anti-map
normalizes Aut(T), so cells it exchanges have equal orbit counts and equal
error contributions. At even order the anti-involution pairs every cell,
making the total a sum of doubled even contributions. This is only a
necessary condition: error divisible by four does not certify self-converseness.
The bad regular-seven class has error14, showing why even order is essential.


## 4. The exact fixed-point count in the two-five-cell sector

Put C5 on each cell with arcs x->x+1,x+2 modulo five. The automorphism group
of C5 is exactly its translations: an automorphism must preserve the unique
source among each vertex's two outneighbors, which identifies successor.
Let B be the crossing matrix from the lower-degree cell to the higher-degree
cell. Every row and column has degree d=1 or 2. The two vertex degrees are
2+d and 7-d, so these are the actual stable cells, intrinsically ordered.
All ordinary isomorphisms act by independent row/column translations C5^2.

A nontrivial stabilizer cannot contain a coordinate-axis translation, since
that would force equal rows or columns and then a uniform matrix. Therefore
it is a nonzero-slope subgroup of order five, or trivial. In the former case
its projections act transitively on both cells and e=0; in the latter case
the tournament is rigid and e=10(5-1)=40.

The complete counts are eight ordinary classes for d=1 (four rigid) and 88
for d=2 (eighty rigid). The labelled matrix counts are 120 and 2040. For d=2,
the bipartite cycle decompositions 10 and 6+4 give respectively 1440 and 600.
These reproduce the inherited structural count, not an extra assumption.

Converse followed by exchanging the cells and negating their coordinates acts
on crossing matrices by

    f(B)_(i,j) = B_(-j,-i).

It normalizes the translation group. The number of translation-orbits fixed
by f is the average, over the 25 translations g, of |Fix(g f)|. This twisted
Burnside statement follows by counting pairs (B,g) with gf(B)=B: every fixed
ordinary orbit contributes exactly 25 such pairs, and every nonfixed orbit
contributes zero.

For the five diagonal choices of translation parameters the fixed condition
is, after changing column coordinates, that B is symmetric. There are S_1=26
symmetric permutation matrices (involutions of five letters) and S_2=112
symmetric degree-two matrices. The latter count decomposes as: five-cycle 12;
path on five vertices with marked endpoints 60; triangle plus marked edge 10;
marked three-path plus marked edge 30. Here a diagonal mark contributes one
to matrix degree, and is not a loop in the original tournament.

For each of the twenty off-diagonal parameter pairs, the squared position
transformation is a nontrivial translation; its fixed matrices are determined
by a subset of the five sums i+j invariant under a reflection. There is one
such degree-one subset and two degree-two subsets. Hence the total ordinary
self-converse counts are

    F_1 = (5*26 + 20*1)/25 = 6,
    F_2 = (5*112 + 20*2)/25 = 24.

These totals include zero-error classes. To remove them correctly, write a
matrix with nontrivial stabilizer as B_(i,j)=1_{j-s i in D}, where s is a
nonzero slope and |D|=d. Modulo translations there are binom(5,d)/5 choices of
D for each slope. The converse sends s to s^(-1). Only slopes +1 and -1
can be fixed. For +1 the subset class is unchanged. For -1 it is negated;
every singleton or two-element subset of F5 is translation-equivalent to
its negative. Thus there are exactly two fixed zero-error classes for d=1
and four for d=2. The positive-error fixed counts are

    6-2=4,    24-4=20.

The general regular-block formulas and the coefficient values are inherited
from the Session40 regular-block packet. The new extraction here retains the
error weight and separates its zero-error fixed classes. `matrix_check.py`
rederives the complete 120/2040 matrix universes and all twisted fixed sets;
`check.py` separately tests unrestricted vertex isomorphisms on each actual
ten-vertex representative.

### 4.1 Error extraction for every prime two-carousel family

The same argument extends, with an explicit family restriction. Let p be an
odd prime and R_p have arcs x->x+j for 1<=j<=(p-1)/2. Take two R_p cells with
a nonuniform d-regular crossing matrix, canonically 1<=d<=(p-1)/2. The internal
arc x->x+j has exactly (p-1)/2-j common outneighbors, so successor is intrinsic
and Aut(R_p)=C_p. As above the two cells are actual degree-distinguished stable
cells and their full automorphisms are the translation stabilizer of the
crossing matrix. Every nontrivial stabilizer is a nonzero-slope line and is
transitive on each cell. Hence every class has error either zero or 2p(p-1),
the latter occurring exactly when the tournament is rigid.

Let A_(p,d) be the number of labelled binary p-by-p d-regular matrices and
S_(p,d) the symmetric count. Put B=binom(p,d) and
c=binom((p-1)/2,floor(d/2)). The inherited Session40 formulas give

    N = (A_(p,d)+(p-1)^2 B)/p^2,
    F = (S_(p,d)+(p-1)c)/p.

For a nontrivial slope s, write the matrix as
B_(i,j)=1_{j-s i in D}. For each slope the d-subsets D have B/p translation
orbits, since a nonempty proper subset of a prime cycle has no nontrivial
period. Distinct slopes cannot represent the same nonuniform matrix: two
independent line stabilizers would generate the full translation group.
Thus the number of zero-error ordinary classes is (p-1)B/p.

Converse sends (s,[D]) to (s^(-1),[s^(-1)D]). Only s=+1 and s=-1 can be fixed.
The +1 slope contributes all B/p subset classes. The -1 slope contributes
only reversible subset necklaces. Their number is c: twisted Burnside for
the p translations averages p reflections, each fixing exactly c d-subsets
(one singleton position and (p-1)/2 paired positions, with the singleton
selected exactly when d is odd). Therefore the exact positive-error ordinary,
self-converse, and merged class counts in this restricted family are

    N_positive = (A_(p,d)-(p-1)B)/p^2,
    F_positive = (S_(p,d)-B-c)/p,
    M_positive = (N_positive+F_positive)/2.

Each contributes the single error value 2p(p-1). The zero-error self-converse
count removed from F is B/p+c, not always 2B/p: not every subset necklace is
reversible. At p=7,d=3 there are five subset necklaces but only three
reversible ones, an explicit hostile to that tempting overgeneralization.

The p=7 coefficient pairs (A,S), independently recomputed by capacity
recursions in `prime_check.py`, are (5040,232), (3110940,6912), and
(68938800,35150) for d=1,2,3. They yield respectively

    ordinary positive:        102, 63486, 1406910;
    self-converse positive:    32,   984,    5016;
    converse-merged positive:  67, 32235,  705963.

Every one of those positive classes has error84. These are exact counts in
the specified two-R7-cell subfamily at order fourteen, not a full order-fourteen
classification. The checker enumerates small subset necklaces and computes
matrix coefficients; it does not enumerate the 68,938,800 degree-three
matrices. The prime hypothesis is load-bearing, and no composite-order
extension of the stabilizer dichotomy is asserted.

## 5. Error-refined SCC series and the order-ten padding exclusion

SCCs of a tournament have a unique linear order, earlier components beating
later components. Their possible outdegree intervals are disjoint, so degree
refinement already separates distinct SCCs. Subsequent stable refinement in
each component agrees with its internal refinement: exterior cells contribute
only uniform constants. Automorphisms preserve each ordered SCC and act
independently inside them. It follows that e is additive over SCCs.

Use formal power series with finite polynomial coefficients in t:

- S(z,t): ordinary nonempty strong classes, weighted z^n t^e;
- A(z,t): all ordinary classes, including the empty class;
- F_s(z,t): self-converse nonempty strong ordinary classes;
- F(z,t): all self-converse ordinary classes, including the empty class;
- M(z,t): converse-merged nodes, including the empty node.

Unique ordered SCC decomposition gives A=1/(1-S). Converse reverses the SCC
list and converses each component. A class is self-converse precisely when
this ordered list is equal, in ordinary component types, to its reversed
converse. Its left half is an arbitrary ordered strong list, its right half
is forced, and its center is either empty or a self-converse strong class.
Every paired component doubles both order and error. Consequently

    F(z,t) = (1+F_s(z,t)) A(z^2,t^2)
           = (1+F_s(z,t))/(1-S(z^2,t^2)),
    M(z,t) = (A(z,t)+F(z,t))/2.

No unordered-half quotient is taken: the ordered SCC positions distinguish
the left half under ordinary isomorphism. The converse-merged identity is
coefficientwise, since converse preserves both n and e.

Let A(z)=A(z,1), C(z)=F(z,1), E_s(z)=partial_t S(z,t)|_(t=1),
E_SC,strong(z)=partial_t F_s(z,t)|_(t=1), and define E and E_SC similarly
from A and F. Differentiation gives the all-order error transfer laws

    E(z) = A(z)^2 E_s(z),
    E_SC(z) = A(z^2) E_SC,strong(z)
              + 2 C(z) A(z^2) E_s(z^2),
    E_merged(z) = (E(z)+E_SC(z))/2.

The factor two is essential and comes from t^2. The first positive-error
strong order is seven, so paired positive-error components cannot appear
before order fourteen. Below fourteen a positive-error self-converse class
must instead have a positive-error self-converse central strong component.
At even order ten such a center has even order at most ten; all smaller
even strong orders have zero error. Therefore every positive-error
self-converse ten-class is strong, proving the padding exclusion.

Do not strengthen this to “all positive-error self-converse classes below
fourteen are strong”: the actual class 1 -> R7 -> 1 is self-converse at order
nine and has error14. The fourteen-vertex class R7 -> R7^op is self-converse
with error28, demonstrating the sharp paired-component boundary and the
doubled error. Both are full-arc-tested in the receipt. `series_check.py`
independently enumerates all ordered SCC words through total order sixteen
on an explicit actual-core alphabet and verifies both formal identities and
the derivative; omission of t^2 first fails at order fourteen in this control.
That bounded alphabet is not asserted to contain all tournaments.

## 6. Evidence and precise scope

`check.py` reads pinned certificates, but imports no source implementation.
It uses its own stable refinement and directed arc backtracking. It:

1. takes the 120 positive strong representatives from the audited structural
   packet and reconstructs all 40 nonstrong representatives from the 17
   positive strong nine-cores and six ordered paddings of the bad seven-core;
2. rebuilds the fifteen regular nine classes from their eight-vertex parent
   representatives and checks their full automorphism errors and distinctness;
3. checks actual stable cells, full automorphism orbits, converse invariance
   of error, and all exact anti-isomorphisms for every ten-representative;
4. checks ordinary pairwise nonisomorphism and finds the unique actual converse
   partner among all 160 representatives, with explicit vertex-map witnesses;
5. powers every anti-map to an anti-involution and verifies its fixed vertices
   and stable-cell action; and tests the same lemma on every labelled
   tournament through order five;
6. records the 5+5 non-self-converse necessity-not-sufficiency hostile, the
   central-padding order-nine control, and the doubled-core order-fourteen
   control.

Normal and optimized Python produce byte-identical receipts. Independent
matrix and bounded SCC-word paths test the two separate quotient mechanisms.
The classification's completeness remains tied to the explicitly inherited
strong-sector proof and complete almost-regular-eight parent slice. This is
not an enumeration of the full order-ten universe, not a total a_10 prediction,
and not a reconstruction of the metagraph's edges. Its output is an exact,
converse-merged, error-weighted node slice and an all-order SCC transfer law.

## 7. Native metagraph incidence, without a wiring inference

Under the repository's line-metagraph convention, a presentation is blue
exactly when its staircase reflection fixes it. THM-644(a),
`THM-644-gridsym-fiber-law-and-anti-redei.md`, proves that this is equivalent
to positional reversal being an anti-automorphism of its tournament. Hence a
blue-incident node is self-converse. Conversely THM-647 proves that for any
anti-involution q, a self-converse tournament has an odd, positive number of
Hamiltonian paths reversed by q. Reading one such path as a presentation
gives a blue port. Only THM-644(a) and the proved THM-647 are used here, not
the older conjectural status of THM-644(c).

Thus among the 92 positive-error merged nodes, exactly 24 are blue-incident
and exactly 68 are black-only. Their error split is:

    e=14:  0 blue-incident,  6 black-only;
    e=18:  0 blue-incident, 23 black-only;
    e=36:  0 blue-incident,  2 black-only;
    e=40: 24 blue-incident, 30 black-only;
    e=72:  0 blue-incident,  7 black-only.

This is an incidence classification, not an identification of the 24 with
all-blue points of symmetry. In fact `port_check.py` directly establishes
that all 24 are **mixed**, and none is all-blue. Each is rigid, so literal
Hamiltonian-path counts equal presentation-port counts without an automorphism
denominator. Indeed two paths give the same normalized base-path presentation
only if the vertex map carrying one path to the other is an automorphism;
rigidity makes that map the identity. Conversely every presentation of this
class supplies such a path. This proves the needed rigid special case directly
and does not promote the historical general fibre-law status of THM-644(b). The script counts all Hamiltonian paths twice, by subset dynamic
programming and by literal word traversal; it counts blue paths independently
by the anti-involution's half-system and compares with a literal staircase-
reflection test on every normalized Hamiltonian word in the full traversal.
It records one legal blue path and one legal black path for every node.

The 24 nodes together have 4200 blue presentation ports and 337304 black
presentation ports. Their blue-port count histogram is

    119:2, 123:1, 131:1, 179:1, 181:2,
    183:4, 185:4, 187:5, 189:2, 191:2.

These are finite port counts, not counts of edges lying inside this node
slice. Their partners may lie outside it; loops and endpoints require the
actual line pairing, which this packet does not reconstruct. In particular
no division of4200 by two is claimed to count induced blue edges, and no
black/blue wiring follows from the refinement error alone.
