# Balanced metagraph ports, a half-size blue polynomial, and an Eulerian cycle lift

24 October 2026. **PROVED** all-size identities and exact algorithms.
**FINITE-EXACT** scalar-input counts through even order twelve, independently
controlled as specified below. No tournament-class atlas is an input.
No new general formula for A000568, full wiring reconstruction, complexity
improvement theorem, or literature-priority claim is made.

## 1. Inheritance, conventions, and what is added

The fixed Hamiltonian-path word carrier, native complement C, reflected-converse
R, and the ordinary word mass h(T)=H(T)/|Aut(T)| are inherited. The precise
sources are pinned in `sources/fixed-path-fibre.txt` (THM-781, Hamiltonian-path
inverse metagraph fibre), `sources/tiling-dictionary.txt` (THM-1430, tiling/class/
merged-metagraph dictionary), and `sources/path-transport.md` (Session40 actual
Hamiltonian-path degree transport). Session51's `sources/score-layer-separation.md`
already proves the only two-score equal-multiplicity layer supporting internal
edges is the balanced layer, with necessary downhill endpoints.

THM-554, `sources/THM-554-score-partition.md`, already gives exactly the
full degree product F_n below (in its reversed vertex gauge). THM-555,
`sources/THM-555-score-multiplicities.md`, already records the balanced total
masses 5, 157, 51,949, and 204,019,829 at n=4,6,8,10. Those totals are
**reproduced**, not first discoveries here. Only the degree product and those
explicit totals are inherited from these mixed-status archives; their other
moment, complement, or maximal-H claims are not dependencies of this packet.

This packet extracts an endpoint- and color-resolved n-only census of that layer and
its internal and boundary ports. It sharpens the endpoint condition to an iff
for the balanced target, folds the blue census to half as many degree variables,
and identifies the entire internal carrier with Eulerian orientations of the
complement of an odd cycle. The latter gives a second independent count, a
regular-parent cycle interpretation, and rotation congruences.

The degree-generating product is inherited both as the standard orientation
construction and as the exact prior repository formula THM-554, not claimed
as a new general enumerative formula. The free
Hamiltonian-path quotient and ordinary/merged mass conventions are not renamed
as discoveries. The regular one-vertex completion and the bijection between balanced classes
and pointed regular classes are already proved in Session51,
`sources/regular-eleven-theorem.md`, Sections 5 and 7, and are used here
without an originality claim. The extra cycle-marked transport in Section 5
below is the connection developed in this packet.

Concept board: full source mass; endpoint score transport; blue reflection
orbits; regular cycle completion; Eulerian orientations; class/port quotient
loss. Canonical hostiles are n=2 (C is not free), n=6 (five ordinary balanced
classes but 157 words), and loops (counted twice by port degree, not removed).

Write n=2p>=2, vertices 0,...,n-1, and retain the directed path
0 -> 1 -> ... -> n-1. Each pair i<j with j>=i+2 is a free bit, one meaning
 i -> j. C flips every free bit. R permutes its coordinate (i,j) to
(n-1-j,n-1-i), so its source is the converse read in reverse path order.
They commute. A blue word is R-fixed; other words are black. A line is a
C-orbit of words, retaining multiplicity even when its endpoint classes coincide.

A balanced tournament has exactly p vertices of degree p and p of degree p-1.
Let A_n count balanced source words; let A_n^{ab}, a,b in {L,H}, distinguish
start/end degree, where L=p-1 and H=p. Let D_n=A_n^{HL}. Superscript blue or
black restricts the word color. These are word/port counts, not class counts.

## 2. The exact internal criterion and complete boundary score types

The inherited actual-path identity is

    d_C(v) = n-1-d(v) + [v=0] - [v=n-1].                 (1)

For a balanced source, the partner is balanced **if and only if** the start
has degree p and the end degree p-1. Indeed each interior degree is simply
interchanged between p and p-1. The start has degree p or p+1, and the end
has degree p-1 or p-2, with the permitted choices occurring exactly at HL.
Consequently the internal word set is closed under both C and R.

For p>=2 the exact possible target degree histograms, conditioned on source
endpoint type, are:

    HL: (p-1)^p, p^p;
    LL: (p-1)^(p+1), p^(p-2), (p+1)^1;
    HH: (p-2)^1, (p-1)^(p-2), p^(p+1);
    LH: (p-2)^1, (p-1)^(p-1), p^(p-1), (p+1)^1.        (2)

A zero multiplicity means the degree is absent. These are consequences for
an existing source word, not existence claims for every conceivable profile.
The LL/HH target types are converse to each other; LH is converse-invariant.
The boundary therefore has only two merged score types. R gives
A_n^{LL}=A_n^{HH}. Blue words have only HL or LH endpoints because their paired
positions have complementary scores. Thus all blue boundary ports are LH.
The respective cyclic-triple losses are 0,1,1,2, by the inherited endpoint law.

For n>=4, C is free. Therefore the balanced induced C-line multigraph has

    internal lines        = D_n/2;
    internal blue lines   = D_n^blue/2;
    internal black lines  = (D_n-D_n^blue)/2.            (3)

There are A_n-D_n boundary lines incident to the layer, each with exactly one
balanced endpoint. Its complete incident-line count is A_n-D_n/2. These
statements include loops and parallel lines. A merged loop contributes two
internal ports. They do not count distinct unordered pairs of nodes, or discard
loops, and they do not additionally quotient the line set by R.
At n=2 there is one balanced blue word and C is the identity; formulas (3)
are deliberately not applied.

For n>=4, internal black words form free four-element orbits under <C,R>.
Indeed C is free, R has no fixed black word by definition, and CR has no fixed
word at all: the reflection-fixed free pair {0,n-1} would require its bit to
equal its complement. Hence D_n-D_n^blue is divisible by four, so the number
of internal black lines is even. Were one additionally to quotient C-lines by
R, the resulting line-orbit count would be (D_n+D_n^blue)/4. That is a
**different** graph convention, not the line count in (3).

## 3. The full degree coefficient polynomial

Let

    F_n(x_0,...,x_(n-1))
      = (prod_(i=0)^(n-2) x_i)
        prod_(0<=i<j<n, j>=i+2) (x_i+x_j).              (4)

The coefficient of x^d is exactly the number of source words with outdegree
vector d: each path edge supplies its tail variable and each free edge chooses
one of its two possible tails. Let [balanced; a,b] mean sum the coefficients
whose degrees all lie in {p-1,p}, with endpoint choices a,b. Then

    A_n^{ab} = [balanced; a,b] F_n.                      (5)

The number of high coordinates is automatically p, since every monomial has
total degree n(n-1)/2. It need not be separately constrained.

The executable `python3 engine.py 12` accepts only the scalar n. It evaluates
all four slices by a row-forgetting DP. A state is the **ordered** vector w
of wins already earned by unprocessed vertices, in their fixed path order.
At r>=2 remaining vertices, the current vertex must beat the next. Choose a
subset S of the other r-2 future vertices to beat, requiring

    w_0 + 1 + |S| in {p-1,p}.

Delete w_0; leave the next vertex's past wins unchanged; add one to every other
future vertex not in S. At r=1 accept w_0=p-1 or p and record the end type.
At the initial step constrain the start to L or H. Memoize identical ordered
vectors; no sorting of the remaining path positions is allowed or performed.

This is exact elimination in (4), not an atlas, sampler, fitted formula, or
constant-width recurrence. A crude state bound is
sum_(i=0)^(n-1) (min(i,p)+1)^(n-i); transition cost remains exponential in the
number of unprocessed free neighbors. No polynomial-time bound is asserted.

## 4. A half-size polynomial for all blue counts

For p variables indexed by the first half of the path, define

    G_p(x) = (prod_(i=0)^(p-1) x_i)
       prod_(0<=i<j<p, j>=i+2) (x_i+x_j)
       prod_(0<=i<j<p) (1+x_i x_j)
       prod_(i=0)^(p-2) (1+x_i).                       (6)

Its exponents are the degrees of those p vertices in a blue word. To prove
this, assign one common bit to each R-orbit of free pairs:

- A pair wholly in the first half is paired with a pair in the second half;
  its two choices contribute x_i or x_j to first-half degrees.
- A nonfixed crossing-pair orbit has two arcs, either both with first-half
  tails i,j or both with second-half tails, contributing x_i x_j or 1.
- A reflection-fixed crossing pair (i,n-1-i) contributes x_i or 1. The
  central adjacent pair is the fixed path arc and is not a free coordinate.

Every first-half vertex has one outgoing path arc, supplying the initial
product. There are p(p-1) independent bits, agreeing with the inherited
blue-cube dimension floor((n-1)^2/4). Positional reflection enforces
 d_i+d_(n-1-i)=2p-1. Therefore

    A_n^blue = sum_[all d_i in {p-1,p}] [x^d] G_p;
    D_n^blue = sum_[all d_i in {p-1,p}, d_0=p] [x^d] G_p. (7)

No extra cardinality condition is needed. Black counts are the differences
from (5). Multiplying (6) while discarding exponents above p is exact because
all coefficients and exponents are nonnegative. There are at most p^p retained
exponent vectors. This is a useful half-dimension reduction, not a claim that
blue counts recover the nonblue carrier or its wiring.

## 5. The regular-cycle/Eulerian orientation bridge

Given a balanced tournament T, add one vertex z, direct z toward every high
vertex, and direct every low vertex toward z. The resulting R(T) is regular
of order N=2p+1 and degree p. It is the **unique** regular extension on that
specified new vertex: each old degree forces its incident arc, and exactly p
arcs leave z. Deleting a vertex from any regular odd tournament recovers a
balanced tournament. The maps commute with isomorphisms and converse.

A Hamiltonian path P in T has HL endpoints iff adjoining z at its beginning
and closing at z gives a Hamilton cycle in R(T). Every directed Hamilton cycle
of R(T), counted up to cyclic rotation but not reversal, has exactly one such
reading beginning at z. Hence

    H_down(T) = number of directed Hamilton cycles of R(T). (8)

At the fixed-word level, the internal slice bijects to regular tournaments
on 0,...,N-1 containing the particular directed cycle
0 -> 1 -> ... -> N-1 -> 0, with z=N-1. Remove its cycle arcs. The remaining
undirected graph is K_N minus C_N, of even degree N-3=2(p-1). Regularity asks
for outdegree p-1 on every vertex, exactly an Eulerian orientation. Thus

    D_(2p) = EO(K_(2p+1) minus C_(2p+1))
       = [prod_(i=0)^(2p) x_i^(p-1)]
         prod_({i,j} not an edge of C_(2p+1)) (x_i+x_j). (9)

This does not enumerate the 2^binom(n-1,2) original words. `eulerian_anticycle_dp`
independently evaluates (9) by fixed-degree elimination on the different
N-vertex graph and agrees with (5) through n=12.

The bridge preserves the operation, not just its count. C on the balanced
word lifts to reversal of every arc outside the specified Hamilton cycle.
Both incident z-arcs belonging to the cycle are retained; the other z-arcs
reverse because the internal vertices' high/low statuses interchange.
Therefore cycle reversal of the complement is the actual transported C,
not an arbitrary graph with matching degree data.

For ordinary classes, let R be a regular parent, and let O be a vertex orbit
of Aut(R), represented by z. Its deletion card T=R-z has
Aut(T)=Aut(R)_z, because every automorphism extends uniquely to fix z. Two
cards are isomorphic iff their deleted vertices lie in one Aut(R)-orbit.
Writing c(R) for the cycle count, their internal word degrees satisfy

    d(T)=H_down(T)/|Aut(T)|
        = c(R)/|Aut(R)_z|
        = |O| c(R)/|Aut(R)|.                           (10)

Thus the internal degrees of all deletion-card classes of a common parent
are proportional to its vertex-orbit sizes, and their sum is
N c(R)/|Aut(R)|. The intermediate ratio c(R)/|Aut(R)| need not be an
integer: a Hamilton cycle marked only up to rotation can have rotational
stabilizers. Only the rooted-path quotients in (10) are asserted integral.
This is a source-structure constraint. It does not say that
C stays among cards of one regular parent: the operation can change R.

At n=6 the five ordinary balanced classes have (H,Aut,H_down) equal to
(41,1,17), (43,1,15), (43,1,15), (45,3,15), (45,3,24).
In particular H=45 and word mass 15 can still have internal degree 5 or 8.
Full source mass alone does not determine internal degree.

## 6. Rotation congruences and explicit blue lower bounds

The rotations of the specified N-cycle act on the Eulerian orientations in
(9). A full one-step rotationally invariant orientation is circulant. There
are p-1 noncycle difference pairs, each independently oriented, so it has
exactly 2^(p-1) choices. All are invariant under reversed-converse reflection,
so they lie in the internal blue slice. Consequently

    D_n^blue >= 2^(p-1),  and D_n >= 2^(p-1).           (11)

When N=n+1 is prime, every nonfixed rotational orbit has size N, yielding

    D_n = 2^(p-1) mod (n+1).                           (12)

For n>=4 and prime n+1, the internal line count therefore satisfies
D_n/2 = 2^(p-2) mod (n+1).

For arbitrary odd N there is an exact divisor refinement. If a rotation has
d=gcd(N,k) vertex orbits of length L=N/d, let f_N(d) count its fixed Eulerian
orientations. For d=1, f_N(1)=2^((N-3)/2). For d>1, put

    r=((d-1)L-2)/2;
    f_N(d)=2^[d(L-1)/2] [prod_i x_i^r]
        prod_(0<=i<j<d) (x_i+x_j)^[L-[{i,j} in E(C_d)]]. (13)

Here C_3 has its usual three distinct edges. Within each vertex orbit, its
(L-1)/2 edge-orbits each give one outneighbor per vertex in either orientation;
this gives the power of two and the fixed degree offset. Between different
vertex orbits there are L matching edge-orbits, with one removed for each
edge of the quotient d-cycle. Choosing each orbit's direction supplies one
factor x_i+x_j. This proves (13). In particular f_N(N)=D_(N-1).

Burnside therefore gives

    cyclic orbits = (1/N) sum_(d|N) phi(N/d) f_N(d).     (14)

All proper-divisor terms have fewer variables, so (14) determines D_(N-1)
modulo N without the largest-dimensional term. This is a congruence and an
orbit-count identity, not a recurrence determining the integer D_(N-1).
Cyclic-orbit counts here are not unmarked tournament-class counts.

## 7. Exact scalar-input counts

All numbers below are computed from n alone; none is fitted or extrapolated.
The A column through n=10 reproduces THM-555, as credited above.

| n | balanced ports A | internal ports D | blue source ports | blue internal ports | internal lines | blue internal lines |
|---|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 1 | 1 | 1 | degenerate | degenerate |
| 4 | 5 | 2 | 3 | 2 | 1 | 1 |
| 6 | 157 | 60 | 21 | 12 | 30 | 6 |
| 8 | 51,949 | 18,152 | 499 | 280 | 9,076 | 140 |
| 10 | 204,019,829 | 67,123,248 | 41,699 | 22,904 | 33,561,624 | 11,452 |
| 12 | 10,188,802,555,033 | 3,213,825,974,560 | 12,569,477 | 6,803,648 | 1,606,912,987,280 | 3,401,824 |

At n=10 the endpoint counts (LL,LH,HL,HH) are
(44,967,565,46,961,451,67,123,248,44,967,565). There are 136,896,581 boundary
ports, of which 18,795 are blue. At n=12 they are
(2,301,770,367,447,2,371,435,845,579,3,213,825,974,560,2,301,770,367,447),
with 6,974,976,580,473 boundary ports, 5,765,829 of them blue.

The ordinary row DP uses 19,291 cached ordered states at n=10 and 426,340 at
n=12. The blue polynomial peaks at 2,309 and 32,877 states, respectively.
The full fixed-path cubes would have 2^36 and 2^55 words. These measurements
state the actual compressed states used, not a claim of uniformly low cost.

## 8. Class, merged-node, edge, and loop loss ledger

For ordinary classes c within the balanced layer,

    A_n = sum_c H(c)/|Aut(c)|;
    D_n = sum_c H_down(c)/|Aut(c)|.                     (15)

For a converse-merged node, the mass and internal degree are respectively
h and d if self-converse, and 2h and 2d if it is a pair of non-self-converse
classes. There is no blue mass on a non-self-converse ordinary class.
Merging class vertices does not identify distinct word ports or C-lines.
An additional R quotient of ports or lines would be a different construction.

At n=6 the literal control finds 5 ordinary classes, 4 merged nodes, 157 source
ports, and 60 internal ports. Its 30 internal lines include 4 ordinary loop
lines and 8 merged loop lines. Thus neither word mass nor internal edge mass
is itself a class count, and merging converse classes can create loops.
All these finite counts are recomputed from actual matrices and permutations.

The aggregate coefficient polynomials forget which words belong to a common
source class, which source class contains their C partner, and which individual
ports form loops. Recovering full induced wiring requires that relational
sidecar. A_n, D_n, their colors, and endpoint partitions do not supply it.
In particular no A000568(n) value is inferred from this degree census.

### Inherited order-ten class counts combined with the new port census

The pinned Session51 receipt `sources/regular-eleven-counts.json` supplies
13,333 ordinary balanced ten-classes and 229 self-converse ones. Hence there
are 6,781 merged balanced nodes, comprising 229 self-converse nodes and 6,552
non-self-converse pairs. These class counts are **inherited inputs**, not
predictions from (4) or (6). The new scalar census gives their aggregate
204,019,829 source ports, 67,123,248 internal ports, and 33,561,624 internal
lines, with 22,904 blue internal ports and 11,452 blue internal lines.
The receipt does not attach those aggregate ports to particular nodes.
No internal-degree distribution, loop count, or connectivity at n=10 follows
from these totals alone.

## 9. Reproduction and stopping boundary

Run from this directory:

    python3 check.py --output normal.json
    python3 -O check.py --output optimized.json
    cmp normal.json optimized.json
    python3 verify_manifest.py

The checker uses explicit exceptions, not assert statements. It evaluates both
all-size DPs through n=12; evaluates the direct n-variable degree product and
exhausts every fixed-path word through even n=8; checks native C, native R,
the exact endpoint criterion, blue restriction, unique regular completion,
and transported cycle operation on every balanced word in those cubes; checks
rotation fixed-count polynomials against actual full-matrix rotational
automorphisms on every internal word through n=8; and
checks actual H/Aut, downhill H/Aut, and parent cycles for every balanced
ordinary class through n=6. It retains the n=2 nonfree hostile and the n=6
class/port/loop distinctions. Results contain deterministic universe digests.

Normal and optimized outputs are required to be byte-identical. A copied
relocation must reproduce them. No target tournament atlas, solver, network,
external library, or known count sequence enters the arithmetic engines.
The proofs establish all-size exactness; finite agreement is only a control.
