# Two rooted five-vertex coordinates resolve the next actual port fibres

23 October 2026. **FINITE-EXACT:** the complete inherited class universes at
orders eight and nine; the stated 48-coordinate minimum; and complete native
port fibres for the 201 ambiguous order-nine target classes. **PROVED:** the
rooted-count projection and converse/parity identities below. **OPEN:** any
fixed small signature that determines every tournament metagraph at all sizes.
This packet does not enumerate the full order-nine black multigraph.

## Inheritance and research board

The closest proved mechanism is actual Hamiltonian-path complement with its
source, path and target-isomorphism sidecars. The inherited successful observer
is one round of rooted strong-four incidence, followed where necessary by a
typed wedge; see the [order-eight frontier](../metagraph_order8_frontier_20261023.md).
The canonical hostile is the same-margin, converse-paired port switch. The
corrected near miss is class separation confused with a full wiring census.
The least-used relevant sidecar was the root of an induced five-tournament.

The board was: source/path ports; rooted versus unrooted motif counts; the
integer capacity transportation problem; converse as an involution; and the
older DRT five-census/forest hierarchy. The anchor is port separation, the
niche is the minimum motif menu, and the wildcard is a converse-closed triple.
The applied research card is **Retain local incidence through counting** in
`00-navigation/META-PATTERNS.md`.

The complete order-nine universe is inherited from the independently audited
[blue quotient census](../blue_quotient_20261020b/census.md): 191,536 pairwise
nonisomorphic classes, orbit mass 2^36, and explicit converse witnesses. The
incoming [Session38 wiring](../session38_incidence_ports_electrical_20261023/README.md)
independently authenticates the order-eight carrier. No external literature
or physical analogy is needed for the mathematical claims here.

## 1. The transitive-five shortcut fails before the full motif observer succeeds

Let v4(v) count strong induced four-sets containing v. Define

    R(v) = (v4(v), multiset{v4(w): v -> w}),
    J(T) = multiset{R(v): v in V(T)}.

The complete menu contains nine specified seed choices and zero, one or two
outneighbor refinement rounds: 27 choices per order. Transitive-five seeds
record the root as source, sink, any member, or all five possible ranks.

| Observer | Order-eight signatures / 6,880 | Order-nine signatures / 191,536 |
|---|---:|---:|
| v4, one round | 6,878 | 191,434 |
| v4, two rounds | 6,880 | 191,530 |
| v4 + source TT5, one round | 6,878 | 191,503 |
| v4 + source and sink TT5, one round | 6,878 | 191,520 |
| v4 + all five TT5 ranks, one round | 6,878 | 191,525 |
| v4 + all five TT5 ranks, two rounds | 6,880 | 191,531 |

Thus the transitive-five separator that worked on the previous order-31 DRT
pair does not resolve either old order-eight pair. This failed transfer is
retained, rather than silently broadening "five incidence" after the test.

There are 48 isomorphism types of a five-vertex tournament with one marked
vertex. For each type tau define c_tau(v) to count induced five-sets whose
marked copy at v has type tau. Full rooted-five vectors at round zero,
appended as a vertex histogram to J(T), separate both order-eight pairs and
all 99 order-nine J-collision groups (201 classes). Only these previously
colliding classes need the extra measurements.

The implemented full-five census has an independent lower-order control:

    sum_tau s4(tau, root) c_tau(v) = (n - 4) v4(v).

Here s4(tau,root) counts the strong four-subsets through the root inside tau.
Every rooted strong four-set extends to exactly n-4 five-sets. This is an
all-size identity for n>=5, and it is checked on every probed class.

## 2. Exactly two coordinates suffice within the declared menu

Rooted types are addressed by the smallest upper-triangle binary word among
permutations fixing vertex zero. Pairs are lexicographic, the first bit is
most significant, and a bit one means i -> j. The exact observer is

    S(T) = ( J(T), multiset{ (c_66(v), c_204(v)) : v in V(T) } ).

It is a **joint two-coordinate histogram at one vertex**, appended to J.
It is not two separate coordinate histograms. It does not join the pair to
R(v) at the same vertex, which would be a different, stronger menu.

S distinguishes every ordinary class at orders eight and nine. The 110
distinct class pairs inside the inherited collision buckets give a cheap
complete minimum test. None of the 48 single-coordinate histograms separates
all 110 pairs; exactly 14 of the 1,128 coordinate pairs do. Those pairs of
rooted words are

    (66,204), (72,204), (73,204), (74,204),
    (194,204), (194,205), (201,204), (201,205),
    (204,206), (204,219), (204,456), (204,457),
    (204,472), (204,475).

Thus cardinality two is minimal **within this 48 raw-coordinate menu with
the fixed observer J**. No claim is made about arbitrary linear combinations,
other motif sets, joint refinement with R, or all mathematical encodings.

The chosen motifs have ordinary structural descriptions:

- **66:** the root r beats one vertex z and loses to a transitive chain
  A -> B -> C. The remaining arcs are z -> A, B -> z, z -> C. Its score
  multiset is (1,1,2,3,3), and it has nine Hamiltonian paths.
- **204:** the root loses to a -> b and beats c -> d. The four cross arcs
  form a -> d -> b -> c -> a. Its score multiset is (1,2,2,2,3), and it has
  thirteen Hamiltonian paths.

These named finite motifs, rather than their incidental integer addresses,
are the predicates being measured.

## 3. Actual ports: 165,570 authenticated source-target presentations

For each of the 201 ambiguous target classes, the C++ producer enumerates
**all** Hamiltonian paths. It retains the open ones, normalizes along each
path, and complements the free nonpath arcs. This gives every closed source
presentation that feeds one of these targets. Every presentation occurs
exactly Aut(target) times before quotienting. The resulting scope is:

| Retained object | Count |
|---|---:|
| Target classes | 201 |
| Source classes | 72,727 |
| Literal closed-source presentations | 165,570 |
| Aggregated source/target endpoint/color ports | 154,070 |
| One-round source buckets | 144,425 |
| Units forced by the old capacity peeling | 38 |
| Units left unresolved by the old capacity peeling | 165,532 |
| Units forced after the two rooted coordinates | 165,570 |

Each 36-byte binary row retains source mask, source class, target class,
source-isomorphism permutation, and the actual target Hamiltonian word.
The producer also computes the two-coordinate profile on **every actual
transported target**, checking it against that target's class record.
The independent Python consumer checks all literal matrix isomorphisms,
path arcs, endpoint transport, colors, and converse-reflection capacities.
An independently organized endpoint-aware Held-Karp count proves that every
selected target fibre is complete; the consumer does not canonicalize.

Source endpoint scores (a,b) transport to target scores (9-a,7-b). Reflected
converse ports carry scores (8-b,8-a). Those indices and actual integer
capacities remain in the certificate. The literal port ledger contains
loops and repeated class pairs without deleting either kind of information.

The old ambiguity is real, not only a stalled peeling algorithm. The saved
`baseline_counterfeit` gives a feasible 2-by-2 switch together with its
converse reflection. It preserves source J-buckets, target endpoint/color
capacities, reflection pairing, and even merged black multiplicities, but
changes merged wiring. The two rooted records reject the false assignment.
The switch is a necessary-constraint counterfeit, not a claim that a second
tournament universe realizes the modified graph.

Since S is injective on the **complete** order-nine class bank, measuring S
on any actual order-nine target identifies its class. The focused replay
grounds this on every port entering an ambiguous target. The remaining
134-million-scale full order-nine line inventory was not generated here;
class lookup, focused port completion and full graph enumeration are distinct.

## 4. A minimal converse-closed triple and its integer parity sidecar

Converse maps rooted type66 to457 and fixes type204. Every successful pair
in section2 fails to be closed under converse. Therefore, in the same raw
coordinate menu, the minimum **converse-closed** separating set has three
coordinates: {66,457,204}. It separates because it contains a successful
pair, and no smaller converse-closed successful set occurs in the exhaustive
single/pair menu.

The existing J observer is itself transportable under converse: its global
v4 histogram minus the root's own entry and its outneighbor histogram is
exactly its inneighbor histogram. Thus the additional closure requirement
does not leave an unmentioned orientation loss in the baseline record.

For every tournament size and each vertex, set

    E = c66 + c457,      O = c66 - c457,      S = c204.

Then converse acts by

    (E,O,S) -> (E,-O,S).

The transformation is lossless on its integer image:

    c66 = (E+O)/2,       c457 = (E-O)/2,
    E >= |O|,           E = O modulo2.

The parity condition is essential if these coordinates are treated as an
integer lattice. The proof is the bijection sending each rooted induced
five-set to its converse; the checker independently tests all 5,120 labelled
rooted five-states. This is a concrete even/odd decomposition with an owned
root. It neither identifies tournament converse with an electrical dyadic
lift nor proves a physical equivalence with quantum Hall or color symmetry.

## 5. An older five-census collision is repaired by its lost root

The earlier [HYP-9371, four-loop law](../../hypotheses/HYP-9371-four-loop-law-doubly-regular-tournaments-see-the-five-vertex-census.md)
already supplied D23_1 and D23_2 with identical full **unrooted** five-census,
TT5=1466, and different Hamiltonian counts. This packet does not claim that
collision as new and does not rerun its large Hamiltonian computation.

Counting transitive five-sets at their unique source gives these new rooted
histograms, computed directly from the inherited matrices:

    D23_1: 55^1, 62^4, 63^1, 64^6, 65^10, 66^1;
    D23_2: 55^1, 61^2, 62^3, 64^5, 65^9, 66^3.

Their totals are both1466 but their root distributions differ. The mechanism
is ownership, not an extra unrooted five-motif. This supplies a useful bridge
to the old loop hierarchy: a globally collapsed five-census can lose a
rooted five-record before a six-census is needed. It does not contradict the
earlier statement that an unrooted six-statistic separates those examples.

## Reproduction and scope

From the repository root, using Python3 and a C++17 g++:

```text
python 05-knowledge/results/metagraph_five_ports_20261023/check.py
python -O 05-knowledge/results/metagraph_five_ports_20261023/probe.py
python -O 05-knowledge/results/metagraph_five_ports_20261023/verify.py
python -O 05-knowledge/results/metagraph_five_ports_20261023/closure_probe.py
```

The first command recompiles both producers in a temporary directory,
recomputes the complete menus and focused port ledger, runs the independent
consumer, and compares the frozen hashes and receipts. `--write` is the
explicit ledger-regeneration option. `manifest.json` pins the core inherited
inputs and producer/consumer files. The later `closure_probe.json` freezes
the addendum separately and binds the native port-byte hash. Binary records use little-endian
`<IQQQQ`; gzip is deterministic. No default command modifies the ledgers.

The all-size portable object is the **path-decorated rooted motif**, not a
bare scalar motif count: the path decides which arcs the complement operation
preserves. Rooted incidence signatures and electrical lifts can independently
reject false rewirings, but neither result here makes one determine the other.
The next live obstruction is beyond this complete order-nine class bank:
find realized collisions of the two-coordinate observer under larger native
families, retaining the path decoration and all source capacities.

## Addendum: converse margins repair some collisions, but no single coordinate at nine

Incoming Session38 commit `e5d18e55f4` supplies the all-size sufficient theorem
in `session38_incidence_ports_electrical_20261023/theory/CONVERSE_MARGINAL_REPAIR.md`.
For a signature X, put one edge per ordinary class c between X(c) on the
left and X(c*) on the right. Exact native reflection supplies both endpoint
margin families. If this bipartite **multigraph** is a forest, leaf peeling
determines each class count. Parallel edges are two-cycles and cannot be
discarded. Merely distinguishing the joint pairs (X(c),X(c*)) is insufficient:
a four-cycle can still support a nonzero signed flow.

We apply that theorem to each of the 48 singleton observers

    X_tau(T) = (J(T), histogram{c_tau(v): v in V(T)}).

The complete J-ambiguous banks are closed under converse. Every other class
is singleton under J, and so contributes an isolated edge to this incidence
graph. Thus the focused bank suffices to decide the whole forest condition.

At order eight, **31** singleton menus give a forest. Eight of them still
have class collisions: rooted words

    65, 73, 198, 203, 206, 219, 459, 475.

Their wiring is nevertheless forced by the reflection-implied margins.
This explicitly separates class-identification minimality from wiring
minimality and gives a positive native application of the incoming theorem.

At order nine, **none** of the 48 singleton menus gives a forest. Each has
parallel class edges, and `closure_probe.json` retains the actual classes,
their converses, and the equal left/right signature histograms. Forest
failure alone would not prove ambiguity, because other capacities or
nonnegativity could still force the answer.

The stronger native test succeeds as well: **all 48 menus admit an explicit
feasible reflected black switch**, grounded in the authenticated 165,570
port presentations. Every saved switch:

1. removes only positive actual port counts;
2. preserves source, endpoint, color and X_tau bucket margins;
3. preserves exact target endpoint/color capacities and converse reflection;
4. preserves even merged black multiplicities; and
5. changes at least one merged edge multiplicity.

The two-coordinate observer rejects every switch. Consequently the minimum
of two is also a **wiring-observer minimum in this specified necessary-constraint
model**, even after the incoming converse repair is used. Individual loop
counts are not imposed as additional fixed margins. This does not construct
a second native tournament metagraph or defeat additional retained constraints
such as higher electrical forms. The cheap closure probe checks
1,952,304 conditions, with normal and optimized receipts agreeing, and needs
no new Hamiltonian traversal or full producer replay.
