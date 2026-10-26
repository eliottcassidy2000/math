# Named attachments, square-sum paths, graceful fibers, and Petersen quotients

Session70, 2026-10-26. Base repository commit:
`fb872272125510c317b8114dd3f697b433c595f9`.

**Result.** The common mechanism is **joint attachment compatibility**, not a
common tournament or a coincidence among the numbers 3, 5, and 6. This packet
gives an exact minimal named-pairing carrier in an explicitly unrestricted
forest-gluing category, an exact guard quotient of an inherited graceful
family, and a port-aware explanation of the Petersen-to-K5 contraction.
The ordinary square-sum interval problem is not reopened as an unsolved
conjecture; neither the Graceful Tree Conjecture nor a minimal all-order
tournament response representation is proved here.

## 1. Bounded inheritance and corrections

The inherited square-sum problem is the simple undirected graph

    V(Q_n)={1,...,n},  ab∈E(Q_n) iff a≠b and a+b is an integer square.

Its targets are Hamiltonian paths, their endpoint incidences, and marked
spanning path covers. This identification follows the literal definitions in
the Session32, square_sum_ports, and Session40 packets, not the word “square.”
The separate Session57 problem of nine distinct square entries in a magic
3×3 grid is a different predicate: eight triple sums with a common total.
That session explicitly rejects replacing its hypergraph/decorated-biclique
conditions by pairwise square-sum adjacency. K3,3 occurs there as a decorated
arithmetic carrier, **not** as evidence for a planarity or tournament theorem.

The closest inherited mechanisms, with repository-relative source paths, are:

1. `05-knowledge/results/session32_square_monster_20261023/ports/PORT_TRANSFER.md`:
   named boundary ports, component-terminal pairing, singleton conventions,
   and weighted forest transfer. Its four-port path-versus-cycle hostile
   already shows why degree marginals fail. The present §2 turns that
   obstruction into an all-size **minimality statement for a specified
   continuation language**, not a novel claim for boundary dynamic programming.
2. `05-knowledge/results/session32_square_monster_20261023/square-sum/README.md`:
   endpoint-aware chain reduction and marked degree-two ear correspondence.
   The September25 correction in `01-canon/MISTAKES.md` forbids treating both
   edges at every graph-degree-two vertex as path-forced: an endpoint may
   have graph degree two. Q23 ending at22 is the internal hostile.
3. `05-knowledge/results/square_sum_ports_20261023/README.md`:
   actual endpoint-conditioned birth recurrence and corrected cycle surgery.
   Its later input correction says the originally supplied `3–12` was a typo
   for `4–12`; this overrides the earlier Session40 “unresolved literal” text.
4. `05-knowledge/results/graceful_upper_tail_two_color_connectivity_20261019/UPPER_TAIL_TWO_COLOR_CONNECTIVITY.md`:
   graphic-matroid closure is the component partition; connected color one
   repairs single-coordinate changes. Upper-tail fibers are connected under
   changes in at most two colors. The exponential bound is superseded by
   `graceful_polynomial_diameter_20261019/POLYNOMIAL_STAIRCASE_DIAMETER.md`:
   diameter at most r²−3r+3 for r≥3, at most one for r=1,2, with its cited
   strict-rank theorem dependency. Intermediate abstract tree shapes change.
5. `05-knowledge/results/graceful_fixed_color_counterexample_20261019/ARITHMETIC_FIXED_COLOR_COUNTEREXAMPLE.md`:
   the earlier open arbitrary-fixed-color width-two question was subsequently
   **refuted**, by five actual arithmetic completions split 3+2. Do not inherit
   packet97's old open status as current truth.
6. `05-knowledge/results/graceful_factorial_fibers_20261019/GRACEFUL_FACTORIAL_FIBERS.md`:
   explicit S_k×{±}^k fibers and color-to-component alignment; §4 below
   determines their replacement-legality quotient, which is smaller.
7. THM-2761,
   `01-canon/theorems/THM-2761-graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge.md`:
   bipartite sign gauge converts sums to oriented differences, but the
   mirror factor is indispensable for absolute differences. Labels 0,1,2
   along a path have distinct oriented differences −1,+1 but duplicate
   absolute difference1. This algebraic bridge does not identify square sums
   with graceful labels, and it does not construct a graceful labeling.

Concept board: boundary pairing / graphic closure / colored attachment
assignment / local port transitions / global connectedness / legal grammar.
Anchor: smallest operation-relative carrier. Niche: graceful guard quotient.
Wildcard: Petersen's explicit K5 and K3,3 maps.

## 2. Exact minimum for unrestricted named forest attachments

**PROVED (self-contained; no priority claim).** Fix r≥1 and 2r distinct named
boundary vertices B. A state is a selected forest consisting of exactly r
nontrivial vertex-disjoint paths, with B precisely its endpoints. All selected
edges are mandatory and remain present. No operation touches an interior
vertex. An allowed connector is a fresh internally disjoint path between
two distinct still-exposed endpoints; choosing length two always avoids
parallel-edge ambiguities in simple graphs. The intermediate guard accepts
exactly when the result remains a linear forest. Final acceptance means that
the union of the selected forest and the chosen connectors is one spanning
path. These are **arbitrary external connectors**, not necessarily Q_n edges.

Let M be the perfect matching on B pairing the two ends of each component.

**Theorem A.** For these operations and Boolean acceptance:

1. M is a sufficient, recursively closed signature.
2. Two states are indistinguishable under every future legal operation and
   final test iff they have the same M.
3. There are exactly (2r−1)!!=(2r)!/(2^r r!) distinguishable states. Any
   deterministic exact signature for this semantics has at least that many
   values, hence at least ceil(log₂((2r−1)!!)) bits in a fixed-length encoding.

**Closure and sufficiency.** A connector ab is legal iff a,b belong to
different matched pairs. If M pairs a with a′ and b with b′, delete those
two pairs and replace them by a′b′. Endpoints a,b are consumed. This exactly
describes merging the two paths; length and shape of their interiors have no
effect on the stated guard. A connector between partners instead closes a
cycle and is rejected. Induction handles all future connector words and
terminal spanning-path tests.

**Necessity, already at one guard.** If M≠N, choose ab∈N\M. The fresh path
a–w–b closes the N-component but joins two M-components, so its legality
distinguishes them. This also works when the local component is a single
edge: the fresh w gives a triangle rather than a duplicated edge.

**Necessity using final tests alone.** After joining the two M-components,
append the other r−2 components in any order by r−2 fresh connectors. The
M-state becomes one spanning path. In the N-state the ab component is already
a closed cycle, untouched by subsequent connectors, so the final union is
not a spanning path. Thus distinct pairings are separated even without
observing an intermediate rejection. The r=1 case has one matching and no
pair to distinguish. At r=0 the empty state is a separate base convention.

The growth is log₂((2r−1)!!)=r log₂ r+O(r). This is an **alignment memory
lower bound**, not a lower bound on every possible algebraic algorithm,
rank-based compression of a family of states, or time complexity. A weighted
family needs its response multiplicities; matching one particular forest
does not reconstruct its hidden interiors or count them.

**Sharp four-port hostile.** M=(01)(23), N=(03)(12), connector1–4–2.
Both states use the same four endpoints and have two one-edge components.
The first union is the path0–1–4–2–3. The second is edge0–3 plus triangle
1–4–2–1. Degree and component-size marginals before gluing do not suffice.

**Boundary of the theorem.** This is not a theorem that every pairing is
realizable within a given Q_n, or that its distinguishing connectors have
square sums. A restricted host/attachment language can identify more states.
Singleton local path components have coincident endpoint vertices and are
excluded from Theorem A's 2r-distinct-endpoint universe. A general carrier
must retain endpoint **occurrences**, their underlying-vertex coincidences,
and vertex capacities; the inherited transfer and Session69 P3 do so.
If cycle closure is an accepted *terminal* operation, it must be distinguished
from legal intermediate forest attachment: closing a pair is accepted only
when that single path exhausts the required current component universe.

## 3. What transfers to the actual square-sum graph

The exact inherited insertion rule, for H with at least two vertices and a
new vertex v with neighbor set N, is

    h(H+v)=Σ_{a∈N} E_H(a)+Σ_{{a,b}⊆N} C_H(a,b).

Here h counts unoriented Hamiltonian paths, E_H(a) counts paths ending at a,
and C_H(a,b) counts spanning **two-path** covers with a and b exposed on
different components. Singleton components are allowed, with their actual
coincidence convention. Deleting v is a bijection: endpoint v yields the
first term; interior v yields the second. No extra factor two occurs.

For growing Q_n the neighbor set is not freely chosen:

    N(n)={s²−n: n<s²≤2n−1},
    |N(n)|=floor(sqrt(2n−1))−floor(sqrt(n)).

Theorem A explains why endpoint/component incidence is a natural sufficient
sidecar, but its unrestricted minimality does not prove minimality for this
arithmetic insertion language.

**FINITE-EXACT, independently replayed here:** path counts through25 are
1 at n=1,15,16,17; 3 at n=23; 10 at n=25; zero at every other n≤25.
The Q23 path ending at22 has graph-degree-two endpoint22 and omits14–22.
Q25 has a path, but contracting the inherited tree blocks

    A={2,7,9,16,18}, B={4,5,11,12}, C={1,3,6,8,10,17}

gives a 13-vertex, 20-edge multigraph without a Hamiltonian path. Deleting
its three block vertices leaves six components, violating c(G−S)≤|S|+1.
The source has 32 edges, 25 vertices, and cycle rank8; tree contraction
preserves that cycle rank but not traceability. Block A is branched at7,
so a Hamiltonian path cannot be modeled as one visit to every whole block.
This is the direct precursor of revisited tournament children in Session69.

**CITED, not independently rebuilt:** Gerbicz's own January21,2018 comment
at [OEIS A090461](https://oeis.org/A090461) states paths for every n≥25 and
cycles for every n≥32. The comment and its scope were rechecked on2026-10-26;
this packet does not audit the full deterministic construction/glue tables.
The inherited cycle-surgery lemma then makes each individual square-sum
edge optional among Hamiltonian paths once n≥max(32,a,b). Simultaneously
prescribed edges, prescribed endpoints, and exact counts are stronger targets.

## 4. Graceful trees: retain colors and distinguish existence from motion

A gracefully labeled m-edge tree on the fixed interval0,...,m is a spanning
tree choosing exactly one edge of every color |a−b|∈{1,...,m}. The Graceful
Tree Conjecture asks for such a labeling **for every prescribed abstract
tree shape**. Connectivity of existing labeled rainbow trees under operations
which change their shapes is a different assertion.

For a fixed forest F and candidate edge set disjoint from F, contract its
components but retain every candidate edge as a distinct colored contact,
including loops and parallel contacts.
Then F∪E is a forest iff the contact multigraph E is a forest; adding E
preserves connectedness exactly through that quotient. The proof contracts
each tree component, preserving the cycle rank; loops and parallel cycles
cannot be silently deleted. Thus the exact component partition is the
graphic-matroid closure sidecar, while the color-to-edge choice is a separate
rainbow constraint. Merely remembering how many components there are loses
which colors can join which components.

### 4.1 The actual arithmetic counterexample to arbitrary frozen colors

On0,...,12, leave colors6,9,11 free and fix

    1:5–6, 2:1–3, 3:8–11, 4:3–7, 5:0–5,
    7:4–11, 8:1–9, 10:0–10, 12:0–12.

The components are A={0,5,6,10,12}, B={1,3,7,9}, C={2}, D={4,8,11}.
Nonloop supports are6:{AD,CD}, 9:{AB,CD}, 11:{AB,AD}. The available contact
tree is B–A–D–C. Any completion must occupy its three links, giving two
assignments, differing in all three free colors. One has3 realizations and
one has2. Direct enumeration of all7·4·2=56 actual free-edge tuples gives
exactly5 completions, whose width-two graph is K3 disjoint union K2.
Therefore a source and target in different components require temporarily
changing an already-correct frozen edge if one uses the global width-two
connectivity theorem. This is not an exponential diameter bound, not a
minimal-m assertion, and not a counterexample to the Graceful Tree Conjecture.

### 4.2 A new exact guard quotient of the inherited factorial family

Fix k≥1, m≥k²+2k, b_i=i(k+1). Freeze every color d>k as edge0–d except
that color b_i uses i–(b_i+i). The fixed forest consists of one tree H and
the k isolated vertices b_i. Each free color d∈{1,...,k} has precisely two
nonloop realizations on each spoke H–b_i: (b_i−d)–b_i and b_i–(b_i+d).
All other color-d edges are loops after contraction.

The inherited state bijection is

    literal completions ↔ (π,ε)∈S_k×{−,+}^k,

where π(d)=i records which isolated vertex color d attaches. There are
k!2^k states; one-color reconfiguration has k! components, each a k-cube;
width-two changes connect them by swapping attachment assignments.

Define a replacement request as: overwrite a specified color coordinate by
a specified actual edge of that color and test whether the result is still
a graceful completion containing F. **Overwriting with the existing edge is
allowed.** For batches, simultaneously overwrite at most two coordinates.
The observation is only validity and future such validity, not the literal
output edge set or the predicate “did this move change the state?”
On a successful request, the state is the actual overwritten selection;
on rejection it remains unchanged. A batch specifies distinct coordinates.

**Theorem B (PROVED).** In this guarded overwrite language the exact
behavioral quotient has k! states, indexed by π; each equivalence class has
2^k different sign realizations. One-color legality already reconstructs π,
but all finite one-/two-color guard words fail to recover ε.

**Proof.** After any overwrite, legality is equivalent to the selected
contact spokes assigning exactly one color to every b_i; left/right sign
does not affect this test. States with the same π have the same validity
for every request and transition to the same next π whenever it is legal.
Induction gives agreement on every finite future request word. Conversely,
delete color d mentally and overwrite it by an edge of color d on H–b_i.
This is valid iff i=π(d): otherwise its previous leaf is unattached and
another leaf is duplicated. Thus single requests determine each π(d).
The state count and fiber size follow. QED.

This makes two distinct points: alignment has factorial memory, **and**
guard reconstruction is operation-relative—these guards do not reconstruct
the full labeled tree. Forbidding no-ops or observing the literal output
would change the language and can expose the signs; that is not this theorem.
For unrestricted tree-labeling problems, vertex labels, permitted colors,
the fixed-forest partition and the desired abstract shape remain distinct
coordinates. No tournament is manufactured from this symmetric relation.

## 5. Petersen, K5, K3,3: exact maps and the lost transition relation

Use vertices u_i=i and v_i=i+5, i modulo5, with outer edges u_i u_(i+1),
inner edges v_i v_(i+2), and spokes u_i v_i.

### 5.1 Contracting the five spokes gives K5, but not Hamiltonicity

Each spoke cluster has two vertices and four external half-edge ports:
two outer ports at u_i and two inner ports at v_i. After contracting the
spokes, outer edges supply the5-cycle and inner edges supply its complement,
so the quotient is K5. The **uncolored** quotient has12 unoriented
Hamiltonian cycles, although Petersen has none.

The exact local spanning-2-factor response is especially small. For a
selected set J of quotient edges let d_O(i),d_I(i) be its selected outer
and inner degrees. Its lift is a2-factor iff

    d_O(i)=d_I(i)∈{1,2} at every i.

The spoke is included iff this common value is1. If the value is1, the
spoke joins the one selected outer port to the one selected inner port:
one run through the cluster. If the value is2, the spoke is omitted and
the two outer ports are paired through singleton u_i, while the two inner
ports are paired through singleton v_i: **two runs through the cluster**.
The original endpoint degree equations d_O+x=d_I+x=2, x∈{0,1}, prove both
directions and the uniqueness of the lift. Global connectedness is still
needed to distinguish a Hamiltonian cycle from a disconnected2-factor.

A simple K5 Hamilton cycle using one visit per cluster could lift only if
outer and inner edges alternated around its odd length5, which is impossible.
This argument excludes only one-run lifts, **not all possible Petersen
Hamiltonian cycles**. The latter can a priori omit spokes and revisit the
same contracted cluster. The full response rule above keeps those cases.

### 5.2 Complete independent non-Hamiltonicity proof

In any cubic graph a spanning2-factor is the complement of a perfect
matching. In Petersen a perfect matching uses an odd number k of spokes:
the remaining5−k outer vertices must be matched among themselves. Thus
k=1,3,5. If k=3, its two unspoked indices would have to be adjacent both
in the outer5-cycle (difference±1) and inner5-cycle (difference±2), impossible.
For k=5 there is one matching. For k=1 at index i, deleting i from each
5-cycle leaves a4-vertex path with a unique perfect matching. Thus there
are exactly six perfect matchings.

The all-spoke complement is the two original5-cycles. By rotational
symmetry it remains to inspect the one-spoke matching at0:

    M={0–5,1–2,3–4,7–9,6–8}.

Its complement has cycles0–1–6–9–4–0 and2–3–8–5–7–2. Hence all six
2-factors are C5⊔C5, not one cycle. This proves non-Hamiltonicity without
assuming one-visit quotient behavior. Petersen nevertheless has a path
0–1–2–3–4–9–6–8–5–7; the exact checker counts120 unoriented paths.

The two-port quotient rule enumerates exactly these same six2-factors:
one uses all10 K5 edges and two runs in every cluster; the other five
use six quotient edges and two runs in one cluster. Degree-correct local
responses still need the global component pairing/connectivity test.

### 5.3 K3,3 as a literal subdivision inside Petersen

Delete vertex0. The degree-two vertices1,4,5 lie on the paths

    2–1–6, 3–4–9, 7–5–8.

Suppressing those paths produces K3,3 with parts {2,8,9} and{3,6,7}.
The three suppressed edges must remain **marked** if asking whether a
quotient Hamiltonian cycle visits every original vertex: it must use them
all. Merely being a K3,3 Hamiltonian cycle need not meet this condition.
This is the same marked-ear rule recovered from Q25, not a new analogy.
The independent audit sharpened this to an exact count: of the six K3,3
Hamilton cycles, exactly two use all three marked edges. Expanding these
gives precisely the two spanning9-cycles of Petersen−0. The producer verifier
now independently checks this equality of expanded cycle edge-set families.

### 5.4 Planarity and tournament predicates must stay separate

K5 is nonplanar because a simple planar graph with n≥3 has m≤3n−6, while
10>9. K3,3 is nonplanar because a simple bipartite planar graph with n≥3
has m≤2n−4, while9>8. These familiar Euler bounds suffice for these exact
graphs; the packet does not rely on a fresh audit of the full Kuratowski
characterization. Petersen's K5 minor and K3,3 subdivision both certify
nonplanarity. Deletion/contraction preserve planarity in the forward
direction; they do not preserve Hamiltonicity.

K5 and K3,3 have Hamiltonian cycles (12 and6 unoriented ones respectively).
A planar star K1,3 has no Hamiltonian path. Nonplanarity therefore is not
a Hamilton-path obstruction and planarity is not a sufficient condition.
Likewise, gluing three planar stars along the same three leaf names produces
K3,3: “each piece planar” is not an adequate gluing signature.

The same two graphs also separate **gracefulness** from those properties.
K3,3 has a graceful labeling with left labels{0,1,2} and right labels{3,6,9}:
its nine edge differences are exactly1,...,9. K5 cannot be graceful. Its
ten differences would have to be1,...,10, hence five odd differences. If a
of its five vertex labels are even, it instead has a(5−a) odd differences,
which is0,4,or6, never5. This elementary parity proof is independently checked
by all84 normalized candidate label sets containing0 and10. These are
non-tree controls, not counterexamples to the Graceful Tree Conjecture.

A tournament has exactly one directed arc on **every** distinct vertex pair.
Q_n, Petersen, K3,3 and difference-colored graceful choices are symmetric
graphs with missing edges; they are not tournaments. Completing and orienting
their missing pairs changes Hamilton-path feasibility. Indeed every finite
tournament has a directed Hamiltonian path, while Q14 does not. The genuine
shared construction is a typed response relation retaining simultaneous
boundary assignments—not an orientation trick.

## 6. Connection contract and research consequence

| Source / map | Preserved predicate | Destroyed information / required sidecar | Decisive test |
|---|---|---|---|
| Selected linear forest → endpoint matching | All future unrestricted path-attachment feasibility | Interior shape is intentionally forgotten; weights needed for counting | Four-port path versus path+cycle; Theorem A |
| Q_n → named block forest responses | Hamiltonian paths, with multiplicities | Whole-block supervertex loses repeated visits, port names and pairing | Q25 tree-block quotient fails cut budget |
| Graceful fixed forest → colored contact multigraph | Acyclic/rainbow completion with frozen edges | Dropping colors, loops, multiplicity or actual component identities breaks predicate | m12 five-state fiber splits3+2 |
| Factorial graceful state → π | All overwrite-legality words | Left/right edge realization is genuinely unobservable for this language | 2^k literal states share every guard; Theorem B |
| Petersen → K5 via spokes | Existence of minor / nonplanarity certificate | Hamiltonicity needs outer/inner port assignment, one/two runs, global connectivity | 12 bare quotient cycles, no alternating odd cycle, six disconnected2-factors |
| Petersen−0 → K3,3 via suppression | Cycle correspondence with expansion | Spanningness needs marked subdivided edges | A cycle omitting a mark skips its subdivision vertex |
| Tournament substitution → Session69 P3 responses | Ordinary literal-target incidence | Ordered endpoint occurrences, weights, child revisits and stabilizer action cannot be discarded | Inherited target codes81/83 from same outer endpoints |

**Positive synthesis.** Composition is performed by summing over jointly
compatible boundary states, not by multiplying unrelated marginals. A small
carrier is exact only after its allowed tests are specified. For unrestricted
forest gluing, legality recovers the whole boundary pairing. For graceful
overwrite gluing it recovers the color/component assignment but not its
sign realization. For the Petersen matching contraction, the two local
vertices and their one-/two-run alternatives are enough to restore exact
2-factor incidence, with global connectedness checked separately. This is
the precise transferable mechanism supporting the tournament investigation;
none of these statements by itself proves tournament reconstruction at all
sizes or identifies the smallest tournament state algebra.

## 7. Verification and primary-source scope

Run from this directory, standard Python3 only:

    python3 verify_graph_alignment.py > /tmp/graph70.json
    python3 -O verify_graph_alignment.py > /tmp/graph70-opt.json
    cmp receipt.json /tmp/graph70.json
    cmp /tmp/graph70.json /tmp/graph70-opt.json

The `--hostile` option intentionally exits nonzero by asserting that all12
bare K5 Hamilton cycles lift. This rejection remains active under `-O`.
There are no Python `assert` statements or external packages. Universe:
all endpoint matchings through r=5, all ordered distinct pairing/context
tests through r=4 (11,136), all Q_n Hamilton paths for1≤n≤25, all56 m12
free-edge tuples, all factorial-family states through k=4, all1,024 K5 edge
subsets for the local lift rule, and all Petersen perfect matchings.
The checker also verifies the explicit graceful K3,3 and all84 normalized
candidate label sets for K5, independent of the parity proof.
These finite controls complement—not replace—the proofs of Theorems A/B.

Primary web material checked on2026-10-26:

- [Gerbicz's own OEIS comment](https://oeis.org/A090461): precise all-n
  path/cycle existence claim, only as **CITED** above. Full glue proof not rerun.
- Letzter–Pokrovskiy–Williams,
  [On the gracesize of trees, v1](https://arxiv.org/html/2511.11331v1),
  Theorem1.3 and Definition1.5/Conjecture1.6: asymptotic gracesize bound
  and exact difference-colored rainbow reformulation. Their result does
  not settle exact gracefulness of every tree; this packet claims none.
- Merino–Mütze–Namrata,
  [Kneser graphs are Hamiltonian, v4](https://arxiv.org/abs/2212.03918),
  revised2025-02-24, identifies Petersen K(5,2) as the exceptional graph.
  The elementary proof here is independent of their theorem.
- Grinberg–Stanley,
  [The Rédei–Berge symmetric function of a directed graph, v1](https://arxiv.org/html/2307.05569v1):
  tournament Hamilton-path existence/parity belongs to the total-orientation
  category. No claim of a new Rédei result or graph-to-tournament equivalence.
- Original Kuratowski scan retrieval returned HTTP403; a full primary proof
  was not inspected. No theorem here depends on importing its complete
  characterization. The explicit obstructions are proved by Euler bounds.

No claim of historical priority is made for the boundary-state or transition
system mechanisms. The next justified mathematical question is **which
restricted allowed attachment languages still separate all pairings or
other joint assignments**, not whether a repeated small numeral proves a
universal bridge.
