# Independent audit: Session66 tournament graph-pair codecs

**Audit date:** 2026-10-26. **Verdict:** the mathematical claims in the sealed
producer report pass, with the decoder contract explicitly restricted to valid
codec images. Root has now integrated and resealed that clarification; its exact
documentation-only delta was independently verified. No counterexample to the claimed faithful encodings, exact metric
transport, no-go scope, or clone interfaces was found. Two documentation actions
are supplied below: a valid-image decoder clarification and a historical THM-209
truth correction. Neither requires changing either codec theorem.

This worker wrote only `session66-audit-graph-pairs/`. The producer packet and
`session66-repo/` were read-only for this audit; no source edit, commit, fetch,
push, or publication was performed. All-order conclusions below are proof
checks, not extrapolations from the finite census. No novelty claim is made.

## 1. Scope, pins, and verification independence

Audited producer files: `REPORT.md`, `check.py`, `normal.json`, `optimized.json`,
`SOURCE_PINS.json`, and `SHA256SUMS` in `session66-graph-pairs/`.
The initial sealed report SHA256 is
`e073dda07857953759b6c09a4a023a23b619d83c77f0a20d3c6e3e4b4ff91e8a`;
checker SHA256 is
`92e4f05295056910508f30f261994ecb0a51e8bd195ca41197a7ef9d27aa2f4f`.
The final clarified producer report SHA256 is
`b58d578afb6e8e3254eb02b7a398b67785d5a9f1f2dbd2f27d1dec71dd0ee7d8`.
Its only change is the agreed valid-image decoder paragraph; original bytes are
preserved in `producer-report-original.md`, and `producer-report-clarification.diff`
records the delta. Code and numerical receipts are unchanged.
Every entry in the original and final producer seals verified. The twelve copied source files were
independently checked three ways: SHA256, Git blob hash calculated from raw
bytes, and `git show` at source commit
`63f52936f8dd8592d4f2fee869748baad2e1997e`.

The assigned fixed read-only pin takes precedence over scripts that would
mutate/fetch the source checkout. `AGENTS.md`, startup navigation, relevant
current frontier and guardrails, research protocol, core papers, corrected
spectrum note, OCF source, homogeneous-module theorem, Session62 clone note,
and Session65 general signed-pair proof were inspected.

The producer checker was rerun normally and under `-O`; both outputs are
byte-identical to each other and to the sealed producer `normal.json`.
The independent checker imports **no producer routines**. Its separate paths are:

- codec construction from edge lists and arm/leaf banks;
- full-image recognition by local recovery followed by exact reconstruction
  of every vertex and edge;
- generic colored graph isomorphism by false-twin quotienting and weighted
  quotient backtracking, without calling codec decoders;
- literal odd-cycle enumeration by minimum-rooted DFS;
- Hamiltonian counts by memoized suffix recursion;
- clone classes compared against an independent exhaustive definition of
  all transitive modules, rather than reconstruction alone;
- strong connectivity by Boolean reachability closure.

Generic graph isomorphism was itself checked against brute permutation search
on every simple graph through four vertices and every four-state symmetric
relation on three vertices: **4,201 graph-pair comparisons and 140 full
isomorphism-group counts**. False-twin compression is exact: identical adjacency
rows define independent interchangeable vertices; an isomorphism permutes these
classes preserving their sizes and quotient colors, and each quotient map lifts
in the product of the class factorials many ways.

## 2. Finite-exact receipts

All gates are explicit exceptions and remain active under `python3 -O`.

| Independent universe | Checks and outcome |
|---|---|
| Every labelled tournament, orders 0–5 | 1,100 tournaments; 4,400 original/scrambled strict codec round trips; all pass |
| Tournament classes, orders 0–5 | Counts `1,1,1,2,4,12`; 188 same-order class-pair comparisons over both codecs; 42 scrambled generic-graph positive controls; all pass |
| Every codec of each small class | 42 full graph automorphism counts, not just class distinctions; exact `Aut(T)` for joint and `12^n |Aut(T)|` for Γ |
| Every labelled tournament, orders 0–5 | 1,100 OCF comparisons using literal-cycle DFS versus Hamiltonian recursion; all pass |
| Every labelled tournament, orders 0–5 | 1,100 clone-component comparisons against all maximal transitive modules and weighted reconstructions; all pass |
| All 21 small tournament classes | Canonical weighted-quotient automorphism counts equal tournament automorphism counts |
| Every labelled tournament, orders 3–6 | 33,864 separate spectrum tests; 7 and 21 absent; complete strong spectra agree with the producer |
| Every tournament and permutation, orders 0–4 | 3,180 codec metric identities; minima over every tournament for all 34 vertex permutations equal the number of even vertex cycles |
| One-edge malformed images, n=1,2, both layers where applicable | Joint: 226 rejections; Γ: 191 rejections, including extra-isolate controls; both valid two-edge arc-reversal controls accepted |
| Empty/small graph boundaries | Every tournament at n=0,1,2 checked explicitly; ten nonempty all-isolate graphs rejected by strict recognition |
| C3 no-go | Exactly masks 0 and 7 among the eight three-vertex simple graphs are C3 invariant |
| Genuine clone-gap hostile | `C3[TT2,TT2,TT2]`: six vertices, three automorphisms, gap polynomial `1+t+t²+t³`, only fixed toggles `0` and `111` descend |

The independent strong spectra are

```
3: 3
4: 5
5: 9,11,13,15
6: 15,17,19,23,25,27,29,31,33,37,41,43,45.
```

At `n=0`, both graphs are empty, both automorphism groups are trivial,
`H=1`, and the empty clone quotient reconstructs the empty tournament.
At `n=1`, the joint code is `(R,empty_7)` with automorphism order one;
Γ has nine vertices and eight edges with automorphism order 12.
At `n=2`, the joint code has 14 vertices and layer edge counts 12 and one;
Γ has 18 vertices and 17 edges with automorphism order 144. Both tournament
orientations are accepted and are isomorphic. The boundary conventions are
therefore compatible with all claims, including the kernel formula.

## 3. All-order proof audit: 7n joint codec

**Accepted for every finite tournament, including n=0.** The seven-vertex
frame is a tree with a unique degree-three center and three arms with distinct
lengths 1, 2, 3. Every automorphism fixes the center, then each arm, then each
vertex. Thus every layer-one isomorphism is exactly a component permutation
with a unique internal extension. It necessarily maps the intrinsic tail and
head roles correctly. Layer two preserves precisely the directed tournament
incidences, so the component permutation is a tournament isomorphism iff it
preserves both layers. Every such tournament map extends uniquely.

Consequently the induced correspondence is a bijection on **all isomorphisms**,
not only existence, and `Aut(G1,G2) ≅ Aut(T)`. Order of layers and a common
vertex identification are indispensable. The carrier and edge counts are
exactly `7n` and `6n + binom(n,2)` across layers. All nonport vertices, five per
block, are isolated in layer two; some ports can be isolated as well.

The natural construction commutes with relabelling by the lifted permutation
`(v,r) -> (σv,r)`. A decoder intrinsically recovers the tournament on the set
of frame components; that recovery is relabelling-equivariant. Numbering those
components by current integer labels produces a representative, not an extra
canonical-labeling algorithm.

**Image recognition strengthening.** A full recognizer can check that every
frame component is exactly the rigid tree, recover every role, require exactly
one of the two possible cross-port edges for each pair of distinct blocks,
and forbid all other layer-two edges. This is necessary and sufficient for
membership in the joint-code image. The independent `strict_decode` implements
these conditions by rebuilding and comparing the entire colored matrix.

## 4. All-order proof audit: 9n single graph and independent pair

**Accepted for every finite tournament.** All six leaves per block have degree
one; every core has degree at least three, including at n=1. Leaf count at a
core is intrinsically one, two, or three, distinguishing center, tail, and head.
Every center's two core neighbors are its unique matched tail and head. Thus
every graph isomorphism recovers the same matching and induces a tournament
isomorphism. Every tournament isomorphism extends. Hence Γ is faithful and
`([Γ(T)],[empty_(9n)])` is faithful even under independent unlabelled graph
isomorphisms and an equal-order requirement.

The automorphism action on recovered centers is surjective onto `Aut(T)`.
Its kernel is exactly independent permutations of the two tail leaves and
three head leaves per block. Choosing the construction's leaf coordinates
splits the extension; the quotient permutes the block factors. Therefore

```
Aut(Γ(T)) ≅ (S2 × S3)^n ⋊ Aut(T).
```

This group has order `12^n |Aut(T)|`; for **every nonempty** T it has even
order. It is not the tournament automorphism group, although the latter is
recovered as its center-action quotient. No contradiction follows from the
odd-order theorem for tournament automorphisms.

Construction and decoding are quadratic in adjacency-table size parameter n.
A graph canonical-form oracle followed by intrinsic decoding gives a canonical
tournament representative. This is a reduction with an explicit oracle, not
an algorithmic breakthrough or asserted complexity equivalence. Γ has `9n`
vertices and `8n+binom(n,2)` edges. A full strict recognizer additionally checks
all node coverage, the matching, the tournament condition, and absence of
unexpected edges; rebuilding the recovered Γ is an exact test.

## 5. Hostiles and scope of the C3 obstruction

**Order loss, minimally n=3.** Relative to `0<1<2`, reversing the adjacent
arc `0->1` in the transitive tournament leaves a transitive tournament; reversing
the long arc `0->2` gives a directed triangle. The backward graph is one edge
plus an isolate in both cases, and its complementary forward graph follows it
under the same isomorphism. Thus even joint forward/backward layers fail when
the total order is discarded. At n<=2 there is only one tournament class, so
this class-level failure cannot occur earlier.

**Independent alignment loss, minimally n=4.** Reversing the tournament swaps
all tails and heads in layer two, which is an ordinary graph isomorphism;
layer one is independent of the tournament. Thus independent layer classes
cannot distinguish a tournament from its converse. `source ⊕ C3` and
`C3 ⊕ sink` have score multisets `{3,1,1,1}` and `{0,2,2,2}`, so are
nonisomorphic. All tournament classes through order three are self-converse.
The global tail/head swap is not a frame isomorphism, which explains why the
joint code survives. The pair `(Ω,J)` also fails here: each Ω is K1 and each J
is the edgeless four-vertex graph.

**Same-carrier no-go.** Assume a map on every labelled tournament is equivariant
under all relabellings and outputs a fixed number of named simple-graph layers
on the original vertex set. On a directed triangle, C3 fixes the tournament and
is transitive on all unordered vertex pairs. Each layer must therefore be
empty or complete, hence is fixed by every transposition. A transposition
reverses the triangle, so equivariance already forces the two different labelled
orientations to have identical outputs. No exact decoder can distinguish them;
in particular no equivariant exact decoder exists. The report's conclusion is
correct (indeed its decoder equivariance assumption is stronger than necessary).

This does **not** rule out an injection of unlabelled tournament classes into
same-size graph-class pairs, arbitrary symmetry breaking, a relabelling action
that also permutes unnamed layers, or any enlarged carrier. None of those
unclaimed statements follows from this proof.

## 6. The actual 7/21 bridge and literal-cycle correction

**CITED + CHECKED:** the external OCF input is Irving–Omar Corollary 20, with
cycle definition in §2.2: cycles are cyclic-shift classes of directed closed
simple walks. A permissible permutation is exactly a collection of disjoint
nontrivial odd directed cycles plus fixed points. Its weight is `2` for each
nontrivial cycle. These collections are independent sets in the literal
cycle-intersection graph, which gives `H(T)=I(Ω(T),2)` by reversal of
Hamiltonian paths. [Primary source](https://arxiv.org/html/2412.10572v3)

The primary digraph complement includes loops; removing them gives the
converse tournament. Loops are irrelevant to Hamiltonian paths. This repairs
a minor literal equality in the inherited THM-002 source without changing any
Session66 conclusion. Its source's second “internal proof” should not be
advertised as independent of OCF; the accepted all-order input here is external
proof 1, not a circular reuse of the identity.

**CITED + CHECKED:** Busch Corollary 1 supplies the exact strong minimum stated
in the report. The piecewise values start `3,5,9,15,25,45` and obey
`f(n+3)=5f(n)`; the displayed residue-class formula is correct. Thus a strong
7-count witness would have at most four vertices, and a strong 21-count witness
at most six. The finite strong spectra exclude these cases. Strong components
are totally ordered by dominance, so Hamiltonian paths concatenate uniquely
and their counts multiply. A product 7 or 21 requires a strong factor 7 or 21;
therefore neither value occurs at any order.
[Primary source, Corollary 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf)

These facts do not classify the full positive-odd spectrum. Completeness
remains **OPEN**; `H=63` refutes a `7·3^k` forbidden tower. The corrected
source's retraction of a novelty claim is retained. The constants 7 and 21
are ordinary Hamiltonian-path counts, not forbidden tournament orders or
forbidden automorphism data.

For general graphs, direct counting gives `I(K3,2)=1+6=7` and
`I(P4,2)=1+8+12=21`. Hence neither graph is an entire Ω-image. No conclusion
about induced subgraphs or arbitrary components follows without a separate
argument. The independent pair `(Γ,Ω)` is faithful because its first
coordinate reconstructs T; Ω separately exposes the statistic. It is not a
polynomial-size claim: for odd n the expected number of Hamiltonian directed
cycles in a random tournament is `(n-1)!/2^n`, so some tournament has at least
that many Ω vertices. This is already superexponential along odd n.

### Minimal support-collapse witness

The independently found five-vertex tournament with lexicographic upper-pair
bitmask 40 has four directed triangles and three distinct directed 5-cycles
on its one full support. Its seven literal cycles all intersect, so Ω is K7
and `H=15`. Collapsing equal supports leaves K5 and gives `I(2)=11`.
See `normal.json` and `PROPOSED_CORRECTIONS.md` for exact arcs and cycle words.
This proves that the still-present “as vertex sets” wording in historical
`THM-209-independence-polynomial.md` is genuinely false, not harmless notation.
The current source HEAD inspected for this file still matched the pinned bytes.

## 7. Exact edge metric and signed-permutation minimum

For **every** finite tournament T and **every** vertex permutation σ, lift σ
while preserving every local gadget coordinate. All frame, matching, and leaf
edges are preserved. Each unordered original pair has exactly two candidate
port edges, of which exactly one is present. A direction disagreement swaps
which one is present, producing precisely two symmetric-difference edges.
Different original pairs have disjoint candidate pairs. Thus for either codec

```
|E(encoded(T)) Δ E(σ encoded(T))| = 2 Dσ(T).
```

For the joint code this concerns layer two; layer one adds zero. Image versus
pullback uses σ or σ^-1, but the counts agree by reindexing. This is not an
optimization over arbitrary graph permutations. In particular, the even
leaf-permutation kernel of Γ cannot invalidate the statement, which fixes a
specific lift of the original vertex permutation.

The inherited minimum also has a checked all-order proof. On an unordered-pair
orbit, orientation transport is a signed cyclic bit system. Positive return
allows zero disagreements. Negative return forces an odd number and permits
one. These orbit systems use disjoint edge variables, so their minima add.
A pair orbit has negative return exactly when its endpoints are antipodes on
one even vertex cycle: different vertex cycles cannot exchange endpoints,
and on a length-d common cycle exchange requires displacement d/2. There is
one such pair orbit per even vertex cycle. Hence

```
min_T Dσ(T) = number of even vertex cycles of σ,
min_T |E(encoded(T)) Δ E(σ encoded(T))| = 2 e(σ).
```

This is a **constrained-image, prescribed-action** theorem. Two transpositions
have even permutation sign but minimum defect two, hence graph difference four.
Arithmetic interpretation still needs the actual arithmetic set and its actual
Frobenius action; none is produced merely by these graph codes.

Every tournament automorphism group has odd order: an even-order finite group
would contain an involution, whose swapped vertex pair would have its arc
reversed. This group-theoretic implication uses Cauchy's theorem (or inverse
pairing to obtain an involution). It applies exactly to the joint group and to
the Γ quotient action, not to Γ's leaf kernel.

## 8. Clone interfaces and preserved action

The inherited module theorem is correct for every finite tournament. A vertex
has at most one outgoing and at most one incoming homogeneous neighbor: two
on either side would distinguish one of the supposed pairs. Hence components
are paths or consistently directed cycles; successive exterior agreement
excludes a cycle and forces all earlier path vertices to dominate all later
ones. Every exterior vertex sees the whole path uniformly, so each component
is a transitive module. Conversely, consecutive vertices of a transitive
module are homogeneous in the full tournament, so every such module lies in
one component. Thus the components are exactly maximal transitive modules.

The quotient directions between components and their positive sizes reconstruct
the tournament up to isomorphism. Every automorphism permutes these canonical
components preserving sizes. Each permitted quotient automorphism lifts
uniquely because a finite transitive tournament is rigid. Neither an undirected
path forest nor sizes alone retain quotient directions. The TT3 hostile shows
that homogeneous-pair adjacency itself is not transitive: its two adjacent
pairs are homogeneous, but its end pair is not.

For a finite odd-order action on actual module gaps, Burnside gives the
cut-subset orbit polynomial as an average of products `∏(1+t^length)`.
Every length is odd. At t=-1, each summand has vanishing order its number of
cycles and positive first nonzero coefficient; therefore the minimum cycle
count is the exact root multiplicity. Fixed symmetric-difference translation
by F descends to cut orbits iff F is invariant: apply it to equivalent F and
gF, whose outputs include the empty set, which alone has an empty-set orbit.
The invariant F form the Boolean group on gap orbits. Its rank need not equal
the root multiplicity, because one element's simultaneous cycle profile is not
the collection of group orbits.

The joint codec preserves all tournament automorphisms. Γ recovers their action
as its action on centers; its leaf kernel acts trivially on those centers and
hence on reconstructed modules/gaps. Both can therefore recover the correct
full gap action if that derived structure is actually retained. A specific
labelled gap/cut is a marked input, not automatically an equivariant operation.
The independent six-vertex example confirms the simplest failure: a single-gap
toggle does not descend, while toggling all three gaps does.

## 9. Decoder API boundary, actions, and blockers

The producer decoder functions are correct **on valid images**, including
relabellings, but are not validators. Concrete tests in
`producer-decoder-boundaries.json` show:

- appending an isolate to the n=1 joint or Γ image still returns the one-vertex
  tournament, despite the graph's wrong order;
- removing the only n=2 port edge returns the all-zero 2×2 matrix, not a
  tournament and not an explicit rejection.

The report never proves arbitrary-input recognition, so this is an API scope
clarification rather than a counterexample to Theorem J or I. The independent
strict recognizers demonstrate that polynomial exact recognition is available
if desired, without modifying the producer API.

| Action | Target | Reason | Done/evidence | Status |
|---|---|---|---|---|
| Verify sealed inputs and reproduce runs | Producer six-file packet and all copied inputs | Reproducibility and integrity | Manifest passed; twelve source pins checked; producer normal/-O identical to seal | **DONE / PASS** |
| Independent proof and hostile checks | Codec/no-go/metric/clone/OCF/spectrum claims | Mathematical validity | Proof audit above; independent normal/-O receipts | **DONE / PASS** |
| State valid-image-only decoder contract | Producer report §8 or publication integration copy | Existing routines accept malformed graphs | Root integrated exact clarification; documentation-only delta and new seal verified | **DONE / PASS** |
| Add current correction above preserved historical body | `THM-209-independence-polynomial.md` | Literal cycles cannot be collapsed to supports; all-order OCF already proved | Exact pin/current-HEAD comparison and minimal n=5 witness supplied | **ROOT INTEGRATION PENDING**; historical truth-surface correction |
| Optional loop-convention clarification | `THM-002-ocf.md`, proof 1 | Primary complement includes loops | Exact harmless repair supplied | **OPTIONAL DOCUMENTATION**, no OCF blocker |
| Commit/push/publication | Parent-controlled repository/publication packet | Worker assigned read-only audit | No worker mutation or publication attempted | **OUTSIDE WORKER SCOPE**, not reported as completed |

There are **no unresolved mathematical blockers** to the Session66 claims with
the scopes retained above. This verdict does not certify arbitrary malformed
input behavior of the producer decoders, spectrum completeness, graph-isomorphism
complexity improvements, new arithmetic realization, or any major open problem.

## 10. Reproduction

From the workspace root:

```bash
(cd session66-graph-pairs && sha256sum -c SHA256SUMS)
python3 session66-graph-pairs/check.py > /tmp/session66-graph-pairs-producer.json
cmp /tmp/session66-graph-pairs-producer.json session66-graph-pairs/normal.json
python3 -O session66-graph-pairs/check.py > /tmp/session66-graph-pairs-producer-O.json
cmp /tmp/session66-graph-pairs-producer.json /tmp/session66-graph-pairs-producer-O.json
python3 session66-audit-graph-pairs/audit.py > /tmp/session66-pair-audit.json
python3 -O session66-audit-graph-pairs/audit.py > /tmp/session66-pair-audit-O.json
cmp /tmp/session66-pair-audit.json /tmp/session66-pair-audit-O.json
cmp /tmp/session66-pair-audit.json session66-audit-graph-pairs/normal.json
python3 session66-audit-graph-pairs/additional.py > /tmp/session66-pair-extra.json
python3 -O session66-audit-graph-pairs/additional.py > /tmp/session66-pair-extra-O.json
cmp /tmp/session66-pair-extra.json /tmp/session66-pair-extra-O.json
cmp /tmp/session66-pair-extra.json session66-audit-graph-pairs/additional-normal.json
(cd session66-audit-graph-pairs && sha256sum -c SHA256SUMS)
```

Both scripts use only Python's standard library. The portable checker requires only the producer packet, not a Git checkout.
It always verifies all copied source SHA256 hashes and computes their Git blob
IDs directly from bytes. By default it locates sibling `session66-graph-pairs/`;
set `S66_PRODUCER` to its actual location after relocation. If desired, set
`S66_REPO` to a Git checkout for an additional `git show` comparison at the
frozen source commit; that opt-in comparison was performed during the original
audit but is not required for offline reproduction. An isolated replay with
only the scripts and producer packet, no source repository, reproduces both
normal and optimized receipts. See `portable-replay.json`. No external package,
network access, or source mutation is needed.
