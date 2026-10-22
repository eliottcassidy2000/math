# Intrinsic seams classify an infinite metagraph boundary

This packet advances the Session29 profile/seam family from edge-weighted statistics to **exact target isomorphism classes and complete family-boundary wiring at every block length**. It does not classify the full tournament metagraph, because target vertices can have other neighbors outside this family.

The input alphabet consists of the five fixed-spine tiles of seed A (canonical tournament mask40, H=15) and the thirteen of seed B (mask41, H=13). Both seeds are self-converse. A source word is their ordered tournament join. Its target flips every nonspine arc. The previous independently audited packet proves that all targets with at least two blocks are strong and gives the exact three-port transfer.

## 1. An intrinsic-cut theorem

For a tournament T and a nonempty proper vertex set U, let d+(U) count arcs from U to its complement. Call U a one-arc cut if d+(U)=1. This is a vertex-set cut; its unique outgoing arc is a strong bridge when T is strong. It is not a bridge of the tournament metagraph.

**Theorem.** In a strong tournament, all one-arc cuts with both sides of size at least two form a chain under inclusion. Consequently there is at most one such cut of each size k, 2<=k<=n-2.

**Proof.** Suppose A and B are incomparable one-arc cuts. Put P=A intersect B, Q=A minus B, R=B minus A, and U=A union B. Direct arc counting gives

 d+(A)+d+(B)=d+(P)+d+(U)+|Q||R|.

If P is nonempty and U is proper, strong connectivity makes the first two terms on the right at least one each, while |Q||R|>=1. This contradicts the left side being two. If P is empty, the arcs between A and B alone number |A||B|>=4 and are all included in the two outgoing-cut budgets, again impossible. If U is all vertices, the arcs between Q and R number |Q||R|>=4 because Q and R are the complements of B and A. These exhaust the cases. QED.

The finite checker independently tests cut-size uniqueness in all568 strong labeled tournaments of orders4 and5. The proof, not this finite test, gives the all-order conclusion.

## 2. Reconstruction for every number of blocks

For r>=2, write the complemented target as r five-vertex blocks X1,...,Xr. All cross-block arcs point backward except the single seam from the exit of Xi to the entry of X(i+1). The prefix sets of sizes5,10,...,5(r-1) are one-arc cuts. The target is strong because its fixed Hamiltonian path closes to a Hamiltonian cycle.

By the theorem, each of these prefix cuts is intrinsic and unique. Every isomorphism preserves them and hence preserves every block. The first block needs only its exit marked, the last only its entry, and each internal block needs both entry and exit marked. Conversely, isomorphisms of these marked blocks assemble into a full target isomorphism because the cross-block arcs and seams are completely prescribed.

For convenient common coordinates, encode the last block by its converse with the entry marked as an exit. The **target address** is then

 (q_left, c2,...,c_(r-1), q_right),

where q's are exit-rooted five-vertex types and c's are doubly rooted types. Under converse the address reverses, and each internal c changes by the operation tau: reverse all arcs and exchange entry and exit.

This yields a necessary AND sufficient isomorphism criterion in this family, not just a separating collection of numerical invariants. Automorphisms factor as the direct product of the marked-block automorphism groups. Additional one-arc cuts inside blocks do not compromise reconstruction; the prescribed cut sizes remain unique.

## 3. The finite rooted alphabet

Complete canonicalization of the five-vertex blocks uses only24 permutations for an exit-rooted type and6 for a doubly rooted type. The18 tile letters give:

- a=13 exit-rooted types;
- b=17 doubly rooted types;
- f=5 tau-fixed doubly rooted types.

Only the two letters755 and949 collide as doubly rooted types. Their class is tau-fixed, but neither letter is mirror-fixed. The other four tau-fixed classes have the individual mirror-fixed letters983,955,757,689.

The exit-rooted multiplicity types (number of A letters, number of B letters) are:

- (1,0): two types;
- (1,1): three types;
- (0,1): six types;
- (0,2): two types.

Thus a collision can cross the source A/B distinction. For example the exit-rooted canonical type76 contains the A tile983 and the B tile955. The complete memberships, marked automorphism orders, transfer matrices, and mirror action are in the JSON files. Canonical integer labels are tied to the documented lexicographic edge-bit convention and are not intrinsic numerical invariants.

## 4. Exact target-class generating functions

Let N_r be the number of ordinary target isomorphism classes reached at block length r>=2. The reconstruction theorem gives

 N_r = 13^2 * 17^(r-2).

The self-converse counts are

 F_(2m)=13*17^(m-1),
 F_(2m+1)=65*17^(m-1),   m>=1.

Burnside's two-element orbit formula gives merged target count (N_r+F_r)/2. Its generating function is

 (1/2) [169 z^2/(1-17z) + 13 z^2(1+5z)/(1-17z^2)].

For r=2,3,4,5,6 these merged counts are91,1469,24531,415701,7059403. Ordinary counts are169,2873,48841,830297,14115049.

Every even-length self-converse target address has a mirror-fixed source tile reaching it. At odd length, the central doubly rooted class {755,949} is an obstruction: it is fixed as a marked isomorphism class but has no fixed tile representative in this alphabet. Therefore exactly

 13*17^((r-3)/2),  r odd>=3,

self-converse target classes receive **no blue edge from this source family**. The generating function for these classes is13 z^3/(1-17z^2). This is a statement about the selected boundary, not an assertion that the target has no blue edge anywhere in the full metagraph.

## 5. Complete multiplicity and colour wiring rules

For each rooted type q and internal type c record its two letter multiplicities m_A,m_B. A prescribed ordinary source word w in {A,B}^r and an ordinary target address s have exactly

 m_(w1)(q_left) * [product_(i=2)^(r-1) m_(wi)(ci)] * m_(wr)(q_right)

lines between them. Zero factors mean no line. The right-root multiplicities agree with the exit-root table because mirror preserves each seed family. This product counts actual tile choices; every choice yields a unique complement line, and reconstruction proves no choices were assigned to a wrong target class.

For converse-merged nodes, sum this formula over the distinct ordinary source words in {w,reverse(w)} and the distinct ordinary target addresses in {s,tau(reverse(s))}. Do not double a fixed orbit.

Blue lines require a mirror-fixed source tile word, so the source word is palindromic and the target address is self-converse. For even r, choose its first half of tiles freely within the required rooted classes; the other half is forced. Thus the blue multiplicity is the product of first-half multiplicities, with the first block using its exit-root class and the remaining first-half blocks using doubly rooted classes. For odd r, multiply by the number of individually mirror-fixed tiles of the required central class and source type. Black multiplicity is total minus blue.

These are sufficient realizability and exact counting rules for **every line from the binary-join family**, including all parallel lines, at every r. They do not supply all incident lines at a target vertex, so they are not a global complete-wiring theorem.

### Explicit r=2 incidence formula

For exit-root classes q and r, write their multiplicities as (a,b) and (c,d). If q<r, the merged target receives

 AA:2ac;   AB/BA:2(ad+bc);   BB:2bd,

all black. On a diagonal q=r it receives

 AA:a^2, including a blue;
 AB/BA:2ab, all black;
 BB:b^2, including b blue.

The JSON table lists all91 target columns and all three source rows. There are324 lines total,18 blue and306 black. It also supplies an explicit vertex permutation, with a converse flag when needed, for every one of the324 target tiles. All witnesses are checked arc by arc.

## 6. Distinct-target Hamiltonian sums, not edge-weighted sums

The three-port transfer descends through the marked-block equivalence. Summing one representative per exit-rooted type gives

 v=(121,17,281)^T.

Summing one transfer per doubly rooted type gives

 B=[[163,-23,383],[23,-4,64],[383,-64,664]],
 J=diag(1,-1,1).

Consequently the sum of H over distinct ordinary target classes has the exact GF

 z^2 v^T J (I-zB)^(-1) v
 = (93313 z^2+3359832 z^3+4057812 z^4)
   /(1-823z-37140z^2-45180z^3).

The r2 value93313 differs from the edge-weighted value195472 because repeated targets are now counted once. The sum over the91 merged r2 targets is50487; their13 self-converse targets contribute7661, and (93313+7661)/2=50487.

For completeness, the merged distinct-target H sum also has an exact rational certificate. Let L be the six-dimensional second-moment transfer summed once over the17 doubly rooted types, and let

 M0=sum_(13 root types q) (Hq^2,Sq^2,Rq^2,Hq Sq,Hq Rq,Sq Rq)^T
   =(1301,39,6399,207,2579,343)^T.

With even cost (1,-1,1,0,0,0), and odd cost obtained from J times the sum of the five tau-fixed matrices, the self-converse H-mass GF is

 z^2 * even_cost (I-z^2 L)^(-1) M0
 +z^3 * odd_cost (I-z^2 L)^(-1) M0.

Half of the sum of this and the ordinary H-mass GF is the merged H-mass GF. All matrices and reduced exact rational expressions are recorded in all-length-results.json. The merged initial H sums for r2..5 are50487,41043219,34890376671,30108479835751. The certificate uses symbolic matrix inversion, not recurrence fitting.

## 7. Bridge certificates and a black-core consequence

The r2 family-boundary multigraph has94 vertices (three sources,91 targets),324 edges, and is connected. Its underlying simple graph has133 edges and55 degree-one target vertices. Its raw multigraph has exactly8 bridges; all are blue edges to singleton diagonal targets:

- AA targets with root keys333 and367, H521 and349;
- BB targets with keys842,335,332,74,77,78, H585,221,849,385,489,401.

These are bridges only in this **restricted boundary graph**. Whether they are bridges in the full metagraph remains unresolved here. For every other line the checker gives an explicit alternate path within this boundary, certifying a cycle. Hence316 of these lines are rigorously nonbridges in the full metagraph as well.

After retaining black edges only, the remaining86 nonisolated vertices form one connected component with306 edges. Each black source-target multiplicity is even. In particular the three source nodes AA, AB/BA and BB lie in one black-connected component. This is an exact local black-core connectivity result, not a global connectivity theorem.

### General black-pair lemma for the raw merged multigraph

Let F complement nonspine arcs and R mirror a tile. They commute; R preserves its converse-merged node. For n>=3, neither F nor FR fixes a tile: its first-last chord is nonspine, fixed as a coordinate by R, and toggled by F. A black F-line therefore has a distinct R-image joining the same merged endpoints. If its R-image were itself, R would either fix an endpoint (blue) or swap the two endpoints (an FR fixed point), both excluded. Thus every black line has a parallel mate, and no black line is a bridge in the raw merged multigraph. Loops are of course never bridges. This elementary observation is given without a novelty claim.

The conclusion must not be transferred to a picture that suppresses parallel edges. Bridge classifications depend on whether that suppression is part of the graph's definition.

## 8. Reproducibility and limits

- classify.py: standard-library rooted canonicalization, all324 arc-by-arc witnesses, all324 balanced-cut checks, exhaustive small hostile tests, and the full r2 incidence table.
- all_lengths.py: finite doubly rooted classification and exact SymPy resolvents for the all-r target and H-mass formulas.
- cycles.py: every nonbridge cycle certificate, the eight restricted bridges, and black-core connectivity.

Normal and optimized-mode outputs agree. Scripts read the frozen Session29 profile data from the sibling session29-profiles directory and never modify it. No global tournament census is used. The proof of intrinsic cuts is what replaces factorial order10 isomorphism search and extends the finite alphabet classification to arbitrarily many blocks.

Remaining obligations include the other neighbors of these strong target nodes, the global status of the eight boundary bridges, and extension from this special seamed family to all strong tournaments. The construction gives a tractable exact subproblem and new generating functions, not a proof of the complete tournament classification or of the global Hamiltonian-count gap.


## 9. General fixed-block reconstruction theorem

The intrinsic-cut argument is not specific to five vertices. Let k>=2 and let a finite mirror-closed alphabet of fixed-Hamiltonian-path k-vertex tiles be given. Form all ordered joins of r>=2 letters and complement the nonspine arcs. Every target is strong, and its prefix cuts of sizes k,2k,...,(r-1)k are uniquely intrinsic. Thus its ordinary isomorphism class factors into the exit-rooted first block, doubly rooted internal blocks, and entry-rooted last block, exactly as above.

If a,b,f are the numbers of outer rooted types, internal doubly rooted types, and mirror-fixed internal types, respectively, the same class-count generating function holds:

 (1/2) [a^2 z^2/(1-bz) + a z^2(1+fz)/(1-bz^2)].

The corresponding product rules give all boundary-edge multiplicities once the alphabet's rooted multiplicity tables are known. This theorem requires fixed block size and the full independent word alphabet. A constrained word language would require its own language-counting step. The claim is not made at k=1, where the cuts need not have both sides of size at least two and directed C3 already has incomparable singleton one-arc cuts.

The source join also has intrinsic size-k segmentation: its strongly connected components have a canonical total order, no component crosses a join boundary, and cumulative component sizes recover each prescribed boundary. Thus nonstrong blocks do not themselves break source reconstruction. The specific three-state Hamiltonian observer used earlier is established here for the finite Session29 alphabet; this general cut theorem does not by itself prove three-state compression for every possible alphabet.

## 10. A polynomial-time seam recognizer

The parent proposed the following constructive speedup, which recognize_cuts.py independently implements and tests. For a strong n-vertex tournament:

1. Choose a root, one spanning out-arborescence and one spanning in-arborescence. Their union contains at most2n-2 arcs.
2. Every strong bridge lies in this union. Otherwise deleting it leaves both arborescences intact, and all vertices still reach one another through the root.
3. For each candidate arc e=u->v, delete e and compute the strongly connected components. Skip the candidate if the graph remains strong.
4. If deletion disconnects strong connectivity, the SCC condensation has a total reachability order. Enumerate its nonempty proper suffix unions. These are exactly the one-arc cuts whose outgoing arc is e.
5. Filter the cuts by both-side size and prescribed block sizes, then canonicalize the small marked blocks.

For step4, all distinct SCC pairs have an arc between them except possibly the pair consisting of singleton u and singleton v, because only e was removed. Even that pair is comparable: a simple path from v to u in the original strong tournament cannot use e, so survives deletion. Thus the condensation reachability order is total. Adding e back makes it strong, so e must go from its last component to its first. The nonempty proper suffix unions are precisely its zero-outgoing-arc cuts, and restoring e gives exactly one outgoing arc. Conversely any one-arc cut becomes such a zero cut on deletion, hence is one of these suffix unions.

With adjacency matrices, each candidate costs O(n^2), so all cuts are found in O(n^3), rather than testing2^n subsets or n! vertex labelings. This is a polynomial recognition step for the intrinsic bottlenecks; canonicalizing arbitrary large blocks remains a separate problem. For fixed small k, the marked-block canonicalization costs a constant per block.

The checker compares this algorithm to all saved exhaustive cuts of the324 order10 targets, and to exhaustive subset testing for every strong labeled tournament at orders3,4,5. It also recovers every prescribed seam in deterministic examples up to20 blocks (100 vertices). Tests use exact arcs, not approximate spectra or heuristic isomorphism hashes.

## 11. The global blue bridge skeleton and parity charge

The fibre parity foundation is already in repository THM-643 (T1-T5), which also points to the earlier black-quotient bridge packet. The present statements are recovered structural consequences, not claimed as new discoveries. The referenced older packet was not locally materialized in this sparse checkout, so no unseen claims from it are imported.

In the full raw converse-merged multigraph, a self-converse node has odd degree H(T)/|Aut(T)|. A nonself-converse merged pair has even degree2H(T)/|Aut(T)|. Loops count twice. Tournament automorphism groups have odd order, and Rédei gives odd H, proving these degree parities.

The black-pair lemma implies every bridge is blue. Deleting a bridge leaves two shores; the degree-sum identity shows that each shore contains an odd number of self-converse nodes. This applies to nonleaf bridges too: it is a parity constraint on their entire shores, not a claim that all bridges are pendant.

More precisely, contract each connected component of the full black-edge subgraph, retaining every blue edge with its multiplicity. A blue edge is a bridge before contraction if and only if it is a bridge in this blue quotient multigraph. Any alternative quotient path expands using black connections within the contracted components, and any original alternative path projects to a quotient walk. Blue edges internal to one black component become loops and are nonbridges.

The degree parity of a quotient vertex is the number of self-converse original nodes in that black component, modulo two. This gives a **bridge-parity skeleton**: blue quotient structure together with a parity charge on its vertices. It is necessary structural information, not a determination of bridges from degree counts alone.

These full-graph charges must not be assigned to the degrees of targets in the restricted boundary graph above, which deliberately omits their other incident lines. Its316 nonbridge certificates remain globally valid because they exhibit actual cycles; its eight restricted bridges have no global verdict here.
