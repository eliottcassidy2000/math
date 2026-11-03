# Complete carriers, canonical ancestry, and a one-vertex native compiler

**PROVED:** the one-temporary-vertex compiler in §5 and the elementary
consequences explicitly proved below. **INHERITED PROVED:** complete fixed-order
reconstruction, guarded native transport, source weights, and path-language
reconstruction. **FINITE-EXACT:** the bounded independent receipts in §7.
No primitive-core inventory, all-target Hamilton-count induction, or new
fixed-order native/black connectivity theorem is claimed.

## 1. Inheritance and what “complete” can mean

The closest complete mechanism is the all-order legal-relocation theorem and
its [count-free reconstruction algorithm](../session67_metagraph_enumeration_information_20261026/session67-root/COUNT_FREE_ENUMERATION.md).
It already reconstructs every ordinary class, every converse-merged class,
and every literal coloured metagraph line at each finite order. It needs no
supplied class count or graph-isomorphism oracle. This is an exact algorithm,
with exponential state space, rather than a small formula for the wiring.

The canonical hostile is the order-four pair of different normalized words
with the same source/target ordinary class pair but different next-order
extension distributions. The [weighted native recursion](../session105_metagraph_recursion_20261029/metagraph/WEIGHTED_NATIVE_RECURSION.md),
§5, supplies literal masks 4 and 18 and their separating strip event. Keeping
the endpoint classes loses the alignment required for the next birth.

The corrected near misses are: a fixed-spine cube is not a full invariant
labelled carrier under arbitrary relabelling; auxiliary path relocations are
not complement lines; converse is not native complement; pure blue does not
mean a leaf; and a degree list does not specify a graph. The current correction
in [THM-643, blue/black line parity](../../../01-canon/theorems/THM-643-blue-black-line-parity-structure.md)
controls the old purity and blue-cap prose. The regular five-vertex class has
three blue presentations and three neighbours. The least-used sidecar restored
here is the actual named path forest after restriction, not its compressed order.

The research board and rounds were:

| Lane | Live object | Consequence of the next hostile test |
|---|---|---|
| Anchor | Full native incidence, including source mass | Recover the existing complete cube/relocation compiler first. |
| Ordinary classification | Strong-component words and adaptive ancestry | Unique peeling gives canonical coordinates, but leaves primitive cores to classify. |
| Symmetry | Relabelling, converse, reflection, stabilizers | Indegree ancestry can change under converse; store both gauges when merging. |
| Niche | Birth and deletion of a named vertex | Internal deletion removes a retained path edge and creates a one-arc defect. |
| Wildcard | Compile an arbitrary arc reversal using this defect | One auxiliary vertex and two legal native moves suffice, with explicit ports. |
| Arithmetic | Hamilton counts and rooted cover profiles | These are observables, not complete metagraph states or freely selectable labels. |

These distinctions give a concrete answer to the completeness question:

* **Complete finite reconstruction:** already proved by legal relocations plus
  complement; the algorithm works at every finite order.
* **Complete object normal form:** ordinary tournaments are ordered words of
  strong components, each a primitive core with a root word modulo its core
  automorphism group (§4). This is a canonical decomposition, not an inventory
  or closed count of primitive cores.
* **Complete named reachability with one temporary vertex:** proved in §5,
  with an explicit operation budget and no hidden relabelling.
* **Complete fixed-order native geometry:** none of these statements by itself
  supplies all components, distances, small sufficient observers, or black-only
  routes. Existing scoped results such as [THM-4819, prescribed black cores](../../../01-canon/theorems/THM-4819-regular-tournaments-under-native-moves-the-transitive-slab-exact-distance-and-the-prescribed-black-core.md)
  retain their stated regular/balanced hypotheses.

## 2. The exact carriers and generators

Fix a named vertex set V of size n≥3. A tournament has exactly one directed
arc on each unordered pair. Write HP(T) for its directed Hamilton vertex words.
For P∈HP(T), the native operation K_P retains the n−1 consecutive path arcs
and reverses all other arcs. The path remains directed, so K_P is its own inverse.

Let X_n consist of tournaments on positions 0,…,n−1 with the spine
0→1→…→n−1. Its q=binom(n−1,2) other arcs are independent bits. Equivalently,
X_n is the set of pairs (T,P) modulo simultaneous relabelling. The permutation
action on pairs is free: a permutation fixing the ordered full path fixes all
vertices. Thus |X_n|=2^q without identifying ordinary tournament classes.

| Carrier | Exact quotient or marking | Data retained / lost |
|---|---|---|
| Named tournament T | No quotient | Every arc and physical vertex name. |
| Named frame (T,P) | A particular legal full word | Every source arc and every path position. |
| Normalized frame x∈X_n | Simultaneous relabelling of (T,P) | Full relation in path order; external names require the transport map. |
| Ordinary class A | Forget P and identify relabellings | Frame fibre has size h_A=H_A/a_A, with a_A=|Aut(A)|. |
| Converse-merged node {A,A^op} | Further quotient by converse | Mass h_A if self-converse, otherwise 2h_A. |
| Literal metagraph line | The pair {x,τx} | Loops and repeated lines remain distinct. |
| Reflection-paired black edge | Further identify a black line and its reflection | Literal black multiplicities, including loops, are divided by two. |
| H or a rooted four-coordinate profile | An observer, not a symmetry quotient | Does not retain classes, paths, incidence alignment, or general continuation. |

Here τ is bit complement on X_n, and σ is reversal of the positions together
with converse; σ only permutes free coordinates. They commute. Blue means
σx=x. Ordinary classes need not be fixed by σ: σ carries their fibres to their
converse fibres. A blue word lies over a self-converse class.

Arbitrary permutations do **not** act on X_n by simply permuting vertices and
leaving its prescribed spine unchanged. The correct re-marking operations are
partial maps: choose a new legal Hamilton word, then reindex the complete
relation in its order. Legal single-vertex relocations generate each ordinary
fibre. The inherited theorem gives a local recognition certificate: a proposed
partition is correct iff each block is relocation-connected and every legal
relocation edge stays inside its block. A spanning forest without edge closure
can split fibres; closure without connectedness can merge them.

The complete reconstruction procedure is therefore exact at every order:

1. Enumerate X_n and its guarded relocation edges; form their components.
2. Those components are ordinary classes. Add σ-relations if merged nodes are wanted.
3. For each τ-pair, emit one literal line with its two component endpoints and
   colour. Keep loops and multiplicity.
4. If wanted, pair reflected black lines. This is an additional quotient,
   not a change in the native operation.

Reflection acts freely on black literal lines: if σ fixes {x,τx}, either
σx=x, making it blue, or σx=τx. The second is impossible because σ has a fixed
free coordinate and complement flips that coordinate. Thus the paired-black
convention has precisely half the literal black multiplicity. If a directed
half-edge matrix has diagonal entry d at a node, it contributes d/2 literal
loops and d/4 paired black loops. The [paired-black source](../merged_metagraph_black_quotient_bridges_20261020/README.md)
fixes these conventions explicitly.

On named objects relabelling obeys πK_P(T)=K_(πP)(πT), and converse obeys
K_reverse(P)(T^op)=K_P(T)^op. For a legal word of native moves,

    T_r = T_0 XOR ((r mod 2) E(K_V)) XOR XOR_j E(P_j).

This inherited sign grammar determines the endpoint after the actual guards
have been checked. It is not permission to apply an arbitrary unguarded word.
We do not assert that involution and covariance relations alone present every
relation of the resulting transformation system.

## 3. Exact source weights and the birth sidecar

For ordinary classes A,B set

    M_AB = #{x∈X_n : c(x)=A, c(τx)=B},
    W_AB = #{P∈HP(T_A) : [K_P(T_A)]=B}.

The inherited exact laws are

    M_AB=M_BA,   W_AB=a_A M_AB,
    Σ_B M_AB=H_A/a_A,   Σ_B W_AB=H_A.

Loops contribute twice to the diagonal of M. Named incidences between all
labelled representatives have mass n!M_AB. A uniform path at a fixed source
induces transition M_AB/h_A; a uniform permutation proposal with invalid words
held fixed uses a different clock and stationary source weights. Neither
clock is a uniform choice of ordinary classes.

The aligned pair (T,K_P(T)) reconstructs P: their agreement arcs are exactly
the spanning directed path. Its simultaneous stabilizer is therefore trivial.
There are exactly 2^q simultaneous-isomorphism types of directed native pairs.
Retaining the full aligned pair is lossless but does not compress the frame.
Independently canonicalizing its endpoints destroys this alignment, as the
order-four strip hostile demonstrates.

Appending a last position gives the inherited exact strip intertwining

    τ_(n+1) E(x,b) = E(τ_n x, complement(b)),

where the last old vertex beats the new one and b specifies all other new arcs.
It is a complete recursion before the final class projection. It does not
close on endpoint class pairs.

For an arbitrary subset S, the [actual-path flag rule](../metagraph_flag_transport_20261023/README.md)
retains F_S=E(P)∩binom(S,2); native complement on S fixes F_S and reverses its
complement. Deleting a path vertex does not add a shortcut. This exact rule
commutes with restriction at every order.

A useful specialization is the **internal-insertion defect**. If P is a legal
old path containing u→v, and Q inserts z between u,v in a one-vertex extension,
then

    K_Q(extension)|V = K_P(T) with the old arc uv reversed.

The old retained edge uv disappeared from Q. If z is appended or prepended,
there is no old-edge defect. When deleting z from Q does not give a legal old
path, use the retained forest itself; the compressed word is not an authorized
Hamilton frame. This defect supplies the compiler below.

## 4. Canonical strong ancestry and the full ordinary carrier

[ANCESTRY.md](ANCESTRY.md) proves the all-order unique inverse. In a strong
T of order >3 there is at most one vertex x of indegree one for which T−x
is strong. Repeated legal deletion reaches a unique irreducible strong core K.
Reversing the history gives a chronological root word: each newborn loses only
to its selected old root and beats every other old vertex.

For fixed core K with G=Aut(K), ordinary descendants at depth m are exactly
these root words modulo G, acting on core letters and fixing chronological
newborn letters. Their number is

    (1/|G|) Σ_(g∈G) ∏_(i=0)^(m−1) (fix_K(g)+i).

The final automorphism group is the subgroup of G fixing every core vertex
that appeared as a root. Newborns are individually characteristic. For C3,
all positive-depth descendants are rigid and the order-n count is (n−1)!/6,
n≥4. General primitive cores can retain nontrivial stabilizers.

Combining this with the unique ordered strong-component decomposition gives
an exact ordinary normal form for **every** tournament: an ordered word whose
letters are singletons or (primitive strong core, root-word orbit) objects.
Automorphism groups multiply across the components; so do H and H/Aut.
This is a complete specification in terms of irreducible objects. It does not
classify the primitive cores themselves or reduce the full class-count problem
to a solved scalar sequence.

It also gives a procedural metagraph compiler. Reconstruct the named component
objects and choose one actual path in each. Concatenate those paths. Native
complement reverses all cross-component arcs except the exit-to-next-entry
bridges, and applies the native rule inside each component. Recompute the SCC
and ancestry normal form of the result. The inherited weighted SCC formula is
exactly the sum over these ported child words, with source factor
∏|Aut(T_i)|. A product of unported child adjacency matrices loses the bridge
exceptions and is insufficient.

**Converse obstruction.** Indegree ancestry is not an invariant of a merged
node. The C3 history (0,1,2) gives

    T    = (10,20,33,6,13,27),
    T^op = (52,41,26,49,34,4),

in outgoing-row-bit notation. T has core C3, but T^op is an irreducible strong
six-vertex core: its indegrees are the outdegrees of T, namely (2,2,2,2,3,4),
so it has no indegree-one vertex. The complete histories through orders four
and five have no such exit; order six is the first for this adaptive family.
Thus a converse-merged normal form must compute both oriented normal forms
and quotient the resulting pair, or specify a canonical converse gauge. Merely
reversing the SCC order and the root word is wrong.

Allowing both indegree-one and outdegree-one legal inverses does not repair
this by a confluent reduction: [TWO_SIDED_PEELING.md](TWO_SIDED_PEELING.md)
gives a minimal six-vertex source reducing to terminal cores of orders three
and five, and a critical family proving nonconfluence at every order n≥6.
Its wrapper W(K)=K→y→b→K embeds every strong class into a two-sided irreducible
core with two additional vertices, preserving its exact automorphism group
and allowing recovery of K by its characteristic ports. Hence the irreducible
core problem itself contains the full strong-class problem; the normal form
must retain that unresolved inventory. The wrapper count is 3H(K)+c_2(K),
with c_2 the ordered nonempty two-path-cover count, not a function of H alone.

**Native nonclosure.** Even the first adaptive class leaves its family under
a native move. The strong source R4=(2,4,9,3), with H=5, has path (0,1,2,3).
Its native partner is the transitive tournament (14,12,8,0), with H=1. The
exact path is retained; strongness and adaptive membership are lost. This is
a failure of closure, not a failure of the adaptive count realizations.

## 5. One temporary vertex compiles any old arc reversal

**Theorem (new synthesis, all orders).** Let T be any named tournament on
n≥2 vertices and e an old unordered pair. There is an explicitly constructible
extension by one temporary vertex z, two legal native moves on n+1 vertices,
and deletion of z whose final restriction is exactly T with e reversed.
Every other old arc and every old name is unchanged. At birth, z has either
at most two in-neighbours or at most two out-neighbours. The construction
retains its complete incident cut, both path words, and the named deletion.

**Traced case, n≥3.** Suppose e is u→v and a Hamilton word P contains u,v
consecutively. Let f and l be its first and last vertices.

* If v≠l, birth z beating only u. Use P_1=Pz. In K_(P_1)(T+z), the word
  Q obtained by inserting z between u and v in P is legal: u→z and z→v
  were reversed, while every other old consecutive P arc was retained.
* If v=l, birth z losing only to v. Use P_1=zP. Since n≥3, u≠f, so the
  same inserted word Q is legal after the first native move.

In both cases E(P_1) and E(Q), restricted to old pairs, differ exactly at uv.
The two-move XOR formula proves the required final old relation. All old
Q edges are explicitly checked, so no unguarded sign calculation is hidden.

The final temporary cut is equally explicit. In the append case z beats
exactly {v,l}; in the prepend case z loses exactly to {u,f}. Thus reversing
a gadget also requires only a two-port exception to a homogeneous cut.

**Untraced case.** By [THM-4807(d), one-guard separation](../../../01-canon/theorems/THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation.md),
the Hamilton languages of T and T^e differ. A distinguishing word must use
the changed pair consecutively. If no path of T uses u→v, some path of T^e
uses v→u. Construct the traced gadget from T^e, then traverse its two
involutions backwards. Its final temporary cut is now the birth cut over T.
Deletion produces exactly T^e. A search over both finite path languages
terminates; the theorem also supplies a structural separating-word decoder.

**Order two.** For u→v birth a sink z. Use P_1=(u,v,z), then Q=(v,z,u).
The first move flips uz, the second flips uv. Delete z.

This proves the theorem. No Hamilton-count hypothesis, strongness assumption,
class atlas, or arbitrary target relabelling is needed.

**Sharp two-step obstruction without a temporary vertex.** Two native moves
at the same order change the symmetric difference of two edge sets of equal
size n−1. Its cardinality is even. Thus they cannot change exactly one named
arc. The auxiliary vertex is necessary for a universal exactly-two-native-step
single-arc compiler. This does not claim a minimal distance among arbitrary
native words or among protocols permitting other old-arc operations.

**Complete stabilized reachability.** Given named T and U on the same n≥2
vertices, run the compiler once for each differing pair, deleting z between
runs. If d pairs differ, the program uses d births, 2d native moves, d named
deletions, and peak order n+1. These are integer operation counts on three
separate clocks; it is not a path of length 2d in the fixed-order metagraph.
The temporary vertex may be reused only with its new explicitly prescribed
birth cut. No claim is made that all these births are adaptive single-in-port
births, that the moves are all black, or that H is monotone or preserved.

The whole family of recipes is equivariant under simultaneous transport of
all names and paths. Choosing a lexicographically first witness adds a naming
gauge; it is not a canonical choice on ordinary classes without that gauge.

## 6. Connection ledger and stopping boundaries

| Source → target / map | Preserved predicate | Lost data and required sidecar | Cheapest decisive hostile |
|---|---|---|---|
| Named frames → ordinary classes | Underlying isomorphism type | Path position and source multiplicity; retain H/Aut and aligned frame | C3 has H=3 but one ordinary frame. |
| Full cube → metagraph via relocations and τ | Every literal incidence, colour and loop | External names; retain normalization map when transporting ports | Relocation is an auxiliary relation, not a native line. |
| Strong tournaments → ancestry coordinates | Exact named reconstruction; ordinary type modulo core group | Chronology if core alone retained; restore root word and stabilizer | Two primitive five-vertex cores have different Burnside laws. |
| Ordinary ancestry → converse merge | Only after recomputing both orientations | Indegree is changed to outdegree; retain both ancestry gauges | Six-vertex pair above has core orders 3 and 6. |
| Adaptive birth → native flag transport | Exact arcs outside the internal insertion defect | Compressed path omits broken edge; retain actual forest | Insert z internally and old uv reverses. |
| Native two-step gadget → single old arc reversal | Every other old arc and name, exact target | Auxiliary incident cut and operation clock | Zero-auxiliary two-step parity obstruction. |
| Rooted profile → arbitrary next-root response | Exact keep/new formulas only | Mixed endpoint covers and ownership by named roots | [ADAPTIVE.md](ADAPTIVE.md), order-six same-profile pair. |

The next structural target is therefore a sufficient state smaller than the
full aligned frame that still disintegrates birth and native incidence exactly.
The endpoint-pair and four-coordinate proposals already fail concrete hostile
checks. The ancestry normal form replaces redundant birth histories by canonical
ones, but does not by itself eliminate primitive-core data or path ownership.
This is a recorded boundary, not a reason to call a numerical spectrum scan a
complete metagraph classification.

The positive bounded repair now appears in [ACTIVE_PORTS.md](ACTIVE_PORTS.md):
restrict roots to current indegree at most d, retain the at most 2d active
ports, and record paired endpoints of path-cover components. Forgotten vertices
are dominated by all future births, so at most one forgotten end can survive.
This gives an exact finite nonnegative integer weighted carrier for each fixed
d. It repairs the counting state within that restricted birth language; it
does not close arbitrary native moves on that carrier or classify its integer
output image. [active_ports_audit.md](active_ports_audit.md) independently checks
the proof and cover-by-cover lumpability using successor-edge sets.

## 7. Reproduction and independent scope

Run from this directory:

    python3 metagraph_probe.py > metagraph_probe.json
    python3 ancestry_audit.py > ancestry_audit.json
    python3 adaptive_independent_audit.py > adaptive_independent_audit.json

`metagraph_probe.py` imports no repository modules. It checks every labelled
tournament and every arc for n=2,3,4,5: respectively 2,24,384,10240 recipes,
10650 total. It replays the two literal guarded native moves, exact old-vertex
restriction, final temporary ports, and the even full-graph difference. It
exercises 1198 inverse gadgets, so untraced arcs are not silently filtered out.
The forward append/prepend counts are 6844/2606; inverse counts are 1029/169;
the two order-two sources exercise the boundary gadget. The R4→TT4 exit is
also literal. This finite bank validates the implementation, not the all-order
one-guard-separation theorem on which the proof relies.

The spectrum lane independently audited the all-order compiler proof and the
two-sided critical family/wrapper embedding, with no corrections requested.

`ancestry_audit.py` independently checks every labelled strong source through
five vertices, exact named histories and automorphism stabilizers; it also
canonicalizes actual depth-one/two descendants of two primitive cores on six
and seven vertices. See [ancestry_audit.md](ancestry_audit.md). This is targeted
structural verification, not a new mass census.

`adaptive_independent_audit.py` counts paths and covers by literal vertex words,
separately from the producer's subset DP. It validates the two rooted updates,
reset, defect, and all fourteen eight-vertex children of the same-profile
hostile. See [adaptive_audit.md](adaptive_audit.md).
