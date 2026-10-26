# Why alignment memory grows: guarded choices, assignments and routed ports

2026-10-26. **PROVED** elementary deductions and explicitly displayed
finite certificates; **INHERITED PROVED** where linked; **FINITE-EXACT**
for enumerations. Historical numeric claims are not silently promoted.
The Graceful Tree Conjecture and LRC(14) remain open. No Collatz convergence
claim is made.

Anchor: all-size tournament guard reconstruction. Niche: dormant interfaces
and address recovery. Wildcard: graceful assignments, square-sum arrivals,
and the distinct six-state matching carriers of Petersen and K3,3.

## 1. The concept board and the recovered mechanism

| Concept | Decisive question | Result of this pass |
|---|---|---|
| Support versus native closure | Does an absent arc remain ambiguous? | Tournament trace support forces it by G1; graph support need not. |
| Legal operation domain | Do raw commuting maps remain legal after each other? | March Vitali example has a three-state guarded component inside a four-state raw orbit. |
| Independently addressed choices | Can a permitted test read each hidden assignment? | C4 reads k bits; graceful replacement guards read a permutation. |
| Component routing | Which endpoints are actually paired, and are all vertices covered? | K3,3/Petersen matching complements separate one cycle from two. |
| Relative frame and symmetry | Which action, rather than just which cardinality, survives? | K3,3's six matchings form regular S3; Petersen's six have a different action. |
| Scale and delayed activation | Can later arrivals expose a currently invisible edge? | G25 and its edge-deleted control agree now but differ17 paths at G27. |

The source ledger pins the passages read rather than claiming a corpus-wide
audit. This pass revisits March's actual reversal code and Walsh support,
July's rooted graceful reconstruction, and October's assignment, boundary,
and square-sum work. The preceding packet retains the June moment and
July continuation bridges. New calculations were compared with every row:
G1 is a positive class-specific recovery theorem, while the three other
lanes show exactly which hypotheses prevent it from transferring.

## 2. The March Vitali square: commuting flips with changing guards

The historical source is
[vitali_holonomy_composition.py](../../../04-computation/vitali_holonomy_composition.py),
March13, and its
[preserved output](../vitali_holonomy_composition.out).
It asked whether overlapping Vitali reversals form a group. A raw reversal
R_S toggles all pairs inside S, so any two such maps commute as XOR masks.
The permitted operation additionally requires the four-vertex score type
(1,1,2,2) and preservation of the current labelled lambda field, where
lambda(u,v) counts cyclic triangles containing that pair.

The first legacy seed42 example is the full upper-triangle mask401896
(lexicographic pairs, least significant bit first), with

    S={0,2,3,4},       Q={1,3,4,6}.

The four raw states are:

| State | Mask | Same labelled lambda as the start? | H | Simple directed c7 |
|---|---:|---|---:|---:|
| T |401896|yes|141|14|
| R_S T |440806|yes|143|15|
| R_Q T |828520|yes|143|15|
| R_Q R_S T=R_S R_Q T |801894|no|95|7|

All allowed four-subset reversals give exactly the first three states as
the start's reachable component. After S, Q is illegal; after Q, S is
illegal. Thus the raw commuting square does not descend to an everywhere
defined C2 x C2 action on that observation fibre. There is no group-theory
paradox in its three-state guarded component.

In finite-difference language, both first differences of lambda vanish at
T, but the mixed difference is lambda(R_S R_Q T)-lambda(T), which is
nonzero. The observer forgets which compatibility conditions will be
broken by another move. This is a precise early version of the present
continuation-memory obstruction, rather than a claim of nonmeasurability.

**Correction recovered during this pass.** The old routine calculated
trace(A^7)/7 and labelled it c7. Those values are227,228,228,154: they count
closed walks with possible repeated vertices, not simple cycles. The source
routine is repaired, the historical output carries a warning, and the
mistakes ledger preserves the lineage. The checker enumerates actual
Hamilton words; the independent audit separately checks all simple cycle
lengths3 through7 and exercises the repaired routine. The raw commutation
and guarded-component conclusion survive the correction.

The other recovered March records,
[overlap hidden dimension](../vitali_overlap_hidden_dim.out) and
[chirality dimension](../vitali_chirality_dimension.out), asked why the
109/111 twins survive many overlap and signed summaries. Their useful
question is about joint realization, not an inference that every finite
observer must fail. G1 supplies the counterindication: a sufficiently
structured support observer can reconstruct everything at every finite n.

## 3. The common amplification lemma and the graceful permutation

**Address lemma (PROVED).** Suppose a carrier contains states x_a indexed
by distinct words a in a subset A of a product of finite alphabets. Suppose
for any a!=b a permitted observation or guard, with its name retained,
distinguishes x_a from x_b by reading a coordinate where they differ.
Every exact continuation quotient has at least |A| cells on these states,
and needs at least ceil(log2 |A|) fixed-length bits. This is the pigeonhole
principle applied to distinguishable futures; it requires the tests to
remain legal in the assembled context.

The preceding tournament C4 has A={0,1}^k. Uniform quotient attachments
keep each local test available, yielding k unavoidable hidden bits even
after both rooted endpoint types and colour are supplied. G2 strengthens
the total literal-state bound to (n-1)(n-2)/2 bits. It does not replace
the different conditional fibrewise bound.

The October [one-sided kernel factorization](../attachment_lattices_square_sum_20261023/KERNELS.md)
is a symmetry version of the same independence test: a motion that fixes
every external owner can act in one cell while other cells stay put.
Breaking all such independent motions costs an additive supply of owners.
Its warning also applies here: quotienting a kernel preserves a literal
path query only if that query is constant on every omitted lift. Killing
a stabilizer and identifying the relative source-target alignment are
different tasks. Neither a group order nor a Frobenius degree alone
replaces the reconstructed relation and requested ancestry labels.

The inherited
[graceful factorial-fibre theorem](../graceful_factorial_fibers_20261019/GRACEFUL_FACTORIAL_FIBERS.md)
provides an orthogonal realization. Freeze its specified upper-tail forest
on labels0,...,m, m>=k^2+2k. Its components are a tree H and isolated
vertices b_i=i(k+1), i=1,...,k. Every free colour d in1,...,k can attach
b_i by either edge {b_i-d,b_i} or {b_i,b_i+d}. A graceful completion is
exactly a permutation pi assigning colours to leaves, plus k signs:

    (pi,epsilon) in S_k x {-,+}^k.

Here is an exact **guard quotient**, beyond the inherited component count.
Take named operations that replace the colour-d edge by a specified
same-colour edge, allowing replacement by itself. Observe only legality,
not the resulting numeric edge list. Replacing it by a spoke at b_i is
legal iff i=pi(d), independently of its sign. Otherwise the old leaf is
isolated and another leaf gets two contacts; an edge internal to H also
fails. Thus these one-step guards recover pi exactly and forget epsilon.
Legal replacements preserve pi, so equality of pi is stable under every
such named operation. Consequently the **smallest guard-only carrier is
exactly S_k**, with k! cells and ceil(log2(k!)) bits. This is an all-k
deduction for that explicitly frozen family and operation menu.

The signs are recovered if the observation exposes the actual chosen
edge, or if an identity replacement is instead declared illegal. The
protocol matters: there is no operation-independent notion of the smallest
state. This example is a positive compression boundary alongside G2's
discrete tournament result.

The July mechanism
[THM-2793 / oriented-ramp-rooted-graceful-reconstruction](../../../01-canon/theorems/THM-2793-oriented-ramp-reference-and-rooted-graceful-gap-reconstruction.md)
explains a different repair: a labelled ramp cross-block restores interval
positions, and a root-zero mark resolves cut complements. One reference
vector supplies one value **per interval**, not one constant-size scalar.
This is compatible with growing information bounds. It reconstructs a
given labeling; it does not supply a graceful labeling for every tree.

## 4. Square sums: delayed activation defeats guard-only graph recovery

The exact prior
[square-sum port census](../square_sum_ports_20261023/README.md)
already proved that no Hamilton path of G25 uses1--15. Therefore G25 and
G25 minus that edge have identical complete named Hamilton-path guards,
despite being different simple graphs. This is a direct failure of G1
outside tournaments. The edge-deleted graph is a hostile control, not
another integer square-sum graph G_n.

Apply identical arrival operations to both: adjoin26 and27 with their
actual square-sum neighbours, retaining every old edge. Exact re-enumeration
gives35 unoriented paths for G27 and18 for its edge-deleted control. The
other17 use1--15. Thus the same current guard profile can conceal different
future arrival responses. The checker stores a literal distinguishing path.

The tournament structural theorem is the missing hypothesis. For arbitrary
graphs, all ordered path-cover guards *do* recover adjacency: use a cover
with one proposed two-vertex path and every other vertex a singleton.
It is legal exactly when that edge exists. This uses arity n-1 and is not
a claim that two outer endpoints suffice. It explains why the prior
square-sum spanning-forest interface and the tournament all-run response
grammar have a common sufficient carrier despite different one-path limits.

The earlier
[graceful/square exchange audit](../graceful_square_exchange_20261019.md)
also supplies the correct operation comparison: graceful midpoint rotations
preserve a numerical colour, whereas square-sum rotations preserve an
arithmetic adjacency predicate. Connectivity and fixed-endpoint restrictions
must be checked separately. A degree-two path-plus-cycle configuration is
not a Hamilton path; the matching example below makes that failure visible.

## 5. K5, Petersen and K3,3 as exact port objects

### Five terminals and the compatibility graph

Take the ten two-terminal subsets of a five-element set. Two are compatible
as **vertex-disjoint** components exactly when the pairs are disjoint.
Their compatibility graph is KG(5,2), the Petersen graph; its complement
is L(K5), recording shared terminals. This is the dictionary already
proved in the [October5 tree/Petersen note](../sumner_t5_petersen_collatz_20261005.md),
now read as an attachment test. Sharing an endpoint is not intrinsically
illegal for path edges; this compatibility predicate specifically concerns
distinct components required to be vertex-disjoint.

The March [THM-071 / Walsh-Fourier diagonalization](../../../01-canon/theorems/THM-071-walsh-fourier-diagonalization.md)
puts L(K5) at the degree-two support of H on five vertices. Thus the same
ten pair addresses have an interaction graph and a disjointness graph.
This is an exact finite incidence identification, not a claim that the
five-tournament metagraph itself is Petersen, or that this Walsh layer
reconstructs all tournaments.

### Two different six-matching carriers

On labelled K3,3, perfect matchings are bijections pi between its two
three-element sides: exactly the six elements of S3. Write both sides
as Z/3. Every pi uniquely has the form pi(j)=i+epsilon*j, epsilon=+/-1.
Hence pi maps to the signed triangle frame (epsilon,i), and then, by the
preceding packet's F2, to e_i or1+e_i in normal-basis F8 coordinates.
Permuting the right side by j->j+1 becomes Frobenius rotation; reflection
j->-j becomes the affine reflection of that six-state model. This is an
exact equivariant map for the specified relabelling action.

If a second matching sigma closes the three routed components, the
alternating cycles of the two-coloured union are the cycles of sigma^-1*pi.
Parallel coincident edges are kept as a two-coloured two-cycle. For three
ports the full6-by6 component-count matrix has six distinct rows. Forgetting
pi while retaining only 'a perfect matching exists' therefore loses a
testable attachment response. More generally, querying sigma=pi gives k
components and is unique, so these closure tests recover a matching on
k+k named terminals.

There is also a different complement map: **K3,3 minus a perfect matching
is a Hamilton six-cycle**, giving a bijection between its six matchings
and six unoriented Hamilton cycles. A simple2-regular bipartite graph on
six vertices cannot split into two cycles, since each would need at least
four vertices.

Petersen likewise has six perfect matchings, but deleting any one gives
**two five-cycles**. The checker exhausts its matching set, and the audit
independently checks all five-edge subsets. Any two of the six matchings
intersect in exactly one edge; every Petersen edge appears in exactly
two matchings. Therefore its15 edges are in bijection with the15 pairs
of the six matchings, an explicit K6 address system. Under the Kneser
five-label transposition(0 1), the six matchings undergo three disjoint
swaps. The full five-label S5 action is faithful and transitive on these
six matchings, as the complete120-permutation check verifies. It is larger
than the specified one-side S3 action on K3,3 matchings. A transposition
has the same triple-swap cycle type in both six-point actions; that single
cycle type does not distinguish the groups. Equal cardinality does not
identify the full actions or the two complement maps.

These are concrete instances of the missing global component data:
both complements have degree two everywhere, but one is a spanning cycle
and the other is disconnected. No actual tournament complement is equated
with any of these graph-complement maps.

### Routed minors, not just nonplanarity labels

Use the standard Petersen edges

    01 04 05 12 16 23 27 34 38 49 57 58 68 69 79.

Contracting the five spokes(i,i+5), i=0,...,4, gives K5: the outer cycle
supplies distance-one pairs and the inner star supplies distance-two pairs.
But Petersen has no K5 subdivision, since a subdivision needs branch
vertices of degree four and Petersen has maximum degree three.

Deleting vertex0 gives an explicit subdivision of K3,3 with branch sets
{2,8,9} and {3,6,7}. Its nine routes are

    2-3, 2-1-6, 2-7;
    8-3, 8-6,   8-5-7;
    9-4-3, 9-6, 9-7.

The interiors1,4,5 are distinct. Retaining these routes certifies the
topological attachment; merely retaining the name of the minor does not.
This matches the earlier
[rooted Kempe reduction](../defect_transitions_20261020/petersen.md):
boundary support, extension multiplicity and actual disjoint wire routes
are different data. The correction in that packet already forbids deriving
a Petersen minor merely from a Petersen-like colour-support table.

Finally the inherited
[square-sum K3,3 bridge](../square_sum_bipartite_bridge_20261019/SQUARE_SUM_BIPARTITE_CLEANUP.md)
gives an actual arithmetic realization, with sides{1,193,673} and{3,168,483}.
The nine sums are squares and there are no within-side edges; these checks
are repeated here. Its proved affine cleanup sends A,B to25A+1,25B-1,
preserving all cross square predicates while deleting within-side ones
modulo5. We do not rerun its separate exhaustive height673 minimality proof.

## 6. What transfers, and what remains separate

### Incoming refinement: exact correlation order and boundary minimum

The concurrently completed [Session70 Vitali report](../session70_attachment_guard_alignment_20261026/session70-vitali/REPORT.md),
merged at `a8fcb1c840`, sharpens the address mechanism. Its A1-A4 are
**INHERITED PROVED** here, not new computations of this packet. The common
strong source of the C4 alignment family is rigid for every k: the lost
relative alignment cannot be blamed on source automorphisms. Moreover,
for every subset A of its k hidden bits there is a single native re-marking
guard with value

    g_A(b) = product_(i in A) b_i.

The full guard matrix is the Boolean zeta matrix, determinant one. Therefore
retaining moments of order at most t has rank sum_(j<=t) binom(k,j), while
exact linear evaluation of all guards on arbitrary signed laws requires
dimension 2^k. Even- and odd-parity laws agree on every proper-block
marginal but differ on the actual all-block guard. This realizes the older
Boolean-moment obstruction inside tournament continuation itself.

The costs must remain distinct: an individual alignment takes k bits;
an arbitrary probability law has 2^k-1 independent affine coordinates.
The March mixed-difference witness above is a concrete two-move precursor
of growing interaction order, not by itself a proof of that rank theorem.
Frobenius phase records cyclic transport, while these guards observe a
joint alignment in a specified frame; field degree alone cannot encode
these growing correlations.

The incoming [graph attachment report, Theorem A](../session70_attachment_guard_alignment_20261026/session70-graph-labeling/REPORT.md)
also supplies an exact minimum for a different alphabet: mandatory
nontrivial disjoint paths, 2r distinct named endpoints, and unrestricted
fresh connector paths. The minimal Boolean feasibility carrier is the
endpoint perfect matching, with (2r-1)!! states and r log2(r)+O(r) bits.
One connector joining an endpoint pair closes a forbidden intermediate
cycle in one matching and merges two paths in the other. This separates
different pairings. A fixed square-sum host may not realize those connectors,
so this minimum does not automatically transfer to its arithmetic arrivals.

Two additional inherited refinements make the graph comparisons operational.
First, the graceful S_k quotient survives simultaneous overwrites of at most
two colours as well as single overwrites, with no-ops allowed; such batches
may change the assignment, but their legality and resulting assignment
still depend only on that assignment. Second, in the displayed subdivision
Petersen-minus-0 -> K3,3, only two of the six quotient Hamilton cycles use
all three marked subdivided edges. Those two expand to the two spanning
nine-cycles. The other four skip a subdivision vertex. Retaining the minor
name or even the quotient cycle without the edge marks loses spanningness.

These refinements and the incoming guard proof are pinned in SOURCES.json.
Their receipts were read as inherited evidence; this packet's replay does
not purport to rerun the separate Session70 computations.

### Connection map

| Source -> target | Preserved predicate | Lost information / required sidecar | Decisive test |
|---|---|---|---|
| Tournament path guards -> trace support -> T | Every labelled arc, then all operations | None within finite tournaments; vertex names retained | All-size separator proof; invisible five-vertex arc |
| Raw Vitali flips -> fixed lambda fibre | Current lambda only | Current operation domains | The fourth corner leaves the fibre |
| Graceful frozen completions -> S_k | Replacement legality and updates | Signs, if actual edges are observed | Literal forest replacements through k=4 |
| G25 -> path guards | Current Hamilton words | Dormant edge1--15 | Identical arrivals yield35 versus18 paths |
| Two-terminal subsets of five labels -> Petersen | Disjointness | Root labels if quotienting by S5 | Explicit Kneser isomorphism |
| K3,3 matching -> signed frame -> six field words | Specified S3 action | Ambient graph operation | Matching-complement cycle versus tournament complement |
| Petersen -> minor name | Existence of a minor | Named branch sets and routes | Displayed K5 contraction and K3,3 subdivision |

The address lemma, the corrected guarded Vitali square, and the exact
graceful quotient explain growing alignment memory without invoking an
infinite-choice obstruction. G1 determines when the native constraints
recover apparently missing data. The useful next compression question is
therefore to specify a *weaker observation or operation menu* and minimize
that protocol, rather than seek a proper quotient that still implements
all named tournament re-markings.

Research cards used: META-PATTERNS' 'Search the statement before the
method', 'Type every analogy and every implication', and 'Retain local
incidence through counting'. Their counterindication here is important:
large fibres alone do not prove lost information, because G1 reconstructs
absent arcs from the native class. No new universal meta-pattern is promoted.
