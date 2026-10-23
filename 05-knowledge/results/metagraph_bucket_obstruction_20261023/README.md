# Incidence buckets force the six- and seven-vertex metagraphs

23 October 2026. **PROVED**: the overlap/disjoint-triangle identity and the
capacity-peeling criterion below. **FINITE-EXACT**: complete class and port
checks through order seven, complete signature census at eight. **OPEN**:
a fixed-size local signature determining the merged metagraph at every size.
No new A000568 enumeration speedup is claimed.

## Inheritance and the question actually answered

Closest proved mechanism: THM-4731, `01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`,
retains the staircase reflection, blue fixed points, and black twin lines.
Closest finite determination result: THM-4739,
`01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md`.
Its target score/curl/second-score signature L3 leaves 48,000 integer
transport tables at order six, three projecting to the native merged graph.
Canonical hostile: the explicit reflected switch in this packet retains
those margins but changes four merged black multiplicities.
Corrected near miss: MISTAKE-648 forbids replacing paired-black multigraph
degree by simple support degree or counting the two converse classes twice.
Least-used sidecar: correlations of cyclic-triangle incidence through EDGES,
then through rooted strong-four subtournaments.

The live board is: (1) source/path/target ports; (2) fixed-point versus free
involution orbits; (3) vertex and edge incidence; (4) transportation capacity;
(5) disjoint odd-cycle packing. The pasted incidence-moments proposal led to
(3), while the incoming session35 squarefree-port work supplied (5).

Our graph is the **converse-merged complement-line multigraph**, with blue
and black multiplicities and loops retained. The paired-black graph divides
every black multiplicity by two. It is not the single-arc-flip metagraph.

## 1. An overlap is exactly the small strong tournament

For a tournament T write c3 for the number of cyclic triples, c_v for the
number through vertex v, and c_e for the number through unordered pair e.
Let S4 be the unique strong tournament on four vertices. Then

    M2(T) := sum_e binom(c_e,2) = number of induced S4 copies.

Proof: two cyclic triangles sharing an edge u->v have their other vertices
w,z satisfying v->{w,z}->u. Either orientation of wz gives S4. Conversely
S4 has exactly two cyclic triples, sharing exactly one edge. This is a
bijection, without automorphism denominators.

Let D33 count UNORDERED pairs of vertex-disjoint cyclic triangles. Classify
pairs of distinct cyclic triangles by intersection size 0,1,2. The vertex
incidence sum counts these with weights 0,1,2, and M2 counts precisely the
intersection-two pairs. Thus, at every order,

    D33 = binom(c3,2) - sum_v binom(c_v,2) + M2.

The identity is elementary; no priority claim is made. The new use here is
that a missing strong-four incidence coordinate is EXACTLY a missing
triangle-disjointness coordinate when vertex incidence has been retained.

This connects to `session35_dimensions_dualities_20261023/recursions/CONTEXT-AND-PORTS.md`
and THM-505's inherited six-vertex formula

    H = 1 + 2(c3+c5) + 4D33.

Its adjacency-cospectral witness has H=29 versus25, entirely accounted for
by D33=1 versus0. Our overlap statistic repairs the same type of lost
compatibility information; a scalar spectrum or a list of separate motif
counts need not retain which motifs overlap.

## 2. Incidence profiles that retain the association between measurements

Put s_v=outdegree(v), and let v4_v be the number of induced S4 copies
containing v. Define multisets

    J(T)  = multiset_v (s_v,c_v, multiset_{v->w}(s_w,c_w)),
    J4(T) = multiset_v (tag_v, multiset_{v->w} tag_w),
    tag_v = (s_v,c_v,v4_v).

These are ordinary isomorphism invariants. They are not asserted to be
complete at arbitrary order. Computing J4 needs only triples, four-sets,
and neighbor-tag multisets; it does not enumerate Hamiltonian paths.

Complete ordinary-class counts:

| n | classes | inherited L3 | J | (J,M2) | J4 |
|---|---:|---:|---:|---:|---:|
| 4 | 4 | 4 | 4 | 4 | 4 |
| 5 | 12 | 12 | 12 | 12 | 12 |
| 6 | 56 | 52 | 55 | 56 | 56 |
| 7 | 456 | 367 | 448 | 454 | 456 |

Entries after the class column count distinct signatures.

At n=6 the sole J collision has fixed-path masks220 and108, outgoing rows

    [48,41,19,5,10,28], [48,49,3,7,12,28].

Both have H=45, |Aut|=3, ordinary fibre h=15, and scores (2,2,2,3,3,3).
Their blue counts are3 and9, so they are NOT a same-blue-fibre collision.
Their M2 values are9 and12. The edge-incidence multisets are
1^6 2^9 and0^3 2^12 respectively. Recording one extra overlap moment
separates them.

At n=7, (J,M2) still collides on masks (1337,1468) and (1430,921).
J4 separates both pairs. The first pair is self-converse with H=133 and137;
the second pair consists of opposite orientations of one non-self-converse
merged node. These are explicit boundaries to the weaker signature,
not claims that the final J4 signature fails at seven.

## 3. Why the pigeonholes force lines, not merely node labels

Use the inherited descending base path n-1 -> ... ->0. A fixed-path tiling
with closed endpoint edge0->n-1 determines a closed source port. If its
endpoint scores are (a,b), complementing all non-path arcs sends it to an
open target port with scores (n-a,n-2-b). Blue/black color is preserved.
The target invariant J, (J,M2), or J4 is computed on that actual complement.

Create a bipartite transportation problem:

* a row retains source ordinary class, endpoint scores, color, and the
  computed target invariant; its capacity is the number of those source
  path-orbits;
* a column retains target ordinary class, target endpoint scores, and
  color; its capacity is the number of its open path-orbits;
* an entry is allowed only when both the target invariant and endpoint
  type agree.

**Capacity-peeling criterion (PROVED).** If a row has only one remaining
column, its entire remaining mass must enter that column. Dually, a column
with one remaining row forces that much mass from the row. Subtract,
remove exhausted vertices, and repeat. If all mass is exhausted, the
integer transport table is unique. Proof is induction: each subtraction
is forced in every feasible table. Native complement transport supplies
feasibility. This is only a sufficient uniqueness algorithm; getting stuck
does not prove another feasible table exists.

Complete fixed-path replay gives:

| n/signature | source rows | allowed entries | unresolved mass after peeling |
|---|---:|---:|---:|
| 6/L3 | 403 | 440 | 67 |
| 6/J | 420 | 421 | 0 |
| 6/(J,M2) | 421 | 421 | 0 |
| 7/J | 13,826 | 14,087 | 425 |
| 7/(J,M2) | 13,848 | 13,971 | 222 |
| 7/J4 | 13,896 | 13,896 | 0 |

Thus **J and the actual endpoint/color capacities already determine every
ordinary source-to-target line count at six**, despite J not identifying
all ordinary target classes. One ambiguous blue source bucket of mass2
must split1+1 between two singleton target buckets. Every other ambiguity
is forced by a zero or singleton residual capacity. Adding M2 removes the
last class ambiguity directly. **J4 determines every aggregated line multiplicity at seven.**
Converse merging then gives the requested colored merged multigraph.
This upgrades the finite L3 determination boundary; it is not an all-order
classification or a replacement for source/path compatibility.

## 4. An explicit alternative arrangement ruled out by the new coordinate

The receipt includes a reflected2x2 switch at n=6. On merged node IDs from
`electric_n6.json`, its black closed-to-open multiplicity changes are

    (16,12):+2, (16,18):-2, (17,12):-2, (17,18):+2.

All blue entries stay unchanged. Every source L3 bucket and every target
endpoint/color bucket has exactly the same total. Black multiplicities
stay even, and the class-level converse-reflection identities hold. Hence
all node tiling amounts, closed/open amounts, blue amounts, black amounts,
and the inherited L3 port constraints survive. Nevertheless this is not
the tournament metagraph: the capacity-peeling theorem makes the J-refined
table unique, so this changed table cannot admit a J-compatible refinement. The receipt retains all eight ordinary-class
edits and literal tiling masks witnessing the native entries.

This is the requested distinction between counting slots and realizing
slots. Parity and margin laws admit the counterfeit; correlations carried
by the actual source path exclude it. The source/path/target map is
preserved in the accepted table; the weak quotient loses neighbor-incidence
association. J is the missing sidecar, and the tiny switch is the decisive
hostile test.

## 5. Order eight: the boundary and the cost of a repair

The complete6880-class census has6878 different J4 signatures. The only
collisions are two pairs of upper-triangle canonical words; converse
exchanges the two pairs, not the two members of either pair:

    (27467784,27467793), (6427781,6427778).

The first pair has literal outgoing rows in `receipt.json`; direct subset
Hamiltonian DP gives H=523 and527, and both automorphism groups are trivial.
Thus **J4 plus h=H/|Aut| separates all6880 ordinary classes**. This is not a
same-fibre counterfeit: it is a sharp failure of J4 alone. The known node
margins already record h, but using h as a transported target invariant
requires knowing the actual complement's h. It is not free information
obtained solely from those margins, and we make no speedup claim for it.

There is a cheaper repair than transporting h. Retain the vertex-level
first-round records before taking their multiset:

    rho1(v) = (tag_v, multiset_{v->w} tag_w),
    J4^(2)(T) = multiset_v (rho1(v), multiset_{v->w} rho1(w)).

The complete n=8 check finds **6880 distinct J4^(2) signatures**. Thus a
second neighborhood-incidence round removes both collision pairs, without
computing Hamiltonian paths or blue multiplicities. It asks which first-
round records are joined, information discarded by J4's outer multiset.
This is an exact success at eight, not an all-order termination claim.

The n=8 statement is a class-invariant census, not a full replay of its
2^21 fixed-path transport matrix. Class injectivity implies that transporting
J4^(2) (or the more expensive (J4,h)) would identify each target, but that conditional reconstruction
statement does not claim the expensive port computations were performed.

## 6. Reproduction and audit boundaries

**Final incoming bridge:** Session36's dyadic electrical law in
`session36_information_curl_contexts_20261023/metagraph/DYADIC-ELECTRICAL-METAGRAPH.md`
sees black wiring at the half-Laplacian layer. The separate
`dyadic_bridge.py` consumer checks our pinned switch: the Laplacian is
unchanged modulo2, yet its half pairing at unit potentials e16,e12 changes
from0 to1. Both vertices are blue-isolated. This is a second exact way to
detect this switch with retained class labels, not an all-order sufficient
wiring invariant. Its output is `dyadic_receipt.json`; normal and optimized
runs agree. The primary incidence checker and its receipt remain frozen.

Run from repository root:

    python -B 05-knowledge/results/metagraph_bucket_obstruction_20261023/check.py
    python -O -B 05-knowledge/results/metagraph_bucket_obstruction_20261023/check.py

Default runs compare against the frozen receipt and never overwrite it.
Use `--write` only for intentional receipt regeneration.
The checker uses no solver and no removable assertions. It independently
expands every permutation of each representative at n=4..7, proves the
ordinary classes disjoint, and covers every labelled tournament:2^6,2^10,
2^15,2^21. It then replays every fixed-path tiling (8,64,1024,32768), checks
H/Aut and blue-fibre counts, and matches every frozen blue/black entry.
At n=8, independent score-canonical minimization checks all6880 classes
using403193 score-preserving permutations. Their distinct canonical forms
and orbit-size sum2^28 establish exhaustive coverage. It independently
checks the four hostile H values by subset DP.

All input hashes, exact class collision masks, port counts, deterministic
peeling trace hashes, counterfeit edits, and positive/hostile controls are
in `receipt.json`. The optimized run must produce the same receipt. The
finite signature successes must not be extrapolated past their stated
orders. No theorem ID was reserved in this subagent packet.

