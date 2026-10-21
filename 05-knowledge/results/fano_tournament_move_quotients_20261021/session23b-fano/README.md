# Tetrahedral Fano incidence, a boundary-orientation obstruction, and exact triangle-move quotients

Status: PROVED structural statements; FINITE-EXACT order-seven censuses; CITED classical Fano/code facts. Session23b research packet, not a claim of literature priority or a scarce theorem-ID reservation. The supplied seven-point-sphere geometry is being analyzed separately. This packet connects actual tetrahedral edge incidence, not an arbitrary seven-point count, to Fano incidence and tournament operations.

## 1. Inheritance and scope

Read current navigation, guardrails, protocol, and core-paper routes. Closest mechanisms:

- `01-canon/theorems/THM-211-paley-design.md`: Paley7's 14 cyclic triples are two block-disjoint Fano planes, each arc in two cyclic triples.
- `05-knowledge/results/arithmetic_braids2_20260917_fano_code.md`: Fano line orientation gauges form the affine binary Hamming coset `96+C`, with `C=[7,4,3]`; 16 orientations have Hamilton-path count 189 and 112 have count 171. This is inherited, not rediscovered.
- `05-knowledge/results/collatz_mod6_20260922_paley_fano_octonion_design.md`: audited 128-orientation census; two Fano systems inside each of the two represented regular tournament classes; all three regular7 classes; octonion sidecar and exact failure boundaries.
- `01-canon/theorems/THM-4529-petersen-and-heawood-families-in-paley-coordinates-anti-circulant-all-odd-tournaments.md`: the two six-vertex H=45 maximizers are different, namely `C3[TT2,TT2,TT2]` and Paley7 minus a vertex. This distinction is exactly the seam encountered below.
- `01-canon/theorems/THM-4553-the-point-stripping-ladder-psl27-borel-torus-and-the-collatz-frobenius.md`: Aut(Paley7) is the order21 affine square-multiplier group.
- `05-knowledge/results/marked_path_tournament_series_20261021/session23-symmetry/README.md`: distinguished ordered partial paths and converse-with-path-reversal; shape-blind scaffold counts do not erase their embedding sidecars.
- `01-canon/theorems/THM-1945-the-invariant-monoid-orbit-dictionary.md`: classical equal-score connectivity by cyclic-triangle reversals. We need no new all-order connectivity theorem.

Canonical hostile: seven cyclic Fano lines do not force Paley or octonions. Corrected near miss: the full Paley cyclic-triple set has fourteen, not seven, blocks. Least-used sidecar: the marked Fano plane, its omitted tetrahedral point, and the residual Hamming codeword after syndrome return.

Concept board: tetrahedral face/opposition incidence; coherent oriented boundary; Fano Hamming code; regular tournament triangle moves; marked subsets/paths; operation-relative generating functions. Every bridge below names its operation and losses. No LRC(14), JC, DC, Collatz, or Hadwiger-Nelson closure follows.

## 2. Exact tetrahedron-to-Fano map

Let the tetrahedron have vertices `0,1,2,3`. Its six edges, together with a formal seventh point `omega`, are the points of a Fano plane. Declare its seven lines to be:

1. The four triples of edges bounding a triangular face.
2. For each of the three pairs of opposite edges, that pair together with `omega`.

Every two adjacent tetrahedral edges occur in exactly one face triangle. Every two opposite edges occur in exactly one line of type2. The seventh point and any edge occur in its unique opposite-pair line. Thus every pair lies on exactly one of the seven triples: this is STS(7), hence the Fano plane.

Explicit binary coordinates use the reduced vertex-edge incidence matrix of K4, deleting vertex3's row:

    edge 03 -> 1=(100), 13 -> 2=(010), 23 -> 4=(001),
    edge 01 -> 3=(110), 02 -> 5=(101), 12 -> 6=(011),
    omega   -> 7=(111).

Here the displayed bit-vector order follows the first/second/third coordinate weights 1,2,4. Fano lines are triples with XOR zero. Equivalently, with ordinary integer bit labels, `L_a={x!=0: parity(a&x)=0}`:

    L1={2,4,6}, L2={1,4,5}, L3={3,4,7}, L4={1,2,3},
    L5={2,5,7}, L6={1,6,7}, L7={3,5,6}.

The four face normals are `{1,2,4,7}`. The three opposite-pair normals are `{3,5,6}`. Forgetting `omega` is exactly the classical matroid identity `F7\omega = M(K4)`. The tetrahedral S4 acts faithfully on its six edges and fixes `omega`; it is the order24 point stabilizer inside GL(3,2) of order168.

Connection contract: source is tetrahedral edge/face/opposition combinatorics; target is a marked Fano plane; map is above; preserved are facial cycles and opposite-edge pairing; lost is every Euclidean metric, sphere radius, center, and orientation. The distinguished formal point is essential. This does not identify a particular sphere's seven metric feature points with these seven combinatorial points.

### Real-representability boundary

A Fano plane cannot be represented by rank-three real linear dependence, or real projective collinearity. A short proof: normalize three independent projective points to `e1,e2,e3` and the points on their first two joining lines to `(1,1,0)` and `(1,0,1)`. The remaining incidences force `(1,1,1)` and `(0,1,1)`. The last required triple has determinant `-2`, impossible over R. It works in characteristic2. This is the classical Fano/non-Fano boundary.

Scope matters: seven generic points on a real sphere can be given seven *selected small circles*, one through each prescribed Fano triple, with no extra points on those circles generically. That realizes incidence restricted to the designated seven points; extra unlabelled circle intersections are not being identified with Fano points. It is an imposed incidence sidecar, not a canonical consequence of cosphericity. The nonrepresentability statement excludes linear/projective or great-circle realizations, not arbitrary selected curves.

Primary background: James Oxley, *What is a matroid?*, Figure8 and Proposition5.3, https://www.cs.cornell.edu/courses/cs6820/2022fa/Handouts/oxley-matroids.pdf . Its searchable source states the deletion identity and characteristic-two obstruction; PDF text extraction is partially garbled. The elementary proofs above are independent of that extraction.

## 3. Coherent tetrahedral boundary orientation is NOT the Paley orientation

Give the tetrahedron an orientation. Its oriented boundary induces a cyclic order on the three edge-vertices of each face. Orient those four Fano lines cyclically. The three lines through `omega` can be oriented independently, giving eight completions to a tournament on the seven Fano points.

**Structural obstruction.** No completion is Paley7.

Proof by twins: partition the six edge-vertices into three opposite-edge pairs. Coherent boundary orientations make the two edges in a pair have identical orientations to the other four edges. Between the three pairs the quotient is a directed triangle. Whichever internal directions are chosen, the six-vertex restriction is therefore `C3[TT2,TT2,TT2]`. The apex `omega` is added so that it and each opposite pair form a directed triangle. In particular, the two members of each opposite pair have exactly two common out-neighbors, both among the other four edges; they cannot both beat the apex. Paley7 instead has exactly one common out-neighbor per pair. Thus the completion is not doubly regular and cannot be Paley7. QED.

All eight completions are isomorphic: interchanging the two vertices of an opposite pair preserves all external arcs and swaps its two possible cyclic completions with the apex. Exact counting gives H=171 for every completion. This produces the third regular7 class, not the rotational H=175 class and not Paley H=189.

**Parity certificate and minimal repair.** Use orientation bit `m(a)=1` for the cycle through the increasingly ordered points of `L_a`, and zero for its converse. The inherited Paley affine code is `m0+C` with `m0=96` (bits6 and7) and

    C={c: XOR over a with c(a)=1 of a = 0}.

A coherent boundary has fixed bits `(m1,m2,m4,m7)=(0,1,0,1)` under the coordinate convention above, or their four-bit complement for opposite boundary orientation. Therefore the parity on the four face lines is even. Every Paley word has odd parity on those four coordinates: take the dot product with `omega=7` of the zero syndrome; the three variable opposite-pair normals have even dot product with7.

Reversing exactly one of the four facial cycles changes this parity to odd. The three remaining line normals `{3,5,6}` span the two-dimensional plane perpendicular to7 and have one relation `3 XOR 5 XOR 6=0`. Consequently precisely TWO of the EIGHT completions are now Paley. This is an exact minimal seam repair in the specified orientation carrier, not a deformation theorem about the original sphere.

Receipt: coherent boundary masks are `66,70,82,86,98,102,114,118`; all H171. Flipping faces1,2,4,7 respectively gives Paley mask pairs `{83,103}`, `{84,96}`, `{74,126}`, `{6,50}`.

## 4. The fixed-Fano line-flip graph is a cube with a K8 syndrome cover

Orient each of the seven Fano lines cyclically and reverse one complete line at a time. Because the lines partition the 21 unordered pairs, this presentation graph is exactly Q7. Every resulting tournament is regular: three lines pass through each point and each contributes one incoming and one outgoing arc.

Define the relative syndrome

    s(m) = XOR over a with (m(a) XOR m0(a))=1 of a.

Flipping coordinate a changes `s` to `s XOR a`. Since a runs over all seven nonzero vectors of F2^3, the syndrome map is a graph covering

    Q7 -> K8

with sixteen vertices over each syndrome and a bijection on every vertex's seven neighbors. Its deck translations are the [7,4,3] Hamming code. This is an operation-preserving bridge, not an isomorphism of the underlying objects.

The inherited/rechecked H census is exactly

    H(m)=189 when s=0; H(m)=171 when s!=0.

Thus the partition into Paley/non-Paley presentations is equitable, with row convention source-to-target and order `(Paley, non-Paley)`:

    Q_F = [[0,7],[1,6]].

Spectrum is `7,-1`. A non-Paley presentation has one and only one line-flip repair to Paley. A Paley presentation has seven line-flip departures, all non-Paley.

An exact Walsh expression is

    H(m) = [693 + 9 sum_(u!=0) (-1)^(u dot s(m))]/4.

This expression's H values rely on the inherited, independently rerun finite census, not on the Hamming-code theorem alone.

### Holonomy: returning to the syndrome is not returning to the presentation

A triangle in K8 with increments `a,b,a XOR b` lifts to three flips whose total codeword is

    e_a + e_b + e_(a XOR b).

That is a NONZERO weight-three Hamming word. The lift ends at a different presentation, even when it ends in the same syndrome and has the same H. Every K8 triangle behaves this way. The seven possible voltage supports are the seven Fano lines on the *line-coordinate* index set.

The parity-enriched quotient `m -> (s(m), wt(m) mod2)` is the graph `K8,8` minus a perfect matching: edges go `(s,p)->(s XOR a,p XOR1)`, with a nonzero. Its fibers have size8. This retains the original cube's bipartiteness, which the syndrome quotient forgets.

Starting at a Paley presentation, the number of length-t line-flip words that end anywhere in the Paley locus is

    R_t = [7^t + 7(-1)^t]/8,
    sum R_t z^t = (1-6z)/[(1-7z)(1+z)].

At t=3 there are42 such words but ZERO returns to the initial presentation. The coarse count must not be called a closed-walk count in Q7. The codeword voltage is the missing sidecar.

For uniformly random line flips, `E[H_t]=693/4 +(H_0-693/4)(-1/7)^t`. From any non-Paley presentation the first hitting time of the Paley locus is geometric with success probability1/7 and mean7. This is an observable/locus statement; the full cube walk remains periodic and does not converge in distribution without laziness.

## 5. The cube sits inside the full regular-tournament move graph

Use cyclic-triangle reversal as the move, NOT single-arc reversal. Let R7 be the graph of all labelled regular tournaments on seven fixed vertices, adjacent when one directed triangle is reversed.

**PROVED:** the fixed-Fano Q7 is an induced, isometric subgraph of R7. A Fano-line reversal remains inside it. Any non-line cyclic triangle uses one pair from each of three distinct Fano lines; reversing it makes those line-triples transitive and exits the cube. If two cube words differ in k coordinates, their tournaments differ in exactly3k arcs. Any triangle reversal changes3 arcs, so at least k moves are necessary, and the k changed Fano lines attain the bound.

Every regular7 tournament has14 cyclic triples, because the transitive triple count is `7*C(3,2)=21` and there are35 triples in all. Therefore every cube vertex has seven internal and seven external neighbors. The Hamming decoder controls only the internal move set.

### FINITE-EXACT full quotient

An exhaustive outdegree-constrained row generator produces all2640 labelled regular7 tournaments. All7! relabelings of one representative in each H channel independently certify exactly three isomorphism classes. In order `(H189,H171,H175)`:

    class sizes: (240,1680,720),
    triangle-neighbor matrix Q_R = [[0,14,0],[2,9,3],[0,7,7]].

Every individual vertex in each class has exactly the displayed neighbor counts. Detailed balance is independently checked from all directed graph edges. Characteristic polynomial:

    det(XI-Q_R)=(X-14)(X^2-2X-14).

Thus the quotient eigenvalues are `14,1+sqrt15,1-sqrt15`. Uniform random cyclic-triangle reversal has stationary class probabilities `(1,7,3)/11`, not equal class weights, with stationary mean H=`1911/11`. Expected first hitting times of the Paley class are10 from H171 and12 from H175, from the two first-step equations. They differ from the restricted Fano cube's mean7 because the move carrier has changed.

For a fixed starting Paley tournament, let W_t be the sum of H over endpoints of all length-t triangle-reversal move words. Then, using the exact three-state quotient,

    sum W_t z^t = (189-630z-1470z^2)/[(1-14z)(1-2z-14z^2)].

All three factors/eigenmodes survive. This is a genuine rational series controlled by the operation quotient, adjacent to the session's marked-path denominator work. It does not identify this move-word series with the session's tournament-order series, path-length series, or renewal model.

Boundary: the three-state quotient is NOT a universal quotient for regular tournaments at larger orders. The special n=7 classification makes it exact here. A future order9 experiment must first test whether a proposed observable partition is equitable; sharing H alone is not a theorem.

## 6. Marked subset and partial-path consequences

For Paley7, `Aut={x->ax+b: a in {1,2,4}, b in F7}` has order21, with21 anti-automorphisms given by nonsquare multipliers. The following counts are independently exhaustive over all5040 permutations and all128 subsets.

A subset is *regular for the group action* when its setwise stabilizer is trivial. This usage is unrelated to a regular tournament's equal outdegrees.

- Paley7 has84 such subsets:21 each of sizes2,3,4,5. The size3 regular subsets are the21 transitive triples; its14 cyclic triples have stabilizer3.
- The bare Fano collineation group of order168 has NO regular subset, since every subset orbit has size at most35, smaller than168.
- The order42 group of automorphisms and anti-automorphisms of Paley7 also has NO regular subset, since35<42. Setwise marking therefore cannot break all converse symmetry at this order even though directed orientation alone has many regular subsets.

This illustrates why marking an ordered directed path is not equivalent to marking its vertex subset. Any Paley path with at least two positions has trivial ordinary stabilizer, because its first ordered arc does. Enumerating directed partial paths gives:

    k:                 1   2   3    4    5    6    7
    labelled paths:    7  21  63  147  273  315  189
    Aut-orbits:        1   1   3    7   13   15    9
    converse-fixed:    1   1   1    3    3    3    3
    converse-merged:   1   1   2    5    8    9    6

Converse acts by reversing the ordered path as in the session23 marked-path packet. These are contributions from ONE underlying tournament, not counts over all tournaments. For k>=2 only the seven involutory anti-actions can fix a path under reversal: its square fixes two vertices and hence is the identity. The fourteen anti-actions of order6 have no such fixed path.

### A single Fano plane cannot carry a tournament anti-symmetry

The root session observed the following useful general obstruction. Every tournament anti-automorphism has at most one fixed point and all other vertex cycles have lengths2 modulo4. For odd n=3 modulo4, the number of its nontrivial cycles is odd, hence its permutation sign is negative. In contrast GL(3,2) on seven nonzero binary vectors lies in A7: elementary transvections generate it, and each transvection has cycle type `1^3 2^2`.

Therefore no anti-automorphism of ANY seven-vertex tournament preserves a single labelled Fano plane. In Paley7 the two Fano planes in the cyclic-triple design must be exchanged by every anti-action. Keeping one plane is a genuine orientation sheet, not disposable decoration. A general diagram of selected small circles need not respect that sheet or that group.

## 7. What can and cannot transfer to the open questions

- **Metagraph series: actual transfer.** A target-preserving finite quotient gives exact endpoint-observable series and recurrence denominators. Here the quotient is explicitly proved/certified equitable. The K8 triangle voltage warns against mistaking quotient return for full-state return.
- **Marked-subset counts: actual transfer.** Retain the acting group and whether converse is included. Paley has84 ordinary regular subsets but zero extended regular subsets; a bare Fano group has zero. A cardinality-only replacement loses exactly the symmetry condition.
- **Tetrahedral geometry: actual combinatorial transfer, metric loss.** Edge-face-opposition data becomes Fano incidence; coherent surface orientation selects the H171 class and obstructs Paley. One facial seam reversal repairs it in exactly two ways. Any seven-point-sphere theorem still needs its own metric equations and identification of features.
- **LRC and coloring: no target map.** A triangle-reversal graph has no clock, separation threshold, speed magnitudes, distances, or forbidden-unit-distance predicate. Until those observables are transported, none of the Fano or Paley results settles LRC or Hadwiger-Nelson. The inherited `THM-640-the-paley-bridge.md` has its own prime/residue scope; modulus14 does not become a Paley tournament.

Precise next hypothesis, not promoted: for larger regular-tournament score fibers, identify the smallest operation-faithful refinement of `(H, cyclic triangle data, marked design)` that yields an equitable quotient. A cheapest hostile is two same-H order9 tournaments with different triangle-neighbor H distributions. A separate coding question is when orientations of a fixed Steiner triple system contain a perfect-code locus for an extremal tournament observable; n=7 alone does not predict it.

## 8. Reproduction and evidence

`check.py`: all128 cyclic orientations of fixed Fano lines; H by subset DP AND all7! orders for each orientation; all Fano and Paley groups by permutation testing; all subsets and marked partial paths; seam completions; code-cover and holonomy tests; return counts through t=12.

`global_triangle_check.py`: all2640 labelled regular7 tournaments by degree-constrained recursion; independent full relabeling orbits; every cyclic-triangle neighbor of every tournament; complete exact quotient and detailed-balance checks.

Run each with Python3 and `-O`; both pairs of JSON receipts are byte-identical. No `assert` gates are used. Runtime is about one second per script in this environment. No random sampling, floating arithmetic, external dependency, or uncommitted repository modification is needed.

Primary octonion background, not a citation for the new move-graph calculation: John Baez, *The Octonions*, Section2.1, https://math.ucr.edu/home/baez/octonions/node4.html . This specifies the signed Fano-line multiplication and binary indexing; arbitrary seven independent line orientations need not give octonions.
