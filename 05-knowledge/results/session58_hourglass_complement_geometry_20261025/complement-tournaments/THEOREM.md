# Complement pairs, ordered observers, and the two different tournament charts

Status: PROVED elementary statements below; FINITE-EXACT n <= 5 controls; RECOVERED Hamilton-marked trijection and Burnside framework, not novelty claims. All graphs are simple, may be disconnected, and their vertex count is part of the object. No LRC, magic-square existence, Keller/JC, or complete-invariant claim follows.

## 1. Same-size chart: the user's pair is exact with an observer

Let < be a fixed linear order on an n-element vertex set V. Define G_<(T) by declaring {u,v}, u<v, an edge precisely when v -> u in T. Conversely orient a graph's edges backward and its nonedges forward. These are inverse bijections between all graphs and all tournaments on the same ordered set, with 2^{n choose 2} objects on each side.

Then G_<(T^op) = complement(G_<(T)). Consequently unordered complementary graph pairs correspond exactly to converse pairs of tournaments WITH THIS FIXED OBSERVER ORDER. For n>=2 there are 2^{binom(n,2)-1} such literal pairs; for n=0,1 the unique graph equals its complement and the count is 1. Swapping the selected member loses precisely the orientation/converse choice in this chart. Reversing the observer order also complements G, so the unordered pair only requires an unoriented order (a linear order modulo reversal).

This is not a bijection between ordinary unlabelled tournament classes and ordinary unlabelled graph-complement pairs. The obstruction is equivariance: relabelling tournament coordinates by a permutation p transforms the backward-edge bit vector x into p*x XOR I(p), where I(p) records pairs whose relative order p reverses. Ordinary graph relabelling uses p*x without this inversion cocycle. Transporting BOTH the graph and the order is equivariant; retaining a fixed ambient order while forgetting its relation to the vertices is not.

### Smallest hostile

On 0<1<2 take G with its only edge 01, versus G' with its only edge 02. They are isomorphic (one edge and an isolate). Their tournaments have arcs respectively (1->0,0->2,1->2), a transitive triple, and (0->1,2->0,1->2), a directed cycle. These remain distinct after converse identification. Thus even replacing G by its unordered pair with its complement does not repair unlabelled descent. At n<=2 every tournament is in the unique class, so n=3 is minimal.

The reverse descent also fails: different order presentations of the transitive triple yield both the empty/full graph pair and the one-edge/two-edge pair. Therefore neither direction of the proposed unlabelled correspondence is a function induced by this chart.

### Sidecar scope and minimality

A linear order is sufficient for the selected-graph/oriented-tournament bijection. An order modulo reversal is sufficient for the unordered-complement-pair/converse-pair bijection. For the universal literal decoder on a fixed named vertex set, two distinct orders disagree on a pair, so decoding the same fixed graph gives a different direction on that pair. Two orders not equal or reverse do not give tournaments equal or converse: the set of directions changed is neither empty nor the complete pair set. Thus these observer states are necessary for this universal literal decoder, up to the stated reversal equivalence. This is NOT a lower bound on arbitrary encodings of isomorphism classes: some particular graphs/tournaments have symmetries, and continuation-equivalent consumers may use smaller sidecars.

## 2. The inherited marked-Hamilton trijection has one fewer graph vertex

Fix P=(v_0 -> ... -> v_{n-1}) in T. Its n-1 path arcs are frozen. The other binom(n-1,2) chords v_i v_j, j>=i+2, map to the pairs {i,j-1} of the ordered set {0,...,n-2}. The corresponding graph bit is 1 iff v_j -> v_i. This gives the recovered trijection:

    (T, marked directed Hamilton path P) up to simultaneous isomorphism
    <-> graph on the ordered set of n-1 positions
    <-> tournament on that ordered set of n-1 positions.

It is not the same-size chart in Section 1. Marking P makes the order intrinsic to the marked object; forgetting P generally loses information. A fixed tournament T contributes H(T)/|Aut(T)| marked types because an automorphism fixing a directed Hamilton path fixes every vertex, hence the action on its paths is free.

In this chart graph complement kappa reverses only non-path arcs:

    K_P(T) = T^op XOR E(P).

It does NOT mean tournament converse. Converse, followed by reading the reversed marked path, instead reflects the (n-1)-vertex ordered graph by i -> n-2-i; call this sigma. Their product sigma*kappa corresponds to reversing only the path arcs and reading the reversed path. This is the inherited Klein dictionary of THM-4733, Proposition A.

The n=3 hostile is immediate: the empty graph on two positions decodes to a transitive triple; its complement decodes to a cyclic triple. Converse cannot do that. At n=4 the eight marked presentations have unmarked-class multiplicities 1,1,1,5, not a 3-1-3 count. Complement K_P also fails to be independent of P already at n=4 (see certificate's chart: masks 2 and 3 both decode to the strong class, while their complemented masks 5 and 4 decode respectively to the strong class and a triangle-plus-sink/source class).

## 3. Exact small counts and the genuine A000568 bridge

The exhaustive universe is every bit vector on all pairs, with no connectivity filters. Graph equivalence uses ordinary S_n relabelling; tournament equivalence uses the inversion-twisted action above.

| n | graphs | graph/complement classes | tournaments | self-converse tournaments | tournament/converse classes |
|---|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 1 | 1 |
| 2 | 2 | 1 | 1 | 1 | 1 |
| 3 | 4 | 2 | 2 | 2 | 2 |
| 4 | 11 | 6 | 4 | 2 | 3 |
| 5 | 34 | 18 | 12 | 8 | 10 |

The equality 2=2 at n=3 does not validate the failed map. The mismatch at n=4 rules out any bijection between those two entire finite sets, not just this chart.

For completeness, Burnside gives t(n)=1/n! sum_p Fix_T(p). If p has an even cycle, a power of p swaps a pair of vertices, forbidding an invariant orientation. If all its cycles are odd, no power swaps a pair, each orbit of unordered pairs has two consistent orientation choices, and Fix_T(p)=2^{c_2(p)}. Thus t(n)=1/n! sum_{p all cycles odd}2^{c_2(p)}. This is the inherited/classical tournament enumeration mechanism, not the ordinary graph Burnside sum, which has no odd-cycle exclusion. The displayed tournament column agrees with A000568, independently checked at https://oeis.org/A000568 on 2026-10-25. That entry identifies the sequence as unlabelled round-robin tournament counts; no further paper theorem is imported here.

Let s(n) count self-converse tournament classes. Since converse is an involution on classes, its quotient has (t(n)+s(n))/2 classes. One must not just divide t(n) by 2: every self-converse class is a fixed point of this class-level involution, although no literal tournament on n>=2 equals its converse. Marked path symmetry and class-level self-converseness likewise must not be conflated.

## 4. A lawful order-free replacement: average the observer, not the vertices

For a tournament T on V, define a histogram A_T indexed by unlabelled complementary graph-pair types q:

    A_T(q) = #{linear orders < of V : {G_<(T), complement(G_<(T))} has type q}.

This is an intrinsic invariant of the converse class: a tournament isomorphism transports orders bijectively; converse swaps the pair members. Its mass is n!, and its normalized law is an order-randomized observer. Coefficients are even for n>=2 because order reversal is a free pairing and does not change q. The measure retains multiplicity of observers rather than pretending a chosen one descends. This is a legitimate orbit-averaging repair, not a graph/tournament bijection.

FINITE-EXACT: these histograms separate all tournament converse classes through n=5 (respectively 1,1,2,3,10 distinct profiles). The actual profiles, not only these totals, are in certificate.json. No all-n completeness is asserted; a future collision would refute completeness without harming invariance. Already the n=3 transitive and cyclic classes differ by whether the empty/full graph-pair type occurs. A000568 counts classes, whereas A_T is an invariant on those classes; their roles differ.

## 5. Shared pins, overlap, and the honest hourglass connection

The exact transferable operation is labelled restriction and compatible gluing, not a shape analogy. On a single globally ordered vertex set, the chart commutes with restriction to subsets, and graph complementation on a subset corresponds to reversing that subset's tournament arcs. Two induced pieces can be glued only after their common vertex identification, common order, and common pair bits agree. Pairs crossing between the two private vertex sets are absent from both pieces and require additional data. Abstract graph types or uncorrelated averages cannot supply these pins or cross-pairs.

Observer choices used by overlapping pieces must be restrictions of the SAME global order. Independently averaging the pieces can invent incompatible order choices, just as independent copying of a shared latent state differs from diagonal copying in the Session56/57 correlated-gluing carrier. Even pairwise compatible partial orders can impose a cyclic global comparison (a<b, b<c, c<a) and admit no global linear order; an explicit global extension condition is needed. Conditional/diagonal averaging over compatible global observers is the safe operation.

This is a precise shared-provenance connection with the hourglass/gluing lane, not an arithmetic realization result. A positional seven-cell 3-1-3 hourglass, a historical Keller inverse-fiber 3,1,3, and tournament path counts are different typed objects. The parent identified the historical FRAGMENT-DECODE-mu3-depth2-tree-klein-S691.md correction; this lane imports none of its Keller claims and provides no predicate-preserving map from that 3-1-3 to tournament classes. A numerical shape match is not evidence of one.

## 6. Reproducibility, inheritance, and stopping boundary

Run python verify.py and python -O verify.py. The standard-library verifier uses explicit exceptions rather than assertions. It exhausts all graphs/tournaments through n=5, enumerates relabelling orbits, checks every permutation's fixed-tournament count against the independent Burnside formula, records the minimal graph-orbit hostile, reproduces the eight n=4 Hamilton-chart rows, and checks converse invariance/mass of every averaged profile. No huge census is needed.

Sources at pinned repository commit 44f15df8c4a24e5624d399862c0cdb1769b2cbb1:
- 01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md, Section 0 and Proposition A: explicit graph mapping and Klein action.
- 05-knowledge/results/trijection_four_nodes_three_edges_20261020.md, inheritance and Theorem 1: marked line/cycle/star carriers, free action, Burnside distinctions.
- 01-canon/MISTAKES.md, 2026-10-20 marked-trijection correction: never transfer path words without relabelling.
- Session57 path-inheritance source S11.md: correlated gluing, diagonal versus independent copying, separator scope, and continuation-equivalent sidecars. This is a local source snapshot, separately hashed in SOURCE_PINS.json.

The trijection-reset subdirectory linked by the extant theorem/restatement is missing in this checkout. We rely only on the extant explicit definitions/proofs, not the missing linked files. This task modifies no repository or remote; all artifacts are confined to its assigned directory.
