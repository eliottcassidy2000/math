# Square-sum paths and magic squares: exact carriers and failed identifications

25 October 2026. Status: **PROVED** elementary carrier equivalences and distinctions; **FINITE-EXACT** standard-library certificates. No integer magic square of nine distinct squares is constructed or excluded. No new finite-field classification is claimed.

## 1. The objects must first be separated

For an integer n≥1, the ordinary square-sum graph G_n has vertices 1,...,n and an edge {x,y} precisely when x≠y and x+y is an integer square. A square-sum ordering is a Hamiltonian path in G_n. Its arithmetic constraints are pair sums, and its global combinatorial constraint is one visit to each vertex. Neither equal triple sums nor a 3×3 line layout is part of the definition.

For a commutative ring R, put Q(R)={r²:r∈R}, including zero. A magic square OF squares is an injective assignment of nine grid cells to Q(R) for which all three rows, all three columns, and both diagonals have the same total T. Over the integers the entries are nonnegative. Over F_p they are field elements, not their integer representatives. A root-decorated version additionally chooses r_v with r_v²=x_v. Root choices are extra data and are not the same as distinct squared entries.

Three separate pair graphs therefore arise:

1. Values x,y with x+y∈Q(R), optionally restricting x,y∈Q(R).
2. Roots r,s with r+s∈Q(R).
3. Roots r,s with r²+s²∈Q(R), the pullback of the value graph through squaring, with collision exclusions specified separately.

They are not interchangeable. For integer roots 3,4 the sum of squares is 25 but the root sum is 7; for roots 1,8 the root sum is 9 but the sum of squares is 65. Over odd fields squaring also identifies r and −r. A graph argument on roots must specify which predicate and what quotient or decoration it keeps.

## 2. The total and opposite-pair identities hold without division

Write c for the central ENTRY, not its root. Let T be the common total. The four lines through the center yield four opposite pairs of sum T−c. Hence the sum of all nine entries is 4T−3c. The three rows also give 3T. Subtraction proves T=3c in every commutative ring. Each opposite pair therefore has sum 2c.

Taking the top corners a,b gives the complete parametrization

    a         3c−a−b     b
    c−a+b     c          c+a−b
    2c−b      a+b−c      2c−a.

Substitution verifies sufficiency; the preceding identities and row/column equations give necessity. Square membership of all entries and their pairwise distinctness remain additional conditions. This is inherited from the Session55 CRT proof, not a new parametrization.

If c=e² is a nonzero integer square, T=3e² is not an integer square: its 3-adic valuation is 1+2v_3(e), which is odd. Similarly 2e² is not an integer square. If c=0, nonnegative entries and T=0 force every entry to vanish. Thus requiring BOTH squared entries and a squared common total makes the distinct-entry integer problem impossible; this does not resolve the classical problem, which never required T to be square.

The separate variant “integer entries with a square common total” is perfectly meaningful: multiplying Lo Shu by 15 yields total 225, without making all entries squares. A third variant allows each line total to be a square but does not require the totals to agree. Neither variant should be silently substituted for the classical one.

Over an odd field, for c=e²≠0 the total is square exactly when 3 is square, and the opposite-pair sum is square exactly when 2 is square (zero counts as square). The zero-center stratum must be retained separately. Cain's F29 example lies in that stratum; division by the center would discard it.

## 3. Exact positive carrier: the fixed-total triple hypergraph

For T∈R define H_T on vertex set Q(R), with a hyperedge for every THREE-ELEMENT set {x,y,z} satisfying x+y+z=T. Let L be the fixed grid hypergraph on cells 0,...,8, with edges

    012, 345, 678, 036, 147, 258, 048, 246.

**Theorem A.** Magic squares of nine distinct squared entries with total T are exactly injective, edge-preserving maps L→H_T. This is a bijection of labelled objects, not of unlabelled isomorphism classes.

Proof. Such a map supplies nine distinct square values and sends each of the eight specified lines to a triple summing to T. Conversely a magic array is exactly such a map. The labels retain the positions, values, and which line constraints share each cell. The image need not be an induced copy: additional equal-sum triples among its nine values are allowed. Root-decorated squares are obtained by independently attaching square-root witnesses to these values; roots are not confused with entries.

The ordinary graph shadow of H_T keeps a pair whenever it extends to an edge, but it is not exact by itself. In F29 at T=0, its vertices 0,1,4 form a triangle: the pair witnesses are respectively {0,1,28}, {0,4,25}, {1,4,24}. All are square residues. But 0+1+4=5≠0. Thus a triangle in the shadow need not be one admissible triple. Retaining the forced third-entry label T−x−y and checking it against the actual third cell repairs this precise failure.

## 4. An exact path formulation, but not Hamiltonicity

Build a four-layer directed graph D_T with vertices (k,s), k=0,1,2,3 and s∈R. For each x∈Q(R), place an edge (k,s)→(k+1,s+x), labelled x. Retain only paths from (0,0) to (3,T). These are layered paths even when R is finite; distinct layers avoid any ambiguity from a modular return to the same residue.

**Theorem B.** Ordered triples of square entries summing to T correspond bijectively to length-three paths in D_T via

    (x,y,z) ↔ (0,0) → (1,x) → (2,x+y) → (3,T),

with edge labels x,y,z. Magic squares are exactly EIGHT such paths, one per ordered grid line, whose edge labels agree at every occurrence of the SAME cell, with all nine cell labels distinct.

Proof. Differences of consecutive cumulative sums recover x,y,z. Conversely their partial sums determine the path. Taking one path per line enforces all eight sums; equality at shared cell identities ensures they describe one array, not eight independently chosen triples. Injectivity gives the required distinctness. Both constructions are inverse. No normalization, reduction to root signs, or loss of magnitudes occurs. Over Z, one can restrict the cumulative sums to 0,...,T when T≥0 because the increments are nonnegative; the layers still handle zero increments.

This is a rigorous version of “intersecting paths” if INTERSECTION means compatibility on shared labelled increments. It is not literal intersection of edge sets in G_n, and these are length-three additive-increment paths, not Hamiltonian paths. Two paths can share a cell value while its different prefix sums place the corresponding edges at different graph vertices. Dropping the shared-cell constraints and merely counting paths would be the independent-copy error identified in Session56's correlated-gluing proof. No universal minimality of this full representation is claimed; certified quotients may suffice for narrower observables.

A second exact graph carrier is the bipartite incidence/factor graph with nine cell variables and eight line factors, an edge for each cell-line incidence, domain Q(R), an all-different constraint on cells, and sum T at every line factor. Its graph alone is insufficient; its arithmetic constraints and value labels are essential.

## 5. In what sense could it look dense?

The line-intersection graph has eight line vertices and an edge when two lines meet at a cell. It has 22 of the possible 28 edges, density 11/14: only pairs of different rows and pairs of different columns are disjoint. This is a fixed property of the 3×3 layout and holds for arbitrary, even nonmagic, cell values. It therefore cannot establish any arithmetic existence result. The full bipartite cell-line incidence graph has 17 vertices, 24 edges, and cycle rank 24−17+1=8. These counts are verified separately.

The square-sum graphs G_n, in contrast, are asymptotically sparse. For fixed x, neighbors arise from squares x+y in [x+1,x+n], hence degree at most floor(sqrt(2n)). Thus their edge density is at most sqrt(2n)/(n−1), tending to zero. This does not preclude Hamiltonian paths. It only shows that the word “dense” needs a stated carrier and parameter; the dense fixed line-intersection graph is not the growing square-sum graph.

Nor does arbitrary graph density force Hamiltonicity or magic compatibility. No impossibility theorem for all imaginable gadget reductions is claimed. The exact positive reductions above clarify what must be retained before a transfer can be asserted.

## 6. Finite hostiles and replay

Run `python verify.py --check` (also verified under `python -O`). No third-party dependencies or network are used. The script exhaustively enumerates all oriented Hamiltonian paths in G_n, 15≤n≤25, and divides by reversal. Counts up to reversal are

    n:      15 16 17 18 19 20 21 22 23 24 25
    count:   1  1  1  0  0  0  0  0  3  0 10.

These are bounded finite computations, not an asymptotic Hamiltonicity theorem. The script retains one witness when present, exact graph edge counts, and the source field examples with all eight additive paths.

Cain's F29 root array [9,11,1;6,0,14;12,16,8] produces nine distinct square residues and center zero. Cain's F59 array [20,12,7;2,1,23;22,25,29] has center 1 and total 3; opposite sums are 2, which is NOT a square mod59. Consequently even its compulsory opposite pairs are absent from the square-sum graph on squared-entry values.

The stronger no-wrap hostile

    4  29 15
    27 16  5
    17  3 28

has every entry a square residue mod59 and all eight sums EQUAL AS INTEGERS to 48. But only entries 4 and16 are actual integer squares. Thus removing carries from line sums alone does not lift square membership to Z. The certificate supplies modular roots for every entry. Lo Shu separately verifies that ordinary equal-sum line incidence need not give square-sum edges between adjacent row entries.

## 7. Source pins and scope

External source: O. M. Cain, “Gaussian Integers, Rings, Finite Fields, and the Magic Square of Squares,” arXiv:1908.03236v2 (12 August 2019), https://arxiv.org/html/1908.03236v2 . Precise pins: Section1 defines distinct squared entries and eight equal sums; Theorem2.1 and Corollary2.1 give the center/total and opposite-pair identities for hourglasses; Section5 displays the F29 witness; Corollary7.1 displays F59. The explicit witnesses are independently replayed here. The source's unproved conjectures are not adopted as proved results. No whole article is copied.

Internal dependency pins: Session55 `crt-collision-wiring/THEOREM.md`, Sections1–2 (universal pivots and fixed cell labels); Session56 `correlated-gluing/THEOREM.md`, Sections1–2 (compatible original choice identities, not independent copies); Session56 `overlap-signatures/THEOREM.md`, Sections1–4 (carrier-dependent idempotence). Local read paths and SHA256 digests are in SOURCE_PINS.json. This packet neither changes these sources nor asserts Hamiltonian-path counting is equivalent to the earlier CRT collision census.

Bottom line: the intersection intuition is salvageable EXACTLY as shared-label compatibility of equal-total triple constraints, or as glued additive paths. Ordinary square-sum Hamiltonian paths and non-Parker magic arrays impose different arithmetic predicates and different incidence requirements.
