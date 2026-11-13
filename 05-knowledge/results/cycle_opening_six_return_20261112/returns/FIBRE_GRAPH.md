# Exact neutral fibers: Johnson graphs and a distributive refinement

Status: **PROVED** for the explicitly specified macro family and rewrite rules. **VERIFIED** exact controls through seven macros. Canonical peeling proves the ordinary isomorphism-class count for a fixed strong source. This is a subfamily classification, not a classification of the full tournament metagraph.

## The objects and observation map

Use the physical eight-vertex macros A=TTTLL and B=TTLTL from [NEUTRAL_LANGUAGE](NEUTRAL_LANGUAGE.md). Fix any strong source T on n labeled vertices, with a marked indegree-one root r. Set q0=2x+2b+d for its exact interface state. For a chronological word w of m>=1 macros, the total count is

    H(w)=51*q0*93^a*97^b,

where a,b count A,B among the first m-1 positions. The last macro E in {A,B} selects the final complete L/T-state q0*93^a*97^b*vE. Each macro adds eight new names.

The word-to-labeled-tournament map is injective. In zero-based macro block i, inspect the internal arc between vertices n+8i+4 and n+8i+5. For A the later vertex beats the earlier; for B the earlier beats the later. Future births never alter this arc. Thus these m private arcs recover every macro choice. More strongly, changing one macro letter toggles exactly the four internal arcs listed in NEUTRAL_LANGUAGE, all within that block. Two words at Hamming distance t therefore differ on exactly 4t tournament arcs. This first argument concerns fully specified labeled graphs, not just formal construction words. The next argument proves that the addresses also survive forgetting all labels.

The two primes dividing 93 are 3 and 31, and 97 is coprime to 93. Consequently different exponent pairs (a,b) give different H for the fixed source q0. At fixed depth m and fixed a+b=m-1, the H fiber therefore consists of exactly

    2*binomial(a+b,a)

canonical labeled representatives. Its two distinct complete L/T-states each have binomial(a+b,a) such representatives. The following peeling theorem shows that all these labeled representatives belong to distinct ordinary isomorphism classes.

## Canonical peeling and ordinary enumeration

This is a local specialization of the earlier general [adaptive ancestry theorem, A1-A3](../../adaptive_metagraph_reformulation_20261103/ANCESTRY.md), including its automorphism formula. The new scope here is the explicit macro word fibers and their count observations, not a new general unique-history theorem. The broader L/K and Fibonacci restrictions are recorded with this provenance in [BINARY_LANGUAGE](BINARY_LANGUAGE.md).

Suppose a tournament has every indegree at least one and p has indegree one. Adding a newborn z beaten only by p leaves exactly two indegree-one vertices: p and z. Every other old vertex gains an incoming arc from z and now has indegree at least two. The sole predecessor of z is p, still of indegree one; the sole predecessor of p is another old vertex and now has indegree at least two. Thus z is intrinsically identifiable as the unique indegree-one vertex whose sole predecessor also has indegree one.

Every birth used by L or T has such a port. T's first birth leaves the old port p at indegree one, so the same argument applies to its second birth. Strongness is preserved at each birth. Therefore an unlabeled output can be canonically peeled, newest vertex first, until the fixed original order n is reached. Every graph isomorphism must preserve this deletion order. In particular it preserves the private arcs encoding each A/B choice, so two different macro words cannot be isomorphic. Stopping at n is essential; the source might itself admit further peeling, but no further deletion is needed.

Let z1 be the earliest chronological newborn, which is the last vertex peeled. Its sole predecessor identifies the original marked root r. More generally, for a fixed source order n and positive depth m, the map from rooted strong source isomorphism classes (T,r), with indegree(r)=1, and macro words to ordinary strong output classes is injective. An isomorphism of outputs would canonically recover both the rooted source class and the word. This statement avoids treating two different root orbits of the same base tournament as one source.

Consequently, for a single fixed rooted strong source, the exact number of ordinary isomorphism classes produced at order n+8m is 2^m. Their exact Hamiltonian-count distribution is

    count value: 51*q0*93^a*97^(m-1-a),
    number of ordinary classes: 2*binomial(m-1,a),  0<=a<=m-1.

There are 2m complete L/T-states at positive depth m; the exact ordinary family of size 2^m compresses to them. The ordinary generating function for this selected family, including its single original source at m=0, is x^n/(1-2x^8). These are certified subfamily counts contributing to A000568, not new evaluations of the whole sequence.

For every positive depth, the output automorphism group is isomorphic to Aut(T,r). Indeed canonical peeling forces every new vertex to be fixed, and recovering the first newborn's predecessor forces r to be fixed. Conversely every rooted automorphism of T extends by fixing every new vertex, because each birth's old incidences depend only on its named port. Thus each output has exactly (n+8m)!/|Aut(T,r)| different labelings. For the directed triangle rooted at one vertex, the stabilizer is trivial, so all 2^m outputs are asymmetric. At order N=3+8m, m>=1, this gives exactly 2^m*N! distinct labeled tournaments in the selected family. Here q0=5, so each H-fiber has 2*binomial(m-1,a) ordinary classes and that number times N! labeled members, at H=255*93^a*97^(m-1-a).

The minimum-indegree premise has a minimal hostile witness: start from the transitive pair with its sink marked, then apply L. The output is the directed triangle, whose three vertices all have indegree one and an indegree-one predecessor. The newest vertex cannot be selected by the canonical rule. The strong-source premise is explicit. The earlier private-arc argument remains valid for any indegree-one source as a labeled assertion; the ordinary-isomorphism promotion here uses minimum indegree at least one at every birth. A135244 is an aliquot inverse sequence and is not being counted by this theorem.

## All neutral moves changing at most two macro choices

Declare a rewrite edge when two words in this fiber differ in at most two macro positions. This is an intrinsic binary relation on the defined labeled objects: the pairwise observable is macro Hamming distance, and the preserved target is the whole Hamiltonian count. No orientation or forced tournament structure is imposed on this graph.

A one-position change is neutral exactly when it changes the final macro. It toggles four graph arcs and exchanges the two interface types while keeping H.

A two-position change is neutral exactly when it exchanges an A and a B among the first m-1 positions. It toggles eight graph arcs, preserves both H and the complete L/T-state, and retains the entire last macro unchanged. Any two-prefix-position exchange is valid; the positions need not be adjacent. There are no other neutral one- or two-position changes: changing one prefix position changes the exponent pair, and also changing the final position cannot repair that numerical change.

Represent a vertex by (S,E), where S is the a-element subset of {1,...,a+b} occupied by A, and E is the final letter. The rewrite relation becomes:

- replace one element of S by one element outside S, holding E fixed; or
- hold S fixed and exchange E.

This is exactly the Cartesian product J(a+b,a) square K2, where J(N,a) is defined as the graph of a-element subsets of an N-set with edges for a single replacement. Hence the selected neutral metagraph family is completely classified:

| quantity | exact value |
|---|---|
| ordinary classes / canonical representatives | 2*binomial(a+b,a) |
| degree at every vertex | ab+1 |
| undirected edges | binomial(a+b,a)*(ab+1) |
| distance between (S,E),(S',E') | size(S minus S') + indicator(E differs from E') |
| diameter | min(a,b)+1 |

From any fixed vertex, the exact distance generating polynomial is (1+z)*sum_j binomial(a,j)*binomial(b,j)*z^j, with 0<=j<=min(a,b). Choose j positions to remove and j to insert, then independently decide whether to flip the final type.

The distance formula follows because each replacement corrects at most one missing chosen position and can do so at every step. The independent final flip contributes exactly one when needed. This proves the diameter and connectivity without a finite extrapolation.

An eight-arc exchange retains the literal final root and its sole predecessor, since both belong to the unchanged final macro. The four-arc final switch also retains the final root and its sole predecessor by the original A/B construction. Thus every edge in this graph preserves that physical root interface, even though only the eight-arc exchanges preserve the full three-coordinate continuation state.

## Adjacent exchanges and distributive structure

Restrict the eight-arc exchanges to adjacent prefix positions. Full count and continuation state are still preserved: locally one uses ABA=BAA or ABB=BAB, with the third macro supplying the needed state collapse. The prefix restriction is essential; swapping a prefix macro with the final macro changes the count in general.

Encode the positions of the a letters A by increasing integers s1,...,sa. The sequence si-i lies between 0 and b and describes a Ferrers diagram in an a-by-b rectangle. An adjacent AB/BA exchange adds or removes one cell. Thus this restricted prefix graph is the cover graph of the distributive lattice of order ideals in the rectangle; the full neutral graph is its Cartesian product with K2 for the final macro.

Its local degree is one plus the number of adjacent unequal prefix pairs. Its diameter is ab+1: moving A^a B^b to B^b A^a requires exactly ab adjacent exchanges, and an opposite final letter requires one further flip. For two general prefix words, the adjacent-exchange distance is sum_i absolute(si-ti), plus the independent final-bit difference. From an extreme prefix word A^a B^b, its distance generating polynomial is (1+z) times the Gaussian binomial [a+b choose a]_z, equivalently the area-generating polynomial of diagrams in the a-by-b rectangle. This follows directly from the adjacent-exchange/one-cell bijection.

This subgraph gives a locality-sensitive refinement of the Johnson fiber. The two diameters must not be conflated: unrestricted exchanges have diameter min(a,b)+1, adjacent exchanges have diameter ab+1.

## A precise five/six connection

With six macros there are five prefix choices and one final interface choice. In the a=2,b=3 fiber, S ranges over the ten two-element subsets of a five-element set. The prefix graph is J(5,2), equivalently the line graph of the complete graph on five vertices: two chosen pairs are adjacent exactly when they share one element. Its complementary relation joins disjoint pairs, giving the usual ten-vertex Petersen presentation. Adding the final K2 produces a 20-vertex, 7-regular neutral fiber of diameter 3; the adjacent-exchange refinement has diameter 7.

The map is explicit: a pair of selected prefix positions maps to an edge of K5. It preserves pair intersection, replacement adjacency, the exact total H, and the marked root interface when translated back to macro rewrites. This subset projection alone discards the detailed old tournament; the fixed rooted source and final interface letter are required sidecars to reconstruct a graph. The full output retains the rooted source intrinsically through canonical peeling. Five here arises from m-1 with m=6. The same construction exists for every m, so this does not identify the five Chan-Pak core tournaments with the vertices of K5 or assert a privileged universal role for five.

## Size growth and scope boundary

At fixed source and depth m there are m H fibers, 2m continuation states, and 2^m pairwise nonisomorphic macro outputs. Increasing depth appends one macro, adds eight vertices, and multiplies the scale by 93 or 97 according to the previous final type. The new final type is freely chosen. This composes exactly with the fiber graphs above.

These are fibers of a defined construction language and its four-/eight-arc rewrites. Additional H-neutral graph changes might connect them to other outputs; that possibility is not excluded. Isomorphisms among the fixed-source addresses have been excluded by canonical peeling. The exact A000568 contribution concerns this selected subfamily only. The theorem supplies a fully classified infinite subfamily and a rigorous test bed for broader metagraph enumeration.

Reproduce:

    python -B 05-knowledge/results/cycle_opening_six_return_20261112/returns/fibre_graph.py

`fibre_graph.json` checks every word through depth seven: exact count matrices, private-arc recovery, canonical peeling before and after a nontrivial relabeling, injectivity, fiber sizes, neutral move classification, graph degrees and diameters. Direct labeled arc-difference counts are additionally checked for all pairs within each H-fiber through depth five. The all-order proof is the subset representation and private-arc argument above.
