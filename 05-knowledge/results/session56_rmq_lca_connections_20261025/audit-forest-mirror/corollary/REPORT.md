# Separate corollary audit

Scope: new parity and weighted component-matching corollaries only. The earlier audit and core producer proof scope are unchanged.

## Independent proof check

Let S be the set of realizable labelled outdegree vectors of a finite loopless undirected graph G, with parallel edges allowed if orientations distinguish edges. Reversing every edge sends r to deg(G)-r. It is an involution on S; its only possible fixed vector is deg(G)/2. This fixed vector exists precisely when every vertex degree is even. Necessity is integrality. For sufficiency, the midpoint of each endpoint segment sums to deg(G)/2 in Z_G, so the proved lattice-rounding theorem realizes it whenever it is integral. Equivalently, orient an Euler circuit in each nontrivial Eulerian component. Therefore |S| is odd iff every degree is even. Since |S| is the forest count, the stated parity follows. 'Eulerian' must mean all degrees even componentwise, allowing disconnected graphs and isolated vertices; no connectedness assumption is needed. This is unrelated to even automorphism permutations or Royle's even-graph nomenclature.

For G_n=K_n-P_n with n>=3, endpoints have degree n-2 and internal vertices degree n-3. If n is odd, endpoints are odd-degree. If n is even (thus at least four), internal vertices are odd-degree. Hence no score is fixed by complementation.

Assuming the audited cycle-move connectivity theorem, contract each cycle-connected component to its score and retain the edges of the free-edge-complement involution. The resulting underlying simple graph is exactly a perfect matching on |S| vertices, with |S|/2 pairs. This wording matters: identifying complement-paired vertices is a further different quotient, not the matching itself. For the fixed-path model, native complement reverses free edges only, preserving the lock; it is not literal reversal of all tournament edges.

The coefficient w_r of product_(uv in E)(x_u+x_v) is the number of orientations of G with outdegree r. Complement gives a bijection of its fibre to that at deg-r. If these fibres differ, there are w_r undirected complement edges between them: counting from one chosen fibre counts each edge once. Thus sum over unordered score pairs of w_r=2^(|E|-1), when no fixed score exists. Dividing by two again within each pair would be wrong. These are labelled orientation counts, not isomorphism/converse orbit weights; no original merged tournament metagraph is identified.

For an all-even graph the sole fixed score has a self-paired fibre. If |E|>0, free-edge complement has no fixed orientation, so this fibre contributes w_center/2 actual undirected internal complement edges, while other score pairs contribute w_r. For |E|=0 there is one fixed orientation, requiring a separate identity/loop convention; the expression 2^(|E|-1) is not a count. The path n>=3 statement automatically avoids this exception. Loops are excluded in this corollary's orientation/edge-involution model; any extension must explicitly set loop orientation and degree conventions rather than importing a loopless orientation count.

## Independent finite checks

`python check.py` and `python -O check.py` enumerate all 1,099 labelled simple graphs with 1..5 vertices and all their orientations. Forest counts are independently enumerated by merging vertex sets. They verify parity, the exact fixed-score count, complementary fibre multiplicities, and matching weight sums. Additional path-complement checks run through n=7:

n=3: 2 fibres, 1 pair, edge mass 1.
n=4: 8 fibres, 4 pairs, edge mass 4.
n=5: 52 fibres, 26 pairs, edge mass 32.
n=6: 478 fibres, 239 pairs, edge mass 512.
n=7: 5706 fibres, 2853 pairs, edge mass 16384.

Acceptance checks use explicit exceptions; an intentionally failing sentinel is caught to verify the guard remains active under optimization. Outputs match exactly. These are finite controls, not the source of the general proof above.

Producer-file review pending at report creation; see ADDENDUM.md for final file pins and disposition.
