# Session 57: square-sum graphs forget magic incidence

## Definitions and status

All field statements use a finite FIELD of odd order, not a residue ring. Let Q={x²:x∈F}, including 0. The simple square-sum graph S⁺(A) on a set A⊆F has an edge xy for x≠y exactly when x+y∈Q. Its nonzero-sum variant S*(A) additionally requires x+y≠0. Vertex values below are squared entries, not roots. A magic grid has nine distinct entries in Q and all eight conventional row, column and diagonal sums equal. Hamiltonian paths are counted modulo reversal, not modulo graph automorphisms.

**PROVED:** magic existence and nine-square square-sum Hamiltonian-path existence do not imply one another. More strongly, in characteristic three the induced graph of any magic entry set is complete or empty according to the quadratic character of −1 (for S⁺). An exact three-uniform hypergraph retains the missing eight-line incidence. **FINITE-EXACT:** all 202 normalized valid grids over the 30 odd prime-power fields q≤101, with both graph conventions, exact HP counts and adjacency characteristic polynomials. No inference about existence of an integer magic square is made.

## 1. The character mismatch

Write a magic grid as

    c+u      c-u-v     c+v
    c-u+v    c         c+u-v
    c-v      c+u+v     c-u.

The total is T=3c and each opposite pair sums to 2c. For a nonzero square c and characteristic not 2, these four pairs are edges of S⁺ precisely when χ(2)=1. If characteristic is not 3, T is square precisely when χ(3)=1. In characteristic 3, T=0: it qualifies only for the zero-allowed convention. At c=0 all opposite-pair sums are 0 and therefore are edges of S⁺ but not S*. Neither a magic line nor its opposite pairs inherently gives ordinary square-sum adjacency. Over positive integer squares, 2e² and 3e² are nonsquares for e≠0, by prime valuations.

## 2. Both implications fail, with exact witnesses

### Hamiltonian path does not imply magic

F17 has nine square values Q={0,1,2,4,8,9,13,15,16} and no valid magic grid. The absence was independently exhaustively checked in both normalized center types. Its square-sum graph nevertheless has 1,152 unoriented HPs when zero sums qualify, and 112 when they do not. A path that works even in S* is

    0, 1, 8, 13, 2, 16, 9, 4, 15.

Its consecutive sums modulo17 are 1,9,4,15,1,8,13,2, all nonzero squares. S⁺ has 20 edges and degree multiset 4^8,8; S* has16 edges and degrees 3^8,8.

### Magic does not imply Hamiltonian path

In F59 the grid

    15 35 12
    57  1  4
    49 26 46

has nine distinct square entries and every line sums to3. Vertex57 is isolated in either square-sum graph: its sums with the other eight entries are 58,2,10,13,24,33,44,47, all nonsquares modulo59. There are14 edges and degrees 0,2,2,3,3,4,4,5,5. Thus neither graph has a Hamiltonian path. All eight normalized valid F59 grids give the same spectral/HP signature. This is a finite-field counterexample to a proposed reduction, not a claim about the unresolved integer problem.

For contrast, the F29 witness [1,23,5;4,0,25;24,6,28] has28 edges and17,688 HPs in S⁺, versus24 edges and3,368 HPs in S*. These counts alone cannot recover its line incidence.

## 3. Characteristic-three collapse theorem

**Theorem.** Let A be the entry set of a valid magic grid over F_(3^n). Then A=c+W, where W is a two-dimensional F3-subspace. For S⁺(A):

- If n is even, S⁺(A)=K9.
- If n is odd, S⁺(A) is the empty nine-vertex graph, and 0∉A.

**Proof.** In characteristic3, the nine offsets su+tv for s,t∈{−1,0,1} are exactly W=span_F3(u,v); distinctness forces dimension2. For x=c+w and y=c+w', put z=−x−y=c−w−w'∈A. Thus x+y=−z. If n is even, −1 is square, so every pair sum is square, including zero. If n is odd, −1 is nonsquare. First 0 cannot belong to A: otherwise A=W and any nonzero w gives both w and −w as square entries, a contradiction. Hence z≠0 and −z is nonsquare, proving there are no edges. ∎

For S*, the odd-n conclusion is unchanged. When n is even, delete exactly the edges {x,−x} whose endpoints occur in A. If 0∈A then A=W, so exactly four disjoint edges are deleted; if 0∉A then A and −A are disjoint cosets, so none are deleted.

An explicit empty-graph witness exists in F2187, represented in the pinned polynomial field as [37,78,5;77,1,42;6,41,73]. Full field modulus and base3 encoding are in witness-2187.json; interpreting these integers modulo2187 would be WRONG. An explicit F729 witness is also supplied. Every normalized F81 magic grid has graph K9 under S⁺, with181,440 unoriented HPs. For the witnessed zero-containing F81 plane, S* is K9 minus four disjoint edges and has70,080 HPs.

**Lost arrangement count.** Any fixed affine plane A of nine square entries has exactly9·8·6=432 labeled magic arrangements: choose c∈A, choose nonzero u∈W, then choose v∈W outside span(u). The universal parametrization proves both sufficiency and necessity. There are9!=362,880 total arrangements, so only1/840 are magic, although every arrangement yields the same unlabeled square-sum graph. Direct enumeration verifies432 in the supplied F81 and F2187 planes. This is stronger than an unexplained correlation: it identifies precisely what the graph forgets.

## 4. A positive incidence-preserving replacement

For each T∈F define the three-uniform hypergraph

    H_T: vertices Q; edges {x,y,z} of three distinct values with x+y+z=T.

Let M be the fixed nine-vertex hypergraph whose eight edges are the three rows, three columns and two diagonals of the3×3 position grid. A valid square-entry magic grid with totalT is EXACTLY an injective edge-preserving map M→H_T. This is a genuine equivalence: mapping an edge means its sum isT, and injectivity supplies distinctness. It is a non-induced embedding, because additional T-summing triples are allowed.

Every H_T is linear: a pair x,y determines at most one third value z=T−x−y. The link at a vertex z is the matching relation x+y=T−z, with invalid or repeated vertices removed. At a magic center c, T=3c and its four prescribed opposite pairs lie in the exact link x+y=2c. This matching records the equality that the square-sum graph replaced with the weaker membership question x+y∈Q.

In characteristic3, H_0 restricted to a valid affine-plane entry set is precisely the12 lines of AG(2,3): distinct x,y determine z=−x−y, and these are exactly affine three-point lines. The eight required grid lines form a distinguished subconfiguration, rather than a Hamiltonian path. Hence affine incidence/hypergraph embedding is a faithful positive framework even when the square-sum graph is K9 or empty.

There is a limited Paley-style observation: the relation χ(x+y)=1 can be written χ(x−(−y))=1 by reflecting ONLY the second coordinate. This gives a relation between two labeled copies of the value set, not an isomorphism of ordinary induced graphs by one simultaneous vertex relabeling. It supplies no magic-to-HP bijection.

## 5. Exact census and verification boundaries

results.json groups grids by the PAIR of exact edge-count, HP-count and adjacency-characteristic-polynomial signatures for S⁺ and S*. Equal signatures are not asserted to imply graph isomorphism. Every group includes an explicit valid matrix, degree list and HP witness when one exists. All normalized grids are enumerated, not merely preselected examples. Nonzero center is normalized to1 by square scaling; zero center has top-left entry normalized to1. Scaling by any nonzero square induces a graph isomorphism for both conventions. Thus this covers all graph signatures for these field orders, while the counts themselves remain normalized labeled counts.

The only orders q≤101 with a valid magic grid but a zero-HP induced signature are59 and83. This bounded statement is not extrapolated. No all-subset graph census is claimed for Parker fields: their reported graphs use only the first nine encoded square values.

analyze.py uses the pinned prior finite-field implementation and exact SymPy characteristic polynomials. check.py independently uses tuple-polynomial field arithmetic, exhaustive normalized-grid enumeration, forward subset DP for HPs, and Newton identities from matrix-power traces for characteristic polynomials. It also brute-forces all9! vertex orders for both F17 conventions and all9! arrangements on the F81/F2187 planes. check-result.json records the successful replay. No floating-point spectrum or probabilistic classification is used. Sources are copied locally and hashed; no live repository state is needed.

## Primary-source context

Cain, *Gaussian Integers, Rings, Finite Fields, and the Magic Square of Squares*, arXiv1908.03236v2, https://arxiv.org/html/1908.03236v2 (actually inspected2026-10-25). Theorem2.1 and Corollary2.1 give the central-total and opposite-pair identities. Sections5–6 analyze Parker fields and square-scaling normalization; Corollary5.4 identifies F17 as Parker, and Corollary7.1 supplies an F59 magic construction. These are inherited context, not new claims. The graph counterexamples, exact path/spectrum census, collapse proof and incidence formulation above are developed and independently checked here. No priority claim is made. No paper text is reproduced.
