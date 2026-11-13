# The exact score-walk recurrence needs memory from order five onward

**Status:** PROVED IN NOTES after independent root audit; see [AUDIT.md](../AUDIT.md) (the weighted memory identity, its sharp iff boundary, and failure of folded-score Markov closure at every order n≥5); FINITE-EXACT (full native matrices at n=3..6, the order-five cubic memory recurrence, and the recorded replay). The general Schur-complement identity is elementary linear algebra, not a novelty claim. Its application here is deliberately to the actual complement-line multigraph. A recurrence generating the full class wiring without class representatives or path information is still OPEN. The Hamiltonian-path value spectrum and Chan–Pak remain OPEN.

## Inheritance and the next interface

The closest proved mechanism is [THM-796, three-sorted recursive incidence](../../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md), especially its literal half-edge matrix and reversible node shadow, Sections 1 and 5. [Session 314's alignment recurrence](../../session314_metagraph_constraint_algebra_20261113/ALIGNMENT_RECURRENCE.md) already evaluates the entire colored matrix from class representatives. Its [all-color score theorem](../../session314_metagraph_constraint_algebra_20261113/ALL_SCORE_SUPPORT.md) already classifies score-quotient support. Those results are inherited, not claimed again.

The canonical hostile to a score-and-degree-only alignment test is the directed path plus triangle in TT2→C3; the determinant connectivity guard in the alignment recurrence is essential. The corrected near miss is the order-four blue loop whose two endpoint-marked outputs split under padding ([marked stabilization](../../session314_metagraph_constraint_algebra_20261113/HISTORICAL_BRIDGES.md)). The least-used relevant sidecar here is the **source mass of each actual dart fibre**. A score block is not a uniformly weighted set of classes.

Portfolio/concept board: anchor—actual repeated complement moves; niche—source-weighted score compression; wildcard—the local triangle twist of the transitive node's unique neighbor. The five live concepts are class identity, score support, path-port multiplicity, reversible conditional expectation, and hidden walk memory. The new triangle twist leaves score support unchanged, separates class identity by local triangle counts, and makes the lost multiplicity visible in an observable two-step return.

## 1. Exact objects, weights, and gauge

Fix n≥3. Nodes are tournament isomorphism classes modulo converse. A dart is an actual fixed-path tiling; opposite darts are related by complementing **all nonpath arcs**, retaining the displayed directed Hamiltonian path. Blue/black is fixed by the path-reflection criterion. Parallel lines are retained, and a loop supplies two darts. Write A=A_B+A_K for the symmetric integer dart adjacency matrix.

For a merged node u represented by T,

    d_u = sum_v A(u,v) = |R(T)| H(T)/|Aut(T)|,

where |R(T)| is 1 for a self-converse ordinary class and 2 otherwise. This is THM-796's fibre formula. Define P=D^(-1)A with D=diag(d). One step chooses a **uniform literal dart from the current merged node** and traverses its line. It is not uniform over neighboring classes. Uniform sampling of physical Hamiltonian paths in any one ordinary representative does induce this same merged-node transition law: the free automorphism-orbit factors and converse duplication cancel after normalization. Those factors remain essential for raw adjacency counts and literal dart identities. P is self-adjoint for <f,g>_D=f^T Dg.

Let s(u) be the sorted score sequence, folded with its degree complement by lexicographic minimum. All following algebra works for any partition, but the sharp failure theorem uses this folded-score partition. Let L lift a function on score blocks to a function constant on each block, let

    m_a = sum_{u:s(u)=a} d_u,
    E(a,u) = d_u/m_a if s(u)=a, and 0 otherwise,
    Q = E P L,
    R = I - L E.

Thus EL=I. E is conditional averaging in the reversible measure, R is the orthogonal projection onto hidden functions ker(E), and

    Q(a,b) = (1/m_a) sum_{u in a, v in b} A(u,v).

Session 314's score product computes the numerator, and its support theorem determines where it is positive. Neither determines the omitted conditional distribution within a score block after a move. The native convention is used throughout; replacing it by path reversal after converse merging gives the same A, but the marked path and label gauge are not silently identified upstairs.

## 2. The exact recurrence and the sharp memory boundary

Define matrix coefficients and their formal generating function by

    F_k = E P^k L,               F(z) = sum_{k≥0} z^k F_k.

F_k(a,b) is the actual k-step probability from a source microstate distributed as d_u/m_a within a, ending anywhere in b. In particular F_0=I and F_1=Q. Put B=RPL and C=(RPR)|ker(E). Define

    K_j = E P R (RPR)^j R P L,
    K(z) = sum_{j≥0} z^j K_j.

Then, as a formal power series over the rationals,

    F(z) = (I - z Q - z² K(z))^(-1),                         (1)
    K(z) = E P R (I-z RPR)^(-1) R P L,

or, coefficientwise,

    F_k = Q F_(k-1) + sum_{j=0}^{k-2} K_j F_(k-2-j), k≥2.  (2)

**Proof.** Split the node-function space orthogonally as im(L)⊕ker(E). The P blocks are Q, the hidden-to-visible map EP, B, and C. Eliminate the hidden block in I−zP. This gives (1); multiplying on the left by its denominator gives (2). All inverses have constant term I, so no analytic convergence assumption is needed. This proves every coefficient, not only a finite walk census.

The first missing coordinate has a particularly sharp diagnostic. With M=diag(m),

    F_2 - Q² = K_0,
    M K_0 = B^T D B.                                       (3)

The second identity follows because E=M^(-1)L^T D and P,R are D-self-adjoint. Consequently the following are equivalent:

1. F_2=Q²;
2. B=0;
3. for any two nodes in the same score block, their P-probabilities of entering each target score block agree;
4. F_k=Q^k for every k, or equivalently K(z)=0.

Indeed (3) makes K_0=0 iff B=0; B=0 is exactly strong lumpability; it kills every K_j. Conversely item 4 includes item 1. Thus a nonzero defect cannot first hide until a later walk length: it is already detected by a positive diagonal entry of the weighted Gram matrix at length two. Individual off-diagonal entries of K_0 or later memory matrices need not be nonnegative.

**Finite recurrence.** If there are N merged nodes and b blocks, C has dimension N−b. Every entry of K(z) has denominator dividing det(I−zC), of degree at most N−b, and hence its coefficients satisfy the corresponding finite linear recurrence. This is an exact algebraic operation on the full incidence matrix; it does not manufacture that matrix from score data alone.

**Colors.** Replace P by P(x,y)=D^(-1)(x A_B+y A_K), keeping the same source weights D and E. Equations (1)–(2) hold over Q[x,y][[z]] and record the blue/black inventory of weighted walks. Real specializations remain D-self-adjoint and obey the Gram identity. Exact closure for all colors requires the hidden maps of both color operators to vanish; a statement only at x=y=1 must not be promoted to colored closure.

## 3. Folded-score closure fails at every n≥5

Let T=TT_n. It has one Hamiltonian path, is rigid and self-converse, and therefore has exactly one dart. Its unique neighbor is the staircase R_n with arcs

    i→i+1, and j→i whenever j≥i+2,

on 0,...,n−1. Its score multiset is

    (1,1,2,3,...,n−3,n−2,n−2).

For n≥4 it is rigid: every score bucket has size at most two, and swapping a two-vertex bucket reverses its internal arc. Reversal of the vertex order is an anti-isomorphism, so R_n is self-converse. Its degree is therefore H(R_n), the inherited tribonacci value from [THM-4731, self-lines/twins/tribonacci partner](../../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md): h_0=h_1=h_2=1 and h_n=h_(n−1)+h_(n−2)+h_(n−3). This cited recurrence is not needed for strictness below.

For every n≥5, reverse the directed triangle on vertices {1,2,3} of R_n, obtaining U_n. Triangle reversal preserves every vertex score. But R_n and U_n are not isomorphic, even after converse.

To prove this, the cyclic triples of R_n are exactly the consecutive triples (i,i+1,i+2): for a<b<c the long arc is c→a, so cyclicity requires a→b and b→c, forcing consecutive indices. The triangle flip retains the cyclic triple 123 in reversed orientation, removes 012 and 234, and adds 013 and 134; there are no other changes. Vertices outside {0,...,4} have uniform directions to {1,2,3}, and cannot add an affected cyclic triple.

At n=5 the vector of cyclic-triangle incidences is

    (1,2,3,2,1) → (1,3,1,3,1).

At every n≥6, the old maximum incidence is 3 and vertex 3 acquires incidence 4. Sorted triangle-incidence vectors are invariant under isomorphism and converse, so the two merged nodes are distinct.

Let a be their common folded-score block and m_a its dart mass. Then

    m_a > d_(R_n) > 0.

The transitive profile is a singleton class. Its sole line reaches R_n, so

    F_2(TT_n,TT_n) = 1/d_(R_n),
    (Q²)(TT_n,TT_n) = 1/m_a,
    K_0(TT_n,TT_n) = 1/d_(R_n) - 1/m_a > 0.               (4)

The second formula holds because Q(TT_n,a)=1 and the entire a-to-TT_n dart count is one. This proves non-lumpability for every n≥5. At n=3,4 every folded-score block is a single merged node (independently checked by the exhaustive replay), so closure holds there. The boundary n=5 is sharp on the nondegenerate complement-line tower n≥3.

**Failure anatomy.** The invalid implication is “one-step support and multiplicity are exact, therefore their quotient iterates are exact.” The missing coordinate is the distribution inside the intermediate score block, with its literal source-dart weights. The strongest survivor is the exact one-step Q. The repair is the hidden recurrence (1), or any refined state partition satisfying item 3 above. This does not refute the all-color support classification.

## 4. Explicit finite consequence and cubic memory

The complete universe is all 2^((n−1)(n−2)/2) fixed-path tilings at n=3..6, with no filters. Fixed-path complement enumeration and a separate enumeration of physical Hamiltonian paths in ordinary representatives give identical weighted native matrices, including loop doubling and converse factors. Score-sorted canonicalization checks permutations inside every score bucket; it is not an assumed classification by scores.

At n=5 there are 10 merged nodes, 7 folded-score blocks and 64 darts. R_5 has degree 9, while its score block (1,1,2,3,3) has mass 18. Thus the exact TT5 return is **1/9**, while the memoryless score quotient gives **1/18**. At n=6 the corresponding values are 1/17 and 1/64.

An equal-weight order-five hostile is stronger: upper-arc masks **1015** and **1013** both have profile (1,1,2,3,3), H=9, |Aut|=1, self-converse status, and dart degree 9. Their probabilities of moving to TT5 are 0 and 1/9. Thus adding H, automorphism order and converse status to the scores does not repair the quotient.

A second transparent order-five witness uses upper-arc masks **983** and **981**, little-endian in lexicographic pair order (01,02,03,04,12,13,14,23,24,34), bit 1 meaning increasing direction. Both have score profile (1,2,2,2,3). Their (H,|Aut|,d) values are (15,3,5) and (13,1,13). Their probabilities of moving to profile (0,2,2,2,4) are respectively 0 and 1/13. Full literal adjacency rows and all code definitions are retained in the receipts.

The order-five hidden space has dimension 3. In the explicit weighted-zero-sum basis constructed by the checker,

    det(I−zC) = (186615 − 6590 z − 16714 z² + 2126 z³)/186615.

Every matrix coefficient of the memory sequence therefore satisfies

    186615 K_(j+3) = 6590 K_(j+2) + 16714 K_(j+1) − 2126 K_j,

and its TT5 diagonal generating function is exactly

    K(z)[TT5,TT5]
      = 5(60 z² − 1068 z + 4147)
        / (2(2126 z³ − 16714 z² − 6590 z + 186615)).

The constant term is 1/18, as required by the actual two-step discrepancy. This is a compact recurrence for the correction to the complete score quotient, not a replacement classification of all metagraph nodes.

## 4.1 A minimal linear completion is already the full class dimension at n=5,6

There is a constructive completion smaller than retaining arbitrary hidden coordinates whenever some of them are invisible to all score observations. Define the increasing rational Krylov spaces

    V_j = span{columns of L, PL, ..., P^j L}.

Stop at the first V_j=V_(j+1). The common space V is P-invariant; self-adjointness makes its orthogonal complement invariant too. Restricting P to V and retaining the score lift therefore reproduces every F_k. This is a known linear-algebra completion, not a claimed new general control theorem.

Its dimension is minimal among exact linear **moment representations** of the form F_k=E' (P')^k L' for all k≥0. To see this, the block matrix with (i,j) block M F_(i+j) is the D-Gram matrix of the columns of L,PL,...,P^jL, and its stabilized rank is dim(V). Any representation of state dimension r factors the same block Hankel matrix through r dimensions, hence r≥dim(V). This statement includes F_0 and assumes no separate direct-feedthrough term; it is not an unqualified claim about rational transfer-function McMillan degree.

Exact fraction elimination gives

| n | successive dimensions dim(V_0),dim(V_1),... | merged nodes |
|---:|---|---:|
| 3 | 2, 2 | 2 |
| 4 | 3, 3 | 3 |
| 5 | 7, 10, 10 | 10 |
| 6 | 14, 26, 34, 34 | 34 |

Thus at n=5,6, preserving all score-observed walk moments linearly already needs all 10 or 34 class dimensions, respectively. The scores recover the missing coordinates after one or two applications of P. This is a finite exact obstruction to further linear compression in this model, not a general theorem that score-walk observability separates every merged class at every n.

## 5. Connection contract and verification

Source: the actual colored complement-line multigraph with dart source weights. Target: score-observed walk series F(z). Map: conditional expectation E and lift L. Preserved predicate: every weighted, optionally color-counted walk coefficient. Destroyed by Q alone: the conditional microstate distribution after an intermediate visit. Restoring sidecar: B and the hidden operator C, equivalently K(z). Cheapest decisive test: TT5's two-step return, equation (4).

[replay.py](replay.py) uses only the Python standard library and always-active `require` gates, so `python -O` retains every check. Positive controls are the singleton score partitions at n=3,4 and independently constructed native matrices. Hostile controls are the explicit n=5 same-score pair and the uniform triangle-twist construction through n=32. The all-order claims rest on the proofs above, not the finite ceiling.

Reproduce from repository root:

```bash
python3 05-knowledge/results/joint_metagraph_dyadic_20261113/metagraph/replay.py
python3 -O 05-knowledge/results/joint_metagraph_dyadic_20261113/metagraph/replay.py
```

[normal.json](normal.json) and [optimized.json](optimized.json) are byte-identical. They contain source SHA-256, native matrix hashes, source block masses, explicit literal rows, the hidden matrix, exact Krylov ranks, and exact walk comparisons through length eight. Their input repository pin is f90cd879d156a37f0027d1719a25b8731918db6a; dependencies named above are read at that pin. This lane makes no changes to canon or scarce theorem namespaces.
