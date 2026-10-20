# Raw Hamilton counts are a Gram diagonal for two-complement transport

**PROVED, all finite orders.** This addendum is separate from the frozen endpoint/icosahedron packet. The matrix factorization is the standard common-neighbour Gram mechanism; the content is identifying the correct labelled score-fibre incidence, its unique path witnesses, and its raw-H diagonal. No novelty claim for Gram positivity itself is intended.

## 1. Exact incidence theorem

Fix a labelled outdegree vector d and let F_d be all tournaments T with that exact vector. Let S_d be the set of tournaments that arise as K_P(T), for T in F_d and a directed Hamiltonian path P in T. Define the binary incidence matrix

    B[T,S] = 1 if S=K_P(T) for some P; 0 otherwise.

The witness P, when it exists, is unique: the unordered edges on which T and S agree are exactly the n-1 edges of P, and their directions in T determine the unique ordered Hamiltonian path. Thus different Hamiltonian paths of T give different labelled intermediate tournaments S, and

    sum_S B[T,S] = H(T).                                 (1)

If S=K_P(T), its score vector satisfies

    d(S)=(n-1)1-d+delta(P).

Therefore the ordered endpoints of P are recoverable from d(S) and the fixed fibre vector d. For n>1 the positive and negative entries of delta(P) identify its start and end. The order-one case is immediate separately.

Let W[T,U] count path-pair presentations T -> S -> U by two K moves, where T,U belong to F_d. Each common S counted by B[T,S] B[U,S] gives exactly one presentation. Indeed the unique U->S witness is also a directed path in S, because K preserves its chosen path. Its endpoints agree with those of the unique T->S witness by the score formula. Conversely every such two-step presentation yields the same common column. Hence

    W = B B^T,               W[T,T]=H(T).                (2)

This is a literal all-order correspondence, not an analogy between equal numerical counts.

## 2. Consequences and what the quotient loses

* W is symmetric, integral, entrywise nonnegative, and positive semidefinite.
* W[T,U] is the intersection size of two explicit sets of intermediate tournaments. In particular 0<=W[T,U]<=min(H(T),H(U)).
* Every principal minor of W is nonnegative. For three states with counts h1,h2,h3 and off-diagonal weights a,b,c, respectively on pairs12,13,23,

      h1 h2 h3 + 2abc >= h1 c^2 + h2 b^2 + h3 a^2.        (3)

* The return presentation mass M(T) is the row sum of W; the simple return degree is the number of positive off-diagonal entries. These are not interchangeable with H(T).
* The construction respects any group of relabellings preserving d. An orthonormal orbit compression X gives X^T W X positive semidefinite, but its diagonal sums within-orbit return weights and generally is not H.

The five-vertex regular example makes the loss explicit. There are 24 labelled states, all in one ordinary isomorphism class. The frozen packet proves that W is the disjoint union of two blocks

    15 I + 2 A_icosahedron.                              (4)

Thus its raw-count diagonal is 15 and every row sum is 25. The one-dimensional orthonormal isomorphism-orbit compression is the scalar25, not15. Meanwhile the path-normalization fibre H/|Aut| is3. The numbers15,25,3 are respectively raw paths, two-step presentation mass, and path-orbit mass. Their differences are forced by three distinct maps.

## 3. The actual intermediate tournaments are the icosahedral faces

The geometric factorization at order five is stronger than the Gram identity alone. Across all 24 regular sources, B has exactly 280 distinct intermediate-tournament columns: 240 columns have degree one, and 40 have degree three. On each 12-vertex sign component there are 120 private columns and 20 face columns. After column ordering,

    B B^T = 10I + F F^T = 15I+2A,                       (4a)

where F is the ordinary vertex-versus-triangular-face incidence matrix of that icosahedron. Thus the faces are represented by actual intermediate tournaments, not supplied by a decorative geometric analogy.

Here is a structural proof. The three path orbits in the frozen packet show that every regular source has ten intermediates with one endpoint-matching return path and five intermediates with three such paths. A degree-three intermediate has the following form, up to relabelling:

    t -> s -> C -> t,                                  (4b)

where C is a directed triangle on the remaining three vertices. Its score is three at s, one at t, and two on C. Its three s-to-t Hamiltonian paths are precisely the three rotations of C between s and t. Complementing along any of them gives a regular tournament, by the endpoint score formula. There are exactly 5*4*2=40 tournaments of form (4b).

Their three regular neighbours are pairwise adjacent, so each column gives a graph triangle. Distinct columns give distinct triples. Indeed, across the three neighbouring regular tournaments, the only edge whose orientation is unanimous is s->t: it is never on a selected path, whereas each s/C and C/t edge occurs on exactly one of the three paths and each triangle edge occurs on two. Thus the triple recovers the ordered pair (s,t); its majority orientation on the other three vertices recovers C and hence the entire intermediate S. The two icosahedra have exactly 40 graph triangles in total, so these 40 distinct triples are exactly all their faces. Every vertex belongs to five faces and every edge to two, proving (4a).

For an explicit representative, intermediate S has outgoing masks (28,1,10,18,6), code654. Its exceptional endpoints are s=0,t=1, its interior cycle is 2->3->4->2, and its regular neighbours have codes947,613,217, the face (N,U_0,U_1) in the displayed certificate.

There are two different natural stochastic kernels. Uniformly choosing one of all 25 path-pair presentations gives W/25=(3/5)I+(2/25)A. Choosing P uniformly among the 15 paths first, then Q uniformly among its one or three matching-endpoint continuations, gives instead (7/9)I+(2/45)A. The latter self-loop probability is 10/15+5/(15*3)=7/9. These are distinct sampling rules, not contradictory normalizations.

## 4. A small exact spectral consequence

For completeness the spectrum of each block (4) can be derived directly from the displayed two-pentagon graph, without numerical eigensolvers. Let A be its adjacency matrix and C its antipodal involution matrix. Its five common-neighbour counts give

    A^2 = 3I + 2J - 2C.

On antipodally even functions the graph quotient is K6, so A has eigenvalues5 once and -1 five times. On antipodally odd functions J=0,C=-I, hence A^2=5I. Its trace there is zero, giving +/-sqrt(5), each three times. Consequently the eigenvalues of (4) are

    25 (multiplicity1), 13 (multiplicity5),
    15+2sqrt(5) (multiplicity3), 15-2sqrt(5) (multiplicity3).

All are positive. This finite spectral observation concerns the derived regular return matrix only, not the spectrum of the original merged metagraph.

## 5. Reproduction

`check_gram_return.py` independently constructs every labelled tournament and every Hamiltonian path through order five, groups tournaments by exact labelled score vectors, and verifies (1)-(2) by a direct enumeration of all second paths from every intermediate state. It additionally verifies the icosahedral matrix relation exactly with integer arithmetic. The test universe is all labelled tournaments, not just regular examples. Normal and optimized Python modes must agree.
