# What one attached vertex can recover

**Status: PROVED for labelled tournaments and the stated query language; hostile tests FINITE-EXACT.**

## Current coloring data does not determine future coloring data

The two twelve-vertex tournaments in the native 81/27 control have exactly the same 67 labelled cyclic triples. Consequently they have exactly the same 547 transitive supports, including empty, and the same answer to *every* current named weighted acyclic-coloring problem. Their automorphism orders remain 81 and 27; their regular binary-mask counts remain 1620 and 1728. Each has dichromatic number 3 and maximum transitive size 6.

Append the same new vertex p=12 to both, with 8→p and p→v for every other old vertex. The outgoing mask of p is 3839. The query on named vertices {8,9,12} has dichromatic number 1 in the first extension and 2 in the second. These are actual extensions of the pinned native pair, not an abstract replacement example.

## A sharp general statement

**A1 (one-port arc recovery).** Let u,v be two named old vertices. Attach p with \(u\to p\to v\); orient its other arcs identically in all comparisons. Then
\[
\chi(T[\{u,p,v\}])=
\begin{cases}
1,&u\to v,\\
2,&v\to u.
\end{cases}
\]
The first case is transitive and the second is a directed triangle. Thus any distinct pair of labelled tournaments is distinguished by one new vertex and a query on three named vertices. Zero added vertices do not suffice for the full palette language, as the 81/27 pair proves. One is therefore sharp.

**A2 (carrier consequence).** A deterministic representation of labelled n-vertex tournaments which supports every such named attachment and query must have at least \(2^{\binom n2}\) distinct states, or at least \(\binom n2\) bits in a fixed-length encoding. Otherwise two distinct inputs would have the same state and receive the same answer under the distinguishing continuation from A1.

This is an information bound for exact recovery under a specified interface. It is not a time lower bound, a circuit lower bound, or an ETH consequence. The much smaller invariant \((\chi,\alpha)\) is sufficient for the restricted substitution recurrences in this packet; allowing arbitrary vertex attachments changes the required information.

## A faithful clique conversion

**A3 (rank-layer clique bijection).** For a tournament T and positive k, form an ordinary graph with vertices \((v,i)\), \(v\in V(T)\), \(i\in\{1,\ldots,k\}\). For i<j, join \((v,i)\) to \((w,j)\) exactly when v≠w and v→w. There are no same-layer edges.

Every k-clique has one vertex in each layer and gives a distinct-vertex forward ordering, hence a transitive k-subtournament. Conversely every transitive k-subset has a unique forward order and gives exactly one clique. This is a bijection, preserving named vertices and product weights assigned to their copies.

This is the standard use of auxiliary position states in a clique encoding; no algorithmic speedup is claimed. The checker compares actual clique-compatible ordered tuples with transitive subsets for P7 and H at k=1,…,4.

On the original vertex set, no ordinary graph can represent all transitive subsets of a tournament containing a cyclic triple as its independent sets: every pair is transitive, so that graph would have no edges and would incorrectly accept the triple. The rank layers repair a type mismatch by keeping order information.

The LRC companion makes a different repair: it replaces whole runner safe sets by their actual interval components. Both constructions use cliques, but the preserved predicates and the necessary added coordinates differ. A graph analogy without those coordinates is insufficient.
