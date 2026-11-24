# The 56-site star versus Hoffman–Singleton

## Decision

Neither is a possible **involution** fixed graph in a degree-57 Moore graph. The established fixed-set classification already picks the star in the involution case; Ishida then rules out the involution itself. A 50-vertex Hoffman–Singleton fixed graph instead belongs to the **order-five** search. This does not establish that the latter search has a solution, or rule out an induced Hoffman–Singleton copy without a pointwise symmetry.

Sources: [Ishida, No involutions in the missing Moore graph](https://arxiv.org/html/2606.29183v1), §§2, 4–5; [Mačaj and Širáň, Search for properties of the missing Moore graph](https://www.sciencedirect.com/science/article/pii/S0024379509003735), Lemmas 4, 12, 16 and Theorem 5. For the earlier trace replay and separate tournament-address results see [the involution correction](../moore_merged_count_addresses_20261124/INVOLUTION-CORRECTION.md).

### One explanation for both involution failures

In either proposed fixed graph \(F\), any outside vertex is adjacent to at most one fixed vertex: fixed pairs are adjacent or already have their unique common neighbor in \(F\). If \(x,tx\) are nonadjacent, their unique common neighbor must itself be fixed. Conversely, if \(x\) has a fixed neighbor then \(x,tx\) cannot be adjacent. Thus the vertices moved to adjacent mates are exactly the outside vertices not touching \(F\). Their number is
\[
a=(3250-|F|)-(57|F|-2e(F)).
\]
The rational \(-8\) spectral trace is \((7|F|+50-a)/15\), and the rank of \(I+J+A_F\) over \(\mathbb F_2\) must equal it by the trace-rank identity.

| Proposed fixed graph | Fixed edges | Outside vertices touching it | \(a\) | Trace | Rank mod 2 |
|---|---:|---:|---:|---:|---:|
| \(K_{1,55}\) | 55 | 3082 | 112 | 22 | 54 |
| Hoffman–Singleton | 175 | 2500 | 700 | \(-20\) | 28 |

The star rank is that of \(I_{55}+J_{55}\) on its leaves. For the second rank, \(I+J+A_F\) is the reduction of the projector \((15I+5A_F-J)/25\) onto the Hoffman–Singleton eigenvalue \(2\) of multiplicity 28. Its integral localized image is a direct summand so reduction preserves its rank. An exact 56-site star fixed set is not an odd-order alternative either: a prime-order element with this fixed set would have to fix 55 neighbors at its center, leaving exactly two nonfixed neighbors; its prime order must divide two. The star can still occur as an induced graph with no such symmetry.

## The order-five continuation: paired matchings must form alternating rings

Assume from here that \(t\) has order five and fixes a Hoffman–Singleton graph \(F\) pointwise. Let \(C_f\) be the 50 outside neighbors of \(f\in F\), \(C=\coprod C_f\), and \(Z\) the 700 sites not attached to \(F\). This yields 500 free orbits in \(C\) and 140 in \(Z\).

* Every \(z\in Z\) has one neighbor in each \(C_f\) and seven neighbors in \(Z\). Each \(C_f\) partitions \(Z\) into 50 sets of size 14.
* The graph on \(C\) is a 50-sheet matching cover of \(\overline F\). For \(x\in C_f\), its common neighbor with \(g\ne f\) is \(f\) if \(fg\) is an edge; otherwise it is the unique matching partner in \(C_g\). There are 42 such \(C\)-neighbors and 14 \(Z\)-neighbors. After taking \(t\)-orbits this is a 10-sheet cover with a phase in \(\mathbb Z/5\) on each edge.
* Each \(Z\)-orbit is independent or an internal pentagon with step class \(\{\pm1\}\) or \(\{\pm2\}\). Two distinct \(Z\)-orbits can be joined by at most two matching voltages. A double matching between \(U,V\) uses one absolute difference class: the two edges through a vertex of \(V\) connect a pair of that difference in \(U\), and vice versa. A repeated difference would give two common neighbors, and three matching voltages would repeat a difference.
* For two distinct vertices of the **same** \(Z\)-orbit, their edge or unique common neighbor must live inside \(Z\): their neighbors in each \(C_f\) are distinct translates of a free vertex. An internal pentagon already accounts for all ten pairs, so has no incident doubled matching. On an independent orbit, all ten pairs must instead be accounted for by doubled matchings, five pairs of each difference class. Consequently there is exactly one doubled edge of each color at every independent orbit.

**Ring lemma.** The doubled edges therefore form vertex-disjoint even alternating rings on the independent orbits, with length at least six. The only shorter candidate is a four-ring. Write a doubled voltage set as \(m+\{-r,r\}\). Two consecutive ring edges are of opposite classes, so their length-two walks give the four distinct offsets \(m_1+m_2+\{\pm1,\pm2\}\). The two ways around a four-ring to its opposite vertex would then give intersecting sets of offsets, impossible by unique common neighbors. This proof does not depend on a choice of phase origins.

**Phase constraints away from the ring.** Along two consecutive doubled edges, *any* additional route of length at most two must use the unique missing offset \(m_1+m_2\); in particular the orbit endpoints can have the same neighbor-orbit in at most one of their 50 attachment coordinates \(C_f/t\). A three-edge ring segment between distinct endpoints uses all five offsets, so its endpoints cannot be joined by an extra simple edge: the joined lift would contain a square. Adjacent cyclic orbits must have opposite internal step classes (equal steps give a square). For an opposite-class adjacent pair the direct matching and the paths using one internal edge already use all five offsets, so they cannot share any attachment coordinate. Each ring can be gauge-normalized to zero center voltage on all but one edge, retaining a single total center flux. This is a normalization of ring data, not a claim that the other voltages of the full graph disappear.

The known five trace cases give a short list of **ring supports**:

| \(a=|\{v:v\sim t(v)\}|\) | Internal pentagons | Independent ring sites | Multiplicity of \(7\) at a nontrivial character |
|---:|---:|---:|---:|
| 25 | 10 | 130 | 341 |
| 100 | 40 | 100 | 340 |
| 175 | 70 | 70 | 339 |
| 250 | 100 | 40 | 338 |
| 325 | 130 | 10 | 337 |

To read the table: the \(7\)-eigenspace has dimension 1729 and trace \((335+a)/15\) at \(t\). Subtracting the trace from 1729 gives five times the nontrivial-character multiplicity. Invariance of the trace under \(t\mapsto t^2\) balances the two internal step classes at \(a/5\) pentagons each. In the last row there is **one alternating decagon**, not a four-ring plus a six-ring. Cyclic vertices have five simple \(Z/t\) neighbors; independent vertices have three simple neighbors and their two doubled edges. Neither the last row nor any other row is proved globally feasible.

### Exact simultaneous test, without flattening the voltages

Use the ring \(\mathbb Z[z]/(z^5-1)\), with \(U=1+z+\cdots+z^4\). Let \(M\) be the \(640\times640\) free-orbit matrix with \(M_{ij}=\sum z^a\) for its matching shifts. Let \(S\) be the \(50\times640\) fixed-to-orbit incidence (one 1 per attached column, zero on residual columns), \(P=A_F\). An intended simple undirected lift uses reciprocal voltages \(M_{ji}(z)=M_{ij}(z^{-1})\), distinct shifts in each entry, no zero shift on the diagonal, ordinary symmetric loopless \(P\), and all lifted degrees 57. It is a Moore graph **if and only if**
\[
\begin{aligned}
P^2+P+5SS^\top &=56I+J,\\
PS+SM(1)+S&=J,\\
M^2+M&=56I+(J-S^\top S)U .
\end{aligned}
\]
These equations are the fixed/fixed, fixed/free and free/free blocks of the common-neighbor identity \(A^2+A=56I+J\), with the last block retaining each displacement separately. In particular \(M(\zeta)^2+M(\zeta)=56I\) at nontrivial fifth roots. The replay verifies the analogous equations on the actual 50-site Hoffman–Singleton graph with a fixed pentagon, and tests all five ring center fluxes at lengths 4, 6, 8 and 10. A ring alone passing its local girth test is not a solution of the full equations. For the 325 case the smallest natural initial branch is the decagon and its attachment-coordinate exclusions, not all 3200 moved sites.

There is also a distinction between the cyclic and reflection uses of a tournament tiling. Since involutions in the host are ruled out, its automorphism group is odd. The normalizer of an order-five subgroup cannot act nontrivially on it (its induced action has order dividing four), so an actual host automorphism cannot exchange its two step classes for a fixed generator. A change of generator can rename both colors; it is a coordinate choice, not a reflection symmetry on the Moore graph. For the tournament encoding, retain the ordered voltage data as a payload of the [guarded +2 chart](../moore_merged_count_addresses_20261124/README.md), rather than erasing it in an unphased count.

## Exact tournament comparison: one class, two different marked exits

Independently of any degree-57 host, take the 56-site star on \(a,L_1,\ldots,L_{55}\) with one helper \(b\). Start with the spine \(a,b,L_1,\ldots,L_{55}\) and reverse exactly the 55 star edges. This is the three-block cycle \(a\to b\to L\to a\), with \(L\) an ordered chain. Generalizing to \(m\) leaves, a Hamiltonian path either uses the consecutive pair \(a,b\), allowing any subset of \(L\) before the pair, or is the single path \(b,L,a\). The count is \(2^m+1\).

The tournament is rigid for \(m>1\): \(a,b\) are the only vertices in \(m\) directed triangles while every leaf is in one, and \(a\to b\) distinguishes the two. Its converse is realized by swapping \(a,b\) and reversing the chain. On Hamiltonian paths this operation takes the subset before \(a,b\) to its reversed complement. For odd \(m\) the central leaf prevents a fixed subset, so \(b,L,a\) gives the **unique blue tiling**, with just one nonspine backward edge. The star-edge spine is black. These are two markings of the *same* node.

The word-mate transfer from the earlier session gives the following exact \(m=55\) addresses. After reflection the black word is \(1^{55}0\), and the blue word is \(10^{55}\). The black mate has nonsymmetric scores (a sink but no source); the blue mate is self-converse by the marked reflection. Every score multiplicity of each mate is at most two, so odd tournament automorphism order makes both rigid.

| Merged metagraph node | Hamiltonian paths of each class | Number of tilings in merged node |
|---|---:|---:|
| Common source | 36,028,797,018,963,969 | 36,028,797,018,963,969 |
| Mate along the marked black star line | 288,053,426,793,161 | 576,106,853,586,322 |
| Mate along the unique blue line | 76,560,404,850,031,802,201,707 | 76,560,404,850,031,802,201,707 |

These are tile-complement exits, not single-arc moves. This distinguishes exactly what the star can do for the metagraph from what it cannot do for the missing graph: the blue exit persists in the tournament, even though it cannot be interpreted as a compatible nontrivial reflection of the degree-57 host. The surviving graph search under this fixed-set assumption is a five-phase Hoffman–Singleton extension, with the two phases of a doubled matching coupled around even rings.
