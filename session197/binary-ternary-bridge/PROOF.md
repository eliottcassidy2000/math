# A physical sign bridge and a returning binary–ternary family

## Result and scope

There is a strong physical two-parameter family with

`H(m,t)=60*2^t*3^m+1`, for all m,t>=1.

It has the same five old ports at every stage. One literal birth returns `2H-1`, and another returns `3H-2`; they may be interleaved arbitrarily. Before the first binary birth, the same construction is the negative clock `60*3^m-1`. The first binary birth therefore acts by `2H+3`, changing the sign of the constant term. Later binary births act by `2H-1`.

This is an explicit returning two-operation arithmetic family, with an exact inverse criterion below. It is not a universal ternary digit grammar: the two operations have different multipliers, and the scalar image remains the stated restricted set. No reverse sign gate is asserted in this packet.

## Graphs, labels, and retained ports

Start with the eight-vertex tournament

`T=[2,12,145,149,99,143,175,19]`.

The marked root is vertex7. Its outneighbors are exactly0,1,4. Keep the ordered five-port list

`[0,1,4,5,6]`

forever. Neither birth replaces an old port.

- R is a right birth with mask7: its newborn beats old0,1,4 and loses to every other old vertex.
- L is a left birth with mask24: its newborn loses to old5,6 and beats every other old vertex.

In the convention that one temporarily exposes the newborn then forgets one of the six exposed vertices, both operations use drop5. Thus every newborn is immediately forgotten as a port while remaining in the graph. The original five physical vertices remain exposed, in the same order, after every operation. No old vertex is re-exposed.

A graph after a word with m-1 occurrences of R and t occurrences of L has three parts:

1. Old vertices0,...,6, with their original incidences.
2. The original root7 and the R newborns, forming a transitive module of length m, ordered by birth time. Each beats old0,1,4 and loses to old2,3,5,6.
3. The L newborns, forming a transitive module of length t, ordered in reverse birth time. Every one beats old0,...,4 and every root-module vertex, and loses to old5,6.

All edges between parts are specified here. Denote this graph by F(m,t), including t=0. A canonical labeling puts old0,...,6 first, then the m root vertices, then the t bridge vertices, with each module in its transitive order.

The relabeling from a physical word to this canonical graph fixes0,...,6. The canonical root list maps to original7 followed by the R newborn labels in chronological order. The canonical bridge list maps to the L newborn labels in reverse chronological order. This proves that all interleavings have the same graph up to the displayed relabeling. Equivalently, opposite-side births commute because the L newborn beats the R newborn in either order.

At m=t=1, the source rows are

`P=[2,12,145,149,99,399,431,19,159]`.

The extra vertex8 loses only to5,6, exactly as specified.

## Strongness at every stage

T has the directed Hamilton cycle `(0,1,2,4,6,5,3,7,0)`, so it is strong. P adds a vertex with both an inneighbor and an outneighbor in T, for example5 and0, and is also strong. Replacing a vertex of a strong tournament by any nonempty transitive module preserves strongness: paths through other quotient vertices connect its internal vertices in either direction. Thus F(m,0) and F(m,t) are strong for every m>=1 and t>=1. The construction never needs a neutral source or sink prefix.

## All-order proof, before any interpolation

For a fixed root replaced by a transitive module of s vertices, every Hamilton path uses j nonempty runs of that module. Removing its vertices gives a literal root skeleton with j root occurrences and no consecutive roots. The exact lift count of one skeleton is `j! S(s,j)`: assign the s ordered vertices surjectively to the j ordered runs, with increasing order inside each run. Thus

`H_s=sum_j a_j j! S(s,j)`.

Each of the j-1 internal fragments starts at an outneighbor and ends at an inneighbor of the root. Choose an out-to-in crossing edge on every such fragment. The fragments are vertex-disjoint, so these edges form a matching in the root's crossing graph. Any vertex cover of that graph bounds the number of internal fragments.

In F(m,t), for the root module the crossing edges are exactly

`1->2, 1->3, 4->5, 4->6`.

The bridge vertices are inneighbors of the root and receive no crossing edge from its outneighbors. The fixed singleton cover `{1,4}` therefore bounds the root run index by3, uniformly in t.

For the bridge module, the inneighbors are exactly5,6. The only crossing edges from its outneighbors into that set come from old4. The fixed singleton cover `{4}` bounds the bridge run index by2, uniformly in m.

Consequently, for positive m,t the count belongs to the six-dimensional tensor span

`{1,2^m,3^m} tensor {1,2^t}`.

To justify the tensor statement rather than assume it, recover the two coefficients in t from t=1,2 using the invertible matrix `[[1,2],[1,4]]`. Each coefficient is a fixed linear combination of physical counts and therefore retains the uniform three-dimensional span in m. The matrix for m=1,2,3 on bases1,2,3 has determinant12, so six grid values determine the function at every positive m,t.

The six physical counts, obtained by the supplied subset counter, are

- t=1, m=1,2,3: `361,1081,3241`;
- t=2, m=1,2,3: `721,2161,6481`.

They uniquely give

`H(F(m,t))=60*2^t*3^m+1` for m,t>=1.

The proof is all-order because the uniform crossing covers establish the tensor span before the finite grid is evaluated. The verifier also computes the independent literal root-skeleton coefficients at the nine-vertex source:

- Root7: `(a0,a1,a2,a3)=(121,361,360,120)`.
- Bridge vertex8: `(a0,a1,a2)=(179,361,180)`.

Here a0 denotes the count after deleting that root; the clone formula itself uses only positive run indices.

## The boundary t=0 and the first sign change

The old root in T has full positive run vector `(a1,a2,a3)=(179,180,60)`, with all higher entries zero. Therefore, for every m>=1,

`H(F(m,0))=179+180(2^m-2)+60(3^m-3*2^m+3)=60*3^m-1`.

The preceding positive-t formula must not be evaluated at t=0: it would give the wrong sign. On the actual boundary the operations are

- R: `3H+2`, staying on t=0;
- first L: `2H+3`, entering t=1.

Once t>=1, they are

- R: `3H-2`, increasing m;
- L: `2H-1`, increasing t.

The first L may occur after any number of R steps. The same five old ports remain available throughout, and any further R/L word stays in the positive family. In particular a fixed interior word with r R births and l L births acts exactly by

`H_new=2^l*3^r*(H_old-1)+1`.

This common physical interface is closed. It does not supply a choice between3H+2 and3H+4 on one interior state.

## Exact inverse criterion inside this family

An integer N is an interior count if and only if

`N-1=2^a*3^b*5`, with a>=3 and b>=2.

For such a target choose `m=b-1`, `t=a-2`, start from T, and perform m-1 R births and t L births in any order. Conversely the proved formula forces precisely this factorization. This gives a target-independent inverse algorithm for the stated family, not for arbitrary tournament counts.

Equivalently, starting from the positive seed count361, repeatedly multiplying N-1 by2 or3 is implemented by the two returning births. The graph has order `7+m+t`.

## Verification and provenance

Run `python3 verify.py` and `python3 -O verify.py`. The deterministic receipts agree. The verifier checks all six tensor-grid values, the full root-run vectors and crossing covers, literal permutation count361, graph identities for every R/L word through length4, direct counts and strongness, additional off-grid pairs, and the inverse criterion on a bounded control interval.

The negative clock source is the previously audited Session197 ternary clock. The new result is the physical sign-changing bridge, the exact two-parameter positive family, its fixed five-port return, and its inverse factorization criterion. The related positive72 clock is a different source and is not used here.
