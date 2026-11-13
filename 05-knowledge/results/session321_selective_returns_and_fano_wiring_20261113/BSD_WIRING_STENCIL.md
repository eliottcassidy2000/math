# A binary path-complement stencil for the descent matrix

This is a separate metagraph constraint using the row-zero tournament
matrix from [BSD parity](BSD_PARITY.md). It applies to a labelled
Hamilton-path alignment. It is an exact binary matrix identity, not
a statement about ranks of elliptic curves or a complete wiring test.

Let A be the row-zero matrix and L be the ordinary undirected
Laplacian of the marked Hamilton path, reduced mod 2. The line partner
is the converse except on the path, so its row-zero matrix is
```
B = A + J + (n mod 2)I + L
```
with J all ones. Because JL=LJ=0 and the descent matrix is unchanged
by the converse, the exact line difference is
```
Phi(B) + Phi(A) = A L + L A + L² + L.
```
Indeed shift by I leaves Phi invariant and adding L to the converse
gives the displayed expansion. Its row sums are zero. It is
symmetric because A^T=A+J+I and L=L^T. Its diagonal is zero:
the i-th diagonal of AL+LA is the path degree of i modulo two,
L² has diagonal zero since every row of L sums to zero, and
the diagonal of L is again that degree. Thus the line changes Phi
by an alternating form with the all-ones vector in its radical.
In fact Phi(A) itself is symmetric, has diagonal all ones and row
sums all ones; these are separate one-node constraints.

Number vertices by the spine order. For two interior positions i,j
with distance at least three, the identity reduces to the
four-neighbor stencil
```
(Phi(B)+Phi(A))[i,j]
 = A[i,j−1]+A[i,j+1]+A[i−1,j]+A[i+1,j].
```
Here there is no L or L² entry at distance at least three, and the
diagonal of L at interior vertices is zero. Near the two ends and
the short diagonals one must use the full matrix formula.

Together with the [wedge inverse](WEDGE_WIRING.md), this suggests a
two-stage metagraph test: first use candidate endpoint/spine-local
constraints (degree and the binary stencil) to prune alignments,
then use the full wedge-difference support to certify a literal
line. The stencil is necessary; equality of descent nullities alone
is neither an alignment certificate nor a preservation theorem.
