# The five-parameter positive ternary clock family

## Scope and physical graph

This appendix gives an exact all-order formula for a physical family surrounding the positive-offset ternary clone clock. All five multiplicities below are positive integers. It proves the complete cancellation classification within this family, and an actual returning four-port construction. It does not classify arbitrary ternary clocks or prove a universal digit grammar.

Let the seven-vertex quotient, in the order `(r,P,U,a,R,b,V)`, have outgoing bit rows

`Q=[14,32,2,22,103,13,47]`.

Replace these vertices by transitive modules of respective lengths

`(m,p,u,1,q,1,v)`.

All edges between modules follow Q; edges within each module point from the earlier to the later vertex. Call the resulting tournament `G(m,p,q,u,v)`. The quotient contains the directed Hamilton cycle `(0,2,1,5,3,4,6,0)`. Thus every graph in the family is strong: between any two vertices, use quotient paths, and between two vertices in one module use a quotient cycle if their internal edge has the wrong direction.

For `m=1,p=q=2,u=v=1`, the rows are

`[30,132,128,6,110,463,399,25,159]`.

The relabeling from these vertices to the positive-clock source

`[30,64,10,66,174,335,21,367,95]`

is `[0,3,1,2,4,7,5,6,8]`, preserving the marked root0. The two length-two modules are the old vertex sets `{1,3}` and `{5,7}`. This is the same positive clock, not a distinct discovery or an isomorphism-class count.

## Exact formula

Put `x=2^p`, `y=2^q`, `z=2^u`, and `w=2^v`. Then, at every positive choice of the five multiplicities,

`H(G)=A*3^m+B*2^m+C`,

where

`A=2zw(x-1)(y-1)`,

`B=w(y-1)(4-x)+z(x-1)(4-y)`,

`C=(x-3)(y-3)`.

In particular, at p=q=2 the second mode cancels exactly and

`H(G(m,2,2,u,v))=18*2^(u+v)*3^m+1`.

This gives all the rays `72*2^t*3^m+1`, with t>=0 and m>=1. The formulas do not extend by simply putting a multiplicity equal to zero: empty-module boundary terms are outside the theorem.

## Why a finite grid proves this formula for all orders

We give uniform physical bounds before using any finite computation. For a vertex replaced by a transitive module of length s, a Hamilton path breaks its module occurrences into j nonempty runs. Delete the module vertices and place one copy of its root in each run. This gives a root skeleton with no consecutive roots. Conversely, a skeleton with j roots lifts in `j! S(s,j)` ways: distribute the s ordered module vertices surjectively among the j ordered runs, keeping each run increasing. Therefore

`H_s=sum_j a_j j! S(s,j)`,

where a_j counts literal skeletons. The j-1 internal fragments start at outneighbors of the root and end at inneighbors. On each such directed fragment, choose an edge crossing from its outneighbor side to its inneighbor side. The fragments are vertex-disjoint, so these chosen edges form a matching in the directed bipartite crossing graph. A vertex cover of size d for that graph gives j<=d+1, uniformly in all the other module lengths.

Here are fixed singleton vertex covers for our quotient:

- For r, all crossing edges are from P to the singleton b, or from the singleton a to R. The cover `{a,b}` has size2, so r has at most3 runs.
- P has only the singleton b as an outneighbor outside P. The cover `{b}` bounds P by2 runs.
- R has only the singleton a as an inneighbor outside R. The cover `{a}` bounds R by2 runs.
- U has only P as outneighbors, and the only crossing edges from P go to b. The cover `{b}` bounds U by2 runs.
- V has only R as inneighbors, and the only crossing edges into R come from a. The cover `{a}` bounds V by2 runs.

The same singleton covers work after every other module is enlarged. Consequently, as a function of each of p,q,u,v separately the count lies in the span of `1,2^s`, and as a function of m it lies in the span of `1,2^m,3^m`.

These separate uniform statements imply membership in the 48-dimensional tensor span. Indeed, recover the two coefficients in p from p=1,2 using the invertible matrix `[[1,2],[1,4]]`. Each recovered coefficient is a fixed linear combination of counts at those p values and retains the stated separate spans in the remaining variables. Repeat for q,u,v, then use m=1,2,3. The final matrix in m is the Vandermonde matrix on bases1,2,3, with nonzero determinant12. No coefficient depends on a chosen sample bound.

The verifier calculates the 48 physical counts on the grid `p,q,u,v in {1,2}`, `m in {1,2,3}` by subset Hamilton-path DP and checks the displayed formula. This proves the all-order identity in the established tensor span. It also calculates all sixteen full root-run vectors independently and checks the conversion

`A=a3`, `B=a2-3a3`, `C=a1-2a2+3a3`.

The source run vector at p=q=2,u=v=1 is `(a0,a1,a2,a3)=(73,217,216,72)`. In particular the constant is1, not a0=73. Here a0 is the Hamilton count after deleting the root, whereas the clone formula uses only positive run indices.

## Complete cancellation classification in this family

The coefficient B vanishes if and only if p=q=2.

Since x,y are powers of two at least2, each is2,4, or at least8. If x=y=4, both summands vanish. If one is4 and the other is not, precisely one summand is nonzero. If both are2, B is positive; if both exceed4, it is negative.

The only remaining case, up to exchanging `(x,z)` with `(y,w)`, has x=2 and y>=8. Vanishing would require

`z/w=2(y-1)/(y-4)`.

The left side is an integral power of two, possibly with negative exponent. The right side lies strictly between2 and4, so equality is impossible. This exhausts every positive integer choice, not merely a bounded search.

As A is positive, linear independence of the bases1,2,3 means that an eventual pure ternary ray `alpha*3^m+beta` in this family is equivalent to B=0. Hence all such rays are exactly the displayed positive-offset family, with beta=1 and alpha=72 times a power of two. In particular this entire family cannot realize the proposed amplitude180 clock or a negative-unit-offset clock. This does not rule out those amplitudes or signs in other quotient families.

## Physical four-port return, with fixed old vertices

Set p=q=2 and u=1; let v>=1 be arbitrary. Start from `G(1,2,2,1,v)`, in the consecutive-module labeling above. Keep the same four physical ports `[1,2,3,4]` forever. These are the two P vertices, the single U vertex, and a.

At each right birth, the newborn beats all four ports and loses to every other old vertex. In mask notation this is `R(mask15,drop4)`, where drop4 forgets only the newborn's temporary exposure. All four old ports remain exposed, with the same labels, after every step. Every added vertex is therefore another transitive clone of the original root0, with every older root clone preceding it.

After m-1 births the graph is isomorphic to `G(m,2,2,1,v)`. If the seed order is N=v+8, the canonical root module maps to physical labels `[0,N,N+1,...,N+m-2]`; the other canonical vertices, in order, map to the unchanged labels1,...,N-1. This describes every incidence, not just the scalar response.

Thus the count at every stage is `72*2^(v-1)*3^m+1`, and the literal same-interface move returns

`H_new=3H_old-2`.

The port count is four; this is not a construction in the earlier two-port or three-port language. All graphs remain strong. The result is a returning single-sign family, not a common-state choice between two ternary digits.

## Verification and provenance

Run `python3 verify.py` and `python3 -O verify.py`. Both use only the Python standard library and print the same deterministic receipt. Besides the complete48-point grid and sixteen full run vectors, the verifier checks source relabeling, literal permutation count217, additional off-grid counts, the matching vertex covers, graph strongness, and the fixed-port return with every edge checked.

The positive source was supplied by the Session197 path-polynomial lane, independently isomorphic to the positive clock from the orbit/symmetry lane. The new content here is the complete five-parameter physical identity, the uniform tensor proof, its all-order cancellation classification, and explicit returning four-port convention. Arbitrary nontransitive module substitutions require a separate whole-module theorem and are not asserted in this appendix.
