# Two nontrivial components at a first birth need at least sixteen vertices

**PROVED by finite witness covers.** Suppose a tournament has the smallest possible order among tournaments with its Hamiltonian path count. Any two of its nontrivial strongly connected components have total order at least 16. At total order 14 their hypothetical joint contribution can always be realized on at most 12 vertices. At total order 15 it can always be realized on at most 14.

Every Hamiltonian path visits the components consecutively in their canonical dominance order, so the count is the product of their individual counts. In a minimum-order witness each component must therefore already have minimum order for its own count. Moreover, if the product of the counts of any two components is attainable on fewer vertices than their combined order, replace those two factors with the smaller tournament. The scalar product allows this even when other components originally lay between them. The result has the same count and smaller order. The [exact product test](birth_cuts.py) covers total orders through 13 using the [order-by-order census](../path_spectrum_two_angles_metagraph_20261106/data). The following two covers extend this criterion.

## What needs a new witness

For total size 14, the possible component order pairs are (3,11), (4,10), (5,9), (6,8), (7,7). The minimum-order heights in every row are known by the census. Form their products and remove every height already attained by order 11. Exactly 29,253 distinct numbers remain, from 80,517 through 285,285. All have a 12-vertex witness.

For total size 15, the pairs with smaller member at least four are again covered by known minimum-order sets. The remaining pair (3,12) has first height 3. If the second height is b, deleting the first vertex of any Hamiltonian path shows that b is at most 12 times the order-11 maximum, or 1,141,140. Indeed for each first vertex the remaining path occurs on its deletion. Use **every** odd b from 3 through that upper bound, omitting the impossible heights 7 and 21. Taking the union of 3b with those other pair products and subtracting heights already attained by order 11 leaves a superset of all remaining possibilities. Its 584,052 distinct numbers run from 80,453 through 3,423,417. Every one has a 13- or 14-vertex witness. In particular, no classification at order 12 is needed to exclude the (3,12) pair.

## One parent supplies a quadratic table of counts

Fix a parent tournament T on n labelled vertices. The added vertex beats the set with indicator bits x. Let e(S,i) and s(S,i) count Hamiltonian paths of a nonempty induced subset S ending or starting at i. For distinct i,j put

\[
 Q_{ij}=\sum_{i\in S,\,j\notin S} e(S,i)s(V-S,j).
\]

Splitting each Hamiltonian path at the new vertex counts each exactly once. At the first position it contributes a path starting at a vertex the apex beats. At the last position it contributes a path ending at a vertex that beats the apex. At an interior position it contributes precisely one of the two-piece paths counted by the corresponding entry of Q. Consequently

\[
 H(x)=\sum_j s(V,j)x_j + \sum_i e(V,i)(1-x_i)
                    +\sum_{i\ne j}Q_{ij}(1-x_i)x_j.
\]

The standard endpoint subset recurrence computes e, and the same recurrence on the converse computes s. If j is adjoined to a bit mask S, the increment in the displayed quadratic is its linear coefficient minus the sum of \(Q_{ij}+Q_{ji}\) for i in S. Thus all mask values can be filled without performing a fresh Hamiltonian path computation for each apex. The formula retains a literal tournament for every value, given by the parent rows and the mask. This specialization supplies finite witnesses; it makes no distributional assumption on the parents.

## A small reproducible description of both covers

The C++ source fixes the integer generator, seed, edge order, and density schedule for each stage. All apex masks of each resulting parent are evaluated by exact integer arithmetic. A marked target is removed only at an actual equality. With both parent rows and a mask given, a direct endpoint recurrence verifies the child independently. The maximum intermediate counts fit 64-bit signed arithmetic: each cross summand counts a concatenated path permutation with one specified split and the sum for fixed endpoints has at most n! possibilities. Here n is at most 13, and the polynomial updates sum at most a few dozen such terms. The checked direct controls include extreme and mixed masks at every stage.

| Component order sum | Residual targets | Parent family | Targets still missing |
|---:|---:|---|---:|
| 14 | 29,253 | 950 parents of order 11 | 0 |
| 15 | 584,052 | 20,000 parents of orders 12 and 13 | 6 |
| 15 | 6 | 1,057 parents of order 13 on a second fixed schedule | 0 |

For example, the largest missing product on the first row is 285,285. One 12-vertex witness has outgoing masks **1712, 3553, 1963, 3155, 2310, 3544, 2325, 3416, 3081, 506, 2640, 517** in vertex order starting at zero.

From this directory, run `python3 -B birth_cut_cover.py 14` or `python3 -B birth_cut_cover.py 15`. The Python driver reads the exhaustive small census, compiles the C++17 oracle in a temporary directory, and requires zero remaining targets before reporting success. Add `--witness 285285` to the first command to print and independently recount the example; any height in the respective residual set is supported. Whether two nontrivial components of total order 16 can occur at a first birth remains open here. [Target generator](birth_cut_cover.py) · [All-mask oracle](birth_cut_cover.cpp)
