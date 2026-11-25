# Which copies of a 13-cycle actually fit together?

Write \(C_\ell\) for an undirected chordless cycle. As the allowed diameter varies, the graphs \(C_{2r+1}\) are the degree-two Moore graphs of diameter \(r\). At diameter two the degree-two example is the pentagon. In degree 57 the problem here always has diameter two. Vertex sets of cycles, their cyclic orders, and orientations chosen on their edges are three different pieces of input.

## The packing polynomial in Hoffman–Singleton

**FINITE-EXACT.** Give a variable \(z\) to each induced 13-cycle chosen, and count unordered collections with disjoint vertex sets in Hoffman–Singleton. The polynomial is
\[
        1+126000z+60984000z^2+8694000z^3.
\]
Fix one induced 13-cycle. Of the 968 choices of a disjoint second cycle, **712 cannot be extended to a third**. Some pairs have eight choices. Fifty sites prohibit four cycles.

To reproduce, compile and run the [intersection census](thirteen_packings.cpp) with C++17 and optimization. The affine construction and fixed cycle use the labelling in the [earlier orbit theorem](../moore_fixed_order_spectra_20261124/THIRTEEN-ODD-ALPHABET.md). The search chooses the least cycle vertex and the smaller of its two cyclic neighbors as its first step. It excludes used sites and all chords at each extension. For intersection sizes 0 through 13 with the fixed cycle, the counts in order are
\[
 (968,8002,25105,36483,30654,16205,6239,1759,468,86,30,0,0,1).
\]
Testing disjointness on the resulting 968 masks gives the extension counts. The previously established transitivity on induced 13-cycles lets us multiply the rooted answers by 126,000 and divide by two or six to forget the ordering.

## Every packing of at most 65 induced 13-cycles extends in degree 57

**PROVED conditional on existence.** Deleting any 857 vertices from a degree-57 diameter-two Moore graph still leaves an induced 13-cycle. In particular **every** disjoint packing of at most 65 induced 13-cycles extends, and every such graph contains at least 66 disjoint induced 13-cycles. The cycles may have edges between them. The hypothesis places no automorphism on the host.

It suffices to treat exactly \(s=857\) deletions by extending any smaller deletion set. For the deleted set \(S\), let \(d_S(v)\) count the neighbors in \(S\) of a uniformly chosen original vertex. Its mean is \(\mu=57s/3250\) and its variance is at most \(V=64s(3250-s)/3250^2\), using the absolute bound 8 for the nonconstant adjacency eigenvalues on the centered indicator of \(S\). By the one-sided second moment inequality at 20, the number of vertices touching at least 20 deleted sites is at most
\[
 \left\lfloor\frac{3250V}{V+(20-\mu)^2}\right\rfloor=1087.
\]
There remain at least 1,306 sites outside \(S\) touching at most 19 deleted neighbors apiece. A set of \(g\) vertices spans average degree at least \(65g/3250-8\), by applying the least eigenvalue bound to its centered indicator. At \(g\ge1306\) this exceeds 18, giving a nonempty induced subgraph of minimum degree at least ten by successively deleting vertices of degree at most nine.

Build an induced path of ten edges in this subgraph. When the current path has \(k\le9\) edges, the predecessor is forbidden, and every older vertex forbids at most one neighbor of the endpoint, by uniqueness of common neighbors and absence of triangles. This excludes at most \(k\) candidates. Every induced ten-edge path in the full graph has at least 40 chordless three-edge closures between its endpoints, with mutually disjoint added vertex pairs. At most \(19+19=38\) of the closures use a deleted vertex, since each of the two added vertices is adjacent to its respective endpoint. This proves the assertion. At most \(65\cdot13=845\) sites are used by the partial packing. [Prescribed-path closure lemma](../moore57_thirteen_20261125/CONTACTS.md)

## What survives when those cycles are directed

**PROVED.** Prescribe a disjoint family of consistently directed cycles of lengths \(\ell_1,\ldots,\ell_k\ge3\) on a labelled set and orient all other pairs as a tournament. Let \(b(m)\) be Busch's exact minimum number of Hamilton paths on a strong \(m\)-vertex tournament. The **exact minimum** over completions is
\[
                         \prod_{i=1}^k b(\ell_i).
\]
Indeed each cycle forces its own induced tournament to be strong. In the directed odd-cycle collection formula one may take an arbitrary collection inside each block independently. The weights multiply, proving the lower bound. Any minimizing strong tournament on \(\ell_i\) sites has a directed Hamilton cycle and so can be labelled to carry the prescribed one. Put the minimizing blocks and unused singleton sites in any ordinal order to attain equality.

Equality has a precise placement criterion. Each block must attain its own minimum, and the blocks and remaining singleton sites must be totally ordered by their cross edges. To see necessity, any directed triangle meeting different blocks or containing an unused site contributes a new positive term absent in the product. An exterior vertex nonuniform on a strong block makes such a triangle because there is an arc from the part it beats to the part that beats it. Thus the blocks are homogeneous, their quotient has no cyclic triangle, and the quotient is transitive. This also proves sufficiency. [Busch's strong minimum](https://doi.org/10.37236/1141) · [Grinberg–Stanley packing formula](https://arxiv.org/abs/2307.05569) · [THM-4744, Moon chains and Busch floors](../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md)

In particular \(b(13)=625\). Orienting the 66 forced disjoint cycles consistently, then completing the pairs arbitrarily to a tournament, gives at least \(625^{66}\) Hamilton paths. One may leave all other edges initially unoriented to attain this minimum exactly by the preceding theorem. This statement uses the Moore graph to select supports. Imposing further edge orientations can raise the minimum. An explicit 13-fold automorphism is a separate condition in the [orbit chart](../moore57_thirteen_20261125/THIRTEEN.md).

The repo’s even-path forest expansion gives an exact averaging formula that retains the positions of individual edges. [THM-5313, switching derivatives and forests](../../../01-canon/theorems/THM-5313-switching-derivatives-of-hamiltonian-paths-and-cycles.md) · [Prime radius balance and the forest coefficients](../moore_fixed_order_spectra_20261124/RADIUS-TWO-PRIME-BALANCE.md) On \(n\) sites fix the disjoint directed cycles and toss a fair independent coin on every remaining pair. Write \(F\) for the union of the prescribed undirected cycle edges. Then
\[
 {\mathbb E}H=2^{1-n}\sum_{J\subseteq F}(n-|J|)!\,2^{c(J)},
\]
where the sum uses precisely those \(J\) whose nonempty components are paths of **even edge length**, and \(c(J)\) counts these components. A permutation contributes the product, along its consecutive pairs, of \(1+\varepsilon\) for prescribed edges and 1 otherwise, times \(2^{1-n}\); \(\varepsilon\) is the traversal sign relative to the prescribed cyclic direction. After expansion, an edge set contributing to a linear order must be a linear forest. Contract its components. Each can be read in its two directions, with sign sum \(1+(-1)^{\text{length}}\), giving the formula. Reversing one fixed edge \(e\) changes only the terms with \(e\in J\), negating each of them. Thus the difference of the conditional means, coherent orientation minus reversed, is exactly twice their original contribution. This is an expectation over completions. An exact one-arc difference in a specified tournament needs its actual exterior orientations.

For example, if the Hoffman–Singleton edges receive independent fair directions, the *expected* packing polynomial of consistently directed induced 13-cycles is the first displayed polynomial at \(2^{-12}z\). Each 13-cycle has two coherent orientations out of \(2^{13}\), and the events are independent for disjoint cycles. The tournament count also includes the other odd directed cycles and their individual compatibilities.
