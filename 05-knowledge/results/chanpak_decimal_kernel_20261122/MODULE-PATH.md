# Which two leaves can support a flip out of a multiple of three?

**PROVED.** Fix any presentation of a tournament by iterated module substitution. Its leaves are the actual vertices. At each internal node \(v\), retain its oriented quotient and the parity of the size of each child. Replace each odd child by one vertex and each even child by a transitive pair. Let \(w_v\) be the Hamilton path count of this small tournament, modulo three. Then
\[
                      H(T)\equiv\prod_v w_v\pmod3 .
\]
If \(3\mid H(T)\) and reversing \(xy\) makes its count nonzero modulo three, **every node with \(w_v=0\) lies on the undirected tree path from leaf \(x\) to leaf \(y\)**. This test applies to every chosen substitution tree, including trees with nonprime quotients. A passing path is only a candidate.

For a nonempty module \(A\), its ordered-slot polynomial on \(\mathbb F_3\) equals \(H(A)j^e\), with \(e=1\) in odd order and \(e=2\) in even order. This follows from its value at 0, its value at 1 and its fixed parity. The monomial \(j^e\) is the slot polynomial of a transitive module of size \(e\). The path count after substitution is linear in the module's ordered-cover readings, so substituting one child at a time proves the product formula. [Finite-field context theorem](FINITE-FIELDS.md)

Let \(a\) be the least common ancestor of the two endpoints. Every strict ancestor of \(a\) keeps its quotient and all child sizes after the reversal, even if the tournament inside one child changes. A zero there still kills the product. A node \(v\) disjoint from both endpoint leaves still defines a module in the edited tournament. Its internal tree is unchanged. A zero at \(v\) makes that entire module divisible by three, so every host containing it also has count divisible by three. Every off-path node belongs to one of these two cases. The proof only uses decompositions that still exist after the flip.

There is also an integer bound. Given \(m\) pairwise disjoint modules of heights divisible by three, every term in the simultaneous skeleton count contains one three-divisible factor per module. Indeed the ordered nonempty run counts \(B_r(A)\) are the forward differences of the slot polynomial at zero, and all the integer evaluations of that polynomial are zero modulo three. Thus \(3^m\mid H(T)\). For odd \(p\ge3\),
\[
 v_3(E_p-1)=v_3(2^p+1)+v_3(2^{p-1}-1)=2+v_3(p(p-1)).
\]
The last step uses lifting the exponent at \(2+1\) and \(4-1\). At 61 there can be at most three pairwise disjoint such modules. If there are three, any two endpoints miss one and the exact exit is impossible. The tree-path test also acts on nested zeros.

## A distinction visible on five vertices

In \(T=C_3[C_3,1,1]\), the inner triangle is the only proper module whose height is divisible by three. Flip an arc inside it. Its own height drops from 3 to 1. The full height drops from 15 to 9: both before and after the edit the outer quotient has three odd children and weight zero. Thus an edit that meets and repairs every bad proper module can still be excluded by the path test.

## Several edits and the exact cost of covering the tree

**PROVED.** For any set of edits whose final count is nonzero modulo three, the union of their endpoint paths contains all the original zero-weight nodes. Indeed if a node is on none of the paths, no edited edge crosses the boundary of its module or joins two distinct children. The module survives in the final host. Its children still form a substitution presentation inside it, with the same quotient and child parities. A zero weight there makes the whole final count divisible by three.

For two or more zero nodes, take their smallest connecting subtree, and let \(\ell\) be its number of leaves. At least \(\lceil\ell/2\rceil\) arcs must be edited. Each leaf of this subtree is a required zero node, and a simple path can cover at most two of them. In fact this number is the exact cost of the path-cover condition alone. Put the leaves in cyclic order around a plane drawing of the subtree and, if necessary, repeat one leaf beside itself to make the length even. Match opposite positions in this list. Each tree edge cuts the list into two nonempty consecutive intervals. A proper consecutive interval cannot be invariant under the half-turn, so at least one matched path crosses that edge. Extend each terminal endpoint to a vertex leaf of the original substitution tree on its outside side. A single zero node requires one covering path. Realizing a nonzero residue after the chosen edits is an additional arithmetic condition.

**A strong family attains the edit bound.** For \(m\ge1\), let \(U\) be the ordinal sum of one singleton and \(2m\) cyclic triangles. Take \(T=C_3[U,1,1]\), of order \(6m+3\). Its outer weight is zero since \(U\) has odd order. The \(2m\) inner triangles also have zero weight. In the presentation using one ordinal node for \(U\), their connecting subtree has \(2m+1\) leaves. Thus reaching any count coprime to three costs at least \(m+1\) reversals. Pair adjacent triangles. Reversing one interblock arc in \(C_3\to C_3\) changes its count from 9 to 17 by direct six-vertex path enumeration. Use one such edit in each pair, then flip the edge between the two exterior singletons. The triangle wrapper becomes an ordinal arrangement and the final count is exactly \(17^m\). The minimum is therefore **\(m+1\)**. This distance is measured in the full tournament graph on the same vertices and allows any final orientation. The starting count is \(3^{2m}+2\cdot12^{2m}\).

## Applying the path condition at the full exponent-61 height

The [concurrent 91- and 93-vertex lower families](../chanpak_e61_compact_volchkov_20261122/STRONG-93.md) have an especially short obstruction. Both are \(C_3[U,1,1]\) with \(U\) odd, so the root weight is zero. Inside \(U\) there are two disjoint three-divisible components: two directed triangles at order 93, and the components of heights 3,963 and 1,599 at order 91. To meet both components a single arc must lie wholly inside \(U\). The root zero then survives. Equivalently the two components and the outer root are three different terminal directions in the connecting subtree. **Every one-flip neighbor of all 3,237,399,360 and 237,918,897,216,000 lower classes, respectively, has height divisible by three.** This gives the entire residue obstruction throughout both families for every permutation of the components. Replacing small components by other classes with identical slot counts enlarges the two families by a factor of 128 each. The mod-three argument applies unchanged, so it covers **414,387,118,080 classes at order 91** and **30,453,618,843,648,000 at order 93**. [Full slot-profile fibers](../chanpak_e61_fibonacci_target_20261122/PROFILE-FIBERS.md) At order 90 the concurrent construction has a different quotient, with weight one; this particular three-direction argument does not apply.

Flatten the six-layer product grammar into one quotient with seven singleton roots and six homogeneous arms. Replacing each arm by the transitive tournament of size one or two according to its parity gives scalar \(\prod_i(1+2^{e_i})\) by the product formula for the rooted construction. The quotient weight therefore vanishes whenever one arm is odd. The 61-arm is odd in all **1,382,400** classes. Immediately the path test excludes **every flip inside a single arm** throughout that family. The four-layer theorem further excludes the cross-arm possibilities. [Exact product fiber and four-layer proof](../chanpak_e61_fibonacci_target_20261122/NEW-LOWER.md)

For an example where the edited arm's zero actually disappears, take the factors \(861,3963,4681,8275,8723,2^{61}+1\) in that word. Use the arms of the product note and the additional rows
\[
 861:(2,8,35,69,15,27,55),\qquad
 8723:(26,28,145,100,104,71,135,59).
\]
The arm heights are \(49,219,257,563,595,1\) in order. Reversing the edge \((4,6)\) in the 3963-arm changes its height from 219 to 215. The full 107-vertex lower host still has a zero quotient strictly above that edit. The example can be replayed by evaluating the at-most-eight-vertex arms separately and multiplying their triangle-return counts. [Integer replay](module_path.py)
