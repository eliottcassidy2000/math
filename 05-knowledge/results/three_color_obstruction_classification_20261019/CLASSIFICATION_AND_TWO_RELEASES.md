# Complete three-color obstruction classification and two-release structure

19 October 2026. Status: proved. The classification is analytic; no finite search is used in its proof. Exact replay checks the explicit arithmetic example and the binary matrix calculation.

Dependency: [the tight-kernel and one-release dichotomy](../graceful_three_color_release_20261019/THREE_COLOR_RELEASE_DICHOTOMY.md).

## 1. Complete classification in a representable rank-three matroid

Let M be a finite rank-three matroid representable over a field, with three color classes A_1,A_2,A_3 and at least one independent transversal. The state graph joins independent transversals differing in at most two coordinates. Distinct actual elements, including parallel elements, remain distinct choices. Loops never occur in states and may be deleted from the classes.

### Theorem 1

The state graph is disconnected if and only if there are three independent rank-one directions P_1,P_2,P_3 such that, after permuting color names and direction names, the nonloop color supports are exactly

    A_1: P_1 and P_2,
    A_2: P_2 and P_3,
    A_3: P_3 and P_1,

with at least one actual element in each indicated color-direction pair.

In that case there are exactly two components. Each consists of all independent choices of actual realizations of one of the two direction assignments. Thus no additional type of three-color obstruction, and no three-color example with more than two components, occurs in this class.

### Proof: first reduce the color ranks

A color of rank one cannot cause disconnection: change its current element to its target parallel element, then change the other two coordinates together to their target values. Both endpoints of each move remain bases.

A color of rank three supplies a full-rank buffer. Consider the independent pairs chosen from the other two classes. Two such pairs with the same closure can be changed directly by a closure-preserving two-coordinate move. If their closures differ, some target element lies outside the current pair's closure; replacing its corresponding current coordinate is valid, and the remaining coordinate can then be changed to its target. This uses at most two singles. Each single lifts to a full basis by simultaneously replacing the buffer if necessary; a closure-preserving double retains the buffer. Finally correct the buffer alone. Hence this case is connected too.

Therefore in a disconnected system every color has rank exactly two.

### Proof: the zero-minor argument

Choose states S=(s_1,s_2,s_3) and T=(t_1,t_2,t_3) from distinct components. Every mixed tuple choosing some but not all coordinates from T and the others from S must be dependent. Otherwise it would share a coordinate with each endpoint and be adjacent to both, connecting the components.

Use S as a vector-space basis, and let Y be the 3-by-3 matrix whose columns are the vectors t_i. It is invertible. The invalid one-coordinate replacements give

    y_11=y_22=y_33=0.

The invalid two-coordinate replacements give zero principal two-by-two minors, hence

    y_ij y_ji=0 for each distinct i,j.

Only the two directed-three-cycle products remain in the determinant:

    det(Y)=y_12 y_23 y_31 + y_13 y_21 y_32.

At least one product is nonzero. If the first is nonzero, the pairwise-product conditions force y_21=y_32=y_13=0; the other case is symmetric. Thus Y is a monomial matrix whose permutation is a three-cycle. Equivalently, every t_i is parallel to a different s_j, and these correspondences form a cyclic permutation.

### Proof: excluding additional directions

Color A_i has rank two and contains the two nonparallel directions s_i,t_i. Any nonloop element z of that color outside these two directions has a representation

    z=alpha s_i + beta t_i, with alpha and beta both nonzero.

Replacing s_i by z in S gives a basis, since the s_i coefficient is nonzero. Replacing t_i by z in T also gives a basis, since the t_i coefficient is nonzero. These two new states share their i-th coordinate z and therefore differ in at most two coordinates. Together with the single changes at their ends they connect S and T, a contradiction.

So each color contains exactly its two cyclic directions. An independent transversal must occupy all three directions once. The cyclic incidence pattern has exactly two such assignments. Within either assignment, every parallel realization can be changed independently; between assignments all three color directions differ. This proves necessity, sufficiency, and the exact count of components. QED.

The proof works over any field, not just the binary field. No claim for arbitrary nonrepresentable matroids is made here.

### Corollary 1.1: graphic and arithmetic fibers

After contracting a fixed forest in a three-free-color completion problem, there are four vertices. Graphic matroids are representable, and three independent directions are the three edges of a spanning tree J. Parallel directions mean precisely the same unordered pair of contracted vertices. Consequently every disconnected three-free-color graphic fiber is exactly a tree-supported chordless-six-cycle matching obstruction: each free color supports two edges of J, omitting a different one.

The earlier arithmetic witnesses on 13 and 15 labels therefore realize the complete obstruction type, not just one possible family. This does not classify fibers with four or more free colors.

### Corollary 1.2: necessary ranks for disconnection

For every disconnected three-color fiber,

    q_F(d)=2 for every free color d,
    q_F({d,e})=3 for every pair of free colors,
    q_F(C)=3.

The pair statement also has a direct proof: if a pair had rank two, its selected edges would always span a fixed rank-two flat; replace both coordinates together, then adjust the remaining coordinate. That would connect the fiber.

## 2. Releasing bridge edges gives an exact product

This statement allows any number r of old free colors and any move width w>=1. Let F be fixed, assume the original fiber is nonempty, and put

    Q=F union all edges of the old free colors.

Thus Q is connected. Let B be any b fixed edges that are bridges of Q, and set G=F minus B. The graph Q minus B has b+1 connected components.

### Theorem 2: bridge-only release product

The enlarged completion state set after releasing the colors of B is a direct product of:

1. the original old-color completion choices;
2. the rainbow spanning-tree choices of the released colors on the b+1 components of Q minus B, with every actual crossing realization retained.

Under moves changing at most w colors, its component set is the Cartesian product of the component sets of these two factors.

Proof. In the graphic matroid after contracting G, the union of the old-color edges has rank

    (m-b) - (m-r-b) = r,

where m+1 is the original number of vertices. Let H be its closure. Every independent old-color selection has r edges and is a basis of H. The original bridge set B is a basis of the quotient by H, so such an old selection is valid precisely when it was a valid original completion choice.

A released-color selection completes any old basis exactly when it is a basis in that same quotient. This proves the state-set product. Any allowed full move projects to allowed moves of width at most w in both factors. Conversely, paths in either factor lift with the other factor fixed. Hence the component product follows. QED.

In particular, no collection consisting only of bridges of the original Q can repair a disconnected old fiber. Even several full-rank released color classes do not remove the old tight-flat obstruction. At least one released edge must be a nonbridge of Q. This is a necessary condition, not a general sufficient one.

## 3. A proved two-release sufficient condition

Return to a disconnected arithmetic fiber with exactly three old free colors C. Color 1 is necessarily fixed: if it were free, its full difference graph would be connected and have contracted rank three, contradicting Corollary 1.2. Let its fixed edge be e_1.

If e_1 is not a bridge of Q, the preceding checkpoint already proves that release 1 alone repairs the fiber. Suppose instead e_1 is a bridge, and let U,V be the two components of Q minus {e_1}.

### Theorem 3: one nonbridge whose color crosses the failed cut

Suppose another fixed edge f of color b satisfies both conditions:

- f is not a bridge of Q;
- the full difference-b class contains at least one edge crossing U,V.

Then releasing the two fixed edges e_1 and f makes the entire enlarged fiber two-color connected. Its diameter is at most 14. The theorem does not require that releasing f alone works.

Proof. Put G=F minus {e_1,f}. Its contracted matroid has rank five. Color 1 supplies a full-rank buffer. Treat the four colors C union {b} as the old classes for the tight-kernel theorem.

Since a bridge belongs to no cycle, f remains a nonbridge after e_1 is removed. Thus Q minus {e_1,f} still has exactly two connected components, and

    q_G(C)=4.

For any pair J of old free colors, Corollary 1.2 gives q_F(J)=3. Therefore F union E_J is connected. Removing e_1 splits it into exactly U,V: it is a connected spanning subgraph of Q containing that bridge. The graph G union E_J union E_b contains this entire graph minus e_1, because E_b restores f. By hypothesis E_b also has an edge crossing U,V. Consequently

    q_G(J union {b})=5.

All triples of the four nonbuffer colors have rank at least four, and the full set has rank five. Hence no subset of size at least three is tight. The maximal tight kernel has at most two colors, and the connected-buffer theorem proves the result. The same path estimate with four nonbuffer colors gives at most 14 moves. QED.

A second elementary restriction is useful for searching for an existential one-release counterexample. Color 2 must also be fixed: the fixed color-1 edge joins the two parity components of E_2, so a free color 2 would have full contracted rank three. Moreover, releasing the fixed color-2 edge always gives a connected color-2 buffer, because the color-1 edge remains fixed. Thus any disconnected three-color example in which every single-edge release fails must have both fixed edges of colors 1 and 2 as bridges of Q.

The next example shows that this necessary bridge pattern is arithmetically possible. It does not show that every single-edge release fails.

## 4. A 17-label example: two full-rank buffers still fail

Use labels 0,...,16, free colors 10,13,15, and fixed edges

    1: 7--8      2: 7--9       3: 3--6       4: 7--11
    5: 10--15    6: 1--7       7: 5--12      8: 2--10
    9: 3--12     11: 3--14     12: 4--16     14: 0--14
    16: 0--16.

The fixed components are

    A={0,3,4,5,6,12,14,16}, B0={1,7,8,9,11},
    C0={2,10,15}, D={13}.

Their nonloop free supports are

    10: AC0 (3 realizations), AD (1)
    13: AB0 (1 realization), AD (1)
    15: AB0 (1 realization), AC0 (1).

The contact tree is the three-leaf star centered at A. The two type assignments give original components of sizes 1 and 3.

Vertices 8 and 9 are incident to no old-free-color edge. They are leaves of Q, attached respectively by 7--8 and 7--9. Thus both edges of colors 1 and 2 are bridges.

Exact outcomes are:

    no release: 4 states, components 1 and 3;
    release 1: 8 states, components 2 and 6;
    release 2: 8 states, components 2 and 6;
    release 1 and 2: 32 states, components 8 and 24;
    release 3: 21 states, one component;
    release 1 and 3: 86 states, one component.

After releasing 1 and 2, there are six contracted vertices and both classes 1 and 2 have full rank five. Color 1 is connected by its path. For color 2, its parity components are joined by the still-fixed odd edge 3--6. Nevertheless the three old classes span a fixed rank-three flat, so Theorem 2 preserves their obstruction.

The multiplicative count can be checked without enumeration. After deleting the two bridges, the quotient sides are {8}, {9}, and the rest O. Color 1 has one realization of each of O--8, 8--9, 9--O. Color 2 has two realizations each of O--8 and O--9. Choosing two edges of different colors that connect the quotient gives 2+4+2=8 buffer states. Because this factor has only two coordinates, its two-color graph is complete. Multiplication by the original component sizes gives 8 and 24.

The positive two-release certificate applies to colors 1 and 3: the fixed color-3 edge is a nonbridge of Q, and E_3 crosses the failed color-1 cut {8}, for example by 5--8. The exact census also shows that color 3 alone already suffices.

## 5. Remaining existential gap

The classification and the bridge product remove ambiguity about the minimal three-color obstruction. They do not settle either of the following existential statements for integer-difference fibers:

- every disconnected nonempty three-free-color fiber has some successful single-edge release;
- every such fiber has some successful release of at most two fixed edges.

The 17-label example disproves only the blanket strategy of releasing colors 1 and 2, or of choosing two connected buffers without checking the old tight flat. A counterexample to the existential statements would need to defeat all permitted choices, including nonbridge edges. Conversely, the two-release theorem reduces one positive route to finding a nonbridge fixed edge whose color crosses the failed color-1 cut. No assertion that arithmetic always supplies that edge is made here.
