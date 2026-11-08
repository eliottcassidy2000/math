# Exact physical membership on deletion levels one and three

Status: PROVED by the structural and normalization arguments below; independent audit requested. These are actual integer-state membership criteria on two fixed deletion levels, not a claim about arbitrary physical states or Chan–Pak completeness. They supplement PROOF.md and its mobile-triangle construction.

## 1. Small Hamilton counts classify the deletion graph

Hamilton-path count is monotone under induced vertex inclusion: extend each old Hamilton path by inserting the new vertex at a deterministically selected legal position; deletion recovers the original path, so the extension is injective. A legal insertion always exists by the usual first/last or consecutive in-to-out transition argument.

A tournament has H=1 iff it is transitive. Indeed every nontransitive tournament contains a directed triangle, whose H=3, while a transitive tournament has a unique Hamilton path.

A tournament has H=3 iff it is an ordered sum

    transitive prefix + directed triangle + transitive suffix,

where either transitive part may be empty. Here ordered sum means all edges point left to right.

Proof: choose a directed triangle C. Any vertex with mixed adjacency to C creates a four-vertex tournament with H=5. This last finite fact follows by listing Hamilton permutations: a root joined to C with k outneighbors has H=3 for k=0,3 and H=5 for k=1,2. Therefore every external vertex belongs either to S, beating all of C, or T, losing to all of C. Neither S nor T contains a directed triangle: such a triangle together with C would be an ordered sum of two triangles, with H=9. Thus S and T are transitive. Finally there cannot be t→s with s∈S,t∈T: together with two triangle vertices c_1→c_2, the edges are s→c_i→t→s, a four-tournament with H=5. Monotonicity excludes each H=5 or9 witness. Thus S→T and the asserted ordered-sum decomposition follows. Conversely such a sum has H=1·3·1=3. No literature-priority claim is made for this elementary classification.

## 2. Deletion one is rigid

If x(A,a) has d=1, the deletion tournament is transitive. In its total order, encode a's outgoing neighborhood by a binary word b. The same path-cover decomposition as PROOF.md gives

    U=1+binary_value(b).

Leading zeros correspond to vertices that are global sources of A. Adding or deleting a global source while retaining an existing root preserves every coordinate (D,E,U,d): each Hamilton path starts at the source; for U, the root-starting distinguished path cannot contain that source, which therefore starts the other path. This is a direct bijection, not a loss of information assumed from a matrix.

Consequently for U≥2 the word bin(U−1), with no leading zeros, supplies the unique full rooted state at deletion level1. For U=1 the state is (1,1,1,1), represented by the singleton root. The scalar and endpoint coordinates can be computed by the six-state transfer of PROOF.md on the indicated singleton blocks. In particular E is the parity of U, and fixed (E,U,d=1) has at most one D. This is exact physical membership: compare the supplied state with this explicit representative.

## 3. Deletion three has an exact finite canonical bank

For a proposed integer state (D,E,U,3), physical membership is equivalent to the following finite test.

1. Require U≥3 and U divisible by3. Set M=U/3−1.
2. If M=0, accept exactly (D,E,U,d)=(3,3,3,3).
3. If M>0, let b=00 concatenated with bin(M). If L is the bit length of M, b has length L+2. Construct T(b,p) from PROOF.md for every p=0,...,L, including the terminal pair. Accept iff the proposed (D,E) equals one of those representatives' coordinates.

The U and d coordinates of every candidate are already (3(1+M),3). For E, if the last block is an ordinary singleton then E=0 when its bit is1 and E=3 when its bit is0. If the last block is the triangle with k outneighbors, then E=3−k, since each triangle vertex terminates exactly one triangle Hamilton path. D is the explicit six-state transfer value, or can be counted directly on this actual finite graph.

Proof of completeness: by Section1 the deletion graph consists of singleton blocks with exactly one triangle. Encode every singleton by its root-neighborhood bit, and replace the triangle by the two-bit binary expansion of k, the number of its vertices beaten by the root. Rooted triangle types of a given k are isomorphic. This produces a word w and an adjacent-pair position p whose physical tournament is exactly T(w,p). Section2 of PROOF.md gives U=3(1+binary_value(w)).

If M=0 all vertices beat the root, so it is a global sink, and the state is (3,3,3,3). Otherwise consider the triangle's position relative to the first1 of w.

- If its two positions are both before the first1, then k=0. All preceding zero blocks and the triangle, as well as subsequent zeros before the first1, form an ordered prefix of H=3 that beats the root and every later block. The full state is exactly3 times the state of the remaining rooted suffix. Indeed Hamilton and ending-root counts multiply by H of the prefix, deletion counts also multiply, and in a root-starting two-path cover all prefix vertices must belong to P and form one Hamilton path. Thus arbitrary such placements have the same full state as the canonical p=0 triangle in 00+bin(M).
- If the triangle straddles the first1, it occupies the last leading0 and that1, so k=1. Remove earlier leading singleton zeros, which are global sources. Add one leading singleton zero to match the canonical p=1 placement. These operations preserve the full state by the global-source bijection.
- If the triangle begins at or after the first1, all leading zeros are singleton global sources. Remove them and prepend exactly two singleton global sources. This gives exactly one of the canonical p=2,...,L placements, preserving the full state.

These cases exhaust the possibilities. Soundness is immediate because every canonical candidate is an actual tournament.

## 4. Consequences and boundaries

For fixed U and d=3, there are at most bitlength(U/3−1)+1 full states when U>3, and at most one at U=3. This also bounds every scalar fiber with fixed (E,U,d=3). The unbounded mobile-triangle family has U=4^r+2, M=(4^r−1)/3 with bit length2r−1, and 2r−2 distinct scalar values at E=0. Thus the logarithmic upper bound is asymptotically sharp, with only a constant additive gap on that sequence.

By Rédei parity, positive deletion counts are odd. The rigid d=1 case and unbounded d=3 fibers show that3 is the first deletion level permitting unbounded scalar fibers. This statement does not assert that the smallest nontrivial fiber has a particular order.

The canonical bank proves state existence by constructing representatives. It does not recognize whether a given graph is a prescribed graft image; Session227's distinction remains intact. Nor does this fixed-d classification supply an all-odd Hamilton-count construction. Multiple-triangle deletion graphs and scalar coverage remain outside this theorem.
