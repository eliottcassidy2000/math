# Unbounded physical scalar fibers from a mobile triangle

Status: PROVED by the combinatorial transfer and recurrence argument below; independent audit requested. Chan–Pak completeness remains OPEN. No scalar coverage cutoff or literature-priority claim is made.

## 1. Inheritance and the new physical move

Use the inherited rooted state x=(D,E,U,d), where D counts Hamilton paths, E those ending at the root, d counts paths after deleting the root, and U counts ordered spanning path pairs (P,R) with P possibly empty and R nonempty starting at the root. The Session227 adaptive-arithmetic packet proves a two-value fiber bound for its particular four-word menu and a physical D-column obstruction for a common graft-containing continuation. Neither result applies to branch-dependent placement of a module in a graph.

Here the deletion graph is an ordered sum of singleton blocks and one directed triangle. The root's neighborhood changes with the block placement. This is an actual graph construction, not a cone or lattice prescription. The live concepts are physical fibers, marked modules, binary positional values, and ordered two-path covers. The hostile is falsely iterating a two-point switch; the retained sidecar is the position of the triangle in the root-deletion graph.

For a binary word b_0...b_(m−1), choose p with 0≤p≤m−2. Replace positions p,p+1 by a directed triangle, with every edge between distinct blocks directed from earlier to later. Add a root a. For an ordinary singleton at position i, direct a to it iff b_i=1. Into the triangle, direct a to exactly k=2b_p+b_(p+1) vertices. A cyclic triangle's rotation group is transitive on subsets of each cardinality, so the choice of which k vertices does not affect the rooted isomorphism type.

Call the result T(b,p). It has m+2 vertices. Its root-deletion graph has H=3, hence d=3.

## 2. Sidecar identity for every word

For every word and allowed position,

    U(T(b,p),a)=3(1+sum_i b_i 2^(m−1−i)).

If the last position is an unreplaced singleton and b_(m−1)=1, then E=0.

Proof of U: remove the initial root from the distinguished path R. The remainder R' and P are two ordered paths covering the deletion graph, with R' empty or its first vertex an outneighbor of a. In an ordered sum, each path's intersection with a block is an induced path, and nonempty block pieces concatenate uniquely. A singleton block has H=1 and ordered two-bin count Q=2; a cyclic triangle has H=3 and Q=12. If the first nonempty R' block is a singleton, its root-permitted contribution is b_i; if it is a triangle, the contribution is 3k. The latter follows because for each possible starting vertex, its nonempty path of length 1, 2, or 3 contributes exactly one complementary path, giving 3. Earlier blocks go wholly into P, while later blocks have arbitrary ordered two-bin covers.

After dividing by the product H=3, a singleton contributes its bit times the product of subsequent Q/H factors. A triangle contributes k times that product. Its ratio Q/H=4 equals that of the two replaced singleton positions, and k=2b_p+b_(p+1) equals their combined binary value. The empty R' adds 1. This proves the formula. Finally a path ending at a must enter it from the last deletion block, since a Hamilton path of an ordered sum finishes there. A final singleton beaten by a cannot do this, proving E=0.

More generally t disjoint adjacent-pair replacements give d=3^t and U=3^t(1+binary value). An unreplaced terminal 1 still gives E=0. This generalization asserts no distinctness or coverage of D.

## 3. Explicit unbounded fiber theorem

For each integer r≥2 let b=(01)^r and allow p=0,...,2r−3, so that the final 1 remains an ordinary singleton. Set

    B_r=(4^r+6r−1)/3.

Then

    x(T(b,0),a)=(B_r,0,4^r+2,3),
    x(T(b,p),a)=(B_r+(-2)^p−2,0,4^r+2,3)  (1≤p≤2r−3).

These 2r−2 scalar coordinates are pairwise distinct: the powers (-2)^p are distinct and none for p≥1 equals 2. Thus actual physical scalar fibers have unbounded cardinality, even with fixed deletion count d=3. This conclusion is not merely an arbitrarily large number of graphs with the same full state.

For r=4, all six rooted graphs have sidecar (0,258,3), and the scalar set is {59,83,89,93,95,107}. The triangle location is a branch-dependent construction choice, not a shared suffix preserving an old scalar difference. The Session227 common-context theorem is therefore respected.

## 4. Exact finite-dimensional proof of the scalar formula

Delete the root from a Hamilton path. This gives an ordered pair (P,R') of paths covering the deletion graph, where P must be empty or end at an inneighbor of a, and R' must be empty or start at an outneighbor of a. Conversely every such pair gives exactly one Hamilton path by concatenating P,a,R'.

Process deletion blocks in order. Track P status 0=empty, 1=last vertex is an inneighbor, 2=last vertex is an outneighbor. Track R' status 0=empty, 1=nonempty with its first vertex an outneighbor. Discard a first R' piece with a bad first vertex. Later pieces do not change the first-vertex test; a new nonempty P piece replaces its last-vertex status. The six states in order are (0,0),(0,1),(1,0),(1,1),(2,0),(2,1). Column vectors are used, with initial v=(1,0,0,0,0,0)^t and terminal functional ell=(1,1,1,1,0,0).

The exact local matrices are:

    S0 = [0 0 0 0 0 0; 0 1 0 0 0 0; 1 0 1 0 1 0;
          0 1 0 2 0 1; 0 0 0 0 0 0; 0 0 0 0 0 1],
    S1 = [0 0 0 0 0 0; 1 1 0 0 0 0; 0 0 0 0 0 0;
          0 0 1 1 0 0; 1 0 1 0 1 0; 0 1 0 1 1 2],
    T1 = [0 0 0 0 0 0; 1 3 0 0 0 0; 2 0 2 0 2 0;
          2 6 3 9 2 6; 1 0 1 0 1 0; 0 3 0 3 1 6],
    T2 = [0 0 0 0 0 0; 2 3 0 0 0 0; 1 0 1 0 1 0;
          2 3 4 6 2 3; 2 0 2 0 2 0; 2 6 2 6 4 9].

Here Sb is a singleton of bit b and Tk a triangle with k outneighbors. Each entry can be computed by the local path-pair rule above, requiring only the six ordered permutations and subsets of a triangle. The supplied checker derives these entries from that rule rather than trusting a table.

Write A=S1 S0. Direct integer multiplication proves

    A(A−4I)(A−I)^2=0,

or A^4−6A^3+9A^2−4A=0. This is a finite matrix identity, not an extrapolation from path counts.

For j,k≥0 put r=j+k+2. For even position p=2j the scalar is

    F(j,k)=ell A^(k+1) T1 A^j v.

For odd position p=2j+1 it is

    G(j,k)=ell A^k S1 T2 S0 A^j v.

Both satisfy the order-four recurrence with characteristic polynomial z(z−4)(z−1)^2 separately in j and k. Their proposed closed forms are

    F(j,k)=B_(j+k+2)+4^j−2+1_(j=0),
    G(j,k)=B_(j+k+2)−2·4^j−2.

These also satisfy that recurrence separately in both variables: constants and affine terms are killed by (z−1)^2, powers of 4 by z−4, and the boundary indicator by z. Equality is therefore proved by the following 4×4 initial grids, rows j=0,...,3 and columns k=0,...,3:

    F: 9   27   93   351
       29  95   353  1379
       107 365  1391 5489
       413 1439 5537 21923

    G: 5   23   89   347
       17  83   341  1367
       59  317  1343 5441
       221 1247 5345 21731

For completeness, first use the recurrence in j to extend each of the four initial columns to all j, and then use the recurrence in k for each such j. Its leading coefficient is 1, so no inversion or divisibility assumption is needed. This proves the formulas at every order. The j=0 boundary indicator is essential; dropping it gives an incorrect D at p=0.

## 5. Scope and next question

This gives actual many-choice sidecar collisions, unlike the two-value word-menu return. It does not give an arithmetic interval, a recursively closed return, or a smaller predecessor for each prescribed odd target. The fiber's spacing is essentially geometric. The source graphs need not be strong: the construction deliberately retains a final singleton sink. Existing O/Z contexts preserve a pure-D fiber, but a common graft does not, by Session227. A useful next question is whether multiple independently mobile triangles can create controlled digit sums or dense scalar fibers without destroying their common sidecar.
