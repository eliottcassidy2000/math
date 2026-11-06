# Ordered root amalgamation and four three-cycle templates

Session 197. Direct counting identities and physical constructive templates. This note does not by itself classify every three-cycle tournament, and makes no literature-priority claim.

## 1. General identities

For rooted nonempty tournaments A and B, identify their roots to one vertex r and orient every other cross arc from A-r to B-r. Denote this tournament A⊕B. Write h,e,s for total, root-ending and root-starting Hamilton-path counts, d=H(T-r), R for the ordered two-cover count with the reserved path starting at r, and L for its ending analogue. The other path may be empty. Write b(T-r) for the two-labelled-bin path-cover count of T-r, allowing either bin to be empty; H(empty)=1.

Then, with subscripts identifying the two inputs,

    h'=(R_A-d_A)(L_B-d_B)+(h_A-e_A)d_B
        +(h_B-s_B)d_A+e_A s_B,                         (1)
    R'=(R_A-d_A)b(B-r)+d_A R_B,                         (2)
    L'=(L_B-d_B)b(A-r)+d_B L_A,                         (3)
    d'=d_A d_B,  e'=d_A e_B,  s'=s_A d_B.              (4)

The e' formula in (4) assumes B is nontrivial, and the s' formula assumes A is nontrivial. If the respective input is the singleton, use the unchanged other input's endpoint count instead. Equations (1)-(3) and the deletion formula hold without that restriction.

For (1), split a Hamilton word at r into a before path and an after path. In each path all A-r vertices precede all B-r vertices. Partition into the four cases according to whether the before path has a B-r vertex and the after path has an A-r vertex. If both occur, the root's outgoing condition is tested only in A and its incoming condition only in B, yielding (R_A-d_A)(L_B-d_B). If only the latter occurs, all B-r vertices lie after r and the A-word does not end at r, yielding (h_A-e_A)d_B. If only the former occurs, the symmetric contribution is (h_B-s_B)d_A. If neither occurs, the word is an A-word ending at r followed by a B-word starting at r, yielding e_A s_B. All choices reconstruct one actual word, so no cancellation or multiplicity assumption is involved.

For (2), partition a two-cover according to whether the reserved r-starting path uses any A-r vertex. If it does, the root incidence is determined in A, and the B-r vertices can occupy its two labelled path slots arbitrarily. Otherwise all A-r vertices occupy the unreserved path, and the B-part is an arbitrary reserved-root two-cover. Formula (3) is the converse. Formula (4) follows from the uniform cross arcs. If either input is the singleton, the corresponding operation is the identity; this gives the endpoint exceptions just stated.

## 2. Four actual C=3 families

Let V and W be nonempty rooted tournaments with exactly one root-starting Hamilton path. Their four-coordinate profiles are v=(D_1,E_1,U_1,d_1) and w=(D_2,E_2,U_2,d_2), where U denotes R. They may be weak tournaments. Use W's converse, so its marked root has exactly one ending path.

Take a directed triangle K with vertices x,y,z in cyclic order. Add the shared root r, with r beating exactly j of its vertices, for j=0,1,2,3. For j=1 use {x}, and for j=2 use {x,y}; each size is a single orbit under the triangle's rotations. Call this rooted four-vertex graph K_j. Form

    T_j(V,W)=V ⊕ K_j ⊕ converse(W).                    (5)

Every Hamilton cycle, started at r, must traverse V-r first, then the triangle, then W-r, and return to r. It cannot skip a nonempty block and later return to it, because all cross arcs go forward. The V segment is its unique r-starting Hamilton path with r removed; the W segment is its unique r-ending path with r removed. The triangle has three Hamilton paths. Consequently (5) is strong and has exactly three Hamilton cycles, for every j. Its order is |V|+|W|+2. Here both sources are initially assumed nontrivial. There are also precisely two allowed singleton-source boundaries: W may be singleton when j=0 and V is nontrivial, or V may be singleton when j=3 and W is nontrivial. In the former case all triangle vertices point to r; in the latter r points to all of them. The same three cycles exist. If just W is singleton, the cycle count is 3-j; if just V is singleton, it is j. With both singleton the result is K_j and never has three cycles. In all allowed cases, the construction also gives s_r(T_j)=3d_2, e_r(T_j)=3d_1 and H(T_j-r)=3d_1d_2.

The exact constants (h,s,e,R,L,d,b(K_j-r)) are, in order j=0,1,2,3,

    (3,0,3,3,12,3,12), (5,1,2,6,9,3,12),
    (5,2,1,9,6,3,12), (3,3,0,12,3,3,12).             (6)

They follow by listing paths in the triangle and its proper subsets. Applying (1)-(4) twice gives H(T_j)=v M_j w^t, with the coordinate order (D,E,U,d) and matrices

    M_0 = [[0,0,0,3], [0,0,0,-3], [0,0,12,-3], [3,0,-12,3]],
    M_1 = [[0,0,0,3], [0,0,0,-2], [0,0,12,-6], [3,-1,-9,5]],
    M_2 = [[0,0,0,3], [0,0,0,-1], [0,0,12,-9], [3,-2,-6,5]],
    M_3 = [[0,0,0,3], [0,0,0,0],  [0,0,12,-12],[3,-3,-3,3]].

Although these coordinate matrices have negative entries, (1) is an explicit sum of disjoint nonnegative physical cases. In particular no arbitrary formal quadruple is asserted realizable. The transpose identity M_j^t=M_(3-j) agrees with the literal converse map T_j(V,W)^op=T_(3-j)(W,V), up to the triangle's cyclic relabelling.

## 3. The canonical triangle-tail branch

This explains a specific part of the complete grammar in the separately reviewed START_THREE.md, Section 2, Type II. Such a rooted start-three input has an initial rooted start-one component U after root deletion, a singleton tail P, one directed triangle, and a singleton tail Q. The root r has arbitrary incidences to P,Q and one of the four triangle incidence patterns, subject to the original strongness condition.

Suppose this input is followed by any ordinary one-entry suffix and the degree-one closure producing a strong three-cycle graph. Bundle that suffix and the closing vertex as a tournament G: all earlier vertices point to G except its first entry g points to r. The bundle is nonempty; when there are no ordinary suffix steps it is just the singleton closing vertex. The usual cycle-cut bijection gives C(final)=3e_g(G), so e_g(G)=1 in either case.

Now put V=r+U+P and put W^op=r+Q+G. Both contain at least two vertices. The root-start count of V is one: a path must enter U at its marked vertex, traverse it and then P. The root-ending count of W^op is one: a path ending at r must traverse Q, then all G and finish at g before r. Ending at r directly from Q would strand G, since no G vertex points back into Q. Thus the complete resulting graph is precisely one of (5). This proves a two-source containment for this triangle-tail branch. It does not establish uniqueness of its presentations, identify overlaps with the Type-I mixed/triple-entry families, or classify the converse-only C3 graphs. Those questions belong to the broader C3 factorization task.

The direct degree-one closure from START_THREE Section 5, with a terminal triangle in its deletion tail, also fits: put the degree-one root, the initial start-one component and all pretriangle singleton components in V, and take W singleton and j=0. If singleton components follow the triangle, include them with the root as W^op instead; every such vertex points to the degree-one root, so its ending count is one. Taking converses gives the opposite boundary.

## 4. Finite controls

The standalone verifier enumerates physical directed words and their proper-subset covers. It checks (1)-(4) on all rooted labelled input pairs of orders at most three, the four literal K_j profiles, and all (5) pairs from rooted labelled start-one sources of orders at most three, subject to the stated singleton restrictions. It also checks the singleton boundaries and the literal converse map. These finite controls test the identities; the proofs above establish their arbitrary-order scope.
