# A complete two-source reduction for three Hamilton cycles

Session 197. All-order physical source reduction and an exact finite inverse for the strong C=3 Hamilton-path spectrum. No claim of universal odd coverage, unique presentation, or a returning port interface is made.

## 1. Sources, conventions, and inherited structural input

H(T) counts Hamilton vertex words; C(T) counts directed Hamilton cycles modulo rotation. H(empty)=1. For a marked vertex c, write s_c and e_c for Hamilton words starting and ending at c. R_c counts ordered path covers (P,Q) of T, with Q nonempty and starting at c and P possibly empty. L_c is the ending analogue. Set d_c=H(T-c). A source is any nonempty rooted tournament (V,c) with s_c(V)=1, and its profile is

    v=(D,E,U,d)=(H(V),e_c(V),R_c(V),H(V-c)).

The singleton has profile o=(1,1,1,1). Write S for the actual source-profile set and S* for the nontrivial-source subset. Profiles are not graph isomorphism classes. Every construction below uses actual rooted sources, not arbitrary integer vectors.

We use the following independently established structural theorems, stated here to make the dependencies and reduction exact. Frozen proof copies are in dependencies/.

1. Every strong C=3 tournament has an outdegree-one or indegree-one vertex (C3_COMPLETE.md). This is an all-order theorem proved by a justified finite cutoff, not census stabilization.
2. Delete the root r of a strong rooted s_r=3 tournament. Its initial SCC U and later SCCs have one of two forms (START_THREE.md, Section 2):
   - Type I: all later SCCs are singletons, and the root enters U at either one start-three vertex, a start-one and a start-two vertex in a C1 graph, or three start-one vertices in a C1 graph.
   - Type II: the root enters U at exactly one start-one vertex, and the later SCCs contain exactly one directed triangle and otherwise singletons. U may itself be a singleton.
   The one-entry Type-I case is recursive. Root incidences to the tail are arbitrary subject to strongness. The empty-tail whole-triangle Type-I source case is excluded.
3. A strong C1 graph's start-one vertices form a consecutive two- or three-vertex block on its unique Hamilton cycle; all except its last have outdegree one. For the triangle all three qualify. Every marked start-two vertex r in a C1 graph is obtained by choosing consecutive start-one vertices a,c with a's sole outneighbor c, then adding r beating exactly a,c. This is transitive cloning of a (START_THREE.md, Section 1; START_TWO.md).
4. The source oracle below generates every actual source, including weak sources, with its actual graph recipe (BILINEAR_C2.md, Section 4, and ROOTED.md). A weak source consists of a strong start-one initial component, or a singleton, followed by singleton SCCs.

For clarity the exact oracle is reproduced. Begin at o. A universal sink O and an indexed one-entry gadget M_L act by

    O(D,E,U,d)=(D,0,2U,d),
    M_L(D,E,U,d)=(g(L)D-(L mod 2)E+LU,
                  D-[L=1]E, Q(L)D+B(L)U, D).

Here g(1)=Q(1)=B(1)=1. For L=2l use (g,Q,B)=(g(l),Q(l),2B(l)); for L=2l-1, l>=2, use (g(l)+l-(l mod 2),Q(l)+B(l),2B(l)). The physical gadget has root0, a transitive chain, and root-incidence bits equal to the bitwise complement of binary L-1 read low-to-high. The singleton is the L=1 gadget. All cross arcs point from the old source to the gadget except that its root beats the old marked root; retain the gadget root.

U strictly increases at every step. All sources with U<=B can be enumerated exactly by retaining only such states, trying O and L<=B/U at each predecessor. Duplicate profiles may be merged because these four coordinates determine every later profile and every scalar formula below; one actual graph recipe is retained. This merging is not a claim that the underlying rooted graphs are isomorphic.

## 2. The arbitrary selected-entry graft

Let T be a nonempty tournament, A a nonempty set of selected old vertices, and G a nonempty tournament with marked vertex r. Orient all cross arcs T->G except r->A. Define

    h=H(T), E_A=sum_(a in A)e_a(T), R_A=sum_(a in A)R_a(T),
    K_A=sum_(a,b in A,a!=b) K_(a,b)(T),

where K_(a,b) counts two nonempty ordered spanning paths, the first ending at a and the second starting at b. Then

    H(graft)=h H(G)-E_A s_r(G)+R_A L_r(G)-K_A H(G-r).       (1)

This formula is valid for arbitrary G, including its singleton boundary.

Proof. A Hamilton word uses at most one backward cross arc, since all such arcs leave r. With no backward seam, concatenating the two Hamilton words counts hH(G), with exactly E_A s_r(G) invalid boundary joins. With one backward seam, cutting it gives an old reserved-start two-cover and a new reserved-end two-cover, counted by R_A L_r(G). The only invalid forward join is when the new path ending at r is the singleton r and the preceding old path ends at another selected entry. The remaining G-r vertices form one Hamilton path, and these choices are exactly K_A H(G-r). These cases are disjoint and deletion of the named seam reconstructs the counted pieces.

If |G|>=2, the cycle law is

    C(graft)=(sum_(a in A)s_a(T)) e_r(G).                  (2)

Indeed each cycle must use exactly one backward seam r->a; its remaining pieces are a Hamilton path of T starting at a and one of G ending at r. The first G vertex is not r, so the other cross join is always present. The hypothesis |G|>=2 is essential: with a singleton G the last old endpoint could itself be selected.

Two useful physical inequalities are

    E_A<=h,   R_A-K_A>=h.                                (3)

For the second, cut each Hamilton word of T just before its first selected vertex. Its preceding endpoint, if any, is not selected, so this is an R_A cover not counted in K_A. Concatenation recovers the word. For a nontrivial ending-one G, its root has an inneighbor, whence L_r(G)>=H(G-r)+1: the singleton reserved-root covers give H(G-r), and an ending-root Hamilton word gives another cover. Also L_r(G)>=H(G), by cutting each Hamilton word after the root. Substitution in (1) gives

    H(graft)>=h L_r(G),   H(graft)>=R_A.                 (4)

Indeed

    H(graft)-R_A = hH(G)-E_A s_r(G)
       +(R_A-K_A)(L_r(G)-1)+K_A(L_r(G)-1-H(G-r)),
    H(graft)-hL_r(G) = hH(G)-E_A s_r(G)
       +K_A(L_r(G)-H(G-r))+(R_A-K_A-h)L_r(G).

Every term on the right in each line is nonnegative by the preceding inequalities.

## 3. Four complete Type-I left templates

For any rooted V,c, let J_c mean adding a vertex beating only c. Every other old vertex beats it. Starting from a source v=(D,E,U,d), construct these left graphs and selected entry sets:

- A: add a=J_c(V), then r beating exactly {a,c}; select {r,a}. Require V nontrivial.
- B: the same graph, select {r,c}. Require V nontrivial.
- C: first add c=J_x(V), where x is V's original root, then a=J_c and r beating {a,c}; select {r,x}. The original V may be singleton.
- T: add b=J_c(V), then u=J_b; select {u,b,c}. V may be singleton.

Every permitted left graph is strong with C=1, and its selected start-count sum is three. A/B/C select a start-two root and a start-one vertex; T selects three start-one vertices. Their complete response rows (h,E_A,R_A,K_A) are

    A: (D-E+3U-2d, D-E+2U-2d, 3D-2E+7U-2d, 2D-2E+2U),
    B: (D-E+3U-2d, U,          2D-E+8U-4d,   U+E),
    C: (D+4U,       U+d,       4D-E+14U-d,   U+E+2d),
    T: (D+2U,       D+U+d,     3D-E+8U-d,    2D+4U).     (5)

These response identities in fact hold for every rooted input tournament, without the start-one assumption. The assumption is used for the C1/start-count interpretation and completeness.

Here is a direct path-cover derivation, rather than a fitted coefficient argument. In the A/B graph the two added vertices are a transitive clone pair r->a, both beating c. Their endpoint counts are

    e_c=d, e_r=U-d, e_a=D-E+U-d.

Their start-cover counts are

    R_c=4U-3d, R_r=2D-E+4U-d, R_a=D-E+3U-d.

For R_r split the reserved path into: [r]; [r,a]; r,a,c,...; r,c,...,a; or r,c,... with a in the other path. Their counts are respectively D-E+U, D, U, U-d, U. R_a follows by viewing a as the final singleton seam on V+r; R_c follows by two singleton births beating only c, or by placing the two clone vertices in the two cover paths. The scalar count is D-E+3U-2d by the same one-seam split twice.

For pair A the two directions of K each equal D-E+U: the reserved-start path is either its singleton new vertex, contributing D-E, or begins with new,c, contributing U. For pair B, K_(r,c)=U-d and K_(c,r)=d+E. In the first, a must be last on the c-starting path, excluding the singleton-old-c case. In the second the r-start path is [r] or [r,a], with contributions d and E.

For C, the added graph is r->a->c->x with r->c, and the remaining incidences follow the stated births. One gets e_r=U, e_x=d, R_x=8U-d, and R_r=4D-E+6U. For the last count the reserved path cases [r], [r,a], [r,a,c], r,a,c,x,..., [r,c], and r,c,x,... with a on either path contribute D+2U, D-E+U, D, U, D, 2U. In K_(r,x), a,c are the final chain of the x-starting path, giving U. The reverse direction has reserved prefixes [r], [r,a], [r,a,c], contributing d,d,E. The scalar count is D+4U, as a one-entry seam with a transitive three-vertex gadget marked at its sink.

For T, direct one-entry splits give

    e_u=D, e_b=U, e_c=d,
    R_u=2D-E+2U, R_b=D+2U, R_c=4U-d,
    H=D+2U.

The nonadjacent pair {u,c} contributes U+E. The other pair contributions are 2D-E+U for {u,b} and 2U for {b,c}. To check the latter claims uniformly, if a has sole outneighbor b, then

    K_(a,b)+K_(b,a)=R_a-e_a.

Put G=T-a. The two K directions are R_b(G) and e_b(G); meanwhile R_a(T)=H(G)+R_b(G) and e_a(T)=H(G)-e_b(G). This proves the identity and yields K_T=2D+4U. Summing the displayed endpoint and cover counts proves (5), including all singleton-input boundaries.

These four templates exhaust the Type-I nonordinary entry configurations. For a mixed pair use the inherited clone description: the old degree-one vertex a can be deleted, leaving an actual nontrivial start-one source V,c. The selected start-one vertex in the cloned graph is either a, c, or the third vertex x of the start-one block. The first two are A/B. The third exists exactly when c has sole outneighbor x, so V=J_x(S) for an actual source S, giving C. The triangle boundary is included. A triple of start-one entries is the entire three-vertex block; deleting its first two vertices leaves an actual source and gives T. Deletion leaves a unique marked starting path even when the remainder is weak.

For v in the indicated source set and w=(g,e,L,f) in S*, define

    I_X(v,w)=h_X(v)g-E_X(v)e+R_X(v)L-K_X(v)f.           (6)

Physically use (5) on V and graft the converse of W at its marked root. Equation (2) gives exactly three Hamilton cycles, hence strongness. No interface return is asserted; the newly marked root can have many starting paths.

## 4. Four shared-root triangle templates

For rooted A,B, identify their roots r and orient all A-r -> B-r; call the result A⊕B. Use two sources V,W, a directed triangle on x,y,z, and a shared root r beating exactly j triangle vertices, j=0,1,2,3. For j=1 choose {x}; for j=2 choose {x,y}; these exhaust the rotation orbits. Let K_j be this rooted four-vertex graph and set

    J_j(V,W)=V ⊕ K_j ⊕ converse(W).                   (7)

Both sources are ordinarily nontrivial. Also allow precisely these singleton boundaries: W singleton and j=0 with V nontrivial, or V singleton and j=3 with W nontrivial. Do not allow both singleton.

Every Hamilton cycle started at r traverses V-r, then the triangle, then W-r, and returns to r. The nonempty source segments have one permissible Hamilton word each; the triangle has three. Thus all allowed graphs are strong C3. With just W singleton the cycle count is 3-j; with just V singleton it is j. This proves exactly the stated boundary restrictions.

The scalar forms are v M_j w^t, in coordinate order (D,E,U,d), where

    M0=[[0,0,0,3],[0,0,0,-3],[0,0,12,-3],[3,0,-12,3]],
    M1=[[0,0,0,3],[0,0,0,-2],[0,0,12,-6],[3,-1,-9,5]],
    M2=[[0,0,0,3],[0,0,0,-1],[0,0,12,-9],[3,-2,-6,5]],
    M3=[[0,0,0,3],[0,0,0,0], [0,0,12,-12],[3,-3,-3,3]]. (8)

For a self-contained derivation use the general ordered-amalgamation identities

    h'=(R_A-d_A)(L_B-d_B)+(h_A-e_A)d_B
          +(h_B-s_B)d_A+e_A s_B,
    R'=(R_A-d_A)B_(B-r)(2)+d_A R_B,
    d'=d_A d_B, e'=d_A e_B.                           (9)

The last endpoint formula requires B nontrivial; if B is singleton the operation is the identity. The scalar and cover formulas have no such exception. Split a Hamilton word before and after r. In each piece the A vertices precede the B vertices. The four cases according to a B vertex before r and an A vertex after r give the four scalar terms in (9). For R', either the reserved r-start path uses A-r, allowing arbitrary two labelled bins of B-r, or it does not and all A-r vertices lie on the other path. This proves the cover formula. Uniform cross arcs prove the deletion and endpoint laws.

The (h,s,e,R,L,d,Bdelete) data of K_j are

    j0: (3,0,3,3,12,3,12), j1: (5,1,2,6,9,3,12),
    j2: (5,2,1,9,6,3,12), j3: (3,3,0,12,3,3,12).

They are obtained by listing the triangle's proper-subset paths. Applying (9) twice gives (8); using the scalar identity when a source is singleton gives the same boundary formulas. Equations (9) and this specialization were independently developed with the separate reviewed ROOT_AMALGAMATION packet; that note's proof is included as a named supporting dependency, but the derivation above is sufficient for (8).

Reversing arcs exchanges the two sources and j with 3-j; the triangle is relabelled by a reversal followed, if needed, by a rotation. Algebraically M_j^t=M_(3-j). We permit converses of the Type-I graphs too. Converse changes no Hamilton-path count.

## 5. Completeness for every strong C3 graph

Theorem. Every strong C3 tournament is isomorphic to a Type-I graph in (6), a Type-II graph in (7) with its stated boundaries, or the converse of one of these. Conversely every graph in this union is strong C3. This is a surjective physical normal form, not an injective parametrization.

Proof. Reverse all arcs if necessary to obtain an outdegree-one vertex v. Let a be its sole outneighbor. Delete v and take the ordered SCC decomposition. Since every v-start Hamilton word closes to v, its starting count is three. The SCC product yields exactly the two cases in START_THREE Section 5: an initial marked start-three component and a singleton tail, or an initial marked start-one component with one triangle and singleton components in its tail.

In the first case peel the recursive one-entry Type-I branches from the initial rooted start-three component. Each peeling strictly decreases its order, so it terminates at either a nonordinary Type-I mixed/triple entry or a Type-II triangle tail.

If it ends at a mixed/triple entry, its initial component is one of the four left templates in Section 3. Bundle the first new entry root, all later vertices removed during peeling, all singleton tails, and the final degree-one closer v into a right tournament G. Across the left component, all arcs point right except that the first entry root beats exactly the chosen entries. This follows at each peeling from the ordered SCC decomposition and its single-entry incidence, so bundling introduces no additional backward cross edge. G includes the entry root and v, hence |G|>=2. The selected start sum is three. Equation (2), together with the target's C=3, forces its root-ending count to be one. Therefore G is the converse of an actual nontrivial source W, as required in (6). This also covers the otherwise invalid empty-tail whole-triangle branch: the final bundle is nontrivial and supplies the return edge.

If peeling ends at Type II, write that rooted graph as r plus an initial start-one component U, a singleton chain P, one triangle, and a singleton chain Q. Bundle all ordinary later suffixes and the final degree-one closer into G. Here G is nonempty, and is just the singleton closer if there were no ordinary suffixes. Every earlier vertex points to G except its first entry g points to r. The single-entry cycle cut gives C(final)=3e_g(G), so e_g(G)=1, including that singleton case.

Take V=r+U+P and take W^op=r+Q+G. Both are nontrivial. A Hamilton path starting at r in V must enter U at its unique selected entry, traverse U, and then P, so s_r(V)=1. A path ending at r in W^op must traverse Q and then G, ending at g before r. Ending directly from Q would strand G, which has no edge back into Q. Thus e_r(W^op)=1. All cross arcs outside r are V-r -> triangle -> W^op-r, so the final graph is literally (7), for the actual root-triangle incidence orbit.

In the second original SCC-product case, v is itself the degree-one closer of a start-one initial component and a triangle tail. Put v, the initial component, and all pretriangle singleton components into V. If singleton components follow the triangle, put these and v into W^op; the ending-v path is unique. If the triangle is terminal, use W singleton. All triangle vertices point to v, giving j=0, precisely the allowed boundary. The converse orientation gives the j=3 boundary. This exhausts both original cases. Soundness was proved in Sections 3 and 4. QED.

The only inherited finite-reduction theorem is the all-order C3 bridge theorem in Section 1. All subsequent reductions here are exact graph identities and terminating decompositions. No order cutoff or observed support agreement is used to infer completeness.

## 6. Exact finite scalar inverse and physical cutoffs

For a nontrivial source, the root has an outneighbor. In each Hamilton word of V-c, insert c immediately before the first vertex it beats. If that vertex is not first, its predecessor points to c. This is an injective insertion into a non-ending-c Hamilton word, proving D-E>=d. Cutting every Hamilton word before its root gives U>=D. The d covers whose reserved path is singleton root contain exactly E covers in that insertion image. Hence

    U>=D+d-E,    U>=2d.                               (10)

Set x=(D-E-d,U-D,E,d). Its coordinates are nonnegative and sum to U. In Type I, the scalar h_X in (5) is at least 2U, and every permitted left graph has h_X>=3. The right source is nontrivial and has U_right>=2. Equations (4) therefore imply

    U_left<=H/4,   U_right<=H/3.                       (11)

The singleton left source also satisfies the common cutoff U_left<=H/3 because all outputs are C3 and have H>=15.

For two nontrivial sources in Type II, substituting x and y into (8) gives these nonnegative matrices:

    P0=[[12,12,12,12],[12,12,12,9],[12,12,12,9],[3,0,3,6]],
    P1=[[12,12,12,9],[12,12,12,6],[12,12,12,7],[6,3,5,8]],
    P2=[[12,12,12,6],[12,12,12,3],[12,12,12,5],[9,6,7,8]],
    P3=[[12,12,12,3],[12,12,12,0],[12,12,12,3],[12,9,9,6]].

For P1/P2 every entry is at least three. For P0 its first three rows have every entry at least nine, while x_4=d<=U_left/2 by (10); thus H>=9(U_left-d)U_right>=9U_left U_right/2. P3 follows by transposition. In particular

    H>=3 U_left U_right.                             (12)

For the singleton boundary J0(v,o), H=3(D-E)+9U-6d>=6U, and the converse boundary is identical. Consequently every source in every allowed representation satisfies the convenient uniform bound

    U<=floor(H/3).                                   (13)

This deliberately conservative bound avoids any claim of optimal inverse complexity.

Algorithm for a target N: generate the exact source-profile bank with U<=floor(N/3), retain an actual predecessor recipe for each state, and evaluate all eight forms on all admissible ordered source pairs. Return a graph recipe for every match. There is a match if and only if N is the Hamilton-path count of a strong C3 tournament, by the complete union and (13). The same algorithm gives the full bounded image through N. Scalar converses need no second evaluation.

All coordinates in the bank lie between zero and its U cap, so there are O(N^4) states. Trying L<=N/U gives a finite numeric-polynomial enumeration, and all ordered pair tests use O(N^8) arithmetic operations as a conservative bound. Integers have O(log N) bits up to fixed-degree intermediate products. This is polynomial in numeric N, not in its bit length. Predecessor recipes terminate because U strictly increases and construct actual tournaments. No inverse closure of a larger formal cone is assumed.

## 7. Verification, reproducibility, and exact scope

Run `python3 verify.py` and `python3 -O verify.py` from a copied packet directory. Both modes use only the Python standard library and write SUPPORT.json and VALIDATION.json beside the verifier. No repository checkout, network, unreviewed graph bank, or third-party package is needed. core.py and source_oracle.py are unchanged copies of the independently reviewed physical helper and complete C2/source-oracle module; their bytes are pinned in RELEASE_FILES.sha256.

The verifier checks all four response rows on every rooted labelled tournament through order four, even outside the source class; the arbitrary selected-entry graft scalar and cycle laws on all labelled input orders at most three; all physical rooted start-one source graph classes through five; literal template graphs and their converse closure through order eight; the entire labelled C3 universe through six; and every inverse witness through 501. Endpoint subset counts and directed-word enumeration are separate counters. A full unpruned scan of the same source bank checks that the safe pair-pruning rules lose no bounded scalar target. Exact totals and the support are recorded in VALIDATION.json, not inferred in this proof.

The all-order result is a complete physical C3 source reduction and a constructive exact membership algorithm. The finite support is only a finite support result. The union is redundant, source profiles need not determine source graphs, and its output roots do not constitute a returning start-one interface. This does not prove that C<=3 realizes every odd number other than 7 and 21, does not prove that arbitrary odd targets admit a scalar descent, and does not settle the unrestricted spectrum.
