# Four generalized pair families and an exact finite inverse

Session197. This packet extends the reviewed fixed-two-port N/M/P/Q scalar formulas from one-apex sources to arbitrary actual rooted tournaments having exactly one Hamilton path from the root. It gives all-order physical constructions and a complete finite membership algorithm for their union, denoted Gamma. It does not identify Gamma with the full C2 class or prove cofinite coverage.

## 1. Source data and elementary inequalities

Let V be nonempty, with root b and exactly one Hamilton path starting at b. Write

    D=H(V), E=e_b(V), U=R_b(V), d=H(V-b).

Here R_b counts ordered spanning covers (P,Q), where P may be empty and Q is nonempty and starts at b. All paths are directed; singleton paths are allowed and H(empty)=1. These four numbers are measured on one physical graph. In particular d cannot be chosen from the triple independently.

Every full path can be cut immediately before b, injectively giving a reserved-start cover. Thus U>=D. If |V|>=2, the unique starting path shows that b has an out-neighbor. Given any path of V-b, insert b immediately before its first out-neighbor on that path. If there is a preceding vertex, it points to b by the choice of the first out-neighbor. The resulting path does not end at b, and deleting b recovers the input. Therefore

    D-E >= d,                 |V|>=2.              (1)

The same insertion argument without an endpoint restriction proves Hamilton-path monotonicity under an arbitrary induced extension: insert one vertex before the first vertex it beats, or at the end if it beats none. Iterating preserves injectivity by deletion. In particular d<=D for nonempty V-b, also true at the singleton with our empty convention.

## 2. Literal gates and all seed responses

The operation J_r appends a vertex whose sole old out-neighbor is r. The new vertex becomes the root if a root is specified. More generally K(T,r;G) orients every cross edge T->G except the edge from the distinguished gadget root to r.

For |V|>=2, append a=J_b(V), then append y beating only b. Call the resulting graph F. Thus a->y, and both a,y point to b and lose to every vertex of V-b. They are a transitive pair replacing the sole-out-neighbor vertex a. The roots used below are exact physical vertices:

- N seed: F, retaining a, the FIRST clone.
- M seed: F with z appended beating only a, retaining y, the LAST clone.
- Q seed: the same F and first clone a as N.
- P seed: attach the transitive three-vertex cell with local rows [0,5,1] to V, entry local0 (last in the transitive order), exit local1 (first). Then append J at that exit and mark its newborn. This construction also permits singleton V.

The full seed profiles (h,e,R,d0) are

    N: (D-E+3U-2d, U-d, 2D-E+4U-d, D-E+U),
    M: (3D-2E+6U-2d, D+U-d, 3D-E+8U-d, D+2U),
    P: (5D-E+9U, D+3U, 5D-E+10U, D+4U).           (2)

Every seed has two Hamilton paths starting at its indicated root. N has C=1; M and P have C=2. All these graphs are strong.

Here are physical derivations, rather than interpolation of the rows. Let T=V+a before cloning. The ordinary J and retained-root Z laws give

    H(T)=D-E+U, e_a=D-E, R_a=D+U,
    e_b=d, R_b=2U-d.

T is strong C1: the one-seam cycle formula gives C=s_b(V)=1, since total order is at least three. Also s_b(T)=1, because in a path starting at b the newborn a must be last. Applying the reviewed adjacent first-clone identities yields N's first three rows in (2). Deleting its first clone leaves a copy of T, yielding d0. Adding the second clone y is a one-seam extension at b, so C(F)=s_b(T)=1. Its first-clone start count is two.

For M, delete y and call the remaining graph W. It is J_a(T), with H(W)=D+2U. A path of W ending at b must finish z,a,b; the earlier vertices form any path of V-b, giving e_b(W)=d. Also e_z(W)=H(T)-e_a(T)=U. Since y's only out-neighbors in W are b,z, appending y to a W-path is valid except for these endpoints, giving e_y=D+U-d. The scalar row follows from J at the first clone in (2).

For R_y, split the reserved path by its beginning. A singleton y contributes H(W)=D+2U. A path beginning y,z either stops there (contribution H(T)), stops after y,z,a (contribution D), or continues through b (contribution U). A path beginning y,b distributes z,a into the two slots of a cover of V reserved at b. There are four placements per cover, except the placement ending in a in the reserved singleton-b slot. The excluded covers are exactly the d paths on V-b. This case contributes 4U-d. Summing gives the stated R_y. The two paths starting at y are y,z,a followed by the unique b-start V path, and y followed by that V path and then z,a. Thus s_y=2. Its cycle count is two by the singleton one-seam law at a.

For P, the reviewed entry-last/exit-first cell acts on the physical triple by

    (D,E,U) -> (D+4U, U, 4D-E+6U).                 (3)

This is the audited distinct-terminal matrix [[1,0,4],[0,0,1],[4,-1,6]]. Its exit has exactly two starting paths: in cell order z->y->r, either z,y,r followed by the unique b-start V path, or z,r followed by that V path and then y. Applying J at z gives (2), with deletion count D+4U and start/cycle counts two. This remains valid for singleton V. Equation(3) is the same physical cell used and proved in the reviewed fixed-two-families packet; the current code also counts every small source directly.

Singleton V is deliberately excluded from N/M/Q: their intermediate T would have order two and the strong C1/first-clone start-two argument would fail. It is not excluded from P.

## 3. Integer suffixes and exact scalar laws

Use the ordinary one-apex gadget G(v) of the rooted calculus. G(1) is a singleton. For v>=2 its transitive base incidence word is the complemented binary digits of v-1, least significant first. Its ending count is one and its deletion is transitive. Put g(v)=H(G(v)), t=v mod2. The exact recursion is

    g(1)=1,
    g(2l)=g(l),
    g(2l-1)=g(l)+l-(l mod2), l>=2.                 (4)

The accompanying constants Q(v),B(v) satisfy the ROOTED recurrences; B(v)=2^bit_length(v-1), with B(1)=1. The ordinary transfer is

    (h,e,R) -> (g(v)h-t e+v R, h-[v=1]e,
                Q(v)h+B(v)R).                    (5)

For N/M/P, retain the seed root and append O/Z according to the binary digits of v-1, least significant first. O is a universal sink and Z beats only the retained root. For v=1 append nothing; otherwise the last digit is1, so the suffix ends in Z. Their exact four-coordinate laws are

    O(h,e,R,d0)=(h,0,2R,d0),
    Z(h,e,R,d0)=(h-e+R,d0,2R-d0,d0).

Induction on the digits, using (4), gives the scalar suffix law

    h_out=h+(g(v)-1)d0+(v-1)(R-d0)+(1-t)(d0-e).    (6)

It also holds at v=1. N requires v>=2; M/P permit v>=1. For Q, instead apply the ordinary gadget G(v) at F's first clone, for v>=1.

Substituting (2) gives the promised formulas:

    N_V(v)=(D-E+U)g(v)+(D+3U-d)v-U-E-t(D-E+d),
    M_V(v)=(D+2U)g(v)+(2D-E+6U-d)v-E-U-t(U+d),
    P_V(v)=(D+4U)g(v)+(4D-E+6U)v-tU,
    Q_V(v)=(D-E+3U-2d)g(v)+(2D-E+4U-d)v-t(U-d).   (7)

Every output is strong with exactly two Hamilton cycles. An O/Z suffix on a graph of at least two vertices preserves its root-start count: any root-start path must visit the newborn last, since the only possible reverse arc points to the already visited root. At the final Z, the one-seam cycle bijection gives exactly the previous root-start count, two. Empty M/P suffixes already have C=2. For Q, (5)'s gadget has ending count one and the old root has start count two, so again C=2. Intermediate O steps may be weak; no cycle-monotonicity assertion about those steps is used.

The operations are physical and append-only once V is chosen. They are compatible with the established changing-two-port compiler for strong C2 graphs, but this packet does not claim two permanently retained physical ports or feedback vertex number at most two for arbitrary V. The previous fixed-port formulas are recovered by the special one-apex data (D,E,U,d)=(g(u),u mod2,u,1). Weak general sources can have several possible d values at one triple.

## 4. Proven source and suffix cutoffs

The smallest relevant N/M/Q graph is the same scalar

    N_V(2)=M_V(1)=Q_V(1)
          =3D-2E+6U-2d
          =7D+6(U-D)+2(D-E-d) >=7D.               (8)

Every N output contains its seed plus one Z vertex induced, every M output contains its seed, and every Q output contains its seed plus the gadget root. Hence induced path monotonicity proves the same lower bound for every permitted parameter. For P, its seed has H=5D-E+9U>=13D, including the singleton, and is induced in every output. Therefore every Gamma presentation of a target h must satisfy

    D <= floor(h/7).                               (9)

For a weak source E=0, d<=D makes the common expression in (8) at least D+6U. P's minimum is larger still. Thus every relevant weak source also satisfies

    U <= floor((h-D)/6).                            (10)

No universal bound on arbitrary rooted R from H is assumed.

Finally g(v)>=1, t<=1 and (1) give the safe explicit bounds

    N: v <= floor((h+E+d)/(D+3U-d)),
    M: v <= floor((h-D+E+d)/(2D-E+6U-d)),
    P: v <= floor((h-D-3U)/(4D-E+6U)),
    Q: v <= floor(h/(2D-E+4U-d)).                  (11)

For example N's omitted constant lower bound is -E-d. Q's omitted constant is D-E+2U-d>=0. All denominators are positive on their permitted domains. These bounds may include impossible parameters, which are rejected by checking (7); they never remove a solution.

## 5. Complete finite source enumeration

The reviewed ROOTED classification obtains every strong start-one source from a singleton by ordinary gadgets (5). The deletion response after a move is the old H: deleting the new root leaves the old graph followed by a transitive gadget base. Consequently the exact four-coordinate recurrence is

    (h,e,R,d) -> (g(L)h-(L mod2)e+L R,
                  h-[L=1]e, Q(L)h+B(L)R, h).       (12)

Every output of order at least three is strong C1 and start-one; the only smaller outputs are the singleton and a two-vertex transient. Formula(12) retains d even when identical triples could otherwise conceal information.

For a source cutoff B, start from (1,1,1,1), try exactly 1<=L<=floor(B/R), and retain each new quadruple with h<=B. This is a complete finite enumeration of the strong/singleton responses up to B. Indeed h_new>=L R and induced H never decreases, so no route to an accepted final source can leave the cutoff. This enumeration is finite without any guessed depth or graph-order limit. The gadget recurrences imply Q(L)<=B(L)<=2L. With h<=R, every accepted successor has R_new<=4B. Also 0<=e,d<=B. Therefore only O(B^4) quadruples can be visited. Each has at most B tested gadgets. The loop stores a physical predecessor recipe for every retained state. The exceptional count-one transient is (1,0,2,1).

For completeness on weak sources, the SCC factorization gives an exact description. A Hamilton path starting at b must begin in the initial SCC. Its count is that component's b-start count times H of all later components. A nontrivial strong component has H>=3 (cut a Hamilton cycle), so start count one forces all later SCCs to be singletons. Thus every weak start-one source is uniquely

    C -> Tr_k, k>=1,
    (D,E,U,d)=(H(C),0,2^k R_b(C),H(C-b)),           (13)

where C is a strong start-one core or singleton. Appending a sink doubles U by placing the sink last in either cover slot, including a formerly empty first slot. It preserves H and d. Conversely every such source is valid. The initial strong-core grammar is complete, and (10) bounds k. For C singleton, this includes every nontrivial transitive source. The two-vertex transient from (12) is omitted as a core and supplied canonically by k=1 instead.

The algorithm therefore computes all relevant quadruples, including every different d at an equal weak triple, from B=floor(h/7), (12), and (13). It then checks every integer in (11) against (7). This is an exact decision algorithm for Gamma and returns an explicit graph recipe for every accepted numerical presentation, using one source witness per quadruple. It need not enumerate every isomorphism class with those responses.

It is polynomial in the numerical target under elementary integer arithmetic, not asserted polynomial in the input bit length. For example the crude state and loop bounds above give polynomial numerical complexity; all generated response integers have polynomially bounded bit length. This result is not an adaptive scalar descent proving all admissible h reachable. The separately reviewed backward-profile inverse and weak-source oracle could replace forward source enumeration, but are not implementation or proof dependencies of this algorithm.

## 6. Bounded support and physical repairs

The exact algorithm at h<=1861 gives 1,316 strong/singleton/transient source quadruples through floor(1861/7)=265 and 2,000 relevant source quadruples after weak tails. Gamma has 883 distinct output values in that range. Independently enumerate the complete rooted C1 grammar through1861 (23,432 response quadruples). The C1 completeness theorem identifies its nontrivial scalar image with the full strong C1 class. Adjoining the singleton value1, the union of C1 and Gamma contains precisely every odd value from1 through1861 except

    7,21,67.

These bounded-output statements hold over all graph orders because the cutoffs were proved. They do not extrapolate to larger numerical targets. The following actual C2 repairs show the generalized gates reaching beyond the fixed-two-port C2 families:

- 533: strong source recipe [7,2], then one sink, data (37,0,74,11), N parameter2.
- 821: strong source recipe [3,2], data (19,5,19,5), P parameter3.
- 1861: strong source recipe [4,1], then one sink, data (9,0,20,5), N parameter21.

A recipe list means ordinary gadget indices applied from the singleton. Their full graph rows and marked roots are stored in VALIDATION.json and counted directly there. Every one of the 883 stored support witnesses is reconstructed and independently checked by subset Hamilton-path/cycle DP, not just scalar recurrence evaluation.

The weak deletion sidecar is necessary, not just a bookkeeping convenience. Strong sources with profiles (9,1,9,d=1) and (9,3,9,d=3) both become the weak triple (9,0,18) after one sink. Their N seeds have different profiles (61,17,89) and (57,15,87), and N parameter2 gives different actual C2 values133 and129. Thus forgetting d would give a false universal three-coordinate gate law. Both quadruples are retained and checked by the verifier.

The support alone does not show that arbitrary sequences of ordinary moves after a pair gate can be absorbed into one of (7). Any such scalar reduction would need a new proof and may require another root or family. Neither gcd one, the strong source grammar, nor this bounded interval supplies such a proof.

## 7. Verification, dependencies and scope

Run python verify.py, save VALIDATION.json, then run python -O verify.py and compare. The verifier writes SUPPORT.json and VALIDATION.json beside itself. It checks literal small graph sources and seeds, all scalar/cycle outputs, complete bounded source/witness recipes, exact inverse membership and an oversized source/parameter rectangle beyond the pruning limits. It uses explicit exception checks, so optimization does not remove validation.

Pinned mathematical dependencies, copied under dependencies/:

- ROOTED.md, SHA256 ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e: exact gadgets and complete strong start-one grammar.
- PAIR_ENTRY.md, SHA256 fcff7ca192c4055243aeab2344582686214aa7aa25504d1233f3b6b33465a856: adjacent-clone responses and J/O/Z laws.
- FIXED_TWO_FAMILIES.md, SHA256 c48c6deca5f5a0fd5b8336ec11a297dd8bf8d97aac83bca726390044bb22d545: the exact physical transitive cell (3) and binary suffix law (6).
- C1_COMPLETE.md, SHA256 d38a53b1c475df3578d58e83d9e9f97785d9455084c0a094aaed752577d8b289: used only when calling the full rooted grammar image the complete C1 scalar class.

The included core.py is the unchanged reviewed pair-entry graph/counter implementation, SHA256 dc32dc1ef3647b694015cb4ff1ef99415c162bfbaf58d3b545175b09f0378042. No external packages, graph census, floating-point logarithms, saturated signed readouts, or unproved scalar coverage hypothesis enter the algorithm.
