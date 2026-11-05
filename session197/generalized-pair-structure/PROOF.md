# An intrinsic graph criterion for the generalized pair families

Session197. Let Gamma be the unrooted graph image of the separately reviewed generalized N/M/P/Q constructions. This note gives an exact graph criterion, a canonical rooted N/M/Q compiler, and explicit obstructions to graph-level absorption. It does not exclude any Hamilton-path scalar from Gamma's image.

For a tournament T, H counts Hamilton paths and C counts Hamilton cycles modulo rotation. At a root r, s_r counts paths starting there. For a strong C2 graph and a root with s_r=2, let U_r be the initial SCC of T-r, in the unique dominance order.

**Theorem.** A strong C2 tournament T belongs to Gamma if and only if it has a root r with

    s_r(T)=2 and C(U_r)=1.                           (1)

For a specified root satisfying (1), its rooted graph has exactly one N, M or Q presentation up to a rooted source isomorphism. The three cases are determined intrinsically by the canonical start-two decomposition. P is redundant in the generalized class, by an exact labelled identity, not merely equality of H.

## 1. The inherited canonical structure

The reviewed start-two grammar says that deleting r leaves one initial nontrivial strong component U, followed by a chain W of singleton SCCs. The entry set inside U consists either of one start-two vertex or of two start-one vertices. In the latter case U is C1. If W is nonempty its last vertex points to r. Every U vertex points to W and every non-entry U vertex points to r. The ordered SCC decomposition and entry set are intrinsic.

In a strong C1 tournament, its start-one vertices form a consecutive block of two or three vertices on the unique cycle; all but the last have outdegree one. The triangle has all three vertices start-one and every unordered pair is cycle-consecutive. A selected pair is thus either an adjacent pair a->b with a's sole out-neighbor b, or the first and last vertices u,c of a genuine three-vertex block u->b->c, where u,b each have outdegree one. The latter case does not occur in the triangle.

Also, a C1 tournament with a distinguished start-two root a is uniquely the first transitive2 clone of an outdegree-one vertex y in a smaller strong C1 graph. Its root a beats exactly y and y's sole out-neighbor b, while a and y lose to every other vertex. Here is the short reason, recalled from START_THREE Lemma1. In the start-two decomposition, a nonempty deletion tail would let both starting paths close and give C2. A single entry likewise gives C2. Therefore the deletion is strong C1 and the two entries must be consecutive on its unique cycle, for exactly one of the two paths to fail to close. The block theorem gives the stated sole-out-neighbor vertex and hence the clone description.

## 2. Every family has the local property

Use the source V and its start-one root b from the generalized-family packet. Its constructions and their marked vertices are fixed as follows.

N has V followed by two transitive clones a->y, both pointing only to b among V, and retains a. Deleting a leaves the strong C1 graph V+y, followed by the suffix's transitive chain. Its initial SCC therefore has C1. The designated root has start count two.

M additionally appends z beating only a and retains the last clone y. Deleting y leaves the two successive singleton one-seam extensions V+a+z, which are strong C1, followed by the transitive suffix. Again (1) holds.

Q attaches an ordinary gadget at the first clone a of the N seed and retains the new gadget root. Deleting that root leaves the C1 clone seed followed by the transitive gadget base. Thus (1) holds. Singleton V is excluded for these three families, exactly as in their definitions.

For P, write its transitive cell in the order c->d->a, where entry a points back to b and exit c is the next singleton attachment's entry. Let J_b(V) append a beating only b and retain a as the source root. Then for every positive suffix parameter v,

    P_V(v) = Q_(J_b(V))(v),                         (2)

as actual labelled rooted graphs under the construction's natural vertex order. In Q, the next two vertices c,d both beat only a and satisfy c->d. Its gadget root then beats only c. These are precisely P's cell, final root and binary tail. No relabeling is needed. The identity includes singleton V: J_b(V) is then the permitted two-vertex transitive source for Q. Thus P also has (1).

## 3. Constructing the inverse in every qualifying case

Suppose (1) holds. Let U=U_r and let W=(w_1,...,w_m) be the singleton tail of T-r. Encode the root-tail incidences by

    v=1+sum_(i:w_i->r) 2^(i-1).                     (3)

For m=0 this is1. For m>0 the final tail vertex points to r, so the highest digit is1 and the binary suffix has exactly length m. This is exactly both the retained O/Z parameter and the ordinary gadget index in the reviewed definitions.

### Two adjacent entries: N

Let the selected entries be a,b with a's only out-neighbor in U equal to b. Delete a from U and retain b in the remaining V. Every b-start path of U must visit a last, because a cannot leave toward its already visited neighbor b. Removing or appending a is a bijection with b-start paths of V. Thus V has exactly one such path, and |V|>=2.

The graph T is literally N_V(v): the new first clone is r and the old last clone is a. The tail is W with incidences (3). Here m cannot be zero, since an empty-tail adjacent pair gives C1. Thus v>=2, exactly N's domain. The rooted source and its parameter are recovered uniquely.

### Two nonadjacent entries: M

The selected vertices are the outer vertices u,c of the three-vertex start-one block u->b->c in U; u and b have outdegree one. Delete u and b, retaining c in V. A c-start path of U must end in the consecutive suffix u,b: b must be last because it points only to the already used c, and u must immediately precede b. Removing or appending this suffix is a bijection with c-start paths of V. Thus V has start count one. The genuine nonadjacent case has |U|>=4, so |V|>=2.

The M graph has, beyond V, its first clone, its retained last clone, and the vertex beating only the first. Map these three vertices respectively to b,r,u. Every edge agrees: b->r, r->u,c, u->b, while all other source incidences follow the two sole-out-neighbor relations. Append W unchanged. This is exactly M_V(v), including v=1. The three-block theorem makes u,b,c and the extracted rooted source intrinsic.

### One start-two entry into C1: Q

Let a be the sole selected vertex in U. It has start count two. By the C1 clone lemma in Section1, U is the first-clone graph of an outdegree-one vertex y with sole neighbor b. Remove a and y, retaining b in V. Deleting a leaves a strong C1 graph in which y points only to b; the preceding path bijection then proves s_b(V)=1. The smaller C1 graph has at least three vertices, so |V|>=2.

The clone seed is exactly F(V), with first clone a and last clone y. The outer root r and tail W form the ordinary gadget G(v). Hence T is Q_V(v), at its prescribed new gadget root. The clone description is unique for the marked root a, so the rooted source and parameter are again uniquely recovered.

These three mutually exclusive cases prove both sufficiency and the canonical rooted inverse. Necessity was proved in Section2. The compiler returns the full vertex map and verifies every edge; it does not infer graph identity from equal numerical responses.

## 4. Exact rooted presentation count

As in the separately reviewed count note, let Delta=1-3z+z^2. All rooted start-one sources, strong or weak, have generating function

    W(z)=z(1-2z)/Delta,

and all retained suffixes, including the empty one, have S(z)=(1-z)/(1-2z). The canonical N/M/Q inverse proves that their presentations now count rooted isomorphism classes satisfying (1) exactly, not just an upper bound. Consequently

    R_Gamma(z)=(W-z)[z^2(S-1)+2z^3 S]
              =z^5(1-z)(3-2z)/[Delta(1-2z)].        (4)

P contributes no additional rooted graph because of (2). Formula(4) is not an unrooted Gamma count: one graph may have several qualifying roots. The prior four-family upper bound remains valid, although redundant P presentations made it less sharp.

## 5. Small explicit graph obstructions and scalar rebasing

Let

    A=[2,4,9,19,103,79,143,63].

This order8 tournament has H=77 and C=2. Its start-two roots are exactly0,1,2. Their initial deletion SCCs have vertex sets respectively1..7,2..7,3..7, all with C=2. Thus A is outside Gamma by the all-order theorem, independently of any presentation census. Its converse does belong to Gamma; the verifier supplies an explicit M compiler witness.

A stronger two-sided example is

    B=[2,4,9,83,39,79,151,319,639,255].

It has order10, H=289 and C=2. Its start-two roots are exactly0,1,2; each initial deletion SCC has C=2. Its converse has start-two roots exactly7,8,9, again all with initial deletion C=2. Hence neither B nor its converse belongs to Gamma. Allowing only a converse operation does not repair graph-level absorption in general.

These are not scalar obstructions. P from the singleton at v=8 gives H77 at order8, and P from the singleton at v=26 gives H289 at order10. The latter actual rows are

    [1020,1009,1002,1010,196,976,896,768,528,16].

Thus even a graph outside Gamma in both orientations can have the same H as a graph inside Gamma at the same order. A scalar-only rebuilding theorem remains a legitimate possibility; none is claimed here.

## 6. Controls and dependencies

The verifier enumerates the complete rooted start-two grammar and complete C2 graph classes through order10, compares the local criterion with every physical generalized-family presentation at each order, and checks the exact compiler's vertex maps for every qualifying root. The complete Gamma graph counts at orders5..10 are1,5,18,57,169,482, compared with complete C2 counts2,5,18,60,191,589. These finite inventories corroborate the theorem, not its completeness proof. It also tests the P=Q(J) identity including singleton inputs, counts the explicit hostile graphs directly, and supplies their scalar repairs. The SCC order is checked structurally; arbitrary vertex intervals are never assumed.

The all-order inputs are START_TWO.md's canonical grammar and C1 block theorem; START_THREE.md Section1's C1/start2 clone lemma; the reviewed generalized-pair-inverse definitions; and the rooted source enumeration in FIBONACCI_CLASS.md. C2_COMPLETE.md is used for the full finite-class inventory controls and the rigid marked-to-unmarked check, not to infer the criterion from data. Exact dependency hashes are in RELEASE_FILES.sha256. Run python verify.py and python -O verify.py on a copy; both write the validation receipt beside the script.
