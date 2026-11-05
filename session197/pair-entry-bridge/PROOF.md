# A physical pair-entry bridge for the complete low-cycle inverse problem

Session 197. This packet gives an exact arithmetic branching event, an actual same-graph root choice, and a finite inverse criterion for the full strong C≤2 class. It does not prove all-target coverage. The old one-seam formulas and their retained-root specializations are used as ingredients, not claimed as new in isolation.

## 1. Responses and the allowed gadgets

Write H(T)=h. At a vertex a let e_a count Hamilton paths ending at a, s_a those starting there, and R_a count ordered spanning covers (P,Q), where P may be empty and Q is a nonempty path starting at a. Let d_a=H(T−a), with H(empty)=1.

For distinct vertices a,b, let K_ab count ordered covers by two nonempty paths, the first ending at a and the second starting at b. Define the actual pair response

    E=e_a+e_b,   R=R_a+R_b,   K=K_ab+K_ba.             (1)

These numbers are measured on one physical graph with both a and b exposed. They are not freely chosen integers.

The allowed gadget G(L), for L≥1, consists of a root r and a transitive tail. G(1) is a singleton. For L≥2, complement the bits of L−1 read least significant first; bit 1 means r beats the corresponding tail vertex. The last bit is 0. Thus G−r is transitive and exactly one Hamilton path of G ends at r. Its constants are

    g=H(G), Q=R_r(G), B=2^(|G|−1),
    L=L_r(G), t=s_r(G)=L mod 2, δ=[L=1].

The exact integer recurrence, inherited from ROOTED.md, is

    (g,Q,B)(1)=(1,1,1),
    (g,Q,B)(2l)=(g(l),Q(l),2B(l)),
    (g,Q,B)(2l−1)=(g(l)+l−(l mod 2),Q(l)+B(l),2B(l))

in the last line for l≥2. The usual single-entry gluing, all old vertices pointing to G except r→a, returns the new root and sends

    M_L(h,e,R)=(g h−t e+L R, h−δ e, Q h+B R).        (2)

## 2. Exact pair-entry transfer

Instead reverse both cross arcs r→a and r→b, with every other cross arc T→G. Mark the new root r. Then

    h' = g h−t E+L R−K,
    e' = h−δ E,
    R' = Q h+B R.                                    (3)

This is an all-order identity for arbitrary T and distinct a,b. Start-one assumptions are needed only for the low-cycle classification below.

For h', paths avoiding the reverse arcs are an old Hamilton path followed by a gadget Hamilton path. The forbidden joins contribute tE. A path using a reverse arc has old runs P,Q, with Q starting at one of a,b, and gadget runs X,Y, with X ending at r. This gives LR possibilities. The only forbidden remaining join is P→X when X is the singleton r and P ends at the other selected entry. The remaining gadget vertices form one transitive path, so the correction is exactly K, with no extra factor or unbounded tail term. No path can use both reverse arcs, because both leave the same vertex r.

A path ending at r cannot use either reverse arc. It contributes h times the unique gadget path ending there, except that a singleton gadget forbids old endpoints a,b. This proves e'.

For R', remove r from its reserved starting path. If the reserved path never enters T, the old vertices form one Hamilton path in the other slot, giving Qh. Otherwise its first old vertex is a or b; the old two-slot cover contributes R_a or R_b, while G−r is distributed between the two ordered slots in B ways. This proves the third row without an inclusion-exclusion assumption.

The gluing is physically realizable with two changing ports. Begin with ports a,b, add r beating both, retain r, and stream the tail in transitive order. Each tail vertex beats only r for bit 0 and beats no old vertex for bit 1. Keep r exposed and forget each newborn tail exposure. No forgotten old vertex is re-exposed.

If T is strong C1 and a,b each have s=1, the new root has s'=2. For nonempty tail, the result has exactly two cycles. For a singleton gadget its cycle count is

    2−p_ab−p_ba,                                    (4)

where p_ab counts Hamilton paths of T starting at a and ending at b. In a C1 graph each such starting path is unique. Formula (4) is 1 for adjacent selected vertices on the unique cycle and 2 for a nonadjacent selected pair. The output is always strong. Every subsequent ordinary move (2) preserves s=2 and has exactly two cycles.

## 3. A same-graph switch between start count one and two

Suppose T is strong C1 and a has sole out-neighbor b. Add f beating exactly a,b. This is literally replacement of a by a transitive pair f→a. Starting from ports [a,b], the move is mask 3, drop 1, leaving [a,f]. Both clones remain exposed.

Put E=e_a+e_b and use the old individual responses. The graph has

    h_* = h−e_b+R_b.                                 (5)

Its first clone f and last clone a have respective profiles

    first: (h_*, h−e_a−e_b, h+R_a+R_b),  start count 2;
    last:  (h_*, h−e_b,     h+R_b),      start count 1. (6)

Thus changing the distinguished root between two already exposed vertices changes (e,R) by (−e_a,+R_a), without changing the graph, vertex count, or H. This is a real interface operation, not an assumed equal-count reconstruction.

To verify (5)–(6), remove a and call the remaining graph V. Since a points only to b,

    R_a(T)=H(V)+R_b(V),
    e_a(T)=H(V)−e_b(V),
    K_ab(T)=R_b(V),  K_ba(T)=e_b(V).

Therefore K_ab+K_ba=R_a−e_a, and (3) with L=1 gives the first profile. A single-entry birth beating only b gives the very same unmarked transitive-pair graph, but with its new root as the last clone; the explicit isomorphism exchanges that newborn with the old vertex a. Formula (2) then gives the last profile. The first clone has two starting paths and the last has one, by placing the other clone either immediately after the first or at the path's end; equivalently apply the exact initial-component start-count formula.

Applying the same subsequent gadget to the two possible roots gives counts differing by

    H_first−H_last = t e_a+L R_a >0.                 (7)

The two outputs need not be the same graph. Equation (7) says exactly what a physical root choice changes; it is not a scalar-only affine digit rule.

## 4. The nonadjacent pair is also a two-birth local construction

The proved C1 structure says its start-one vertices form a block of two or three consecutive vertices, except that all three triangle vertices are start-one. The only nonadjacent pair consists of the first and last vertices u,c in a three-vertex block u→b→c, where u and b have outdegree one.

Deleting u leaves a strong C1 graph S, with b still pointing only to c. Indeed the unique cycle shortcuts past u: its predecessor points to b because b's sole out-neighbor is c and the order is at least four. The shorter graph is strong and cannot gain cycles under this strong-induced deletion.

Conversely, start with such S and ports [b,c]. Perform exactly

    mask 3, drop 1: add f beating b,c; ports [b,f];
    mask 2, drop 1: add z beating only f; ports [b,z].    (8)

The first birth clones b. The second uses the first clone as its sole entry but retains the last clone b as the distinguished root.

Compare with the canonical nonadjacent-pair graph: first add u beating only b to S, then add r beating u,c. The explicit graph isomorphism sends old b to f, u to z, r to old b, and fixes every other S vertex. This checks every incidence, including the directed triangle on the three local vertices. In particular the desired new root r is the LAST clone b in (8), not the first clone f. Confusing these roots gives a false response identity.

The scalar count of (8) is

    H = h+e_b+R_b+2R_c,                              (9)

measured on S. This follows by applying (2) with L=1 to the first-clone profile in (6); H does not depend on which surviving vertex is then marked. For further gluing at the desired root, its full response must instead be taken at the last clone b, or computed by (3) on the canonical source S+u.

For S=C3 with rows [2,4,1] and (b,c)=(0,1), (8) gives

    rows [18,20,25,3,8], desired root 0,
    (H,e,R)=(13,2,18).

The other old clone, vertex 3, has profile (13,3,15). Both have start count 2. This explicit difference is included as a hostile control against accidental rerooting.

## 5. Returning tail operations and the formal inverse boundary

For any rooted graph with state (h,e,R,d), where d=H(T−a), the following actual births are closed on these four numbers:

- J: a newborn beats only the current root, and becomes the root:

      (h−e+R, h−e, h+R, h).

- Z: the same newborn beats only the root, but the OLD root is retained:

      (h−e+R, d, 2R−d, d).

- O: a newborn loses to every old vertex, and the old root is retained:

      (h, 0, 2R, d).

The Z and O identities are retained-root specializations of path-cover gluing. Directly, in Z a path ending at the old root ends with the newborn and then that root, giving d. A reserved root-start path cover places the newborn last in either slot; the forbidden singleton-root case removes d from the two copies of R. O simply appends its sink to either slot. The J identity is (2) at L=1 with the deletion response added.

These operations can keep both old and new candidate roots exposed when desired. However the four numerical responses do not certify that an arbitrary formal integer state is physically reachable.

For example, in (a,b,c)=(e,h−e,R−h), J acts by the unimodular map

    (a,b,c) -> (b,a+b+c,a).

Its nonnegative formal inverse is (c,a,b−a−c). The exact physical C1 state (5,1,7), obtained from a singleton using L=3, formally inverts to (3,2,4). This is not a rooted tournament profile at all. A graph with H=3 consists of one directed-triangle strong component and singleton components: every nontrivial strong component has H≥3, and Busch's bound excludes a strong H=3 component of order above three. Its endpoint counts are consequently only 0,1, or 3, never 2. The inverse is therefore an arithmetic aid, not a proved physical descent.

## 6. Complete low-cycle reduction and a finite inverse criterion

The structural dependencies prove that every strong C2 graph except the regular five-tournament has a start-two root. Deleting that root gives either one start-two entry or two start-one entries in an initial strong component, followed by a transitive tail. Following the single-entry cases backwards must eventually reach the two-entry case. Its source is a strong C1 graph.

Define ℒ={1} union {H(T): T is strong and has one or two Hamilton cycles}. Consequently ℒ is exactly the union of:

1. the complete C1 one-seam image;
2. the regular-five value 15;
3. all states obtained by (3) from an actual C1 source with a pair of start-one vertices, followed by any finite word of ordinary moves (2), taking their H coordinates.

The pair-source alternatives are exhausted by Sections 3–4. The graph compiler supplies the required exposed pair before the event. Every operation has an explicit physical meaning and preserves the asserted root; this is not merely a formal state semigroup inclusion.

There is also an exact finite test for any specified target bound N. The initial strong C1 source has h≤N and order at most 1+floor(log_5 N^3), by induced Hamilton-path monotonicity and Busch's strong-order bound. Enumerate those sources using the complete C1 grammar and take their actual pair responses.

For a pair gadget with L≥2, K≤R and

    g h−t E≥h, hence h'≥h+(L−1)R.                   (10)

Here E≤2h; if t=1 and L≥3 then g≥3, while t=0 is immediate. Thus only L=1 and

    2≤L≤1+floor((N−h)/R)

can yield a target at most N. Record every such pair output whose H≤N.

For its subsequent ordinary moves, h'≥L R and R'≥h+R. Therefore a preterminal state on a route to H≤N must have R≤N, and its next gadget satisfies L≤floor(N/R). Keep only states with h≤R≤N and 0≤e≤h, but record terminal H values before discarding a successor with R>N. These finitely many states give an exact arbitrary-depth closure. This criterion can reconstruct physical witnesses by retaining predecessor words.

The source enumeration, pair response and positive-depth continuation all have proven finite bounds. What is missing is a uniform theorem that every desired odd N admits such a physical predecessor. Neither the cone, the root switch, nor a finite interval of coverage proves that assertion.

## 7. Dependencies, verification and scope

The structural dependencies are session197-lowcycle/{C1_COMPLETE.md,START_TWO.md,C2_COMPLETE.md}, and the exact physical pair-exposure compiler is session197/lowcycle-two-port/PROOF.md. The one-seam recurrence is inherited from ROOTED.md. Their separately frozen snapshots are not changed here.

The standalone verifier checks the pair transfer on actual source graphs and all tested gadgets; independently enumerates Hamilton words for small outputs; tests both clone-root profiles and the two-birth nonadjacent graph isomorphism; and checks all three four-coordinate birth identities. The full all-order statements follow from the displayed bijections and structural induction, not from those controls.

A separate finite-coverage appendix may use these formulas to provide explicit recipes at specified targets. No universal ternary or binary digit grammar, arbitrary scalar rebasing permission, or full {7,21} theorem is asserted here.

### Frozen dependency hashes

- `session197-lowcycle/C1_COMPLETE.md`: `d38a53b1c475df3578d58e83d9e9f97785d9455084c0a094aaed752577d8b289`.
- `session197-lowcycle/START_TWO.md`: `a746d2e33b0fb78a98252afa602dc8985c25f718225c2f9ce37923fd20b0ada3`.
- `session197-lowcycle/C2_COMPLETE.md`: `4cc0bf11efc435d7a2746daf6d8f1a58c49c94e57e608911b24464facdd1057d`.
- `session197-lowcycle/ROOTED.md`: `ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e`.
- `session197/lowcycle-two-port/PROOF.md`: `c4a16105ebc9f9aa9cb3cd61887ec4c437db58c361b2469af630899d536c5c73`.
