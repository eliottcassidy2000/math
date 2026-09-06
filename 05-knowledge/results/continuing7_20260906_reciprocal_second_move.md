# A second reciprocal transposition: packet obstruction and an infinite failure family

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** After the prescribed exchange of outputs0 and1 in the completed reciprocal permutation, two or more native P-line packets prevent every second ordinary output transposition from achieving integer no-three-in-line. Even one packet need not be repairable: every prime p=37 mod360 admits no successful second transposition. At p43, in contrast, exactly two second moves work. These statements concern a p-point permutation graph, with no intervening affine relabeling, and do not claim a 2p-point construction or a classification of all first moves.

[Independent proof and exact audit](continuing7_20260906_reciprocal_second_move_audit.md) passes.

## 1. Inheritance, target and retained coordinates

The closest proved mechanism is [continuing6_20260906_reciprocal_native.md](continuing6_20260906_reciprocal_native.md), Sections2--5, and its independent audit. Its integer-line packet theorem gives the complete actual triple set after the prescribed first move. The earlier [continuing5_20260906_wildcard_swap_compiler.md](continuing5_20260906_wildcard_swap_compiler.md) distinguishes completed affine moves from bare transpositions; its successful completed moves do not imply that their intermediate bare swaps improve the board.

The canonical hostile is p37, the first prime with character(5/p)=-1 whose zero-one repair fails. The corrected near miss here is the proposed second move involving the retained fixed point p-1: that point belongs to no old bad triple and cannot help a two-label move cover a transposed pair. The least-used sidecar is the actual output support of every old triangle, together with the involution pairing between two lines' endpoints.

Anchor: the unchanged LRC14 frontier remains open; root and another agent handle its native location and compatible-owner calculations. Niche: the moments lane retains original evaluated signs on shared-root fibers. This wildcard tests a deterministic integer-board operation. The live concept board is:

| Source | Map / operation | Preserved target information | Information needing a sidecar |
|---|---|---|---|
| Native reciprocal line packet | List its triangle's output support | Which old triangles survive a swap | Incidence between distinct packets |
| Involutive board | Transpose coordinates | Actual integer collinearity | The conjugated output move g tau g |
| Two candidate moved labels | Hit every old support | Necessary removal of old triangles | Newly created triangles |
| Retained reciprocal points | Separate off-conic anchors | Complete native triple set | Primitive integer directions, not modular slopes |
| Fixed p37 packet | Vary p in one congruence family | Explicit obstructing coordinates | Native lift and unchanged-output inequalities |

No random conditional-column law is used for these selected moves. The fixed-row and restart results require their stated uniform sampling hypotheses.

## 2. Board and inherited packet geometry

Let p be an odd prime. Write g(0)=1, g(1)=0, and g(x)=[x^(-1)]_p for2<=x<=p-1. The initial board is

    B_p={(x,g(x)):0<=x<p},
    P=(0,1), Q=(1,0), F=(p-1,p-1).

It is an involution graph. A second move exchanges two distinct output labels u,v; the resulting graph has function tau g, where tau=(u v). All coordinates are the standard integers0,...,p-1. Every collinearity below means determinant zero over the integers.

By the inherited theorem, all bad triples occur in transposed pairs

    {P,R_i,S_i}, {Q,R_i^T,S_i^T},   1<=i<=N_p.       (1)

Here R_i,S_i are the two untouched reciprocal points on a primitive line ell_i throughP. Put X_i={x(R_i),x(S_i)} and Y_i={y(R_i),y(S_i)}. Thus their output supports are

    A_i={1} union Y_i,     B_i={0} union X_i.          (2)

The X_i are pairwise disjoint, as are the Y_i: distinct lines throughP share no other point, and the graph is a permutation. Also X_i and Y_i are disjoint for the same i. Indeed a common coordinate would either produce a reciprocal point fixed by transposition, or make the two points R_i,S_i transposes of one another. The latter would force ell_i to equal its transpose, hence containP andQ, which is impossible for positive reciprocal points on the line x+y=1. The only fixed reciprocal points modulo an odd prime are(1,1) andF; the first is absent and the second cannot lie in (1).

For clarity, the last exclusion follows directly from the packet coordinates. They are (bk,ak+1),(bl,al+1), where ah=p-1,2<=k<l,k+l=h. Consequently their y coordinates are at most p-2a<=p-2. A point with x=p-1 would also have y=p-1, so neither coordinate can equal p-1.

For i!=j there is at most one label in X_i intersect Y_j. Such a label v specifies the point(v,g(v)) in ell_i intersect ell_j^T. These two integer lines are distinct: equality would putP,Q and a positive reciprocal point on x+y=1. Finally inversion gives the bijection

    g : X_i intersect Y_j -> Y_i intersect X_j.       (3)

These are actual line-incidence facts, stronger than treating the triangle supports as arbitrary sets.

## 3. Universal two-packet obstruction

**Theorem.** If N_p>=2, no second output transposition makes B_p integer no-three-in-line.

First, a successful move must meet every support in (2). A disjoint old support leaves all three of its points unchanged. If the moved pair contains both0 and1, it undoes the first swap and restores the actual diagonal triple(0,0),(1,1),F. If it contains exactly one of0,1, its other label must hit every one of the disjoint two-element sets belonging to the opposite anchor. This is impossible when N_p>=2.

It remains that neither0 nor1 moves. With N_p>=3, two labels cannot meet all pairwise-disjoint Y_i. With N_p=2, the moved labels must lie one in each Y_i and one in each X_i. The within-packet disjointness therefore forces them into Y_1 intersect X_2 and Y_2 intersect X_1. Each intersection has at most one point, and (3) pairs them by g. Thus, if any such two-label cover exists, it is a reciprocal pair {v,g(v)}.

Exchanging these outputs replaces the two transposed reciprocal points by the fixed points(v,v) and(g(v),g(v)). They are distinct and neither isF. The unchangedF then creates a new actual diagonal triple. Every possible second move has been excluded. This proves the theorem.

The same support argument gives a separate useful fact: whenever N_p>=1, no second transposition involving output p-1 can repair B_p. This label belongs to neither member of any pair of disjoint supports A_i,B_i, while the one other moved label cannot meet both. The plausible fixed-point repair direction is therefore blocked before any search.

These results make N_p=1 a necessary condition for repairing a failed B_p by this operation. It is not sufficient, as the following unbounded family proves.

## 4. An infinite single-packet-based obstruction

**Theorem.** For every prime p=37 mod360, no second ordinary output transposition repairs B_p. This claim allows additional native packets; it does not assume N_p=1 throughout the progression.

Set

    X1=3(p+3)/5, Y1=(4p+5)/9,
    X2=3(p+3)/4, Y2=(5p+4)/9.

The inherited native template supplies the two actual triples P,(X1,Y1),(X2,Y2) and its transpose. For p=37+360r,r>=0 their coordinates are integers and

    0<1<Y1<Y2<X1<X2<p-1.                             (4)

Any successful second move must join a label of {1,Y1,Y2} to a label of {0,X1,X2}. These are exactly nine possible unordered pairs. The move(1,0) restores the original diagonal; the moves(Y1,X1),(Y2,X2) create the two reciprocal fixed points and retainF, so they also fail.

The other six moves reduce to three by transposition. The transpose of the graph tau g has function g tau=(g tau g)g, because g is an involution. Accordingly a move(u,v) is carried to the move(g(u),g(v)); this is a bijection preserving the actual no-three predicate. The following table gives a newly created integer triangle for one representative of each remaining pair:

| Swapped output labels | Moved point | Two unchanged reciprocal points |
|---|---|---|
| (1,X1) | (Y1,1) | ((p+1)/2,2), ((8p+1)/9,9) |
| (1,X2) | (Y2,1) | ((2p+1)/3,3), ((5p+1)/6,6) |
| (Y1,X2) | (X1,X2) | ((3p+4)/5,(3p+5)/4), ((p-5)/4,(2p-4)/5) |

Each displayed determinant is identically zero. Each unchanged point has xy=1 modp, distinct native coordinates in2,...,p-1, and output different from both moved labels. All these assertions follow by substituting p=37+360r and checking affine inequalities and polynomial identities; their coefficients and exact quotients of xy-1 by p are in the certificate.

For the last row all three points lie on the integer line y-x=3(p+3)/20. The moved point is the first unchanged point plus(1,1). The six nontrivial moved pairs consist of these three and their conjugates(0,Y1),(0,Y2),(X1,Y2). Thus all nine necessary covers fail, proving the theorem.

Since gcd(37,360)=1, the progression contains infinitely many primes by **CITED** Dirichlet's theorem. The primary source and exact use are inherited from [Dirichlet's1837 theorem, English translation, arXiv0808.1408v2](https://arxiv.org/abs/0808.1408v2). This external input concerns only infinitude, not the coordinate proof. These primes have character(5/p)=-1 as recorded in the inherited template theorem; the negative character does not supply this stronger second-move obstruction by itself.

## 5. Positive control and the precise stopping boundary

At p37 there is exactly one native packet, with points(24,17),(30,21). The complete set of666 second transpositions contains no success. This is a finite control for the unbounded theorem, not its proof.

At p43 there is exactly one packet, with points(10,13),(25,31). Among all903 second transpositions, precisely(0,13) and(1,10) produce no actual integer triple. They are conjugate under the displayed transpose rule. Thus a second move can work, and N_p=1 is not itself an obstruction. The complete safe permutations are retained in the certificate.

The p43 success does not extend unchanged to its inherited packet progression43 mod70. At the next prime113 there are two native packets, so the universal theorem excludes every second move. The exact packet endpoint pairs are (30,49),(40,65) and(24,33),(60,81). A three-packet control at p197 tests the other branch of the proof. These are prescribed positive/hostile controls, not an expanding prime census.

The missing coordinate is the incidence and creation geometry of moved endpoints. Packet count already closes the N>=2 branch, but N=1 still has nine potential covers and arithmetic constraints on their newly created lines. A uniform successful family in that residual class is open here. No all-prime two-swap construction, arbitrary-first-swap impossibility, or impossibility after further moves follows.

## 6. Complete native decoder and reproducibility

The source independently verifies each finite moved board by two paths. The first enumerates every triple and its literal integer determinant. The second separates points on the modular conic xy=1 from off-conic anchors. No three retained conic points can be collinear. For each one-anchor case it groups the conic points by primitive unoriented integer direction from the anchor; every pair in one group gives a triangle. It then checks the two-anchor/one-conic and three-anchor cases separately. After the two moves there are at most four off-conic points, so this is a complete bounded-anchor decoder, not merely a test of surviving old lines.

Every second transposition is checked by both paths at p37 andp43. At p113 andp197 the source enumerates all candidate label pairs, rejects those retaining an old triple by their exact support, and evaluates every remaining cover by both geometric paths. It separately verifies all parameter identities and inequalities for the infinite progression, including the complete nine-cover list. The certificate records every finite move count, complete triangle lists for the covers, safe permutations, and the affine family witnesses.

The source has no repository imports. It writes its certificate adjacent outside the repository, or into05-knowledge/results when filed under04-computation. Both `python continuing7_20260906_reciprocal_second_move.py` and `python -O continuing7_20260906_reciprocal_second_move.py` pass **30,428 always-active gates**, with byte-identical raw LF output and unchanged JSON.

- Source SHA256: `1e1571440d0ea569b7e3f65ffdc1ebe9ae61b2ae3d1f11cc39a7a6667d4b708c`.
- Output SHA256: `060493fb23ff746151420891fd9cce0df5eb3e0a7ebd7b4e0a56d3a84f2d328c`.
- Certificate SHA256: `f5117590d335c3df7579b8bb483d2a0a9d0906ab67c06bbb1cfd01ff7a579c3f`.

## 7. Cross-thread consequence

In the LRC lane, hitting each bad marginal event is likewise only a necessary removal test: the newly created diagonal triple is the exact hostile to equating coverage with successful repair. Here the retained sidecar is a paired integer-line incidence, rather than a compatible divisor-word owner. In the moments lane, the two surviving reciprocal endpoints play the role of a shared carrier whose original evaluated sign must still be checked after elimination. Neither comparison transfers a theorem between different predicates.

The useful operation lesson is precise: first compute the complete surviving-obstruction cover, then exploit the source symmetry to classify all covers, and finally certify the new target violations. Its counterindication is an unstructured collection of supports without the reciprocal incidence (3); a general two-element hitting set need not force a diagonal trap. The next decisive question is whether an explicit unbounded one-packet class admits a cover whose complete off-conic-anchor decoder is empty.
