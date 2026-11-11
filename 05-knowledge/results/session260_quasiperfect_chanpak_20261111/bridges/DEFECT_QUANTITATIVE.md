# Quantitative defect: a positive original-arc cut and the completed injection

**Status: PROVED marked-cover decomposition, original-arc cut law, reverse-capacity circulation, selective fibre bound, and spectator transport. The two formerly open global quantitative bounds are now PROVED by the subsequent exchange theorem linked below; the original finite evidence remains separately labelled.**

**Forward status repair, 11 November 2026:**
[the closed-tree exchange theorem](../../multiplicand_saturation_chanpak_20261111/exchange/STARS.md)
proves (Q) for every tournament and cut by a simultaneous injection of
all old Hamilton paths and all but at most one dominating cut. Consequently
`D≥H(U)−1≥2t_3(U)`. The historical calculations and former missing step
below are retained as provenance; neither bound remains conjectural.

## 1. The sharper target, now proved

Use BIPARTITE_DEFECT.md's notation: U=T-z, A=N_T^+(z), B=N_T^-(z), and G_z has cross edge ab iff a in A, b in B and a->b in U. Write

    D=H(T)-d(G_z),  Delta=H(T)-H(U).

The stronger inequality, now proved by the linked exchange theorem, is

    D >= H(U)-1,

or equivalently

    d(G_z)-1 <= Delta.                                    (Q)

The second proved bound is D>=2t_3(U). It follows from (Q), because OCF gives H(U)-1>=2t_3(U). Neither follows merely from the earlier nonnegative cut fibres: those fibres prove D>=0, not (Q).

The original receipt records exhaustive rooted tournaments through order6 and seeded random rooted tournaments of orders7 through10, without a counterexample. That receipt is finite evidence; the subsequent exchange argument supplies the all-order proof. Both bounds are sharp already on a cyclic example: the2K2 input with compatible choices of its two shore orders has H(T)=13, d(G_z)=9, H(U)=5, and t_3(U)=2, so both bounds are equalities.

## 2. An exact subtraction-free formula for Delta

Fix an old arc u->v of U. Let

- L_uv be the number of Hamilton paths of U using u,v consecutively;
- R_uv be the number of Hamilton paths gained when that single old arc is reversed, namely paths of the flipped tournament using v,u consecutively;
- K_uv=L_uv+R_uv.

These are counts of literal objects. R_uv is equivalently the number of ordered spanning two-path covers whose first path ends at v and second path starts at u. Then for every A subset V(U),

    H(U+z_A)-H(U)
      = sum_{u in A,v outside A,u->v in U} K_uv.            (3)

Thus the extension increment is a positive directed cut supported entirely on the original arcs. The old raw OCF kernel has positive entries on some nonarcs; (3) does not drop those contributions, but redistributes them exactly as explained in section3.

### Proof by erasing the newborn

A Hamilton path of U+z falls into two classes.

**Deletion remains a Hamilton path.** Fix an old Hamilton path P and read its A/B colour word. The allowed insertion slots for z are its B-to-A transitions, its beginning if P starts in A, and its end if P ends in B. For every nonempty binary word,

    #BA + [starts A] + [ends B] = 1 + #AB.

Consequently these new paths number H(U) plus the number of old Hamilton paths with a marked A-to-B adjacency. The latter is sum L_uv over the cut's forward old arcs.

**Deletion is not a Hamilton path.** Then z is internal, its predecessor is v in B, its successor u in A, and u->v is an old arc. Replace the substring v,z,u by v,u and reverse just u->v in U. This gives a Hamilton path counted by R_uv. The operation is bijective, including its marked old arc. Such paths contribute sum R_uv.

Adding the two classes proves (3). For empty U the statement is trivially0=0.

### Even positive weights

The arc flip satisfies

    H(U^uv)-H(U)=R_uv-L_uv,

which is even by Redei parity. Hence K_uv is even. In fact K_uv>=2, as the next section identifies its half with a sum of positive path capacities containing the two-vertex path u->v.

## 3. OCF capacities and a previously hidden circulation

Recover Session258's exact raw endpoint kernel. For distinct u,v let

    c_uv = sum_P H(U-V(P)),

where P runs over directed simple u-to-v paths having an even number of vertices; H(empty)=1. Splitting the unique odd cycle through a newborn z gives

    Delta(A)=2 sum_{u in A,v outside A} c_uv.               (4)

This is inherited from `session258_endpoint_seed_20261111/endpoints/PROOF.md`, not asserted as a new OCF identity. Compare (4) with (3), as multilinear quadratic polynomials in the indicator variables x_u=[u in A]. Their squarefree x_u x_v coefficients give, for every old arc u->v,

    K_uv=2(c_uv+c_vu).                                     (5)

The degree-one coefficient at u then gives

    sum_{v:u->v} c_vu = sum_{v:v->u} c_uv.                  (6)

Define a nonnegative flow on each original arc u->v by

    r(u,v)=c_vu.

Equation (6) says this flow is **a circulation**: its outgoing and incoming masses agree at every vertex. Thus the raw capacity of paths going opposite the original arc orientation can be transported across any cut without changing its total. This gives a direct explanation of the resummation:

    raw forward-crossing c
      = forward-arc c + backward-arc reverse-flow
      = forward-arc (c + reverse-flow).

This is a precise positive quotient, not permission to delete nonarc capacities individually. In particular, the raw c_uv and the resummed K_uv/2 are different endpoint statistics.

Because c_uv>=H(U-{u,v})>=1 for an old arc u->v, (5) proves K_uv>=2. The circulation can be decomposed into directed cycles by the elementary finite-flow decomposition. It vanishes on transitive U, where all even paths follow the total order.

This structural formula offers a useful new interface for exact extension searches: use the genuine old-arc capacities K/2 rather than a dense raw kernel, while preserving every newborn count. It does not by itself show which cut weight a prescribed target realizes.

### Transitive specialization recovers the chain/apex word law

If U is the transitive tournament on positions1,...,n and u<v, then c_vu=0. Put t=v-u. For t=1 the only even-vertex path from u to v is the two-vertex path. For t>=2 the admissible paths choose an even number of the t-1 intermediate vertices, giving2^(t-2) choices. Every complementary subtournament is transitive and has H=1. Consequently

    K_uv=2 for t=1, and K_uv=2^(t-1) for t>=2.

Since H(U)=1, formula (3) is exactly the weighted-pair formula of DOMINATION_TRANSFER.md. Thus the domination adapter is the zero-circulation, transitive specialization of the universal positive arc-cut law, not an unrelated numerical coincidence.

## 4. A proved quantitative bound from selected dominating fibres

The positive cut-fibre formula of BIPARTITE_DEFECT.md yields a smaller but explicit bound without assuming (Q). Retain the fibre of S=V and the fibres of S=V minus one nonisolated graph vertex. Their contributions are disjoint. With endpoint counts e_B and s_A as in that note,

    D >= H(U[A])H(U[B])-1
         + sum_{b in B,deg_G(b)>0}
             (H(U[B-b]) s_A(A union {b})-1)
         + sum_{a in A,deg_G(a)>0}
             (e_B(B union {a}) H(U[A-a])-1).                (7)

Every displayed summand is nonnegative. The all-selected set is always dominating; a set missing one vertex is dominating exactly when that vertex is nonisolated. Formula (7) therefore follows by discarding other nonnegative fibres, not by estimating a signed expression.

This bound records internal-shore complexity and multiple mixed-triangle witnesses simultaneously. It is often stronger than the qualitative gap2, while not being claimed to dominate H(U)-1. Computing it needs the indicated endpoint sidecars rather than another complete tournament census.

## 5. Physical spectator transport

There is an exact lawful way to amplify, rather than cancel, defect. Put an arbitrary tournament W before the entire rooted T in an ordinal sum (the dual operation after T is identical). With k=H(W), the new marked graph adds |W| isolated vertices to G_z, so its domination count stays d(G_z). Meanwhile

    H(T_new)=kH(T),   H(U_new)=kH(U),
    D_new=kD+(k-1)d(G_z),
    Delta_new=kDelta.                                     (8)

In particular the proposed slack Q=Delta-(d-1) transforms as

    Q_new=kQ+(k-1)(d-1).

Thus (Q) is stable under these spectator products whenever it holds initially. This is a real graph operation with a measured change, but it moves the defect in the wrong direction for target realization. Transitive spectators have k=1 and change nothing.

## 6. The former injection problem, now closed

The marked-cover proof turns (Q) into the concrete question

    |Dom(G_z) minus {one distinguished set}|
       <= #marked AB old paths + #reverse-seam repair paths. (9)

Equivalently, seek a simultaneous injection of all old Hamilton paths and all but one dominating set into Hamilton paths of T. A plausible route is to reserve one legal z-insertion for each old path and reroute the dominating-cut representatives around those reservations via alternating exchanges.

The obstruction to simply asserting this is precise: different dominating cuts can use insertions of different old paths, and a naive canonical choice does not rule out consuming all legal slots of more than one old path. The nonnegative fibre proof does not supply the requisite global matching or exchange lemma. Neither a dominating-set injection into (9), nor a Hall condition controlling the reserved paths, was obtained in this focused pass.

This was the stopping point of the original pass. The subsequent exchange
proof now closes it: the literal path/cut incidence graph has at most one
closed tree component, every such tree is a star, and componentwise matching
gives exactly the simultaneous injection required by (9). Thus (Q) and
its triangle-count consequence are **PROVED**. The positive original-arc
formula, circulation explanation, selected-fibre bound, and transport law
remain valid independently; the original finite receipt is unchanged.

## 7. Reproduction

Run `python3 05-knowledge/results/session260_quasiperfect_chanpak_20261111/bridges/defect_quantitative_check.py`. The default exhaustive range is orders1 through6. `--max-order 5` shortens it. The script checks both inequalities, separately checks every old-arc weight and every newborn cut on bases through order4, and adds200 seeded random tournaments at each order7 through10. The unchanged historical receipt labels the inequalities as conjectural; the forward proof status above supersedes that label.
