# The marked mixed-cusp covering has mapping degree at least sixteen

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This is a necessary degree bound for the actual marked mixed-cusp Keller-support problem. It does not assert a representation exists at degree sixteen, or exclude the entire curve family or JC(2).

## 1. Actual supplier and the stronger consequence of small-cycle closure

Retain the actual geometric hypotheses of [odd_retention, Section 1](planar_jc48_sep08_odd_retention.md): an irreducible whole nonproper support with normalization A1, three finite unibranch cusps of types (2,3),(2,3),(2,5), at least three nodes, and the proved simultaneous marked H4 meridians. The [actual mixed-cusp supplier](planar_jc48_sep08_mixed_cusp_braid.md) and its [good-family transport](planar_jc48_sep08_three_cusp_boundary.md) pay these hypotheses for their declared family. The meridians a,b,c,d generate the actual transitive D-sheet monodromy action and satisfy

    aba=bab,  bcb=cbc,  cdcdc=dcdcd,
    ac=ca,  ad=da,  bd=db.

There is an actual positive retained-sheet count k along the smooth support. All four meridians have the same moved count t, and 1<=k<=D-t. This is the actual retained set; no full-fixed-point assumption is imposed initially.

**Theorem.** Under these hypotheses, D>=16.

The closest proved mechanism is [odd_retention](planar_jc48_sep08_odd_retention.md), which pays the arbitrary-subset inequalities and an earlier floor twelve. The new input is the complete closure of moved supports through six, using the proved [single-cycle](planar_jc48_sep08_h4_single_cycle.md), [involution](planar_jc48_sep08_h4_involution.md), [mixed(3)(2)](planar_jc48_sep08_h4_mixed32.md), [mixed(3)(3)](planar_jc48_sep08_h4_mixed33.md), and [mixed(4)(2)](planar_jc48_sep08_h4_mixed42.md) results. The last of these has passed independent audit and is promoted; no RESERVED dependency enters this candidate.

Indeed all nontrivial cycle partitions moving at most six letters are

    (2), (3), (4), (2,2), (5), (3,2),
    (6), (4,2), (3,3), (2,2,2).

The listed dependencies exclude each in the actual covering. A nonidentity permutation cannot move exactly one letter. If all meridians were identity they could not generate a transitive action of degree D>=2. Hence t>=7.

Live concepts are actual retained subsets, typed moved supports, partial versus full retention, complementary supports, and the Euler total. The map from a monodromy action to scalar counts loses the common incidence structure of its four moved supports. The decisive sidecar is precisely that structure, restored at the three marked commuting nodes. The canonical hostile is the nonempty scalar frontier retained at the end; scalar feasibility is not a realized monodromy passport. The corrected near miss is to replace a partially retained set by the entire fixed set. The D=15,k=7 case below explicitly pays its two one-sheet deficits.

## 2. All necessary scalar inequalities, with their direction retained

Let S_a,S_b,S_c,S_d be the four moved supports, each of size t. Let n_ab,n_bc,n_cd denote the three actual retained cusp overlaps, and let W be the sum of all actual node overlaps, including any additional nodes. The proved Euler and odd-retention statements give

    1=-2k+n_ab+n_bc+n_cd+W,
    n_ab,n_bc >= max(0,k-floor(t/2),ceil((3k-D)/2)),
    n_cd >= max(0,k-floor(2t/3),ceil((5k-2D)/3)),
    W >= 3 max(0,D-2k).                              (1)

The last bound uses only the three marked nodes; extra nodes add nonnegative terms. Every marked node overlap is also at least the intersection of its two moved supports. If two meridians commute, that support intersection is invariant under them and a union of nontrivial cycles. Thus a positive such intersection has at least two letters; an intersection of size one is impossible.

The two ordinary edges have

    |S_a intersect S_b|, |S_b intersect S_c| >= ceil(t/2). (2)

These are the proved ordinary support inequalities, not an extrapolation from fifth-cusp pairs.

If k=D-t, retention is full. Then each cusp overlap is
D-2t plus its corresponding moved-support intersection, and every marked node overlap is exactly its moved-support intersection. All six pairs occur among the three cusps and three marked nodes. The inherited incidence argument gives

    1 >= D-4t+sum_(i<j)|S_i intersect S_j| >= 4t-2D,

hence D>=2t. To recall the last step, a letter occurring in m of the four supports contributes binomial(m,2) to the pair sum, and binomial(m,2)>=2m-3 for m=0,1,2,3,4. Summing uses total support incidence 4t. This proof does not need transitivity to remove unused letters.

For partial retention, write delta=(D-t)-k. The two actual retained sets at an ordinary cusp are subsets of their fixed sets with at most delta omitted letters each. Therefore, separately from (1),

    n_ab >= D-2t+|S_a intersect S_b|-2delta.          (3)

This bound is true for any two such retained subsets in the paid local marking. It does not require choosing globally consistent retained subsets at all nodes, nor imposing full retention.

## 3. Complete finite scalar reduction below sixteen

Enumerate every integer D=2..15, t=7..D-1 and k=1..D-t using the necessary inequalities (1). There are no surviving rows below D=12. The complete remaining list is:

| D | t | k | Minimum n_ab,n_bc | Minimum n_cd | Maximum W |
|---:|---:|---:|---:|---:|---:|
|12|7|5|2|1|6|
|13|7|6|3|2|5|
|13|8|5|1|0|9|
|14|7|7|4|3|4|
|14|8|6|2|1|8|
|15|7|7|4|3|4|
|15|7|8|5|4|3|
|15|8|7|3|2|7|
|15|9|6|2|0|9|

The source gives an exact reproducible enumeration, not a numerical test of infinitely many degrees. Every degree being excluded is explicitly in this finite universe. These are necessary scalar rows only.

All rows except the D=15,t=7,k=7 row have full retention. Applying D>=2t eliminates every row except

    (D,t,k)=(14,7,7), (15,7,7), (15,7,8).           (4)

The following actual support arguments eliminate all three. No later small-cycle classification is assumed.

## 4. Degree fourteen: the complementary supports are incompatible

Take D=14,t=k=7. By (2), S_a and S_c each meet the seven-element set S_b in at least four points. They therefore meet each other. Since a,c commute, their intersection has size at least two. Under full retention this is also the actual marked node overlap omega_ac.

The table gives W<=4. The other two marked commuting intersections cannot both be positive: each positive one has size at least two, which together with omega_ac>=2 would make W>=6. Hence at least one of S_a intersect S_d or S_b intersect S_d is empty.

Write ab,bc,cd,ac,ad,bd for the six support intersection sizes. If ad=0, the two seven-sets S_a,S_d in fourteen letters are complementary. Thus bd=7-ab and cd=7-ac. The contribution from the three cusps and three marked nodes to the Euler expression is exactly

    -14+ab+bc+cd+ac+ad+bd = bc >=4.

Additional nodes only increase it, contradicting Euler value one.

If instead bd=0, S_b,S_d are complementary. Then ad=7-ab and cd=7-bc. The same expression equals ac>=2, again a contradiction. This exhausts the degree-fourteen row without inferring actual permutations from a support diagram.

## 5. Degree fifteen with one-sheet deficits

Take D=15,t=k=7. Each fixed set has eight elements, so delta=1. As before, the two ordinary overlaps force |S_a intersect S_c|>=2. Hence omega_ac>=2. Each other marked node has omega>=D-2k=1, while W<=4. The only possible allocation for the three marked nodes is

    (omega_ac,omega_ad,omega_bd)=(2,1,1).

In particular their moved-support intersections ad and bd are at most one. Commutation forbids a positive singleton, so both are zero. S_d is therefore disjoint from S_a union S_b. Since S_d has seven letters, their union has size at most eight in the fifteen-letter universe, yielding ab>=6.

The common fixed set of a,b has size D-2t+ab>=7. Each of the two actual retained sets omits only one letter from its full fixed set. Formula (3) gives n_ab>=5. The other cusp bounds in the row remain n_bc>=4 and n_cd>=3. The three marked node overlaps already sum to four. Thus

    -2k+n_ab+n_bc+n_cd+W >= -14+5+4+3+4=2,

contradicting Euler value one. This step retains the missing sheets explicitly; it does not identify the actual retained set with all eight fixed sheets.

## 6. Degree fifteen with full retention

Take D=15,t=7,k=8. Again ac>=2, and now actual marked node overlaps equal the moved-support intersections. Since W<=3, neither ad nor bd can be positive: either would contribute at least two in addition to ac. Consequently ad=bd=0, so the same union bound gives ab>=6.

Full retention gives n_ab=D-2t+ab>=7. The other cusp bounds are n_bc>=5 and n_cd>=4, and W>=ac>=2. Therefore

    -2k+n_ab+n_bc+n_cd+W >= -16+7+5+4+2=2,

the final contradiction. Equations (4) are exhausted, proving D>=16.

The use of at least three nodes is sufficient throughout. Every argument only lower-bounds the sum of their overlaps, and additional nodes contribute nonnegatively. No equality between the actual curve complement group and the abstract H4 presentation is needed.

## 7. Scope, controls and the next frontier

The degree bound is a necessary property of the actual marked covering with positive retained count and the stated Euler ledger. It is not a property of arbitrary abstract H4 actions without those sidecars. The proof does not claim a degree-sixteen realization.

The finite source preserves four scalar rows at degree sixteen:

    (t,k,n3,n5,Wmax)=(7,8,5,4,3),(7,9,6,5,2),
                     (8,7,3,2,7),(8,8,4,3,6).

For example the last row satisfies the formal Euler value -16+4+4+3+6=1. This is only a scalar hostile to a claim that the inequalities alone eliminate the next degree. It is not an actual set configuration, permutation quadruple, or Keller map. The remaining sidecar is the full cycle/incidence structure at these rows or a new global restriction.

The standalone source imports no inherited implementation. It lists every nonidentity moved cycle type through six, retains all integer triples in the declared degree range, checks every table row and full-retention deletion, verifies the support-complement identities on their entire stated integer intervals, and retains the partial-retention deficit arithmetic and the next-degree scalar controls. Its always-active checks use explicit exceptions.

Reproduce from the worktree root:

    python3 04-computation/planar_jc48_sep08_h4_degree16.py
    python3 -O 04-computation/planar_jc48_sep08_h4_degree16.py

Both modes pass 561 gates and reproduce the frozen output. The [independent complete audit](planar_jc48_sep08_h4_degree16_audit.md)
accepts the actual marked maps, every retention case, both561-gate
replays and an independently enumerated full scalar head. The analytic support arguments in Sections 4–6, rather than the existence of a scalar table, prove that the remaining rows are impossible.
