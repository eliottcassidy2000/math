# Gold partitions from an arbitrary chain of extra successors

Date: 2026-10-19. Status: proved for arbitrary finite ambient posets, with an invariant-weight extension. A separate six-point example shows why incomparable extra successors cannot simply be substituted while retaining the same endpoint restriction. No literature priority claim is made; the general Gold Partition Conjecture remains open.

## 1. Statement and exact convention

Use the repository's [weak Gold convention](gold_terminal_context_20261019.md): if E=e(P), a strategy is Gold when every realized two-answer path satisfies E≥t1+t2. Determined comparisons are allowed, and a chain branch has surviving count1. The first and second comparisons may be chosen from the known poset and its extension counts before any answers arrive; no efficient counting algorithm is asserted.

Let u,v be incomparable vertices of a finite poset P. Suppose the full ambient strict neighborhoods satisfy

    D(u)=D(v),
    U(u)=U(v)∪{t1,...,tk},
    t1<...<tk, k≥1,

where none of the ti belongs to U(v).

**Chain-surplus theorem.** P has a Gold strategy consisting of comparisons between v and two adjacent vertices of the marked chain

    c0=u<c1=t1<...<ck=tk.

The exact choice and order of those two comparisons are given in Section4. They depend on the extension counts. The theorem holds for arbitrary finite ambient size and arbitrary other relations compatible with the displayed neighborhoods.

More strongly, the conclusion holds under any nonnegative weighting determined solely by the outside skeleton obtained by replacing the k+2 marked labels u,v,t1,...,tk by identical blanks while retaining every occupied position and outside label. For each such positive-total measure the rule chooses a Gold strategy. This quantifier does not claim that one fixed pair and ordering of comparisons works for every possible reweighting.

If k=0, u,v are exact twins and the half-split strategy from [the one-relation twin theorem](one_relation_twin_gold_theorem_20261019.md) applies. Order duality gives the analogous theorem for a chain of extra predecessors and equal successor sets.

## 2. Monotone insertion bins

Each ti is incomparable with v: v<ti would put ti in U(v), and ti<v would imply u<v. Therefore every extension restricts to the chain c0<...<ck with v inserted in one of its k+2 gaps.

Let q_i be the number of extensions having exactly i marked chain vertices before v, for0≤i≤k+1. Thus

    q0: v<u<t1<...<tk,
    q1: u<v<t1<...<tk,
    q_(i+1): u<t1<...<ti<v<t_(i+1)<...<tk.

Then

    q0=q1≥q2≥...≥q_(k+1)≥0.                    (1)

**First equality.** Swap u,v between the first two bins. Their downsets agree; their common successors follow both occupied positions; every exceptional successor ti also follows both. This is a bijection.

**Each later inequality.** For i≥1, take an extension in which ti is the last marked chain vertex before v. Swap the labels ti,v at their occupied positions. The new v is valid earlier because D(v)=D(u), and all these predecessors occur before u. The vertex ti moves later, so its predecessor constraints remain valid. Every successor of ti is either a later member of the extra chain or a member of U(v). Later marked chain vertices were already after v; members of U(v) also follow v. Hence the new ti is valid in the old v position. The map is an injection from bin i+1 to bin i.

All these maps preserve every outside label and the entire set of marked occupied positions. Consequently(1) holds separately on every outside skeleton and under any nonnegative skeleton-invariant weighting. The proof never assumes independent insertion weights.

## 3. A monotone-bin Gold lemma

Let q0,...,q_(m−1) be a nonnegative sequence with

    q0=q1≥q2≥...≥q_(m−1), E=sum(q_i)>0, m≥3.

Suppose a comparison is available for every nontrivial prefix cut, so its increasing answer has mass q0+...+q_j for0≤j≤m−2. Then two adjacent prefix cuts, in a count-dependent order, form a Gold strategy.

If2q0=E, the equal second bin exhausts all remaining mass. Use the first prefix cut, which splits E in half, and any available second cut. Every path is Gold since t2≤t1=E/2. In particular the second adjacent cut is available when m≥3.

Otherwise2q0<E. Let ℓ be the first index satisfying

    2(q0+...+q_ℓ)≥E.

Set

    a=q0+...+q_(ℓ−1), x=q_ℓ, b=E−a−x.

Then a<E/2, b≤E/2, and x≤a. The last inequality follows from monotonicity and the fact that ℓ≥1. Both cuts are genuine nontrivial cuts: ℓ cannot equal m−1, because the prefix before the last bin has mass E−q_(m−1)≥E−q0>E/2. Thus1≤ℓ≤m−2.

- If2a+x≥E, equivalently b≤a, first use the prefix-a cut, then the prefix-(a+x) cut. The heavy first branch has mass x+b; its second cells are x,b, both at most a. The other first branch has mass a≤E/2. The three realized deficits are E−2a, a−x, a−b, all nonnegative.
- If2a+x≤E, equivalently a≤b, reverse those comparisons. The heavy first branch has mass a+x; its second cells are a,x, both at most b. The other first branch has mass b≤E/2. The three deficits are b−a, b−x, E−2b, all nonnegative.

Equality between the two cases is harmless; either ordering works. This proves the lemma with the original weak Gold threshold, including zero bins and half-mass boundary cases.

## 4. Exact comparison indices in the poset

For the insertion law(1), the prefix cut after q_j is precisely the comparison

    v<c_j, 0≤j≤k.

If2q0=E, compare v with c0=u first and v with c1=t1 second.

Otherwise compute ℓ,a,x as in Section3. The two chain vertices are c_(ℓ−1) and c_ℓ, adjacent in the marked chain (not necessarily a Hasse cover in the full poset).

- If2a+x≥E: first compare v with c_(ℓ−1), then v with c_ℓ.
- If2a+x≤E: first compare v with c_ℓ, then v with c_(ℓ−1).

The index bound1≤ℓ≤k guarantees both endpoints exist. All comparisons are fixed before receiving the first answer. Since each c_j is incomparable with v in P, each is a legitimate uncertain comparison initially, although its answer may be determined on one later branch.

This proves the theorem. The first comparison is necessarily1/3-balanced, and every two-answer cell has mass at most E/2. These also follow directly from the Gold inequalities: each first branch is the union of at most two second cells, each no larger than the opposite first branch, and t1≥t2.

## 5. Two-extra-chain specialization

For k=2 the four bins are A,A,B,C with A≥B≥C. Put S=B+C≤2A.

- If S≤A, compare u with v, then v with t1. The joint table, with suitable answer order, is[[A,S],[A,0]], giving deficits0,A−S,S.
- If S≥A, compare v with t1, then u with v. Its table is[[A,A],[S,0]], giving deficits S−A,S−A,2A−S.

No extra assumptions on D(t1) or D(t2) are required. This is the precise extension of the earlier one-extra-successor argument; applying that earlier fixed strategy without changing its first comparison would fail in part of the parameter range.

## 6. Relation to the earlier prefix-side theorem

The [prefix-side-poset theorem](gold_prefix_side_poset_20261019.md) concerns two independent backbone chains of lengths m,n≥2, arbitrary side attachments through one prefix, and a prescribed pair of terminal comparisons. Its universal parameter range is a quadratic inequality. The present theorem instead assumes an exact equality of two full downsets and a chain difference of their full upsets, and chooses adjacent thresholds using the actual insertion counts.

These are complementary hypotheses, not a reformulation or a general extension of that earlier theorem:

- A chain of length five together with one isolated vertex satisfies the present theorem. It has no induced2+2, whereas every member of the earlier construction contains the two independent backbone chains of length at least two. Thus this example is outside that literal construction.
- Three disjoint two-element chains satisfy the earlier theorem with m=n=2,j=1, an empty attachment upset, and the third chain as the side poset. They have no pair with equal downsets and nested upsets: the only distinct vertices with equal downsets are the three minima, whose nonempty singleton upsets are disjoint. Thus this example does not satisfy the present neighborhood hypothesis.

This comparison is only between the two stated constructions. It does not establish novelty relative to the wider Gold literature or exclude other proofs of these examples.

The outside-skeleton extension follows the same rigorous invariant-map principle as [the disjoint-gadget composition theorem](disjoint_twin_gold_composition_20261019.md). Here the selected thresholds may change after conditioning; no transport of an unchanged strategy is claimed.

## 7. Incomparable extras: a genuine endpoint-local obstruction

The chain hypothesis cannot simply be replaced by two incomparable extras if the strategy must remain within the four marked vertices.

Take the actual six-point poset with full predecessor masks

    (0,0,1,7,3,1).

Here u=0,v=1,s=2,t=5. The other vertices are3 and4. The conditions are

    D(u)=D(v)=empty,
    U(u)={2,3,4,5}, U(v)={3,4}, s∥t.

There are37 extensions. The eight marked-order counts are

    (u,v,s,t):8; (v,u,s,t):8;
    (u,v,t,s):4; (v,u,t,s):4;
    (u,s,v,t):6; (u,s,t,v):2;
    (u,t,v,s):3; (u,t,s,v):2.

Optimizing the second comparison separately after each possible first answer, the best worst-path deficits for the four uncertain first pairs are

    u:v → −3; v:s → −5; v:t → −11; s:t → −3.

The two other marked pairs are already comparable and cannot give a successful first step. The checker includes them and allows all six marked pairs, including determined comparisons, as possible second comparisons. Hence no adaptive two-comparison Gold strategy restricted to{u,v,s,t} exists. The best restricted deficit is exactly−3.

The full poset does have a strict fixed Gold strategy: compare4 with5, then2 with5. With false/true row and column order, its table is

    [[11,11],[2,13]],

whose minimum Gold deficit is4. Thus this is only an endpoint-local obstruction, not a Gold or balance conjecture counterexample. No minimality or literature priority claim is made.

## 8. Exact controls

[The standalone checker](../experiments/verify_chain_surplus_gold_20261019.py) independently checks the monotone-bin rule on bounded integer laws including zero tails, every upper- or lower-chain-surplus marking in every naturally labeled poset through six vertices, and every individual marked skeleton in those posets. It also enumerates all37 extensions of the incomparable-extra obstruction, checks the total by subset DP, optimizes every restricted adaptive strategy, and verifies the unrestricted strict Gold control.

The [normal](chain_surplus_gold_controls_20261019.json) and [optimized](chain_surplus_gold_controls_optimized_20261019.json) receipts agree. Arbitrary chain length and ambient size follow from the injections and monotone-bin proof, not enumeration. Run from the repository root with `python3 05-knowledge/experiments/verify_chain_surplus_gold_20261019.py`; the checker prints its receipt without modifying files.
