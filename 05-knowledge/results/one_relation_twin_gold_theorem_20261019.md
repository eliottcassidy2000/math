# A two-comparison Gold strategy for twins differing by one relation

Date: 2026-10-19. Status: proved for arbitrary finite ambient posets. This is a small sufficient class theorem, with no novelty claim. The general Gold Partition Conjecture remains open.

## 1. Exact Gold convention

Follow [the terminal-context convention](gold_terminal_context_20261019.md) and [the prefix-side-poset theorem](gold_prefix_side_poset_20261019.md) in the repository. For a finite nonchain poset P, let t0=e(P). A Gold strategy chooses two comparisons so that every realized answer path satisfies

    t0 ≥ t1+t2,

where ti is the number of linear extensions consistent with the first i answers. The inequality is weak. Comparisons whose answers are already determined are allowed. If a first answer leaves a chain, t2=1, consistent with performing a determined second comparison. Two comparisons fixed in advance suffice for the result below; no adaptive endpoint choice is needed after receiving the first answer.

## 2. Ambient neighborhood theorem

Let u,v be incomparable vertices in a finite poset P. Suppose their full strict predecessor and successor sets satisfy

    D(u)=D(v),
    U(u)=U(v)∪{t}, with t∉U(v).

These are conditions in the entire ambient poset, including any exterior vertices.

**Theorem.** The fixed strategy

    first compare u with v;
    then compare v with t

is Gold. In addition, 1/2≤P(u<v)≤2/3.

The theorem imposes no degree bound and no restriction on the number or internal structure of the other vertices. It does not assert robustness under arbitrary reweightings of linear extensions.

## 3. Exactly three marked orders

We have u<t. Also v is incomparable with t: v<t would put t in U(v), while t<v would imply u<v. Thus every extension has exactly one of the following marked relative orders:

    u<v<t, with count A;
    u<t<v, with count B;
    v<u<t, with count C.

These three orders exhaust all extensions, because u<t is mandatory. In particular E=e(P)=A+B+C. Cells corresponding to any other marked order have count zero and are not realized paths.

## 4. Exact equality and injection

**Equality A=C.** In an extension of the first type, exchange the labels u,v at their occupied positions, leaving all other labels fixed. Their predecessor sets agree. Their common successors were after both original positions. The one exceptional successor t remains after both positions. The exchanged word is therefore a valid extension of the third type. The same swap is its inverse, proving A=C.

**Injection B into C.** In an extension with u<t<v, replace the labels at these three occupied positions by v,u,t respectively, leaving all other labels fixed.

- Moving v to the old u position is valid because D(v)=D(u).
- Moving u to the old t position does not cross any predecessor of u. Its successors other than t belong to U(v) and hence follow the old v position; t is now in that old v position.
- Moving t later does not violate a predecessor constraint. Every successor s of t is a successor of u distinct from t, hence belongs to U(v); it therefore follows the old v position.

The output has v<u<t and is valid. Its three occupied positions determine the input uniquely, so this is an injection. Therefore B≤C=A.

Consequently

    E=2A+B, 0≤B≤A.

This already gives P(u<v)=(A+B)/(2A+B) between1/2 and2/3. Positivity of the relevant cells is not needed in the following inequalities; all realized cells are covered even if a degeneracy is allowed.

## 5. The Gold path calculation

Order the first answers as u<v, v<u, and the second answers as v<t, t<v. The exact joint table is

    [[A,B],
     [A,0]].

The three realized path deficits E−t1−t2 are:

    u<v<t: E−(A+B)−A = 0;
    u<t<v: E−(A+B)−B = A−B ≥ 0;
    v<u<t: E−A−A = B ≥ 0.

On the third path the second comparison is already determined: v<u<t implies v<t. It leaves t2=A, exactly as used above. Thus the strategy is valid under the repository's determined-comparison convention, without an implicit change of threshold or stopping rule.

The zero first deficit shows why this is an exact weak Gold certificate rather than a strictly positive uniform margin claim.

## 6. Dual and exact twins

By order duality, if

    U(u)=U(v), D(v)=D(u)∪{s},

then compare u with v first and s with u second. The same three-cell argument in the dual proves Gold.

If u,v are exact twins, their label swap gives P(u<v)=1/2. Either first answer leaves E/2 extensions. Any valid second comparison, including a determined one, leaves at most E/2, so every path is Gold. This includes the two-element antichain with E=2 and t1=t2=1.

Equivalently, any incomparable pair whose full neighborhoods are identical or differ in just one predecessor/successor membership supplies such a strategy. The differing relation identifies the second endpoint when there is one.

## 7. Consequences for the seven-thin mixed-(1,2), c=4 core

Use roots x,y, lower slave b, upper slaves a1,a2, four common neighbors, and exterior A below both roots and B above both. A common vertex with k core comparable vertices has core incomparability degree8−k and hence at most k−1 exterior incomparable vertices. A core vertex cannot lie below any lower exterior or above any upper exterior, because it is incomparable with at least one root. Thus k=1 forces every lower exterior below the vertex and every upper exterior above it; k=2 permits at most one missing exterior relation. These observations prove the ambient neighborhood equalities used below.

1. If a common z has only the core comparability b<z, it is an exact ambient twin of y. The twin strategy is Gold.
2. If a common z has only the core comparability z<a_i, saturation gives D(x)=D(z) and U(x)=U(z)∪{a_j}, j≠i. Compare x with z, then z with a_j. This is Gold for every actual permitted exterior.
3. If a common z has exactly the two core comparabilities z<a1,a2, it differs from x in at most one exterior relation. If none is missing, use the twin strategy. If an upper exterior t is incomparable with z, compare x with z and then z with t. If a lower exterior e is incomparable with z, compare z with x and then e with z. These are Gold strategies, and the latter two explicitly permit an exterior second endpoint.

These are new consequences within this proof packet, not a claim of priority over known low-slave or twin sufficient criteria.

## 8. Distinction from the general directed-defect lemma

The weaker condition

    |D(v)\D(u)|+|U(u)\U(v)|≤1

only yields a one-sided probability bound by two injections. It does not assert the opposite neighborhood inclusions or equality D(u)=D(v), and it need not give A=C. This note does not infer a Gold strategy from that weaker directed-defect condition alone, nor from an arbitrary balanced pair.

## 9. Independent exact controls

[The standalone checker](../experiments/verify_one_relation_twin_gold_20261019.py) generates all naturally labeled posets through six vertices, finds every marked pair satisfying the one-extra-successor condition and its dual, enumerates extensions independently, checks all three counts and Gold path deficits, and directly validates the swap and three-cycle images and injectivity. The [normal receipt](one_relation_twin_gold_controls_20261019.json) and [optimized receipt](one_relation_twin_gold_controls_optimized_20261019.json) agree. These checks control the argument; arbitrary ambient size follows from the explicit bijection and injection above.

Reproduce from the repository root with `python3 05-knowledge/experiments/verify_one_relation_twin_gold_20261019.py`; the checker prints its receipt without modifying files.
