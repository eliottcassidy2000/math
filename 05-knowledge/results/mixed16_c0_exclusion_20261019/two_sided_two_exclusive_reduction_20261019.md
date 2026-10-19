# Both exclusive-neighbor sides must contain at least two vertices

Date: 2026-10-19. Established dependency corollary after complete mixed(1,6), c=0 certification and independent replay. No full seven-thin balance theorem or novelty claim.

## Statement

In a hypothetical finite seven-thin poset P with no balanced pair, choose an incomparable pair x before y of maximum index span in the strict >2/3 majority order. Then

    |D(y)\D(x)|≥2 and |U(x)\U(y)|≥2.

This is a reduction of possible counterexample neighborhoods, not a full seven-thin balance theorem. Morphologies with at least two exclusive neighbors on both sides remain to be treated.

## Critical pair and exact degree budgets

The strict >2/3 comparison tournament is transitive: a directed triangle would have cyclic probabilities summing to more than two, while any linear extension satisfies at most two cyclic indicators. Its order M extends P.

The selected pair is critical. If z<x but z not<y, then z∥y (the alternative y<z would imply y<x), and M has z before x before y. The incomparable pair z,y would have greater span, a contradiction. Thus D(x)⊆D(y). Dually U(y)⊆U(x).

Write

    l=|D(y)\D(x)|, u=|U(x)\U(y)|,
    c=|I(x)∩I(y)|.

Every lower exclusive vertex is below y and incomparable with x, and every upper exclusive vertex is above x and incomparable with y. Therefore

    d_P(x)=1+l+c≤7, d_P(y)=1+u+c≤7.

In particular l,u≤6 and c≤6−max(l,u).

## Exhaustion of zero- and one-exclusive-side cases

If l=0 or u=0, the [combined one-sided exclusion](../five_six_upper_exclusion_20261019/one_sided_corollary_20261019.md) contradicts the hypotheses. It includes the zero-slave endpoint and order duality.

It remains to exclude min(l,u)=1. By order duality, take l=1 and 1≤u≤6. The following established results cover every allowed c:

1. u=1: the [two-slave all-common theorem](../two_slave_all_common_balance_20261019.md), valid for all c, covers the required c≤5.
2. u=2: the [mixed(1,2), c≤3 theorem](../mixed12_witness_migration_20261019/mixed12_witness_migration_20261019.md) and [c=4 theorem](../mixed12_c4_exact_exclusion_20261019/mixed12_c4_exact_exclusion_20261019.md) cover c≤4.
3. u=3: the [full mixed(1,3) corollary](../mixed13_c3_exclusion_20261019/mixed13_c3_exclusion_20261019.md) covers c≤3.
4. u=4: the [full mixed(1,4) corollary](../mixed14_c2_exclusion_20261019/mixed14_c2_exclusion_20261019.md) covers c≤2.
5. u=5: the [full mixed(1,5) corollary](../mixed15_c1_exclusion_20261019/mixed15_c1_exclusion_20261019.md) covers c≤1.
6. u=6: the [mixed(1,6) theorem](mixed16_c0_exclusion_20261019.md) covers the sole possibility c=0.

Every zero- or one-exclusive-side configuration is impossible. Hence l,u≥2 as stated.

## Scope and attribution

This corollary is a direct exhaustion using the linked precisely scoped theorems and their exact certificates; it is not a new independent computational scan. The constituent files include their proofs, replay scope and classical-source attribution. No novelty or priority claim is made. In particular, this does not assert that every seven-thin poset has a balanced pair: the genuinely two-sided l,u≥2 frontier remains.
