# Rigidity of the minimum-amplitude ternary clock

Separate Session 197 equality-case appendix to `MINIMUM_TERNARY.md`; its frozen source is unchanged.

**Theorem.** Every rooted tournament whose transitive root clones satisfy

    H_m=60*3^m-1

for all positive `m`, or eventually, has one of exactly two rooted strong cores:

    A=[2,20,41,35,45,195,31,95], root 3;
    B=[4,5,56,67,43,75,151,63], root 4.

The two are converses, and each full automorphism group is trivial. The complete arbitrary-order sources are obtained by placing any transitive prefix before the core and any transitive suffix after it, with every cross arc forward. No other core, marked root, or nontrivial outside component is possible.

Consequently, at every order `n>=8`, there are exactly

    2(n-7) unlabelled sources,
    2(n-7)n! labelled sources.

Each source has a unique root that gives this clock, so these counts are the same whether or not that root is included in the notation. At orders below eight there are none. These counts classify the amplitude-60 stratum only, not all larger-amplitude ternary clocks.

## 1. Complete finite reduction

The preceding packet proves the exact matching/run correspondence, eventual-to-exact implication, SCC gcd reduction and safe strong-core generation. Apply those results at `A=60`. All non-root strong components have product dividing `gcd(179,539)=1`, so they are singleton components. The rooted strong core has Hamilton count exactly 179. Busch's bound and

    5^10=9765625 > 179^3=5735339

force its order to be at most ten.

The strong deletion lemma permits generation from the triangle by every mixed one-vertex extension. Induced monotonicity of Hamilton counts makes the cap 179 safe at all parent orders. Exact outdegree-preserving canonicalization is the same proved isomorphism procedure as in the minimum packet.

## 2. Exhaustive equality inventory

The complete capped strong banks are:

| Order | Extension candidates | Accepted before deduplication | Classes with H<=179 |
|---|---:|---:|---:|
| 3 | initial triangle | initial triangle | 1 |
| 4 | 6 | 6 | 1 |
| 5 | 14 | 14 | 6 |
| 6 | 180 | 180 | 35 |
| 7 | 2170 | 2169 | 352 |
| 8 | 44352 | 4248 | 736 |
| 9 | 186944 | 2017 | 380 |
| 10 | 193800 | 442 | 84 |

In every graph of exact count 179, test every root of matching number two. There are **152** such candidates: 114 at order eight, 26 at order nine and 12 at order ten. The exact equations

    H_2=539,  H_3=1619

leave only the two marked eight-vertex cores displayed above. By the proved run bound, both have exact coefficients `(a_1,a_2,a_3)=(179,180,60)` and therefore give the stated all-multiplicity family. No order-nine or order-ten core survives.

The complete bank contains 1595 strong graphs. The Python checker independently checks their counts and canonical forms and retests every candidate root, using a different Hamilton recurrence and a different matching implementation from the C++ generator.

## 3. Rooted converse and automorphism certificates

The permutation

    [7,6,5,4,3,2,1,0]

maps the root of `A` to the root of `B` and reverses every arc. The original hybrid construction `[2,12,145,149,99,143,175,19]`, rooted at 7, maps to `A` by

    [0,1,2,4,5,6,7,3].

The two canonical adjacency tuples differ, so the cores are not isomorphic. Direct full permutation checks give a trivial automorphism group for each. The root inventory already checks every vertex and leaves exactly one valid marked root in each underlying graph.

## 4. Padding and the class formulas

All other strong components must be singletons by the gcd argument. The ordered SCC chain therefore has a transitive prefix, one of the two nontrivial cores, and a transitive suffix. Conversely every such padding preserves every clone count, so this description is both necessary and sufficient.

For total order `n`, the prefix can have any of `0,1,...,n-8` vertices; the suffix size is then determined. The two distinct cores give `2(n-7)` classes. SCC order is invariant under isomorphism, so different prefix sizes cannot become isomorphic. Each core is rigid, as is each transitive padding chain, and no automorphism can permute the ordered SCCs. Thus each full source is rigid and contributes exactly `n!` labelled graphs.

A padding vertex is a singleton SCC and has constant clone count. The unique clock root in the core remains the only one in the full source. This justifies the same formulas for marked and unmarked source graphs.

## 5. Reproduction and scope

Run `verify_ternary_equality.py` normally and under `python -O`; it compiles `verify_ternary_equality.cpp` in a temporary directory. Its Python helper dependency is `verify_minimum_ternary.py`. A C++17 compiler is required. The deterministic receipt includes the full bank hash, all finite inventory counts, independent root checks, rooted isomorphisms and automorphism orders.

This equality classification is specific to the minimum-amplitude transitive-clone family. It gives no universal arithmetic coverage theorem or classification of unrelated ternary graph-changing operations.
