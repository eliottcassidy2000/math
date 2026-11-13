# Exact-target repairs and an indefinitely reusable marked grammar

2026-11-12, America/Denver. **PROVED IN NOTES / FINITE-EXACT**, independently
reviewed. **The complete Chan–Pak spectrum remains OPEN.**

This session generalizes the earlier H17 preparation, finds repairs that
also preserve the next exact target, and turns an underused old port move
into a closed physical recursion. The new recursion needs no fresh source
lookup after entry. Its target image is nevertheless still sparse.

## 1. The two H9 five-cores have different next operations

[ONE_ARC_H9](history/ONE_ARC_H9.md) uses the two five-vertex cores with the
same entire unrooted polynomial `x^5+6x^3+2x`. Their directed-triangle
arc-overlap graphs are P3 and K3. Reverse one off-root arc, keeping every
incident edge at the selected root, including its unique incoming arc.
The rooted response changes

    (9,3,3,9) → (9,4,5,10).

Thus readiness q−s changes0→1 at the same count9. A first N birth keeps
the exact next count15 while repairing the named-pair/root divisor3→1:

    (15,6,9,18) → (15,5,9,19).

One more R gives primitive balanced `(29,9,9,29)`. The repaired intermediate
H15 source also enters the new recursive chart at H147. This uses the
operational distinction between the two equal-count cores, rather than
declaring their triangle-overlap drawings equivalent to their responses.
The scalar counts themselves and the core catalogue are inherited.

## 2. A closed actual recursion preserves primitive continuation forever

[CHAIN_REUSE](constructors/CHAIN_REUSE.md) retains the newest vertex in
the R-chain instead of forgetting it. Two live ports suffice. The resulting
actual rooted sources close on

    C(q,s)=(H,s,q,F)=(q+3s,s,q,4q+2s),
    q odd>0, s even>0.

For every k≥2 the literal macro returns to the same chart with

    s'=2^(k−1)(3q+2s),       q'=q+3s',
    gcd(q',s')=gcd(q,s).

The underlying small matrix was recovered from the old28-map grammar;
the closed chart, primitive source entries and exact inverse make it
usable as an indefinitely repeatable construction. Each finite stage
is an actual tournament with the named port available for the next step.
The macro has fixed-slot words `1^k02` and `3^k45`, alternating the active
slot; it never retrieves a forgotten mark.

Examples: H17 preparation enters C(135,38), count249. The one-arc H9
repair enters the smaller C(87,20), count147. Arbitrary admissible exact-
prefix repairs below also enter for every k≥2. Primitive balanced states
and this chart have explicit short physical entry/exit bridges using the
inherited R and C3 moves, with all residue guards retained.

[Exact addresses](CHAIN_LANGUAGE.md) recover the last k and parent from
(q,s), give the constructed rooted-family generating function
`(1−z)/(1−z−z^4)`, and prove the fixed-seed target count is O(sqrt X).
Thus indefinite reuse closes an interface gap; it does not yet supply
arbitrary odd targets. For constant k=2, actual counts satisfy the
recurrence `H_(j+2)=23H_(j+1)−4H_j`.

## 3. Exact-next-target repairs with arbitrary admissible left factors

[EXACT_REPAIR](EXACT_REPAIR.md) starts with a different five-vertex
two-arc change `(15,3,3,21)→(15,3,5,21)`. Put any nonempty tournament L
before that block in an ordinal join. With h=H(L), f=F_L(2), the first
root-child birth has the SAME exact count on both sides:

    M=15h+21f.

At the same physical root/child pair its exact continuation divisors are

    3 gcd(h,7f) → gcd(5h,21f).

Consequently `gcd(5h,21f)=1` gives an exact-target3→1 repair. For L=TT_a,
the targets are `21·2^a+15`, beginning57,99,183,351. Every admissible L
then enters the closed primitive chart, with no special clock choice.
The first birth makes the initially weak ordinal source strong.

The [H17 wrapper theorem](constructors/PROOF.md) separately gives a
one-arc preparation at every count2^n+1. Incoming Session298 independently
proved this all-n preparation with generic period clocks. Our refinements
include a factor-free explicit clock, all-n controller entry, a native
implementation, and the distinction between a primitive retained pair
and a primitive chosen root. Their scopes are credited explicitly.

## 4. The old unrelated work changed which repair was searched

The [history lane](history/HISTORY.md) follows July's exact planar-Jacobian
target-shear fibres, August's principal-part transgression and moving-
numerator continuation, and June's Hodge lift sidecars. The concrete
tournament transfer is the exact current/next-count fibre

    Delta(H,s,q,F)=(0,v,w,v).

Searching that fibre found physical repairs preserving the next target;
the old monotone H17 position shift cannot do this. Its complete R/N
response and the pure-deletion repair's exact future equality language
are proved in [CONTINUATION_FIBRES](history/CONTINUATION_FIBRES.md).

The [Hadwiger–Nelson lane](ports/PROOF.md) keeps joint boundary responses.
It gives marked slot-factorization and cyclic-wrapper formulas, then an
actual compensated block swap preserving both H27,775,899 and its next
target213,127,587. The remaining rooted F coordinate changes by147,456;
the divisor stays9. This supplies a different physical direction in a
response fibre, with its obstruction retained. No Jacobian, Hodge or
Hadwiger–Nelson conjecture is inferred from these transfers.

## 5. Direct ready handoffs and their failure boundary

[READY29](READY29.md) continues the H9 unit-ready branch: three off-root
arc changes repair the exact target29, preserve the entire root star and
predecessor, and create a neutral sink that can be deleted. The resulting
strong ready source exits at89. Three edits are minimal in the specified
seven-vertex root-star-preserving local ball. This is a finite handoff,
not indefinite exact H-to-3H+2 rearming.

[H9_INFLATION](constructors/H9_INFLATION.md) proves the all-parameter
inflation laws. For the two H9 cores, common
homogeneous substitution preserves both count equality and unit
readiness only when all three modules are singletons. The repaired
side nevertheless gives an explicit infinite unit-ready family.
The analogous H29 repair has target drift exactly 3 C_2(M) for an
arbitrary inserted module M, so it preserves that target only for a
singleton. Target and readiness must both be retained in a global proof.

The remaining decisive obligation is a target-directed source supply or
rearming move whose actual positive outputs cover the missing targets.
Neither primitive gcd, congruence coverage, finite handoffs nor the new
fixed-seed recursion proves that only7 and21 are omitted.

See [REVIEW.md](REVIEW.md) for proof reviews, controls and reproduction,
and [INCOMING.md](INCOMING.md) for concurrent-result attribution.
