# Rooting collapses on the self-converse almost-regular boundary

Status: PROVED finite-class bijections, plus FINITE-EXACT order8 control.
No claim of literature novelty.

Throughout the group-theoretic proof assume even n≥2. If n=0 is included,
both sides of the counting identity equal1 directly; the one-vertex completion
does not have disjoint ordinary/anti permutation sets.

Let a_n be the number of almost-regular isomorphism classes on even n vertices,
and let r_q range over regular tournament classes on odd q=n+1. Unique regular
completion gives

    a_n = sum_[R regular q] number of Aut(R)-orbits on vertices.

Labelled counts satisfy A_n=R_(n+1): deleting the prescribed label * and using
the unique completion are inverse maps. Thus the labelled exponential generating
function of almost-regular tournaments is the derivative of the regular EGF,
with the empty object included consistently. Ordinary class counts are NOT
obtained by simply multiplying each regular class by its order.

A stronger equality holds on the self-converse subfamilies:

    # self-converse almost-regular classes of order n
      = # self-converse regular classes of order n+1.                 (1)

## Proof of(1)

For any finite tournament R its automorphism group G has odd order. An
involution in G would swap a pair of vertices and send their arc to its reverse,
which is impossible; Cauchy's theorem therefore excludes even |G|.

If R is self-converse, its automorphisms and anti-automorphisms form a group
G± of order2|G|. Every Sylow2 subgroup has order2, so an anti-involution exists.
All anti-involutions are conjugate under G: Sylow conjugacy initially supplies
h∈G±; if h is anti, replacing h by h times the first anti-involution makes the
conjugator an automorphism without changing its conjugation result.

For odd-order R an anti-involution fixes exactly one vertex: at least one by
parity, and at most one because an arc between two fixed vertices could not
reverse. Hence the fixed vertices of all anti-involutions form exactly one
G-orbit.

A rooted class(R,v) is self-converse exactly when an anti-automorphism fixesv.
Such an anti-automorphism has order2m with m odd (by Lagrange in G±), so its
m-th power is an anti-involution still fixingv. Consequently exactly one vertex
orbit of each self-converse odd R gives a self-converse rooted class. A
non-self-converse R gives none. Unique regular completion transports these
rooted classes to almost-regular classes, proving(1).

This is a fixed-point/orbit mechanism, not an assertion that all vertex deletions
of a self-converse regular tournament are self-converse. Most need not be.

## Complete small control and metagraph consequence

The separate mass-certified census contains15 regular classes at order9 and85
almost-regular classes at order8, both with labelled mass3230080. Independent
converse-isomorphism testing finds all15 regular classes self-converse but only15
of the85 almost-regular classes self-converse. The remaining70 form35 converse
pairs. Therefore the almost-regular merged-node count is

    (85+15)/2 = 50.

This is an exact count of nodes restricted to almost-regular tournaments, NOT a
claim that this restricted set is closed under nonspine-complement wiring.
Regular completion retains the root as a necessary sidecar for the full85-class
count; anti-involution conjugacy removes that multiplicity only for the
self-converse count.

Reproduce: python3 check_selfconverse_completion.py, using the accompanying
session30-almost-regular/census8_bank.json. The class completeness comes from
that packet's independent degree-sequence count and automorphism-mass proof;
this script independently checks converse existence and rooted completion.
