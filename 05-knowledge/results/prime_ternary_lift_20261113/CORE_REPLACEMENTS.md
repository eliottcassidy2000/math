# Exact ternary children by replacing the fixed core

2026-11-13, Codex. **PROVED IN NOTES**, using finite integer certificates
and the universal module expansion below. Chan--Pak remains **OPEN**.
These are constructions and interface identities, not new spectrum records.

## 1. The interface that is retained

Let H(T) count directed Hamiltonian paths. For a nonempty tournament M,
let K_j(M) count ordered j-tuples of nonempty directed paths partitioning
its vertices. A rooted frame (A,r) is inflated by replacing r by M, with
every vertex of M inheriting all exterior incidences of r.

Let a_j(A,r) count words containing every exterior vertex exactly once
and r exactly j times, with no adjacent copies of r and every consecutive
arc legal. Splitting a Hamiltonian path at its visits to M gives a bijection:

    H(A[r <- M]) = sum_j a_j(A,r) K_j(M).                 (1)

This is the inherited run expansion, used for example in
[THM-5133, flip-and-spin kernel signature](../../../01-canon/theorems/THM-5133-flip-and-spin-kernel-signature-sign-tournaments-of-twisted-algebras-and-the-pg32-extremal-theorems.md)
and [the clone-clock notes, section 0](../clone_clocks_negative_collatz_20261113/CLOCKS.md).
We retain the **entire response vector**, rather than just H or a paired
endpoint direction. Equation (1) does not assert that either endpoint
profile is preserved.

For M=TT_m, K_j(M)=j! S(m,j), the number of onto maps from m vertices
to j ordered bins. In particular K_1=1 and K_2=2^m-2. If the bipartite
arcs from the out-neighbours to the in-neighbours of r contain no matching
of size two, then a_j=0 for j>=3: each of two disjoint internal excursions
would contain a distinct out-to-in crossing arc, giving that matching.

## 2. One source, both digits, every m>=1

Rows below are out-neighbour bitmasks, vertices numbered from zero.
The tables are certified by enumerating **all permutations of the exterior
vertices and all sets of root-insertion slots**, independently of the path DP.

| frame | rows | r | a_1 | a_2 | a_j, j>=3 |
|---|---|---:|---:|---:|---:|
| S | (2,20,9,35,45,7) | 0 | 29 | 15 | 0 |
| B2 | (2,12,9,49,135,215,159,15) | 4 | 89 | 45 | 0 |
| B4 | W(10010110) | 4 | 91 | 45 | 0 |

Here W(w) is a transitive spine, in the order of the letters of w, plus
one apex which beats precisely the letters marked 1. All three contracted
frames are strong and satisfy the two-run guard above.

For every nonempty tournament M, equation (1) gives

    H(B2[M]) = 3 H(S[M]) + 2 H(M),
    H(B4[M]) = 3 H(S[M]) + 4 H(M).                      (2)

Consequently, for **every integer m>=1**, the same source S[TT_m] has count

    h_m = 15*2^m - 1,

and the two targets have counts 3h_m+2 and 3h_m+4. The source has m+5
vertices; the targets have m+7 and m+8. Each is strong: inflation of a
vertex in a strong frame remains strong when the inflated set is nonempty,
since exterior paths connect any two of its vertices.

| m | source | +2 target | +4 target |
|---:|---:|---:|---:|
| 1 | 29 | 89 | 91 = 7*13 |
| 2 | 59 | 179 | 181 |
| 3 | 119 | 359 | 361 |
| 4 | 239 | 719 | 721 = 7*103 |

These supply genuine prime-lane instances, without requiring the source
itself to be prime. No infinitude of prime values in any displayed family
is claimed.

## 3. A particularly small physical +4 move

Retain all old vertices. Map the five fixed source vertices to B4 as

    1->8, 2->3, 3->1, 4->0, 5->2;  source module -> target module 4.

Append target vertices 5,6,7. Every old arc incident to the transitive
module is retained. The **only old arc reversed is 2->3** in source labels.
Thus the +4 replacement uses one old core reversal and three new vertices,
independently of m. Internal module arcs are untouched.

Equivalently, S[TT_m] can be drawn as W(1001 0^m) with the single spine arc
3->1 pointing backwards. Straighten that arc, then append the spine
letters 110. The result is W(1001 0^m 110). The intermediate count is
9*2^m-3; the reversal is **not** a count-neutral preparation.

For +2 use the embedding into B2

    1->0, 2->1, 3->3, 4->6, 5->5;  source module -> target module 4.

Append target vertices 2,7. This changes m+2 old arcs: one module
incidence per clone and two fixed core arcs. `core_moves.py` minimizes
these two costs lexicographically over all core injections and verifies
the resulting labelled graphs. This is not a claim of optimality among
all possible source and target tournaments.

## 4. What has and has not escaped the earlier obstruction

The [factor-five theorem](../strategy_reframe_20261113/GROWTH.md) applies
to returns preserving a particular paired direction in identified-root
contexts. These replacements change the old core and retain the module
response instead. There is no conflict with that theorem.

Equation (2) also records the first failed implication in a tempting
recursion: a generated tournament cannot replace TT_m without checking
its Hamilton count. With M=C3, K_1=3 and K_2=6, the source count is177;
the targets are537 and543, namely 3*177+6 and 3*177+12. The desired fixed
digits 2 and4 have become 2H(M) and4H(M).

The target frame is not automatically a source frame for the next digit.
The next obligation is to find an explicit response-vector transition
out of that target, or a proved replacement into another applicable
source class. Attainment of its scalar count alone does not supply that
transition. Conversely, the existential prime-lane reduction does not
require retaining this particular witness; a target-specific alternative
source is allowed.

## 5. Certificates and boundaries

`core_moves.py` / `moves.json`: independent skeleton counts; every old
vertex retained; literal Hamilton-path DP, strongness, tournament validity,
and arc-edit budgets for m=1..7. The all-m proof is (1) plus the certified
finite coefficient table, not extrapolation from those seven tests.

The minus-one rays themselves overlap established results, including
[THM-5145, dyadic neighbourhoods and two-block families](../../../01-canon/theorems/THM-5145-dyadic-neighbourhoods-and-two-block-families-of-hamiltonian-path-counts-complete-census-of-root-types-of-order-ten.md).
The contribution here is a shared source and explicit core edits for the
two exact ternary targets. No global absence statement follows from the
order-eight catalog, and no endpoint-neutrality claim is made.
