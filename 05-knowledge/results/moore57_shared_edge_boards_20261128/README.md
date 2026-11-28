# Shared-edge pentagon boards on the fifteen-neighborhood chart

2026-11-28. This extends the [fifteen-star witness](../moore57_fifteenth_complete_star_20261128) and the [pentagon board tests](../finite_zeta_and_moore_pentagon_boards_20261128). It freezes the existing stars and inserts actual edges, not separate hypothetical matchings on overlapping boards.

## 1. Complete one entire branch-pair matching at a time

Root at the existing fully filled radius-two tree. Its 57 branches each have 56 leaves. Within a branch there can be no edge, because the two leaves already have their common branch hub. Between each two branches the leaf-to-leaf edges must be a matching: two edges with the same endpoint would make a square through the opposite hub. In a completion each branch pair is a *perfect* matching, since the unique route from a leaf to a different branch hub must go through one leaf of that branch.

Delete the endpoints of already present matches in the nominated branch pair. Make an edge available precisely when its endpoints have current distance at least four. **Any perfect matching of this residual bipartite graph can be installed all at once, preserving girth at least five.** A triangle cannot use two edges of a matching. A square with two new edges would either need an edge internal to a branch, or old cross edges between vertices that were unmatched in the branch pair. Neither is possible. One-new-edge triangles or squares are excluded by the distance test. Each leaf still has at most one contact per other branch, hence the degree cap stays 57. This is a construction for a single batch, not a claim that arbitrary batches commute. With every branch pair filled, the graph would be 57-regular of girth at least five on 3,250 vertices, so the Moore bound forces diameter two.

A single installed edge xy occurs as a fixed match in exactly s(x)s(y) pentagon boards, where s(x) counts the saturated neighbors of x (in the frozen chart). Choose one saturated neighbor u of x and one v of y. They must be distinct and nonadjacent, or xy closes a triangle or square respectively. Their common neighbor is known in the frozen chart and is neither x nor y, again by absence of triangles. Conversely the two centers of any board using xy give just such a choice. There is no independence between those boards: they receive the *same* literal edge.

## 2. The common equation behind the boards

Let S be the saturated set in the frozen chart, A0 its adjacency, X the remaining adjacency, and E the coordinate selector for S. Put K=A0²+A0−(d−1)I, with d=57. Since EᵀX=0, the required frozen-row system is linear:

    Eᵀ A0 X = Eᵀ (J − K).

In a simple degree-capped supergraph of this chart, the exact frozen-row system is equivalent to solving all the pentagon boards with shared edge variables. Here the chart hypotheses are that the rooted Moore tree is present, the frozen graph has girth at least five, and every pair of frozen saturated sites already has distance at most two. For one direction, the system gives the unique two-step route from each left slot to the opposite center, hence the board matching. For the other, take a still-uncovered pair (v,x), v saturated. The site x is a leaf and has a frozen saturated parent u. The centers u and v cannot coincide or be adjacent, since otherwise x is already within two of v. Their frozen common neighbor w cannot be x, so x lies on the u-side of the u-w-v board. The board gives x its route to v. Now every off-diagonal entry in the row at v is at least one, and its row sum before subtracting the diagonal is bounded by the degree cap; specifically, there are at most d² nontrivial walks from v of length one or two with last vertex different from v. On d²+1 vertices this forces all those entries to equal one. The equivalence does not test holes between two unsaturated vertices.

This compression suggests a complete search (not run here): choose branch-pair matchings, enter their variables globally, propagate unique supported edges, and backtrack on either a failed board or a physical collision between forced edges. Every prefix of any actual completion survives these tests. A simple example of why row feasibility is weaker than shared feasibility is three exact-one rows x+y=y+z=z+x=1: every row individually supports both edges, while adding the rows modulo 2 gives 0=1.

## 3. Replayed actual extension

Starting from 13,284 edges, seeded matchings installed **10,565** edges in **190** previously incomplete root branch pairs. Alternating-cycle support filtering across all overlapping boards then forced **16 physical edges**. After inserting them globally and recomputing, it forced **one more**, (942,2081). The final graph has:

| Quantity | Frozen chart | Extension |
|---|---:|---:|
| Edges | 13,284 | 23,866 |
| Completed pentagon boards | 729 | 919 |
| Fixed board incidences | 154,070 | 346,408 |
| Uncovered ordered pairs from saturated sites | 676,068 | 585,840 |
| Pentagons | 60,756 | 302,883 |

All **28,467** boards pass the final Hall and alternating-support screen. The original 241 saturated sites stay saturated, no new sites do, and a separate enumeration checks no triangle or square. Exactly 26,506 boards have at least one newly installed match; each physical edge serves 16, 20, or 25 boards. The next strip in this particular unguarded seeded order already fails five marginal Hall tests even before installing the forced edges. That rejects this continuation, not the 190-strip chart. There are still 1,247,744 unfilled board incidences. We have not found a simultaneous full assignment or settled existence of the degree-57 Moore graph. Next construction target: vary the matching or branch order at this failure, keeping global forced-edge feedback as an acceptance guard.

## 4. Supersingular primes: which overlaps are structural?

[Ogg's criterion, as stated and developed by Duncan–Ono](https://uva.theopenscholar.com/files/ken-ono/files/1411.5354.pdf), selects the 15 primes at which all supersingular j-invariants are in the prime field, equivalently the prime-level Fricke quotient has genus zero. They are the prime divisors of the Monster. The shared lesson is that descending many local objects to the same ground data is the test, not merely counting the objects. Here the exact arithmetic footprint is different:

* The three certified center-to-center transport groups in the original chart are S56. Their prime divisors are all primes ≤53. They contain **13 of the 15 Ogg primes**, but not 59 or 71; in the other direction they contain 37, 43 and 53, which are not Ogg primes. The chart itself is rigid; transport is not an automorphism action on the whole graph.
* A hypothetical Moore adjacency has integral eigenvalue labels 57, 7, −8. Their pairwise differences are 50, 65, 15, so the exact spectral collision primes are **2, 3, 5, 13**, all Ogg primes. The intrinsic 15 here is the root separation √(4·57−3), not the number of completed exterior sites. At 3 and 5 the two centered labels coalesce; at 2 and 13 the principal label collides with one centered label. This supplies concrete primes at which a separated-eigenvalue argument needs different handling, not an explanation of the Ogg list.
* **One small prime suffices for an exact final certificate.** For any symmetric simple adjacency on d²+1 vertices of maximum degree d, it suffices that every off-diagonal entry of A²+A be 1 modulo any integer m≥2. Each is a nonnegative integer and therefore at least one. The row's d² nontrivial walk-budget bound forces every entry to be exactly one, and equality forces degree d. So a Boolean adjacency with the degree cap satisfying the Moore matrix identity modulo 2 is a genuine Moore graph, not a mere local residue match. The distinction is critical: solving for unrestricted field elements or forgetting the Boolean degree cap loses the implication. This is also the exact reason that an Ogg-prime CRT bank is unnecessary for verification. Compare the [Chan–Pak finite readout example](../finite_zeta_and_moore_pentagon_boards_20261128): there, equal marginal counts lose disjoint-cycle data; here the degree budget, not a parity analogy, recovers the missing exact integer information.

## Replay

From this directory, with C++17 and Python 3:

~~~sh
python3 ../moore57_fifteenth_complete_star_20261128/replay.py --adjacency >/tmp/moore-frozen
g++ -O3 -std=c++17 feedback.cpp -o /tmp/moore-feedback
python3 certify.py /tmp/moore-frozen matching_witness.txt /tmp/moore-feedback /tmp/moore-shared --matchings
~~~

The witness is one line per branch pair: two zero-based branch indices and a 56-character permutation, encoded by the first 56 characters of Python's ascii_letters followed by digits. Thus replay does not depend on a seeded shuffle. To regenerate the matchings from the search, compile strips.cpp with C++17 and run it with arguments 20261128 190, reading the frozen adjacency on stdin. It emits the augmented adjacency before feedback. Running the same search with 191 reproduces the reported failed next step with g++ on Linux.
