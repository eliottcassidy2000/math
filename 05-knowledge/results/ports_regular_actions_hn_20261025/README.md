# Nonconsecutive repairs, smaller response states, and the six/seven/eight bridge

2026-10-25. **PROVED** statements are scoped below and in their proof files;
**FINITE-EXACT** receipts are supporting controls. Full black-core connectivity
remains **OPEN**. The ordinary A000568 recurrence is unrestricted and has a
strictly smaller response graph; it is still derived from Burnside and has
no proved polynomial bound. No new Hadwiger–Nelson bound is claimed.

## The strongest advances

1. [Protected nonconsecutive ports](ports/NEXT_BOUNDARY.md): for a literal
   closed blue word P=(0,...,n-1), write N^-(0)={n-1} union I. If n>=9,
   min I>=5 and {n-4,n-2} is contained in I, there is a three-black detour
   to the exact labelled K_P target. The other ports and all middle arcs
   are unrestricted. A local arc contradiction removes the former size
   bound; the only surviving cyclic case is the already repaired tail.
   Consequently the one-hole family
   I={n-k-1} union {n-k+1,...,n-2} is covered for k>=5,n>=k+6.
   The distinct [first endpoint-four hole](ports/README.md),
   I={n-5,n-3,n-2}, is covered for every n>=12.
2. [A smaller unrestricted counting state](counting/README.md): at fixed
   target degree N, replace the past gcd response f(d) and implicit cycle
   count m=f(1) by D(d)=dm-f(d), and normalize weight by 2^(Rm), where R
   is the actual remaining degree. This preserves every individual future
   contribution, and strictly merges reachable histories. At N=100 the
   comparable response graphs have 16,230 versus13,444 states and62,766
   versus50,044 contributions. These are graph sizes, not runtime claims.
3. [The clone-substitution boundary is now exact](counting/EXCESS_EXTENSION.md):
   a second, signed excess grading restores the same gauge for pair-free
   core coefficients. Without separating actual and excess degree, every
   old response can be recovered from individual futures; scalar
   normalization cannot produce the same merger. This gives both a valid
   recurrence and a precise obstruction, while retaining the distinction
   between aggregate core counts and literal clone owners. The corrected
   graded coefficients count signed vertex-weighted tournaments. A
   [direct clone-cut involution](counting/CLONE_INVOLUTION.md) now pairs their
   contributions inside each fixed inflated tournament: toggle the first gap
   throughout one automorphism orbit of nontrivial clone blocks. Its odd
   size reverses the sign; exactly the pair-free tournaments survive.
4. [The regular conic action has a geometric boundary](geometry/README.md):
   the inherited free-hourglass/regular-eight-conics criterion is retained
   with its marked direction map. The three divisor-relation modes of
   M=-2I+2A(Q3) are exactly the centered coordinate-function space for
   planar unit cube embeddings, subject to explicit norm and collision
   conditions. More generally Cay(F2^r,S) embeds injectively as a planar
   unit-distance graph iff S is independent. A dependent translation
   model cannot supply a planar six-colour forcing graph by itself.

## Inheritance and concept board

The anchor was nonconsecutive port repair; the niche was fixed-target
enumeration; the wildcard was the marked Fano/Hadwiger–Nelson bridge.
The six live concepts were literal endpoints, future responses, conserved
degree, common exterior assignments, translation squares, and affordable
arithmetic queries. The [review](AUDIT.md) records their proof boundaries.

Closest proved mechanisms: [the consecutive-tail repair and affine compiler](../trijection_hourglass_ports_20261025/README.md),
[the future-response transfer](../defect_boundary_20261020/boundary.md),
[the circle-port decomposition](../circle_ports_polygon_palettes_20261021/README.md),
and [the audited Fano-code criterion](../six_seven_geometry_20261025/fano_codes/RESULT.md).
The hostile examples are a removed insertion predecessor, substitution
excess pretending to be actual degree, independently satisfiable local
tests with no common exterior assignment, and a noninjective unit drawing.
The least-used sidecars were a guard's own arc, the conserved final degree,
the cycle anchor, and the full square-face equations.

The old 6/4/6 incidence dimensions retain their own
[common-kernel carrier](../trijection_hourglass_ports_20261025/incidence/common_kernel.md).
They are not an orbit partition of the sixteen free codes. The newly arrived
[norm-seven torus construction](../session61_hadwiger_six_seven_20261025/README.md)
instead selects their actual four-element geometric orbit. Its two ideal
choices and two face choices continue explicitly to four invertible maps
on the conic cube; the direction vectors are recorded in the geometry proof.

## Connections with explicit maps and losses

| Source and target | Map and preserved question | Lost data; needed sidecar | Cheapest decisive check |
|---|---|---|---|
| Closed blue path to a black detour | Literal words with exact path-edge XOR | No quotient is taken; every intermediate join and reversal is retained | The removed-predecessor h=4 hostile |
| Past permutation cycles to counting states | f(d) to D(d), weight to 2^(Rm)w; all individual completions | Cycle count is absorbed only on a fixed actual-degree fibre | Histories (23) and (5,9,9) at future lengths1,3 |
| Signed clone substitution to graded states | Separate actual remainder R and excess Q | Individual owners and ports were already absent; retain both grades | One versus two actual singleton continuations |
| Fano code to eight branch conics | a to (a.d1,a.d2,a.d3); freeness to regularity | Component-set action is not freeness on the whole K3 | Rank3/2/1 and orbit8/4/2 controls |
| K3 cube relations to planar cube coordinates | Singleton Walsh coefficient vectors | Curve geometry and intersection multiplicities; retain norms and injectivity | Dependent K4,4 directions and the R3 tetrahedron |
| Circle boundary lists to a path language | Exact attainable set to dead/free/forced(c) | Multiplicities and closure; retain count vector or cycle anchor | Same list sizes, opposite C6 outcomes |
| Local tests to universal forcing or menu coverage | Join tests on the same source assignment | Independent marginals lose shared choices; retain boundary relation | The unit pentagon and odd signed-cycle certificates |
| 30030 offsets to affordable gcd queries | Endpoint integers become cycle lengths | Degree and accumulated weight must still be kept | First new query59 separates equal-degree histories |

## The six/seven/eight and 30030 statements, precisely

The [minimal path-list automaton](hn/README.md) has q+2 states: dead, free,
and one forced state for each named colour. Thus four, five and six colours
give exactly six, seven and eight Boolean states. Its eight-state palette
action fixes dead and free, so it is not equivariantly the regular eight-conic
action. Weighted singleton continuations recover the full colour-count vector.
The [relational lift](hn/RELATIONAL_LIFT.md) retains the exterior assignment
shared by all path/C6 tests and yields an exact finite certificate language.
Applied to detours, it also allows a menu to mix legal and illegal ledgers:
each joint failure has a signed parity certificate. Uniformly small menus
and boundary relations remain separate open targets.

The arithmetic comparison is exact. Both cycle histories

    (30030-17,30030+17) = (30013,30047),
    (30030-1,30030+1)  = (30029,30031)

have degree60060 and the same affordable responses for every remaining
degree R<59. At d=59 the responses are2 and60, because30031=59*509.
The first three endpoints are prime. This is a finite-horizon merger already
visible in the old response representation, distinct from the new gauge
merger. The [prime-offset recovery](hn/README.md) also retains the sharp
composite-offset counterexample d=17^2; a primorial wheel is not a prime test.

## Next targets and actual remaining scope

At endpoint four and n>=15, the former h=n-5 remainder is removed. The
inherited residual bank can now be restricted to I={h,n-3,n-2},
3<=h<=n-6, closed middle and failure of all applicable prefix guards.
These conditions describe what is outside this bank; they do not assert
that such a source is disconnected. The next useful structural repair
should remove or relocate the need for the protected port n-4, keeping
the insertion predecessor and exact target visible.

For counting, compare an excess-graded engine against the old core engine
before claiming a net improvement there. For HN, seek bounded shared
boundary relations in an actual geometric construction. The small local
automaton alone does not bound their size. All these proposals are **OPEN**.

## Reproduction and publication

The packet uses Python's standard library. `python replay.py` reruns the
maintained checks and the independent fenced-code reviews in ordinary and
optimized modes; it checks the results against the saved receipts without
rewriting them. `python verify.py` checks source pins, local links and receipt
agreement. The [audit](AUDIT.md) separates structural proofs, finite universes,
inherited inputs and external citations. Incoming geometry at47eea7cad5 was
read as mathematical evidence and integrated after a clean checkpoint.
