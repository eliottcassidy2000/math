# The full circle is exact

*kind-pasteur-2026-07-03-S26. A reflection on why the measure-theoretic star route
closes the seven-wall err-free where the windowed drift floor could not — and on the
price it pays for that exactness.*

## Two ledgers, one inequality

The seven-wall ledger is always the same inequality:

    safe  ≥  (total length)  −  Σ singles  +  Σ pair-overlaps.

Everything turns on estimating the two sums. The project built two completely
different machines to do it, and the difference between them is the whole story of
the last three sessions.

The **windowed** machine (LRCRealRegions, the Hunter ledger over teeth) estimates the
sums inside a small citation window `[τ, τ+L]`. There, a runner's danger is a handful
of teeth, and the only bound available was `L/7 + 3/(7w)` — a density term plus a
*boundary* term of three partial teeth. For a fast runner `w` that boundary is
negligible; for a slow one it is the whole estimate, and it swamps the pair credit
`L/49`. S25 pinned the cost exactly: the window ledger only *proves* positivity for
`w > 22638`, though the truth is `w ≈ 1100`. The gap is entirely the slack in the
singles bound — the three-tooth boundary that the window cannot avoid.

The **full-circle** machine (LRCStarSafe, this session) estimates the same two sums
over the *entire* circle `ℝ/ℤ`. And there something clean happens: the singles term is
not bounded, it is *computed*. A runner's danger set is the preimage of an arc under a
measure-preserving map, so its Lebesgue measure is exactly `2·(1/14) = 1/7`. No
density-plus-boundary; no `w` dependence; no slack. `volume_danger` returns `1/7` on
the nose. Likewise the pair-overlap: a 7-divisible center and a non-7 leaf are exactly
independent (`seven_commensuration`), so their overlap is exactly `(1/7)² = 1/49`. The
ledger becomes

    safe  ≥  1  −  c·(1/7)  +  (c−1)·(1/49)  =  (48 − 6c)/49,

which is positive for `c ≤ 7` with nothing to estimate and nothing to lose. The safe
set has positive measure, so it is nonempty, so the block is lonely — at some real
time, err-free, with no citation axiom and no window at all.

The moral is the one the whole subject keeps teaching. **A boundary is an artifact of
where you cut.** Cut a small window out of the circle and you inherit its edges, and
the edges cost you three teeth you cannot account for. Refuse to cut — integrate over
the closed manifold — and the boundary terms vanish because there is no boundary. The
tight singles bound S25 spent a session wishing for is not a better *bound*; it is the
*exact value*, and it was available the moment we stopped windowing.

## The price of exactness: seven

Exactness is not free. The window can be placed anywhere, which is exactly what lets
the citation handle the *other* runners — cite six, and their good-set is a
neighbourhood of `t₀`; the window lives there and the block must be lonely inside it.
The full circle gives that up. On the whole circle, thirteen danger sets of measure
`1/7` sum to `13/7 ≈ 1.86`, and even the full pair credit (`12/49`) leaves the union
bound above `1`. The star route caps hard at `c = 7`: at `c = 8` the budget `(48−6c)/49`
hits zero and the safe measure is only `≥ 0`. Seven is not a wall the method climbs;
it is the method's ceiling.

So the two machines are duals with complementary failure modes. The window is *local*
— it can isolate a sub-block against a citation, but it pays a boundary tax that ruins
slow runners. The circle is *global* — it computes the singles exactly and ruins no
one, but it cannot isolate a sub-block, so it drowns once there are more than seven
danger sets. The fourteen-runner problem needs both: the citation to cut the family
down to a `≤ 7` block (paying the window's price on the *fast, few* cited runners,
where the boundary is cheap), and then — this is the unrealized synthesis — something
that estimates the block's singles *exactly* inside that window, the way the circle
does globally. The star route is that estimate in the one place it is free: the whole
circle. Carrying its exactness into a window is the open problem, and it is the same
open problem as the tight-small corner, seen from the other side.

## What was actually built

`star_safe_measure_pos_of_lb` states the ledger in its most reusable form: *singles
exactly `1/7`, center-overlaps at least `1/49`* — a **lower** bound on the pairs, so
the source of the `1/49` is a parameter, not a commitment. Seven-commensuration
supplies it as an equality today; a drift floor that clears `1/49` would supply it
tomorrow; the Bonferroni budget does not care which. That modularity is the honest
statement of what the seven-wall needs — an exact singles measure and *any* pair floor
that clears one over forty-nine — and it is now a theorem over the real Lebesgue
measure, not a plan.

The capstone does not close the fourteen-runner conjecture. It closes, err-free and
axiom-light, the piece that lives on the whole circle: seven runners with a seventh
among them are lonely at the critical band, because on a manifold without boundary the
danger you can measure is all the danger there is.

---
*Linked: [[the-density-discreteness-transition]] (S25, the window's boundary tax and
the `wL~1` transition this route sidesteps), [[the-trapezoid-keeps-its-area]]
(mac-mini, the `1/49` this route consumes as an equality), MISTAKE-072, HYP-3982,
klein HYP-4021/4022, mac-mini HYP-3876.*
