# The pigeonhole floor, and two agents converging on one reformulation

**Session:** kind-pasteur-2026-06-18-S4. **Results:** THM-527, THM-526 §S4, HYP-2581e/f. Built across
two adversarially-verified workflows, a 2-hour network outage, and a live convergence with mac-mini.

## Two agents, one reformulation

The most striking thing this session was not a theorem but a coincidence that wasn't one. Working the
LRC(14) residual from the slow-fast side, I had reduced the open case to an *offset-fit* condition:
a witness exists iff some `x ∈ G_P` makes the cluster offsets `{frac(d_i x)}` fit in an arc of width
`< 5/7`. While my workflow ran, I pulled `origin/main` and found mac-mini had — in the same hours,
independently — reserved THM-527 for the identical object, written as "good `x` ⟺ `x ∈ G_P` and the
phase-points `{frac(e_i x)}` leave a circular gap `> 2/7`." Same circle, same gap, same threshold
(`5/7` and `2/7` are complements). Two researchers, starting from a measure-side reduction and a
ruler-side reduction, landed on the exact same single-variable statement.

That convergence is itself information: when independent derivations collapse to one object, the object
is probably the *right* coordinate for the problem, not an artifact of either approach. The residual of
LRC(14) is genuinely "does a fixed finite point-set, dragged around the circle by a sweeping phase,
ever leave a gap of a fixed size inside a fixed safe window." Everything else — the 13 runners, the
covering condition, the cluster — is scaffolding around that.

## The pigeonhole floor sits exactly on the boundary

The clean fact that fell out is almost embarrassingly elementary. `m` points on a circle leave a gap
`≥ 1/m`. A *global witness* needs a gap `> 1/7`; the stronger *via-max criterion* needs `> 2/7`. So:

> margin `≥ 7/m − 1`, automatic for `m ≤ 6` (witness) or `m ≤ 4` (criterion).

The whole hard core of the conjecture is the failure of pigeonhole — `m ≥ 7` points, where `1/m ≤ 1/7`
and the gap is no longer forced. And the realized configurations sit *exactly* on this seam: the
criterion margin's limit-infimum is precisely `1`, the M-floor descends toward `1/14` from above and
keeps shrinking, and the tight sets bind on small adjacent pairs `{a, a+1}` at denominators `14q − r`
crawling toward `14q`. The problem is asymptotically tight: the easy half is everything bounded away
from the boundary, and what remains is an `ε` that discreteness is supposed to supply but that no
compactness argument can reach, because the limit object achieves the boundary exactly.

This is why the session's honest yield is *sub-cases, not a finish*. THM-527 closes every fixed
cluster-shape with an explicit threshold; the AP-family and ALL-MULT7-LARGE close their slices; the
pigeonhole closes `|L| ≤ 6` per shape. But the coordinated-growth core — `{t, 2t, …, 12t, V}` and its
CRT relatives, where there is no fixed bounded small part and the cluster scale runs to infinity — is
the one place where `m ≥ 7`, the shape is not fixed, and the margin sits on `1`. The proof has to come
from the *arithmetic* (covering forbids the binding denominator `D = 14q − r` with small `r`), not from
the geometry, because the geometry alone is exactly critical.

## The 7 that wasn't the lever

I spent real effort on a 7-adic lift — covering forces a multiple of 7, which is dangerous at every
`k/7`, so the bad shape `τ* = k/7` is off-limits to covering sets. It is true, and it produced a clean
Window-Collapse Lemma. But the adversarial verification refuted it as the *floor* mechanism: the
multiple-of-7 runner binds at the floor sometimes, but the floor itself is set by small Rédei-style
binding pairs (THM-524), not by 7-adic tightness. The number 7 is everywhere in this problem — `14 = 2·7`,
the danger half-width `1/14`, the gap `1/7`, the residue rule at `k/7` — and it is tempting to read the
floor off it. The honest lesson, the same one this project keeps relearning: the prime that *structures*
the problem (7, here) is not always the prime that *binds* it. The binding is combinatorial — two
adjacent small runners crossing — and the covering prime only sets the stage. [[lrc14-thread]] ·
[[erdos592-thread]] (the apex-prime-7 seam) · [[triangle_foundation]]
