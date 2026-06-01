---
id: HYP-2028
status: PARTIALLY-TRUE
source: codex-2026-06-01-S540
related:
  - HYP-2022
  - HYP-2023
  - HYP-2024
  - HYP-2020
  - HYP-2026
  - HYP-2008
  - HYP-2003
---

# HYP-2028: all sector-vectors are existentially realizable; LRC is forced hitting of the observer-empty face

**Claim.** Every composition `c=(c_0,...,c_{n-1})` of `n-1` into `n` sector
counts is realized by some primitive distinct speed set at an open time.  Thus
the meaningful sector-vector problem is not global existence.  It is whether
every fixed primitive clock menu intersects the observer-empty face
`c_0=c_{n-1}=0`.

**Constructive proof idea.** Choose `q=nL` with `L` large.  Sector `k` contains
open residues `kL+1,...,(k+1)L-1`.  Pick `c_k` distinct residues in each sector
with total gcd `1`; use those residues as speeds and set `t=1/q`.  Then each
runner lands strictly inside its assigned sector, and the sector-vector is
exactly `c`.

**Evidence.** `lrc_sector_vector_realizability_s540.py` constructs witnesses
for every sector-vector through `n=7`:

```text
n=3: 6/6
n=4: 20/20
n=5: 70/70
n=6: 252/252
n=7: 924/924
```

The good face size is `C(2n-4,n-3)`: `1,4,15,56,210` for `n=3..7`.
Bounded exact searches see all raw and good vectors through `n=6`; at `n=7`,
`B<=14` sees `897/924` raw and `202/210` good vectors.

**Forced-vector finding.** In bounded boxes, the vectors common to every
primitive speed set are bad boundary fans, not LRC witnesses.  For `n>=5` the
intersection stabilizes to:

```text
(n-1,0,0,...,0)
(n-2,1,0,...,0)
(0,...,0,1,n-2)
(0,...,0,n-1)
```

These come from times near `0` and `1`, where almost all runners are adjacent
to the observer.  The forced intersection has no good vector.

**Interpretation.** HYP-2022's sector restriction is not a restriction of
global sector-vector existence.  It is a restriction of fixed-clock menus,
sector-tournament quotients, low-complexity realization, and forced target
hitting.  LRC should be restated as:

```text
For every primitive speed set V, the sector-vector menu M(V) intersects
the good face c_0=c_{n-1}=0, possibly only on the compactified boundary
for AP-like tight rows.
```

**Concurrent integration.** Upstream HYP-2024 studies section/boundary functors
as pure quotients, and upstream HYP-2026 reframes LRC as support-flow
cancellation versus zero-cut cover flow.  This HYP supplies the ambient
sector-vector simplex for those views: all sector-vectors exist somewhere, so
section/boundary and flow/cut arguments must prove forced intersection with
the good face for each fixed clock, not mere existence of good vectors.

**Predictions.**

1. The constructive proof works for all `n`.
2. For `n>=5`, the only universally forced open sector-vectors are the four
   bad boundary fans above.
3. The AP initial segment is the unique family whose sector-vector menu hits
   the good face only on the boundary.
4. The right tournament target is a good-face ideal in anchored sector/hole
   classes, not a singleton sector-vector class.

**Files.** `04-computation/lrc_sector_vector_realizability_s540.py`;
`05-knowledge/results/lrc_sector_vector_realizability_s540.out`;
`07-reflections/lrc-sector-vector-realizability-s540.md`.
