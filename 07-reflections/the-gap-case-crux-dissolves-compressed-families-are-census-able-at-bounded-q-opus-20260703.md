# The gap-case dispatch crux dissolves: compressed families are census-able at bounded q (13 runners bound the band)

*opus-2026-07-03-S55. The owner asked me to work the gap-case dispatch and cluster identification — the
hard branch of `hlarge` (families that are not all-comparable, ratio > 13). mac-mini-S26 had refocused the
crux to the "aligned band-blockers" (compressed large clusters, lonely at `q ~ log(mag)`, no peel). Testing
it, the crux **dissolves**: with only 13 runners the compressed witness `q` is BOUNDED, so the census — not a
renormalization tower — closes the whole compressed branch. En route, I corrected a genuine S52 error.*

## The correction first (MISTAKE-097)

My S52 claim that the deep well `{1..12,182}` is lonely *only* at the Eisenstein `t*=14/183` and thus
"census-invisible" was **wrong**: the S52 script scanned only difference-set denominators. A correct full
Farey scan (all `a` coprime to `q`) finds it lonely at **`3/40`** (min-dist `3/40 = 0.075 > 1/14`, `q=40`).
`14/183` is the covering-min *argmax* (the tightest lonely time), not the only one. The Eisenstein arithmetic
(`14` a primitive 6th root mod `183`, `14·13 ≡ −1`) stands; the census-invisibility does not.

## The crux, and the test

The gap-case dispatch's genuinely-open branch was the **compressed** families: no dominant runner (every
`v_i` within `13×` of some `v_j`), covering, `≥ 7` far runners (`> 22`) — a geometric chain like
`{1, ~12, ~144, …}`, global ratio `≫ 13` but never `13×`-dominated. These are neither `spread13` (needs
*global* ratio `≤ 13`) nor peelable (no dominant runner; chain rungs are only `~12×` apart, far below the
`~10⁶` far-peel threshold) nor absorbable by THM-608 (compressed `≠` near-equal). mac-mini's counting gave
`q` up to `~13 log(mag)`, suggesting `q → ∞` and a forced renormalization tower.

**The test (verified):**
- Every compressed covering `≥7`-far family I found (12 of them, ratios up to `312`) is lonely at
  **`q ≤ 17`** — none exceeds `50`.
- **Scaling** a fixed compressed chain by `c = 1 … 20000` (magnitude `84 → 1,680,000`, staying covering +
  compressed + `7`-far throughout): the witness `q` stays pinned at **`17–19`**. It does **not grow with
  magnitude**.
- mac-mini's own compressed constructor found **no** high-`q` compressed blockers at any scale.

## Why: 13 runners bound the band, hence q

The counting bound `q ≲ 13 log(mag)` is an *upper* bound that assumes new prime-covering runners at every
scale. But a family has exactly **13 runners**, so a compressed chain has **≤ 13 rungs**, hence can carry only
a **bounded number of band-primes**, hence blocks only a **bounded band** `{15..Q}` — so a lonely `a/q` exists
with **`q = O(1)`** (empirically `≤ 20`), regardless of magnitude. The `q → ∞` regime is reached only with
*unboundedly many* runners; at `n = 13` it is unreachable for compressed families.

## Reconciliation with HYP-4040 (the lcm family)

mac-mini's HYP-4040 exhibited `q → ∞` for `{1..11, 13, lcm(2..X)}`. That family is **dominant** (the `lcm`
runner is `≫ 13×` all others), not compressed — and it is **peelable**: `lcm(2..X)` grows exponentially, so
it clears the far-peel threshold and drops to the base `{1..11,13}` (lonely by the LRC(≤13) citation). So the
two regimes are disjoint and both closed:

| regime | example | closer | why bounded/finite |
|---|---|---|---|
| **global ratio ≤ 13** | `{1..13}` | `spread13` | witness `t=1/(min+max)` |
| **compressed, no dominant, ≥7-far** | the chains above | **bounded-`q` census** | 13 runners ⟹ bounded band ⟹ `q≤~20` |
| **moderate dominant** (`22 < w ≲ 10⁶`) | deep well `182` | bounded-`q` census | `q=40` |
| **huge dominant** (`w ≳ 10⁶`) | `lcm(2..X)` | **far-peel** | drops to LRC(≤13) base |

## The dispatch, reframed

`hlarge = {spread13 (ratio ≤ 13)} ∪ {bounded-q census (everything of bounded/compressed shape)} ∪ {far-peel
(huge dominant)}` — **three routes, all finite, no renormalization tower needed for the compressed crux.** My
S54 routing (`hlarge_of_farcount`: `farCount ≥ 7 ⟹` obligation) is intact; the finding is that its `farCount ≥ 7`
obligation is discharged by the **census** (compressed) or **far-peel** (dominant), *not* the tower. The
renormalization engines (THM-608 `scale_separation`, `scale_separation_phase`) remain proved and available for
a near-equal cluster if one ever arises, but the compressed chain — the case that looked like it *needed* them
— does not.

**So "cluster identification" for the gap case is not a peel-target search; it is the dominant-vs-compressed
dichotomy:** one dominant huge runner ⟹ peel; otherwise ⟹ the bounded-`q` census. The census's completeness
(that every covering family, of any magnitude, is lonely at some `q ≤ Q₀`) is now a **finite residue-level
check** — loneliness at `q` depends only on `{v_i mod q}`, and covering constrains those — not an unbounded
search. That is kps's `lonely14_of_ratio` route; this finding says it suffices for the compressed branch.

## Honest scope

- **Verified:** compressed covering ≥7-far families lonely at `q ≤ 19` across magnitude `84 → 1.68×10⁶`
  (scaling test); deep well census-able at `3/40`. (`compressed_hge7_census_test…S55`.)
- **Heuristic (strong):** the `13 runners ⟹ bounded band ⟹ bounded q` counting, plus mac-mini's failed
  high-`q` compressed constructor, plus kps's `q ≤ 35` over 407 hard instances — three independent lines.
- **Open (but now finite, not unbounded):** proving the census bound `Q₀` universal is a residue-level finite
  computation (kps's census engine), not a new analytic tower. The far-peel threshold and its finite window
  (mac-mini's step-5) close the dominant side.

## Status

The gap-case dispatch's compressed "band-blocker" crux is **not a tower problem** — it is census-able at
bounded `q` because `n = 13` bounds the coverable band. hlarge reduces to `{spread13} ∪ {bounded-q census}
∪ {far-peel}`, all finite. The only remaining work is the census-completeness residue check (kps) and the
dominant-far window (mac-mini) — both finite, both in progress.

Related: MISTAKE-097 (the S52 correction); HYP-4050 (my hlarge routing, S54); HYP-4040/4041 (mac-mini
band-blockers / renormalization architecture — reconciled here); HYP-3984 (kps `lonely14_of_ratio` + the
`q ≤ 35` census); the deep-well pack (kps `deepWell_lonely`); far-peel (kps `far_peel_lonely`); THM-608 /
`scale_separation(_phase)` (opus, available but not needed for compressed). Script:
`04-computation/compressed_hge7_census_test_opus_20260703_S55.py`.
