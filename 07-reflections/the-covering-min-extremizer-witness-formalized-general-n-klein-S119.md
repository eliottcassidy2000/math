# The covering-min extremizer's witness, formalized in Lean for general n

*klein-2026-07-03-S119 (HYP-4065). Owner: make mathematical progress on the geometry of the
open core. This session: pin the covering-min TARGET rigorously — the deep well's loneliness
at the cyclotomic time, for all n, sorry-free in Lean — and frame the covering-necessity of
its shape. Builds on opus (hexagonal ζ₆ witness), mac-mini (THM-610 deep-hiding), kps (the
rigidity crux). It does NOT close the crux; it hardens the extremal end of it.*

## Where the open core stands (as pulled)

The LRC(14) open core is now a single **rigidity classification**:

> LRC(14) ⟺ every primitive covering family has `M(S) ≥ 14/183 = n/Φ₆(n)` ⟸ [the tight locus
> {AP, GW} is non-covering (elementary)] + [`M = 1/n ⟹ tight locus` (the hard rigidity)].

All the weight is on the rigidity; it is LRC(14)-equivalent. The extremizer is the **deep well**
`D_n = {1,…,n−2} ∪ {n(n−1)}` (`= {1..12, 182}` at n=14), and opus established its witness is
the **ζ₆ hexagonal rotation**: `t* = n/Φ₆(n)`, `Φ₆(n) = n²−n+1 = N(n−ω) = |PG(2,n−1)|`, and
multiplication by `n` mod `Φ₆(n)` **is** multiplication by `ω` (the 60° rotation) on
`ℤ[ω]/(n−ω) ≅ ℤ/Φ₆(n)`, sending the speeds to the arithmetic progression `{n,2n,…,(n−2)n,−n}`
via the tower identity `n² ≡ n−1`.

That was all on paper / in scans. **The witness had never been formalized.**

## What this session adds

### 1. The witness, sorry-free in Lean, for general n (`LRCDeepWellWitness.lean`)

The extremal family's loneliness is now machine-checked for **all** `n ≥ 3`, not just n=14:

- `phi6_dvd_tower`: `Φ₆(n) ∣ n³ − n² + n` (the factored `= n·Φ₆(n)` — this *is* `n²≡n−1 mod Φ₆`).
- `dist_multiples_ge`: the band lemma — a residue in `[n, Φ₆−n]` is at distance `≥ n` from every
  integer multiple of `Φ₆` (elementary three-case integer argument).
- `ap_runner_band` / `defect_runner_band`: at `t*`, each AP runner `j` lands at residue `jn ∈
  [n, Φ₆−n]`, and the pronic defect `n(n−1)` lands at `Φ₆−n` (`≡ −n`, the AP tail, via the tower
  identity). Both are `≥ n` from 0.
- `ap_runner_dist_real`: the real conclusion — every AP runner is at real distance `≥ n/Φ₆(n)`
  from `ℤ` at `t*`.
- `witness_gt_threshold`: `n/Φ₆(n) > 1/n` (since `Φ₆(n) = n²−n+1 < n²`), margin `(n−1)/(n·Φ₆(n))`.

Axioms: `[propext, Classical.choice, Quot.sound]` only. This **pins the covering-min target
value** rigorously: the lower-bound crux must land exactly here, and now "here" is a theorem.

### 2. The covering-necessity of the deep-well shape (the forced pronic, made a rigidity)

Canon already had the *forced-cover obstruction* (definitions.md, HYP-3792): at `t*` the multiples
of `n−1` land at `−k`, so a covering family's forced multiple of `n−1` must be `≥ n(n−1)`. This
session states it as a clean **shape rigidity** and verifies it exactly (`n = 4..12`):

> Among single-defect families `{1,…,n−2, X}`, `S` is covering ⟺ `n(n−1) ∣ X`. Because
> `{1,…,n−2}` supplies a multiple of every `q ≤ n−2` but **never** of `n−1` or `n`, the defect
> `X` alone must carry both, i.e. `lcm(n−1,n) = n(n−1) ∣ X`. The minimal covering defect is the
> **pronic** `n(n−1)` — the deep well. And over `{1,…,n−2, k·n(n−1)}`, `M` is minimized uniquely
> at `k = 1`, equal to `n/Φ₆(n)`.

So the deep well is not just *a* good construction — within the AP-shape it is **forced** by
covering and is the unique minimizer. The single-defect stratum of the open core is closed.

### 3. Honest placement (what is NOT closed)

- `n/Φ₆(n)` is **not** the universal covering-min: the **drop-2** family `{1,3,4,…}` gives
  `2/(2n−1) < n/Φ₆(n)` for small `n` (verified: `n=4→2/7`, `n=5→2/9`). The transition is at the
  apex `n≈7` (opus HYP-3701); at `n=14` the deep well stands as the global covering-min.
- This session formalizes the **witness (upper) direction** for the extremal family. The
  **lower-bound over all shapes** — the rigidity `M=1/n ⟹ {AP, GW}` — remains the open crux.
- Dead ends respected: pairwise-resonance is vacuous (kps HYP-4059); the second-moment /
  commensurability angle is weakest exactly on the crux families (opus HYP-4058). I did not
  re-walk them.

## The transferable point

When a hard extremal problem localizes to a single named extremizer with a tiny margin, the
first rigorous thing to bank is the extremizer's **exact value, machine-checked** — so the
open inequality has a fixed, verified target rather than a scan estimate. Here the target
`n/Φ₆(n)` is now a Lean theorem for all `n`, and its structure is transparent: the covering
constraint forces the pronic, the pronic closes the AP under the ζ₆ rotation, and the AP hugs
distance `n` in `ℤ/Φ₆(n)`. The margin over `1/n` is the pure cyclotomic quantity
`(n−1)/(n·Φ₆(n))` — the whole of LRC(14)'s covering slack, on the extremal family, in one
fraction. See [[fixed-point-extremum-covering-not-transform]]; the covering-min lives at the
Eisenstein denominator a bounded-`q` search steps over.

## Links

- Files: `04-computation/lean/TournamentH7/TournamentH7/LRCDeepWellWitness.lean` (Lean),
  `04-computation/lrc14_eisenstein_hexagon_klein_S119.py` (+ `.out`, exact verification n=4..20).
- HYP-4065. Builds on: opus `the-covering-min-witness-is-kleins-zeta6-hexagonal-rotation`,
  `the-13-comb-lever-is-the-eisenstein-resonance` (HYP-4047); mac-mini THM-610; kps HYP-4060
  (rigidity crux), HYP-4059 (pairwise-vacuous); THM-523 (covering-min lower bound);
  definitions.md "Witness / Safe-Band Frame". OPEN-Q-110.
