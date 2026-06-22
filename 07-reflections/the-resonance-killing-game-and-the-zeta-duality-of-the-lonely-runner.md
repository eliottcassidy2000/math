# The resonance-killing game: a proof↔disproof dialectic, and the ζ(2)/ζ(−1) duality

*kind-pasteur-2026-06-22-S31p. The owner asked to run a proof↔disproof dialectic over finite/infinite
families, inspired by `M({1..11,13}) = 1/12` and `1+2+3+… = ζ(−1) = −1/12`. The dialectic converges on
the project's own consec-extremality + ζ(2) floor, seen through a new "resonance-killing" lens.*

## The resonance hierarchy (the board the game is played on)
For a speed set `S`, `f(t) = min_{s∈S} ‖st‖` and `M(S) = max_t f(t)`. Computed exactly
(`lrc_resonance_hierarchy_kps.py`), the local maxima of `f` (the **witnesses**) sit at rationals
`t = a/b`, and the lonely VALUE there is `1/b` whenever some speed `≡ ±1 (mod b)`. So the witnesses form
a **harmonic ladder** `1/b`, `b = 14, 15, 16, …` for the AP, and

> **`M(S) = 1/(smallest resonance b that SURVIVES)`**, where `b` is *killed* iff some `s∈S` has `b | s`
> (that runner sits at the origin at `t = a/b`).

The AP `{1..13}` kills `b = 1,…,13` — each by the single speed `s = b` — leaving `b = 14` as the first
survivor: `M = 1/14`. Verified for the whole infinite family (`lrc_infinite_ap_family_kps.py`):
**`M({1..n}) = 1/(n+1)` exactly, n = 2..16.** The `1/12`-core `{1..11,13}` skips the killer `12`, so
`b = 12` survives: `M = 1/12`. Adding `24` (`= 0 mod 12`) re-kills `b = 12` → drops to `b = 14` → the GW
tiler at `M = 1/14`.

## The disproof, and why it is self-defeating
A speed `s` kills EVERY resonance `b | s` at once, so the obvious attack on LRC(14) is: use
highly-divisible speeds to kill `b = 1,…,14` simultaneously, leaving only `b ≥ 15` → `M ≤ 1/15 < 1/14`,
a counterexample. **This fails decisively** (`lrc_disproof_killer_attack_kps.py`): every set that kills
all of `b = 1..14` has `M ≥ 1/8 ≫ 1/14` — *larger*, not smaller (`{2..14}`: `1/8`; divisible spreads:
`1/8`–`1/5`; big killers: `16/83`). The mechanism is a **tension built into the board**:

> To kill resonance `b` you need a speed `≡ 0 (mod b)` — i.e. a *spread-out / large* speed. But
> spread-out speeds open *bigger gaps* in the orbit `{st}`, which *raises* `M`. Killing and compactness
> pull in opposite directions; you cannot kill the small resonances without spreading, and spreading
> re-lonelies the set.

The AP `{1..13}` is the **unique least-spread killer of `{1,…,13}`** (the speeds ARE `1,…,13`, the
minimal multiples). Any attempt to also kill `14` forces extra spread and pushes `M` back above `1/14`.
This is exactly the project's *consec-extremality* (the sector route's consec-max, THM-557/HYP-2607),
now read as a min–max of the killing game: **the AP minimizes `M` because it is the tightest packing of
resonance-killers.** The disproof's failure IS the proof's mechanism.

## Why the floor never reaches 0 — the ζ(2) half
Killing a resonance at the *exact* point `a/b` (a divisible speed) does **not** kill its lonely
*neighborhood* (width `~1/(bV)`): a large killer has a *thin* danger arc `1/(7s)` that misses the
neighborhood. The surviving neighborhoods near the small-denominator Farey points `{0,1/2,1/3,2/3,…}`
carry a mass bounded below by `Σ_{b<q} 2δ_b → 3/π² = 1/(2ζ(2))` (HYP-2856,
`zeta2-governs-the-lonely-runner-floor.md`). So a lonely point always survives: **`ζ(2)` is the density
of resonances the runners cannot all dodge.** The disproof, pushed to its limit (covering sets), runs
straight into this floor — which is precisely where HYP-2895's equidistribution lives (large committed
speeds equidistribute, their thin arcs cover only `~1/7` of the lonely set).

## The ζ(−1) = −1/12 half (the owner's hint, placed)
The lonely values are the reciprocals `1/b` of the killed/surviving resonances, and the AP kills the
initial segment `{1,…,n}` with killer-sum `1+2+…+n = n(n+1)/2`. The **infinite** AP `{1,2,3,…}` kills
*every* resonance, `M → 0`, and its regularized killer-sum is `1+2+3+… = ζ(−1) = −1/12`. So:

> `ζ(−1) = −1/12` is the *infinite-family limit* of the killing game — the regularized "total killing
> capacity" at which loneliness vanishes — and `M({1..11,13}) = 1/12` is its **finite avatar** (one
> killer removed, the first resonance `b = 12` re-opens with lonely value `1/12`).

And the two zeta values are **functional-equation duals**: `ζ(−1) = 2(2π)^{−2}cos(π)Γ(2)ζ(2) = −1/12`
sends `ζ(2) = π²/6 ↦ −1/12`. This is not decoration — the repo already found the discrepancy/cap side
is a **second-Bernoulli / E₂** object (HYP-2745, `B₂ = 1/6`, `−ζ(−1) = B₂/2 = 1/12`) and the Dedekind
ladder (THM-563), while the floor side is the `ζ(2)` Farey density. The LRC bound `floor ≥ cap` is the
two sides of one zeta seen across `s ↔ 1−s`: **`ζ(2)` counts the resonances (floor), `ζ(−1)` weights the
lonely values (cap), and the functional equation is why a Diophantine-avoidance density and a
Bernoulli–Dedekind discrepancy describe the same `1/14`.**

## Net (proof and disproof, each strengthening the other)
- **Disproof → Proof:** every way to kill the small resonances spreads the speeds and *raises* `M`; the
  AP is the unique tight balance at `1/14`. The disproof is not just refuted — its failure mode *is* the
  consec-extremality proof.
- **Proof → Disproof:** the only sets that can even *try* for `M < 1/14` are the covering sets (kill all
  `b ≤ 14`, THM-523); the disproof correctly identifies them as the sole battlefield, and they are
  exactly where the `ζ(2)` neighborhood floor (and HYP-2895 equidistribution) holds the line.
- **Not a new proof** — it re-narrates consec-extremality + the ζ(2) floor through the resonance-killing
  lens, and places the owner's `ζ(−1) = −1/12` as the infinite-AP regularized killer-sum, functional-
  equation dual to the `ζ(2)` floor. The residual is unchanged: the rigorous `ζ(2)` neighborhood-width
  floor over all covering clusters (OPEN-Q-108 / the witness route).

→ `zeta2-governs-the-lonely-runner-floor.md`, HYP-2856, HYP-2895, THM-523, THM-560, THM-563, HYP-2745,
THM-557/HYP-2607 (consec-extremality), [[lrc14-thread]].
