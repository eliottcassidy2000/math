# HYP-2180 — The iterated-log depth of a bound = the altitude in the log-tower where the dynamics linearizes

**Session:** claudebox-2026-06-03-S615. **Lineage:** HYP-2175 (Collatz↔LRC), Lemma A (the "almost-all / randomness" face).
**Prompt:** understand the deeper abstraction behind Tao-style loglog / logloglog inequalities and produce my own.

## The abstraction
Iterated logarithms are not magic constants; they count **how many times you must take a logarithm before the
process you are iterating becomes a geometric contraction.** Equivalently (the probabilistic reading): the depth =
the number of nested "for almost all / on average" reductions, because each averaging over a dyadic range of scales
turns "all scales" into "log-many scales."

Define the **log-tower** `L₀(N)=N, L₁=log N, L₂=log log N, …`. Let a renormalization step `R` act on scale.
- If `R` contracts `L₀` by a constant factor (`R(N) ≤ ρN`, 0<ρ<1) → reaching O(1) costs `Θ(L₁(N)) = log N` steps.
- If `R` contracts `L₁` (`log R(N) ≤ ρ log N + C`, e.g. `R=√N`) → costs `Θ(L₂(N)) = log log N` steps.
- If `R` contracts `Lⱼ` → costs `Θ(L_{j+1}(N))`.

**Altitude principle.** The iterated-log depth of the iteration bound is `j+1`, where `j` is the tower-level at which
`R` is a contraction. The *coefficient* of the leading iterated log is `1 / log(1/ρ)`, ρ = the contraction ratio at
that altitude. (Proof core formalized: `Math/IteratedLog/Altitude.lean`, the geometric-with-offset descent
`aᵢ ≤ ρⁱa₀ + C/(1−ρ)` — the altitude falls geometrically, so it reaches O(1) in `log(a₀)/log(1/ρ)` steps.)

## My own instance — a genuine DOUBLE LOG hiding inside Collatz (verified)
The raw shortcut Collatz step-count is a **single** log: the value drifts by `√3/2` per step (level-0 contraction),
so an orbit reaches O(1) in `~log₂ n` steps (measured steps/log₂n ≈ 4.8 ≈ 1/0.2075). Now coarse-grain. Define an
**epoch** = run the map for (current bit-length)-many steps. Over one epoch the *bit-length* contracts:
`log₂(value)` drops by factor `ρ ≈ 1 − 0.2075 = 0.79` per epoch (level-1 contraction). Hence:

> **Epoch double-log (heuristic, verified).** For almost every n, the number of epochs to reach O(1) is
> `~ (1/log₂(1/ρ))·log₂log₂ n ≈ 2.98 · log₂log₂ n`.

Measured: epochs = `2.818·log₂log₂ n − 3.6`, **R²=0.9987** over n of 16…1024 bits (slope vs predicted 2.98 — the
small deficit is the shrinking final epochs / threshold). The single log counts *map steps*; the double log counts
*epochs*. **The two logs of Collatz are step-altitude and epoch-altitude — successive floors of the same tower.**

## The altitude ledger for Tao's four logs (loglogloglog n)
Each log in "almost all orbits attain almost bounded values" is one accounting layer, one floor up the tower:
1. `log` — bit-length of the value (the address of n).
2. `log log` — epoch count: the bit-length itself contracts geometrically (this note's mechanism).
3. `log log log` — union bound over the ~log-many dyadic scales on which the "almost all" estimate is run.
4. `log log log log` — the slack from letting the target threshold f(n) be *any* function → ∞.
Four nested averagings ⇒ four logs. The depth is the count of "for-almost-all" quantifiers, nothing more exotic.

## Collatz↔LRC (extends HYP-2175): same altitude, both faces
LRC's "almost all configs are very lonely" (random gap `≫ 1/n`, the Lemma-A face) carries the identical signature:
the gap-deficit is controlled by a union bound over dyadic frequency scales — a loglog in the exceptional-set
measure. Resonance energy `E(v)=Σ∏|ĝ(mᵢ)|` summed over the relation lattice's successive minima (log-spaced) is the
LRC altitude tower. **Conjecture (LRC altitude):** the energy-bound failure depth (HYP-2165, dies at n=8) is a
level-1 phenomenon, so the honest almost-all LRC gap improvement is `loglog`-deep, matching Collatz's epoch log.

## Verified / formalized
- `04-computation/iterated_log_altitude_s615.py`: altitude sanity (√N → loglog exactly); Collatz epoch double-log
  fit R²=0.9987; step-count single-log contrast.
- `Math/IteratedLog/Altitude.lean` (math-lean): the geometric-with-offset descent lemma — the provable heart of the
  altitude principle (altitude falls geometrically ⇒ iterated-log iteration count).

## Open
- Make the epoch double-log a theorem (it needs Tao's almost-all drift, currently heuristic + strong numerics).
- The LRC altitude conjecture: is the almost-all gap improvement provably loglog-deep via the relation-lattice minima?
- log* / Ackermann floor: a process contracting at level `j(N)` growing with N (a "moving altitude") yields `log* N`.
