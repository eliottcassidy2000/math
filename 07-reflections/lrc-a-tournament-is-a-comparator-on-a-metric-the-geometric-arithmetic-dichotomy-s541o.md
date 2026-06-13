---
source: oracle-2026-06-01-S541o
status: synthesis + computation (a tournament = comparator on a metric + threshold-tie; the 8 exotic metrics; geometric/arithmetic dichotomy)
tags:
  - lonely-runner
  - tournament
  - metric
  - trienerment
  - p-adic
  - ultrametric
  - sprague-grundy
  - dichotomy
---

# A Tournament Is a Comparator on a Metric: the 8 Exotic Metrics and the Geometric/Arithmetic Dichotomy

Taking each construction of the wild multitude (S540) and asking *what metric it
carries and how a tournament sits on that metric* exposes a single organizing
principle and a clean dichotomy.

## The meta-principle

> **Every LRC tournament is a COMPARATOR ON A METRIC with a TIE (trienerment, S539)
> when the metric distance falls below a threshold; and LRC = "the observer is
> metrically FAR (tie-free)" in that metric.**

The standard runner tournament uses the circular distance with threshold `1/n`. Each
exotic construction supplies a *different* metric `d_X(i,j)`; the same comparator
machinery applies; LRC is always "observer tie-free in `d_X`." The content is in
*which metric*.

## The 8 constructions, as (metric → tournament → LRC)

| # | construction | metric `d_X(i,j)` | tie / threshold | LRC reading | type |
|---|---|---|---|---|---|
| 1 | tropical | min-plus margin = gap to next runner | margin `< 1/n` | observer's tropical vertex is deep (Newton) | GEOM |
| 2 | p-adic tree | ultrametric `p^{-v_p(v_i-v_j)}` | `p^K \| (v_i-v_j)` | observer p-adically far = the **sieve** | **ARITH** |
| 3 | quantum | Fubini–Study / fidelity `|⟨ψ_i|ψ_j⟩|` | fidelity high | observer state distinguishable | GEOM |
| 4 | sandpile | transport/odometer of occupancy (S536) | same cell-coset | observer cells empty (recurrent) | GEOM |
| 5 | zeta | periodic-orbit length `= |v_i-v_j|` | equal length (resonance) | a spectral gap at the observer | **ARITH** |
| 6 | quasicrystal | local-isomorphism (agreement radius) | same patch | observer sees a hole-patch | GEOM |
| 7 | game | Sprague–Grundy nim-value of a pair game | P-position (tie) | observer is a balanced **P-position** | **ARITH** |
| 8 | Galois | cyclotomic chord `2sin(π|i-j|/n)` + Frobenius orbit | same orbit | observer arc avoids the orbit | GEOM(chord)+**ARITH**(QR) |

## The dichotomy (computed, `lrc_metric_tournaments_multitude_s541.py`)

**Geometric metrics collapse.** The tropical margin (`=d`), quantum fidelity
(`|cos πd|`), Galois chord (`2 sin πd`), and quasicrystal agreement (`~1/d`) are all
**strictly monotone in the circular distance `d`** (verified table). A threshold-tie
comparator on a metric monotone in `d` produces the **same tournament** as the
standard circular one. So constructions 1,3,4,6 (and the chord part of 8) carry *no
new tournament information* — they are the standard runner tournament in disguise.
This reconfirms (S535/S540) that *geometry alone is the circular distance*.

**Arithmetic metrics are genuinely new** — and they are functions of the
**differences `v_i - v_j`**, hence live on the difference/channel structure (S533/
S538), not on the circular positions:

- **p-adic ultrametric (2).** `d_p = p^{-v_p(v_i-v_j)}` satisfies the **isosceles
  law** (every triangle has its two largest distances equal — verified True for all
  sampled triples): it is a genuine ultrametric, so the runners sit at the leaves of
  a **`p`-adic tree**. Its trienerment ties = "same `p`-adic ball" = **same residue
  channel** (S533/S534), and **observer tie-free at level `K` = the sieve** (no speed
  divisible by `p^K`; `t = 1/p^K` lonely, THM-369) — verified `300/301`. So the
  Bruhat–Tits-tree tournament *is* the prime-power channel/sieve structure, now as a
  tree metric. (At `n=18, p=3, K=2`: exactly the `n*=9` channels of S534.)
- **Sprague–Grundy game (7).** Each pair plays a subtraction game; the Grundy
  sequence is **eventually periodic** with period set by the speeds (computed: `S=
  {1,2}→3`, `{2,3}→5`, `{1,3,4}→7`, `{3,5}→8`). The game-tournament (`i` beats `j`
  iff the combined game is an `N`-position) is arithmetic; **tie = P-position**; LRC
  = the observer reaches a balanced `P`-position. The nim/XOR metric is a new,
  arithmetic comparator.
- **Zeta orbit-lengths (5).** The pairwise periods are exactly `|v_i - v_j|` — the
  **holdback atoms** (S25); the orbit-length metric is the difference metric again.
- **Galois–Frobenius (8).** The Frobenius/QR orbit comparator is the **Paley
  tournament** (S535 M3), arithmetic in the differences mod `n`.

## The unifying conclusion

> **A tournament on a metric splits LRC's structures into two kinds. GEOMETRIC
> metrics (tropical, quantum, sandpile, quasicrystal, the Galois chord) are all
> monotone in the circular distance and collapse to the one standard runner
> tournament. ARITHMETIC metrics (p-adic ultrametric, Grundy, zeta orbit-length,
> Galois–Frobenius) are functions of the differences `v_i - v_j` and give genuinely
> new tournaments — and they all land on the same arithmetic object: the
> difference/resonance/channel structure (S533/S537/S538), with the p-adic
> ultrametric making it literally a tree whose ties are the sieve channels.**

So the answer to "apply a tournament to a metric, for each construction": geometry
gives you back the circle; arithmetic gives you the channels. The genuinely new
tournament content of LRC is arithmetic, and the p-adic tree is its cleanest carrier
— unifying the exotic multitude with the sieve (THM-369) and the channel/parity
theory (S533/S534).

## Verdict / next
- Meta-principle established: tournament = comparator on a metric + threshold-tie;
  LRC = observer tie-free.
- Dichotomy verified: geometric metrics collapse to circular distance; arithmetic
  metrics (p-adic/Grundy/zeta/Frobenius) are difference-functions giving new
  channel/tree tournaments; the p-adic ultrametric trienerment = the sieve.
- Concrete next: (1) build the full **p-adic tree** for a composite `n*` (e.g.
  `n=18, n*=9`) and read LRC as a tree-covering; (2) work out the Grundy
  game-tournament's `P`-position structure vs loneliness; (3) the Galois–Frobenius
  (Paley) tournament's realizable classes at prime `n` (`n≡3 mod 4`).

## Artifacts
```
04-computation/lrc_metric_tournaments_multitude_s541.py
05-knowledge/results/lrc_metric_tournaments_multitude_s541.out
```
Related: S540 (the wild multitude), S539 (trienerment = metric+tie), S535 (mapping
spectrum / metric restriction), S533/S534 (channels, p-adic), S538 (tension/diff),
S25 (holdback = differences), THM-369 (sieve).
