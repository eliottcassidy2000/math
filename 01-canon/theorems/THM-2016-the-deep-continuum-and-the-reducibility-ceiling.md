---
id: THM-2016
title: "THE DEEP CONTINUUM AND THE REDUCIBILITY CEILING — three worked threads on the continuum coordinates (THM-2013). (A) The L3 coordinate is LOCAL SUBTOURNAMENT DENSITY (the k-profile / flag-limit coordinate): global traces (cycle spectrum=char_A) and the global signed |R| collapse in the hot center, and the k=4,5 induced-subtournament census resolves most of the residue — BUT the deepest hot-shell residue (3/47 at the n=7 τ=6/7 shell) survives even (score, char_A, |R|, 4-profile, 5-profile): the continuum's center is invariant-resistant. (B) REDUCIBILITY CEILING (PROVED + verified): max c₃ over reducible tournaments = c₃_max(n−1), because c₃=Σc₃(SCC) and the largest strong component has ≤ n−1 vertices; hence the condensation temperature τ_c = c₃_max(n−1)/c₃_max(n) (= 1/2,2/5,5/8,4/7 for n=4..7) above which every class is strong. (C) H IS A THERMOMETER: mean/min/max H increase monotonically with cyclic temperature (mean 1→178 as τ:0→1 at n=7), fine spread carried by the free coordinates — locating death-star's H≥disc binding case at the hot center τ=1."
status: >
  VERIFIED (boxeph-2026-07-21-S200), exhaustive n≤7. (A) In the n=7 hot shells: (char_A,|R|) resolves
  36/47 (τ=6/7), 41/52 (τ=11/14); adding the 4-profile → 36/47, 50/52; adding the 5-profile → 44/47,
  50/52. So local density is the L3 coordinate, but 3/47 remain unresolved by all of
  (score,char_A,|R|,4-profile,5-profile) — invariant-twins in the deep center. (B) max-reducible-c₃
  = c₃_max(n−1) exactly for n=4..7 (1,2,5,8 = c₃_max(3..6)); c₃_max(n)=1,2,5,8,14 (n=3..7,
  = n(n²−1)/24 odd, n(n²−4)/24 even); τ_c=c₃_max(n−1)/c₃_max(n)→1. Proof: c₃(T)=Σ_SCC c₃(SCC) ≤
  c₃_max(largest SCC) ≤ c₃_max(n−1). (C) H strictly increases in mean with τ across all 15 shells at
  n=7 (1,3,6.7,10.5,14.4,23.9,33.4,44,60.9,79.7,96,114.5,134.7,153.9,178.3); the spread widens toward
  mid-hot. Uses THM-2013 (temperature), THM-1926 (cycle spectrum), THM-1966 (|R|), THM-1862 (c₃ over
  SCCs), death-star-S84 (H≥disc quasirandom).
source: boxeph-2026-07-21-S200 (owner: work the S199 open threads + one or two picked up from other agents)
depends_on: []
related:
  - THM-2013  # continuum coordinates (this resolves its L3 residue + condensation threshold)
  - THM-1979  # the spectrum (single point -> continuum)
  - THM-1966  # mac-mini |R| (the L2 coordinate that collapses in the deep continuum)
  - THM-1926  # my zeta / cycle spectrum (=char_A, the L1 coordinate)
  - THM-1862  # order-join: c₃ = Σ over strong components (the reducibility ceiling)
  - "07-reflections/coordinates-for-the-continuum-cyclic-temperature-and-the-cycle-spectrum-boxeph-S199.md"
script: 04-computation/continuum_threads_boxeph_S200.py (+ .out)
---

# THM-2016 — the deep continuum and the reducibility ceiling

## Thread A — the L3 coordinate is local subtournament density (and the center resists even that)

Extending the coordinate budget of THM-2013 (`L0` temperature, `L1` cycle spectrum = char_A, `L2`
`|R|`), the natural `L3` is **local structure** — the `k`-profile (census of induced `k`-vertex
subtournaments), the flag/tournament-limit coordinate. In the n=7 hot shells:

| shell | (char_A,\|R\|) | +4-profile | +4&5-profile |
|---|---|---|---|
| τ=6/7 (c₃=12, 47 cls) | 36 | 36 | **44** |
| τ=11/14 (c₃=11, 52) | 41 | **50** | 50 |
| τ=13/14 (c₃=13, 15) | 15 | 15 | 15 |

Local density genuinely resolves — the 4-profile cracks the τ=11/14 shell, the 5-profile the τ=6/7
shell. **But 3/47 survive even (score, char_A, |R|, 4-profile, 5-profile):** the center of the
continuum holds *invariant-twins* — non-isomorphic tournaments agreeing on temperature, the full cycle
spectrum, the signed Rédei count, and all 4- and 5-vertex densities. The deep continuum is
**invariant-resistant**: no cheap global (spectral/signed) or local (flag) coordinate separates it; it
needs essentially the whole object. This is the sharp form of "enumeration is forced at the hot core."

## Thread B — the reducibility ceiling (proved)

> **Reducibility ceiling.** The maximum `c₃` over **reducible** (non-strongly-connected) tournaments
> on n vertices is `c₃_max(n−1)`.

*Proof.* `c₃(T) = Σ_{strong components S} c₃(S)` (a 3-cycle lies inside one strong component). A
reducible T has ≥2 components, so its largest has ≤ n−1 vertices, giving
`c₃(T) ≤ c₃_max(n−1)`; equality by taking the max-cyclic (n−1)-tournament joined to a single vertex. ∎

Verified exactly: max-reducible-c₃ = `c₃_max(n−1)` = 1,2,5,8 for n=4..7. So the **condensation
temperature**
```
        τ_c = c₃_max(n−1)/c₃_max(n)   (= 1/2, 2/5, 5/8, 4/7 for n=4..7)
```
is exact — for τ > τ_c every class is strongly connected. Since `c₃_max(n) = n(n²−1)/24` (odd),
`n(n²−4)/24` (even), `τ_c → 1`: reducible tournaments persist to nearly the hot edge (as a vanishing
fraction, THM-1978), and only the sliver `(τ_c, 1]` is purely strong.

## Thread C — H is a thermometer (locating the H≥disc binding case)

The Hamiltonian-path count H, the flagship #P-hard invariant, is **macroscopically a function of
cyclic temperature**: mean H over the n=7 shells rises monotonically
`1, 3, 6.7, 10.5, 14.4, 23.9, 33.4, 44, 60.9, 79.7, 96, 114.5, 134.7, 153.9, 178.3` as τ goes 0→1,
with the *spread* (the free-coordinate contribution) widening through the mid-hot shells. So `τ`
predicts H's scale and the cycle-spectrum/|R| coordinates predict its fine value — consistent with
"scores explain ~97% of H-variance." This locates death-star-S84's result: `H ≥ disc` is saturated at
the **hot center τ=1** (regular/quasirandom), where H is largest and the inequality is tightest —
exactly where enumeration fails and the coordinate description takes over.

## Net

The continuum coordinates (THM-2013) now have a sharp top and bottom: **H is the thermometer** (the
temperature is even readable off the flagship invariant), the **reducibility ceiling** `c₃_max(n−1)`
pins the cold edge exactly, and the **deep hot center is invariant-resistant** — the one place where
local densities, the spectrum, and `|R|` all fail together, and there is nothing left to do but hold
the whole tournament.
