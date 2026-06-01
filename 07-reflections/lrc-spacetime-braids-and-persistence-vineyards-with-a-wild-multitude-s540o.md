---
source: oracle-2026-06-01-S540o
status: synthesis + computation (exotic LRC lifts: the spacetime pure braid; persistence vineyards; a wild posed multitude)
tags:
  - lonely-runner
  - braid-group
  - pure-braid
  - persistent-homology
  - vineyards
  - configuration-space
  - exotic
  - multitude
---

# Spacetime Braids and Persistence Vineyards: LRC Lifted Far Out of the Box

Two genuinely exotic homes for LRC — the **pure braid group** and **persistent
homology** — each reading the conjecture as a single structural feature, and each
*independently recovering the same invariants* we built by elementary means
(tension, holdback, apex). Plus a wild posed multitude.

## Star 1 — the spacetime pure braid

The runner worldlines `x_i(t) = frac(v_i t)` on the cylinder `S^1 × [0,1]` (time ×
circle) are `n` strands; the observer is the vertical strand at `0`. Over `[0,1)`
they form a **pure braid** on `n` strands — at `t=1` every `frac(v_i) = 0`, so all
strands return to their start (verified: all sampled systems are pure braids). Three
facts, all verified (`lrc_braid_persistence_exotic_s540.py`):

- **Linking = tension.** The signed linking number is `lk(i,j) = v_i - v_j` (`i`
  laps `j` exactly `v_i - v_j` times). The **linking matrix** `L = v⊗1 - 1⊗v` is a
  coboundary — it obeys the cocycle `L_{ij}+L_{jk}+L_{ki}=0` (verified). So the
  braid's linking data is **exactly the S538 tension / difference-speed structure**.
- **Word length = holdback.** The braid word (the time-ordered sequence of adjacent
  transpositions at crossing times `t = m/(v_i-v_j)`) has length
  `Σ_{i<j} |v_i - v_j|` — the **holdback** of S25 (avg `121, 200, 346` for `n=5,6,7`).
- **LRC = a fat channel.** Loneliness is a *time-slice where the observer strand is
  at circular distance `≥ 1/n` from every other strand* — a **fat tube** around
  strand `0` in the braid. Reached `60/60` for the sampled sets.

So LRC lives in the **configuration space of `n` points on the annulus** (`= B_n` of
the cylinder), and the runner system is the *homogeneous / torus* pure braid with
linking matrix `v_i - v_j`. The realizable braids are a thin (torus-type) slice of
`P_n`, and LRC asks every such braid to have a fat observer-tube. The braid is the
spacetime form of the wiring diagram (S535 MAP-wire), now with full linking data.

## Star 2 — persistent homology and the observer vineyard

At each time `t` the `n` points are a cloud on the circle; the Rips **`H_0` barcode
is exactly the gap structure** — a gap `g` merges its two sides at radius `g/2`, so
the death-radii are `{g_i/2}` and the **longest bar = the largest gap = the apex**
(S530). The observer's component persists until

```
r_obs(t) = (min flanking gap of the observer)/2 .
```

Loneliness is `r_obs(t) ≥ 1/(2n)`, so:

> **LRC@n ⟺ `max_t r_obs(t) ≥ 1/(2n)`** — the observer's `H_0` bar outlives radius
> `1/(2n)` at some time.

As `t` varies the diagram moves: a **vineyard**. The observer's "vine" oscillates,
and LRC = *the vine crosses the line `r = 1/(2n)`*. Verified: avg `max_t r_obs =
0.179, 0.162, 0.152` against thresholds `0.100, 0.083, 0.071` for `n=5,6,7`; reached
`60/60`. Example barcode at the optimal time (`n=5`): `[0.017,0.052,0.052,0.19,0.19]`
— the two long bars are the observer's two flanking gaps (the fat channel again).

## Why both matter: the same invariants from alien mathematics

The braid (a *topological/group-theoretic* object) and persistence (a
*homological/metric* object) independently reproduce the trio we found elementarily:
**linking = tension (S538)**, **word length = holdback (S25)**, **longest bar = apex
= largest gap (S530)**. Two completely different lifts converging on the same data is
strong evidence these are the *intrinsic* invariants of the runner system — and that
LRC is the single statement "the observer-tube is fat / the observer-vine is tall."

## The wild multitude (posed)

Each lifts LRC into a different exotic category, reading loneliness as one feature:

- **Tropical.** Newton polygon / tropical curve of the speed polynomial; min-plus
  comparisons; LRC = a tropical cell containing the observer.
- **p-adic Bruhat–Tits tree.** Speeds' `p`-adic valuations (`p | n*`; `p=3` for
  `n=18`, `n*=9`, S534); the dynamics = a walk on the tree; loneliness = reaching a
  far vertex — the tree form of the prime-power channels.
- **Quantum / operator.** Unitaries `U_i = e^{2πi v_i t}`; a tournament from
  operator ordering / non-commutativity defect; LRC = a near-classical (commuting)
  window.
- **Abelian sandpile.** Sector occupancy (S536) = chips on the cycle `C_n`; toppling;
  sandpile group `ℤ_n`; LRC = a recurrent configuration with the observer cells empty.
- **Dynamical zeta.** The wall-crossing flow on the torus; periodic orbits at
  rational `t`; the Ruelle zeta / resonances; LRC = a spectral gap.
- **Quasicrystal hull.** Irrational speeds → cut-and-project set; the hull's local
  patches; LRC at the rational boundary = a forbidden patch.
- **Sprague–Grundy game.** Each runner-pair = a subtraction game of period
  `|v_i-v_j|`; Grundy values → a game-tournament; LRC = a `P`-position of the observer.
- **Galois / Frobenius.** The `n`-th roots of unity with `Gal(ℚ(ζ_n)/ℚ)`; the speed
  action; Frobenius orbits; LRC = an orbit avoiding the observer arc.

## The organizing meta-principle

> Across every exotic lift, LRC becomes **"a single feature is always attainable"** —
> a fat tube (braid), a tall vine (persistence), a tropical cell, a far tree vertex,
> a commuting window, a recurrent empty-observer sandpile, a spectral gap, a `P`-
> position, an avoiding orbit. The braid and persistence lifts both recover
> tension/holdback/apex, certifying these as the runner system's true invariants. The
> conjecture is robust to the language; its difficulty (S519/S529) survives every
> translation, but each gives a different *tool* (braid groups, vineyard stability,
> tropical convexity, building theory, operator norms, sandpile groups, Ruelle
> spectra, CGT, Galois) to attack the *same* fat-tube/tall-vine statement.

## Verdict / next
- Two computed exotic lifts: the **spacetime pure braid** (linking = tension, word =
  holdback, LRC = fat observer-tube) and **persistence vineyards** (`H_0` = gaps,
  LRC = observer vine crosses `1/(2n)`), both verified and both recovering the
  intrinsic invariants.
- A wild multitude posed (tropical, p-adic tree, quantum, sandpile, zeta,
  quasicrystal, Grundy, Galois).
- Concrete next: (1) compute the **annular braid group** image of the runner systems
  (which torus pure braids are realizable) and test fat-tube via braid invariants
  (Burau/Lawrence–Krammer); (2) use **vineyard stability** (bottleneck distance) to
  bound how the observer vine moves with the speeds; (3) build the **sandpile**
  model on `C_n` and test the recurrent-empty-observer characterization.

## Artifacts
```
04-computation/lrc_braid_persistence_exotic_s540.py
05-knowledge/results/lrc_braid_persistence_exotic_s540.out
```
Related: S538 (tension/pairs), S25 (holdback), S530 (apex/largest gap), S535
(wiring/mapping spectrum), S536 (sectors/sandpile substrate), S537 (flows),
S539 (trienerment/Gabor), S534 (p-adic channels).
