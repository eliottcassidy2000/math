---
id: THM-2013
title: "COORDINATES FOR THE CONTINUUM — describing near-regular tournament structure without enumeration. CYCLIC TEMPERATURE τ=c₃/c₃,max=(σ²_tr−σ²)/(σ²_tr−σ²_min)∈[0,1] is the parity-uniform macroscopic score coordinate; σ²_min=0 for odd n and 1/4 for even n. ISO-CYCLIC SHELLS stack tournament space and STRUCTURAL ENTROPY peaks near τ≈0.7 at n=7. The CYCLE SPECTRUM resolves within a shell, with N₃ frozen and N₄ first free. CORRECTED ABSOLUTE-R BUDGET: L2 is the magnitude |R| from THM-1966, not signed R; at the n=7 c₃=12 shell char_A resolves 21/47 and (char_A,|R|) resolves 28/47, not 36/47. Local 4&5-profiles then reach 41/47 (THM-2016)."
status: >
  VERIFIED (boxeph-2026-07-21-S199), exhaustive over all iso classes n≤7. (1) N₁=N₂=0, N₃=3c₃ never
  varies within a score sequence; N₄ is the first free moment (varies from n=5). (2) Cyclic
  temperature τ=c₃/c₃,max, with c₃,max=n(n²−1)/24 for odd n and n(n²−4)/24 for even n;
  shell entropy S(τ)=log₂|shell| peaks at τ*≈0.7 (n=7: 5/7). The first all-strong shell
  τ_all has c₃≥9 (n=7, τ=9/14), c₃≥6 (n=6),
  c₃≥3 (n=5). (3) Corrected coordinate budget in the hot shells: cycle spectrum (=char_A) resolves 21/47 of the
  τ=6/7 n=7 shell, +|R| resolves 28/47 (and 13/15 of the τ=13/14 shell). The old 36/47 and 15/15
  rows used signed R while labeling it |R|. The affine c₃(σ²) law is
  THM-1979; N_k are the zeta moments THM-1926; |R| is mac-mini THM-1966. The contribution is the
  COORDINATE SYSTEM + terms (cyclic temperature, iso-cyclic shell, cycle spectrum, structural entropy,
  condensation threshold) that describe the continuum without enumeration.
source: boxeph-2026-07-21-S199 (owner: find reframes / invent terms & lenses for the continuum so we stop enumerating and focus on near-regular behavior)
depends_on: []
related:
  - THM-1979  # tournament space is a spectrum (single point -> continuum) — this coordinatizes the continuum
  - THM-1926  # my zeta / cycle moments N_k (the cycle spectrum)
  - THM-1966  # |R| refines spectrum from n=6 and becomes independent of (spectrum,H) at n=7
  - THM-1960  # opus modular seeds (the cold/reducible rim, through τ_red)
  - THM-2010  # kps sequence-invariant catalog (|cyc|,|R|,|disc| = free-coordinate sequences)
  - MISTAKE-214  # signed-R carrier was mislabeled as the invariant magnitude |R|
  - "07-reflections/coordinates-for-the-continuum-cyclic-temperature-and-the-cycle-spectrum-boxeph-S199.md"
script: 04-computation/continuum_coordinates_boxeph_S199.py (+ .out)
---

# THM-2013 — coordinates for the continuum

The continuum (near-regular interior, THM-1979) is too big to enumerate (47 classes in one n=7 shell,
6880 tournaments at n=8). These coordinates describe it intrinsically.

## The coordinates

- **Cyclic temperature** `τ = c₃/c₃,max = (σ²_tr−σ²)/(σ²_tr−σ²_min) ∈ [0,1]` — one
  macroscopic coordinate from the scores alone (THM-1979), where
  `σ²_tr=(n²−1)/12` and `σ²_min` is `0` for odd `n`, `1/4` for even `n`.
  Thus `τ=0` is transitive and `τ=1` is the regular (odd `n`) or near-regular
  (even `n`) maximum-cyclic edge.
- **Iso-cyclic shell** `𝒮_τ` = classes at fixed `τ` (fixed `c₃`). Tournament space is a stack of shells.
- **Structural entropy** `S(τ) = log₂|𝒮_τ|`.
- **Cycle spectrum** `(N₄,…,N_n)`, `N_k = tr(Aᵏ)` (zeta moments, THM-1926): `N₁=N₂=0`, `N₃=3c₃`
  **frozen** by `τ`, overtones **free** (first free = `N₄`).

## The two revealed features (verified n≤7)

1. **Diversity maximum at intermediate temperature.** `S(τ)` peaks *inside* the hot edge, not at the
   center: at n=7 the peak is `τ=5/7` (`c₃=10`, 79 classes, `S=6.30`), while `τ=1` holds only 3
   classes. The fattest part of the continuum is `τ*≈0.7`.
2. **First all-strong shell.** Strong-fraction becomes 1 at `τ_all`: every class
   with `c₃ ≥ 9` (n=7, `τ≥9/14≈0.64`) is strongly connected. Reducible
   classes reach the attained ceiling `τ_red=4/7`, one integer shell below.
   (`τ_all`: 3/5 at n=5, 3/4 at n=6, 9/14 at n=7.)

## The coordinate budget (how few numbers pin a near-regular tournament)

```
   L0  cyclic temperature τ        1 real, from the scores  — places the shell
   L1  cycle spectrum N₄…N_n       = char_A                 — resolves 21/47 of the τ=6/7 n=7 shell
   L2  beyond-spectral |R|         mac-mini THM-1966        — (τ,char_A,|R|) → 28/47, and 13/15 at τ=13/14
```

`(τ, char_A, |R|)` improves the spectral address but leaves a substantial residue at the very center;
THM-2016 shows that local profiles refine it further. So near-regular tournaments have a layered address —
temperature + a short cycle spectrum + `|R|` — and the continuum is a coordinate cloud with a
temperature axis, an entropy profile peaking at `τ*≈0.7`, and a first
all-strong shell at `τ_all=9/14≈0.64` for n=7.

## Two lenses

**Thermodynamic:** `τ` = temperature, `S` = entropy, transitive = `T=0` ground state, regular =
quasirandom hot phase, score-spread = order parameter, the n=7 perfection-break = the phase
transition, `τ_all` = first all-strong shell, `τ*` = the specific-heat-like diversity peak. Reduction principles
are the low-`τ` expansion; the continuum is beyond its radius (death-star-S84: H≥disc saturates at the
quasirandom `τ=1`). **Harmonic:** `N_k = Σλ_jᵏ`, the cycle spectrum is the char-poly, `N₃` the frozen
fundamental, `N₄⁺` the overtones/timbre; where it collapses (cospectral) `|R|` is the beyond-harmonic
coordinate.
