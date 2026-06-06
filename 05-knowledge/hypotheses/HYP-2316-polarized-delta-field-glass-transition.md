---
id: HYP-2316
title: The polarized delta field = gradient of a frustrated (antiferromagnetic) energy; a glass transition at the even-n 2-adic seam
status: OPEN (synthesis); glass transition n=5→6 + forbidden-H growth VERIFIED; antiferro/Walsh = canon
source: claudebox-2026-06-03-S637
related:
  - HYP-2311  # H,delta,dichromatic = one conflict graph Ω(T)
  - THM-290   # H is antiferromagnetic (log-submodular, frustrated)
  - THM-163   # H Walsh even-degree; THM-260 bandlimited; THM-259 degree
  - THM-338   # impossible H values; THM-115 (H=21)
  - THM-254   # instant-MCMC (the H energy landscape, sampling)
---

# HYP-2316 — the polarized delta field and its glass transition

The sharpest open question from S636, pursued. The **delta field** `δ = dH` is the **discrete
gradient of `H` on the tournament hypercube** `Q_m` (`m=C(n,2)`; vertices = tournaments, edges = arc
flips). `H = #Hamiltonian paths = I(Ω,2) = OCF` (Rédei: odd ⟹ δ even). By THM-290 `H` is
**antiferromagnetic** (external field `λ_k=log(1+2^{s-1})`, all couplings `λ_{ij}<0`, frustrated),
Walsh **even-degree** & **bandlimited** (THM-163/259/260). So:

> **δ is the gradient of a frustrated Ising/Potts energy `log H` on the spin cube.** "Polarized"
> = the gradient is signed; the polarization points toward the H-maximum (the regular tournament).
> The **second-order propagation operator (S636) = the discrete Hessian = the coupling matrix**
> (antiferromagnetic).

## What it does at larger n (VERIFIED)

**1. Glass transition at n=6 (the headline).** Local maxima of `H` and their values:
- **n=5:** all 64 local maxima have `H=15` = the global max (the 64 regular tournaments). **Single
  basin, no metastability.**
- **n=6:** **720 metastable local maxima at `H=37 < 45`** (the global max), plus 480 global. So the
  landscape becomes **rugged** — the spin-glass signature of frustration appears exactly at `n=6`.

`n=6` is **even** = the **2-adic seam**. So *the onset of metastability is the even-`n` hardness*: the
polarized delta field has a single basin for `n=5` and a glassy, multi-basin landscape from `n=6` on.
The LRC frontier `n=14` (even) sits deep in the rugged phase.

**2. The forbidden-H holes proliferate.** Non-realizable odd `H`: `{7}` (n=5) → `{7,21,35,39}` (n=6)
→ `{7,21,63,107,119,149,…}` (n=7). These are the unrealizable independence vectors of `Ω` (THM-338/115);
the metastable trap `H=37` is flanked by forbidden `35, 39` — the holes shape the basins.

**3. Defects = level edges.** Arc flips with `δ=0` (the gradient cancels) average ≈ 2 per config
(n=5,6) — the **frustration loci** where polarization vanishes; the field's "vortices".

## Where polarized delta fields appear (the abstraction)

The same object — *the signed discrete gradient of a partition function / energy on a configuration
cube* — recurs:
- **Boolean-function analysis:** `δ_k = ` the discrete derivative / **influence** of `H`; band-
  limitedness (THM-260) = `H` is low-degree ⟹ δ is sparse/structured; polarization = the bias
  (degree-1 mass; here zero in the ±1 model since `H` is even, so δ is purely odd-degree).
- **Spin glasses / frustration:** `H` antiferromagnetic (THM-290), rugged landscape (metastable
  states from n=6) — a genuine glass. The instant-MCMC (THM-254) samples this energy.
- **Discrete Morse theory:** `H` a Morse function on `Q_m`; critical cells = extremal/metastable
  tournaments; the forbidden-H = gaps in the Morse spectrum; δ = the Morse matching/gradient.
- **Partition-function gradients:** the LRC covering-depth `Z` (HYP-2245) has the same gradient field;
  loneliness = its ground state, δ = its arc-derivative.

## The unification (extends HYP-2311)
`Ω(T)` and `H=I(Ω,2)` give a *scalar field* on the spin cube; `δ=dH` is its *gradient* (a frustrated
antiferromagnetic field); the Hessian = the coupling = S636's propagation operator; the **landscape
is glassy for even `n` (≥6)** — the 2-adic seam as a frustration/glass transition; the character-ratio
(Walsh) spectrum diagonalizes it and is bandlimited.

## To do
1. Pin the glass transition: count metastable basins vs `n` (even vs odd); is the metastable value
   (37 at n=6) and its forbidden flanks (35,39) a general pattern? Does ruggedness grow toward n=14?
2. Walsh/influence spectrum of δ across n; relate the bandlimit `2⌊(n-1)/2⌋` to the field's structure.
3. Discrete-Morse the H-landscape: Morse complex, the metastable basins as the LRC "worry" sectors.
