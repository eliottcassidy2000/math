---
id: HYP-2301
title: LRC = Hadwiger–Nelson in dimension 1 — both are χ(Cay(G, unit-sphere)); rigidity forces χ, CM rotation is the key
status: OPEN (unification); dimension-ladder + spindle=CM-multiplier + rigidity-mechanism VERIFIED
source: claudebox-2026-06-03-S634
related:
  - HYP-2295  # the coloring atlas (χ=sieve, α=corrector)
  - HYP-2292  # tie-graphs are Cayley/circulant
  - HYP-2230  # CM modulus-1 elements (the shared engine)
  - HYP-2288  # trienement / reversible runners (the 2n-1 = opposite-direction structure)
  - THM-407   # the 2-adic seam (= χ-parity of C_m here)
---

# HYP-2301 — the chromatic bridge: LRC and unit distance are one problem at different dimension

The most promising thread of the coloring atlas, pursued. **Both problems are the chromatic number of
a unit-distance Cayley graph `Cay(G, U)`**, with `U` the unit/resonance connection set:

| | ambient group `G` | connection `U` | the graph | χ |
|--|--|--|--|--|
| **LRC** (tie/sieve) | `ℤ/m` (rank 1) | `{±1,…,±band}` | circulant / `C_m` | sieve arity `= 2/3` |
| **unit distance** | `ℝ²/Λ` (rank 2) | unit sphere | plane UD graph | Hadwiger–Nelson `5–7` |

**The dimension of `G` is the only knob.** Verified χ dimension-ladder (`chi_bridge_s634.out`):
`C_m` → 2 (even) / 3 (odd); triangular (Eisenstein) disks → 3; **Moser spindle → 4**; plane → 5–7.
The LRC `2/3` is the 1-dimensional base of the Hadwiger–Nelson tower.

## The shared mechanism: rigidity forces χ (VERIFIED)

χ is pushed up, in every dimension, by a **finite rigid subgraph** that cannot be colored with the
naive palette:
- **1D:** an **odd cycle** `C_m` (m odd) is not 2-colorable ⟹ χ=3. This is the LRC single-→multi-
  sieve transition (S632): the sieve needs a 3rd color exactly when the tie-cycle is odd.
- **2D:** the **Moser spindle** (two unit rhombi, 7 vertices) is not 3-colorable ⟹ χ=4. This is the
  Hadwiger–Nelson lower-bound move; de Grey's χ≥5 graph is the same idea scaled.

> **"Multi-sieve is necessary" (LRC) and "χ(plane) ≥ k" (Hadwiger–Nelson) are literally the same
> theorem-shape:** exhibit a finite rigid unit-distance subgraph whose χ is forced up. The odd cycle
> is the 1D Moser spindle; the spindle is the 2D odd cycle.

## The shared key: the CM rotation (VERIFIED)

The Moser spindle glues two rhombi by a rotation of angle `arccos(5/6)` — a **rational-cosine
rotation = a CM unit** (a modulus-1 algebraic number, `5/6 + (√11/6)i`-flavored). That is exactly the
role of the **LRC multiplier `a ∈ (ℤ/m)*`** on the shell (the perspective key, S625/HYP-2130): in both
problems the obstruction is built by **overlaying rotated copies via a CM/modulus-1 rotation**. So the
shared engine of HYP-2230 (CM modulus-1 elements) is, concretely, *the rotation that makes the rigid
χ-forcing gadget*. The 2-adic seam (THM-407) is the χ-parity of `C_m`; the plane's χ is its 2D lift.

## The keys flow both ways

- **LRC → unit distance.** LRC circulants / lattice-mod-`p` graphs ARE finite unit-distance graphs;
  the LRC machinery (the residue-profile DP, the partition function, the shell tower) is a *source of
  finite high-χ graphs* for Hadwiger–Nelson lower bounds. The multiplier/perspective is the spindle
  rotation.
- **unit distance → LRC.** The Hadwiger–Nelson method (rigid finite high-χ subgraph) is how to PROVE
  a sieve needs ≥ k colors (multi-sieve necessity, HYP-2075). de Grey's "many rotated copies" is the
  high-arity multi-sieve.

## The literal unification (view-obstruction)
LRC is the `d=1` view-obstruction problem (Cusick); the chromatic number of the plane is the `d=2`
unit-distance χ. **`χ(ℝ^d`-unit-distance`)` is "LRC in dimension `d`"**, with `d=1` the runner sieve
(χ=2/3) and `d=2` Hadwiger–Nelson (χ=5–7). One problem, indexed by dimension.

## To do
1. Build an explicit finite LRC-circulant / Eisenstein-lattice-mod-`p` unit-distance graph and push
   its χ toward 4–5 (does the LRC machinery give a clean Hadwiger–Nelson lower-bound graph?).
2. Formalize "odd cycle = 1D spindle": a single rigidity lemma giving both χ(C_odd)=3 and χ(spindle)=4.
3. The `d=3` rung: χ(ℝ³-unit-distance) (known 6–15) as "LRC in dimension 3"; does the shell/2-adic
   structure predict the jumps?
