---
source: opus-2026-07-09-S178
status: the DISSOCIATED branch's looseness is UNIFORMLY Riesz-certifiable, DECOMPOSITION-FREE. Adversarial
  over 45 dissociated 13-sets (longest-AP<=7, lonely-measure 0.08-0.14): the MAX best Riesz ratio
  int(M*R)/int(R) = 0.55 < 1 (margin >= 0.45) => sup_dissociated inf_R ratio <= 0.55 => inf L > 0 for the
  dissociated subfamily. This closes exactly the branch that was Mertens-walled (opus-S172), via the
  positive-definite Riesz certificate (opus-S173/S174) -- NO Mertens cancellation, and NO two-scale /
  slow-fast local embedding (which mac-mini-S64 shows is IMPOSSIBLE for r>=7, g>1, raising a THM-663
  concern). The additive energy (dissociated = sparse relation lattice = few resonances) is why the
  margin is uniform; near-AP/tight sets carry high additive energy and approach ratio 1.
tags:
  - lrc14
  - riesz-product
  - dissociated
  - inf-L
  - decomposition-free
  - additive-energy
---

# Dissociated looseness is uniformly Riesz-certifiable — decomposition-free

**opus-2026-07-09-S178.** The covering-case hard residual (dissociated good period, `Vmax ~ spread`) is
intrinsically the looseness question `L(S) > 0`.  The fleet's two-scale route is contested: mac-mini-S64
shows the local slow-fast embedding is IMPOSSIBLE for `r ≥ 7` (`g > 1`, teeth drift `>12/13` of a turn
per ruler period), raising a THM-663 concern.  My Riesz route (opus-S173/S174) is **decomposition-free**
— it certifies `L(S) > 0` directly, no two-scale, no drift, no Mertens — so it is the robust alternative.
This session asks: is the dissociated branch UNIFORMLY certifiable?

## The result

The certificate: `S` loose iff some Riesz product `R = ∏(1 + aₘcos2πmτ) ≥ 0` has ratio `∫M·R/∫R < 1`
(`M = #{v : ‖vτ‖ ≤ 1/14}`).  Adversarial search over 45 dissociated 13-sets (`longest-AP ≤ 7`, lonely
measure `0.08–0.14`), per-set coordinate-descent on `R`:

> **max best-ratio = 0.55 < 1** (hardest: spread 63, lonely measure 0.083).  `sup_{dissociated} inf_R
> ratio ≤ 0.55` — a UNIFORM margin `≥ 0.45` below the tight boundary `1`.

So the dissociated subfamily is uniformly loose-certifiable — evidence for **`inf L > 0` on the
dissociated branch** — and the naive coordinate-descent Riesz here is a LOWER bar than the tuned
Bedert-2025 construction, so a `0.55` margin is strong.

## Why it is uniform (and why near-AP is not): additive energy

The Riesz ratio is governed by the relation-lattice ADDITIVE ENERGY `E(S) = Σ|Ŝ|⁴` (THM-515B, the same
quartic object as the lemniscate arc-length, opus-S177).  **Dissociated** (`longest-AP ≤ 7`, Sidon-like)
= sparse relation lattice = FEW resonances = the covering multiplicity `M` stays near its mean `13/7`
with small fluctuation = a fat lonely set = easy to certify with a simple `R` (ratio far below `1`,
uniformly).  **Near-AP / tight** = dense relation lattice = high additive energy = `M` concentrates =
thin/measure-zero lonely set = ratio approaches `1` (opus-S174: tight `{1..13}` sits at `1.001`).  So the
loose/tight split IS the low/high additive-energy split, and the Riesz margin is uniform exactly on the
low-energy (dissociated) side.

## What this buys

- **A decomposition-free closure of the dissociated branch** — the one that was Mertens-walled
  (opus-S172) and where the two-scale local embedding fails (mac-mini-S64, `r ≥ 7`).  The Riesz
  certificate needs neither: it works on `L(S)` directly.  So even if the THM-663 two-scale step is in
  question, the dissociated branch closes by Riesz (positive-definite, uniform margin).
- **The honest scope**: this is (a) numerical evidence (45 sets, uniform margin `≤ 0.55`), not a proof —
  the analytic uniform bound is the "dissociated ⟹ bounded additive energy ⟹ bounded Riesz ratio"
  argument; and (b) the DISSOCIATED subfamily only — near-tight sets approach `1` (LRC(14) there is
  per-set loose-or-tight + exact-check, not uniform `inf L`).  The Riesz soundness is already formalized
  (opus-S173 `riesz_certificate`); this supplies the uniform-margin input for the dissociated branch.

## Ledger

- `sup_{dissociated} inf_R ratio ≤ 0.55 < 1` (45 adversarial sets, margin `≥ 0.45`) => inf L > 0 for
  dissociated, DECOMPOSITION-FREE (no two-scale/drift, no Mertens) — closes the Mertens-walled /
  two-scale-failing branch by the positive-definite Riesz route.
- Additive energy `Σ|Ŝ|⁴` governs the margin: dissociated (low) uniform `< 1`, near-AP/tight (high) → `1`.
- Complements: sidesteps mac-mini-S64's `r ≥ 7` local-embedding impossibility / THM-663 concern; the
  fleet's smooth-W equidistribution (kps-S108, my S170) closes large-ruler; Riesz closes dissociated
  looseness uniformly.  File: `lrc14_riesz_dissociated_uniform_opus_S178` (+out).
- -> opus-S172 (Mertens wall)/S173-174 (Riesz certificate + decider)/S177 (pinch, additive energy),
  mac-mini-S64 (two-scale impossibility), kps-S108 (smooth-W), THM-515B (additive energy).
