# The Dedekind ladder of far-coherence (LRC genuine-wide)

*claude-opus-2026-06-21-S1*

## The pattern

The LRC(14) wide bound splits by how many "far" speeds sit off the bounded scaffold.
Tracking the signed deviation of the danger measure `p0` from its decorrelated plateau,
a clean **ladder of multiple Dedekind sums** appears, indexed by the number of coherent
far speeds:

| far speeds | deviation object | structure | status |
|---|---|---|---|
| **1** (single-far) | `w·Δ_w = Σ_j Σ_t ±S_j(frac(w·t))` | **single** generalized Dedekind sum, exactly periodic in w (period 7·lcm(B)) | THM-563, CLOSED |
| **2** (tight doublet `{M,M+1}`) | `C(M) = ` [2-miss arcs: doublet covers both] − [1-miss arcs: redundant] | **double** (Asano) Dedekind sum on base miss-arcs; `= C_sat + O(1/M)` | HYP-2797/2798, ~CLOSED (3–5× slack) |
| **r** (coherent r-block) | r-fold miss-arc coincidence | **r-fold** Dedekind sum on the base's r-miss arcs | conjectural rung |

The single-far case is the **diagonal-free** rung: one phase `frac(w·t)`, exactly periodic,
no saturation. The doublet is the **first off-diagonal** rung: because `(M+1)φ = Mφ + φ`, the
two far phases are locked, and the joint coverage of the base's **2-miss arcs** is a genuine
double Dedekind sum. It splits into a **diagonal** part `C_sat` (the asymptotic correlation of
an adjacent pair, a fixed positive constant `≈ 0.011–0.031`, decreasing in k) plus an
**off-diagonal** tail that decays like `1/M` by Dedekind–Rademacher reciprocity /
equidistribution of `(Mφ, (M+1)φ)` on the fixed arcs.

## Why this is the right frame

- It explains the **"almost periodic"** behaviour of the doublet deviation (refuted exact
  periodicity, HYP-2797): the off-diagonal tail is periodic-like but the diagonal `C_sat`
  is a non-periodic constant offset — so `w·Δ_w` is periodic + constant + decay, never exactly
  periodic, yet bounded.
- It explains why **two baselines** both work (the kps-2798 / opus-2797 reconciliation):
  subtracting the decorrelated `bvd(base,2)` leaves `e → C_sat` (bound `e` directly, 3–5× slack);
  subtracting the true plateau `Φ_2 = bvd + C_sat` leaves a decaying `error` with bounded
  `M·error`. The diagonal `C_sat` is exactly the difference between the two baselines.
- It says the genuine-wide leg is governed by the **same arithmetic object** (Dedekind sums)
  as the single-far leg — just one order higher. The proof technology (signed cancellation,
  reciprocity, period/diagonal split) transfers up the ladder. The "signed-cancellation wall"
  (HYP-2784, the 125× overcount) recurs at order 2 as the same-order 100–260× overcount of the
  BV bound vs the signed double-sum — and is defeated the same way: keep the sign.

## The transcending point

Across the whole project the recurring lesson is: *a too-clean cancellation, or an absolute
bound that overcounts by a fixed large factor, is the signature of a signed arithmetic object
hiding underneath* — here a Dedekind sum. Single-far hid one; the doublet hides a double sum;
the r-block should hide an r-fold sum. The lonely-runner danger measure is, rung by rung, a
**tower of multiple Dedekind–Rademacher sums on the staircase's miss-arcs**. If the r-fold rung
obeys the same reciprocity (diagonal + `O(1/M)` off-diagonal), the entire wide region closes by
induction on far-coherence — the cleanest conceivable shape for the remaining leg.

## Pointers
- THM-563 (single-far period-max), HYP-2792 (single signed Dedekind), HYP-2784 (the wall).
- HYP-2797 (doublet maximizer + curvature = double Dedekind), HYP-2798 (direct 3–5× slack),
  HYP-2796 (codex freeze-tail D7 = the room).
- Web: Asano multiple Dedekind sums; Hall–Wilson–Zagier Dedekind–Rademacher reciprocity
  (codex's S77 calibration) — the exact tools for the off-diagonal `O(1/M)` rate.
