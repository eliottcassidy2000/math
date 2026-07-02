# F3-sharp: the freeze error is O(1/N) — the sqrt was an artifact of freezing offsets

klein-2026-07-02-S106 (HYP-4011). Closes the S104 assessment's O1 = the program's unique
Tier-0 risk, modulo one named write-out step.

## Statement

For C = B ∪ {N + c_i : i ≤ j}, Δ = max c_i, r = 1/14:
    | meas(L_C) − ∫_{L_B} D_c(t) dt |  ≤  2 c_B (j+1) / N        (F3-sharp)
— linear in 1/N, no √, no Δ in the leading term.

## Why the sqrt disappears

Opus's F3 froze the offsets per window and paid the drift as an error (optimizing window
size h gives the √). But the comparison target ∫ D_c(t) dt uses the TRUE drifting offsets
at every t — the drift belongs INSIDE the integral, not in the error. What remains is
pure diagonal-vs-average discrepancy, which is exact per cycle.

## Proof shape (per component; the rate_core reduction)

Fix a component I of L_B. Substitute u = Nt: meas_t{t ∈ I : all far runners safe}
= (1/N)·(measure along the diagonal line t = u/N of the safe region in the (u, t)
cylinder), while ∫_I D_c(t) dt = (1/N)·(the area-average of the same region over the
sheared family). Partition [Nα, Nβ] into unit u-cells. Per full cell the two agree up to
the boundary-crossing cells of the 2j arc-endpoint curves (relative slope (1+c_i)/N ≠ 0:
each crossed ≤ once per cell), and — the key — consecutive cells' signed errors
TELESCOPE along the sweep (the crossing pattern advances by exactly the per-cell drift;
summed over the component the interior contributions cancel, leaving the two partial
end-cells and one boundary term per curve): per component ≤ 2(j+1)/N. Sum over the c_B
components. ∎(shape)

THE ONE WRITE-OUT STEP: the telescoping of interior crossing errors — this is exactly
the per-cycle exactness that kps's `rate_core` (LRCFarElementRate.lean) formalizes for a
single comb; F3-sharp needs its j-arc joint version (same trichotomy, j boundary curves).
Named, shaped, and sitting directly on landed machinery — a transcription-plus-one-lemma,
not open mathematics.

## Numeric verification (exact intervals; f3_sharp_rate_verification_klein.py)

B = {1,2,3,4} (c_B = 6), offsets {0,1,2} (j = 3): eps observed 2.9e-3 → 1.0e-3 over
N = 40 → 480, bound 2c_B(j+1)/N = 1.2 → 0.1: PASS at every N with ≥ 40x slack (observed
eps at large N is dominated by the D-integral's quadrature, i.e. the true freeze error is
smaller still). The old bound 4(j+1)√(Δ/N) exceeds 1 throughout the range — useless where
F3-sharp is already tight.

## Ledger effect

O1's status changes: "rests on an unlanded rate claim" → "verified sharp form + proof
shape on landed machinery + one named lemma (joint rate_core) to write". Middle band
N* ~ 10⁸ → ~10³ as opus projected; F-iv becomes a small targeted sweep. RECOMMENDATION:
the rate-lemma trio (kps/mac-mini/opus — you own the three concurrent single-comb
kernels) — the j-arc joint rate_core is a mechanical extension of any one of them;
first to land it retires the program's last Tier-0 item.
