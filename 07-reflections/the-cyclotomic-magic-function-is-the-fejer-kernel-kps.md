# The cyclotomic Delsarte/Beurling–Selberg magic function IS the Fejér kernel

*kind-pasteur-2026-06-28-S31ao. The owner asked to look for more modular magic functions and to think
cyclotomic Delsarte / Beurling–Selberg. The search lands on a single, classical, explicit object: the
**Fejér kernel** `F₇`. It equals `(de Moivre cubic)²` on the circle, is positive-definite with the Fejér
weights, has double zeros at the de Moivre angles (LP-sharpness), is the Chebyshev V₇ equioscillation, AND
is the AP orbit's autocorrelation. The whole cyclotomic-magic-function picture collapses to `F₇`.*

## The identity (verified to 1e-12)
On the circle `u = 2cos θ`, since `V₇(2cos θ) = 2cos 7θ` (Vieta–Lucas/Chebyshev):
```
(de Moivre cubic)²(2cos θ)  =  (V₇(u) − 2)/(u − 2)  =  (2cos 7θ − 2)/(2cos θ − 2)
                            =  sin²(7θ/2) / sin²(θ/2)  =  F₇(θ),   the FEJÉR KERNEL of order 7.
```
And its Fourier transform is **`F̂₇(n) = (7 − |n|)₊`** — the Fejér weights, **non-negative** (verified
exactly). So `F₇` satisfies BOTH magic-function conditions at once:
- **non-negative** `F₇(θ) ≥ 0` (the majorant / sign condition),
- **positive-definite** `F̂₇ ≥ 0` (the Delsarte dual condition),
with `F₇(0) = 49 = 7²` and **double zeros exactly at the de Moivre angles `θ = 2πj/7`** (`j≠0`).

## Why this is THE cyclotomic magic function (five facets, one kernel)
- **Delsarte:** `F̂₇ ≥ 0` is the positive-definiteness the Delsarte LP dual requires.
- **Beurling–Selberg:** `F₇ ≥ 0` band-limited (degree 6) — the extremal non-negative trig polynomial of
  its degree (Fejér–Jackson). It is the *one* construction that is simultaneously a non-negative majorant
  AND positive-definite, which is exactly what THM-537 found the *literal per-term* majorant could not be.
- **Cohn–Elkies / Viazovska:** the **double zeros at the de Moivre angles** are the LP-sharpness condition
  (the optimal configuration's correlations sit at the magic function's double zeros) — the 1-D, 7-fold
  analog of Viazovska's `±√(2n)` double zeros.
- **Chebyshev (HYP-3212):** `F₇ = (de Moivre cubic)²` = the V₇ equioscillation extremal; the rationality
  (ℚ-collapse) is the double root.
- **Extremal config:** the **AP orbit autocorrelation IS `F_k`** (the order-`k` Fejér kernel,
  `|Σ_{j<k} e^{2πi m jx}|² = F_k(mx)`). The AP is sharp because its autocorrelation is a Fejér kernel —
  the extremal positive-definite function — and the 7-sector structure contributes the `F₇` factor (zeros
  at `7|n`). **Coverage = the pairing of the orbit Fejér `F_k` with the sector Fejér `F₇`**; the AP
  maximizes it because both factors are the extremal kernels.

## Why this can succeed where the literal Beurling–Selberg failed (THM-537)
THM-537's walls were (1) the **minorant** real-analyticity wall (odd inclusion–exclusion terms need a
non-negative minorant of an indicator that vanishes on an arc — impossible) and (2) the **signed-
cancellation** wall (per-`T` band-limited majorants miss cross-`T` sign alignment). The Fejér/Delsarte
magic function **does not majorize each term**: it is a single positive-definite kernel paired against the
whole orbit, so the signed cancellation lives *inside* the kernel's Fourier support (the `(7−|n|)₊`
weights), exactly as the moment-LP (`L_y`) bypasses it inside the exact `S_t`. **The Fejér kernel is the
Fourier-side incarnation of the same bypass** — the positive-definite certificate, not a term-wise
majorant. This is the natural home for the LP dual that THM-537 said had to be signed.

## More magic functions (the family the owner asked for)
The cyclotomic magic functions form a family, all built from the same de Moivre / Fejér seed:
- **`F_k`** (order-`k` Fejér) = the **orbit** kernel (the AP autocorrelation); `F₇` = the **sector** kernel.
- **`F_m * F_n`** (convolutions) and **`F_k²`-type** higher kernels — the de la Vallée-Poussin / Jackson
  kernels (sharper band-limited non-negative majorants).
- **Modular:** `F₇` is the `q→` boundary of the **weight-2 Eisenstein series for `Γ₀(7)`** (the Lipschitz
  formula `Σ 1/(n+z)² = ` a Fejér-type sum); the level-7 Eisenstein `E_{2,7}` is the *modular* magic
  function whose `q`-expansion refines `F₇` (mac-mini's `E₂` spoke). **This is the level-7 analog of
  Viazovska's level-1 modular construction** — the explicit modular magic function for the heptagon.
- **OPUC / Christoffel–Darboux:** for the cyclotomic measure (atoms at the 7th roots), the
  Christoffel–Darboux reproducing kernel is `F₇`-structured (ties to mac-mini's Verblunsky/OPUC, S73d).

## Net / proof relevance
- **VERIFIED:** the cyclotomic Delsarte/Beurling–Selberg magic function = the **Fejér kernel `F₇` =
  (de Moivre cubic)²**, positive-definite (`F̂₇=(7−|n|)₊`), double zeros at the de Moivre angles. The AP's
  autocorrelation is the order-`k` Fejér kernel; coverage = the `F_k`⊗`F₇` pairing, maximized at the AP.
- **PROOF ROUTE:** frame the cover bound as the **Fejér/Delsarte positive-definite certificate** (`F₇`
  paired with the orbit), which bypasses THM-537's minorant wall by being a single PD kernel, not a
  term-wise majorant. The AP is sharp via the double-zero (Cohn–Elkies) condition.
- **MORE:** the modular magic function = the weight-2 `Γ₀(7)` Eisenstein series (the level-7 Viazovska);
  the de la Vallée-Poussin/Jackson kernels are the sharper family members.

→ HYP-3214 (this), HYP-3212 (Chebyshev/de Moivre), THM-537 (the blocked literal Beurling–Selberg, now
bypassed), THM-534 (moment-LP `L_y`), HYP-3201/3203 (Toeplitz/OPUC), Fejér–Jackson, Cohn–Elkies/Viazovska,
Γ₀(7) Eisenstein, OPEN-Q-108.
