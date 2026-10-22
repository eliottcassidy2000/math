---
id: HYP-9364
title: "Self-blow-up principle: for a regular tournament T on n vertices, H(T) 2^(n-1)/n! = rho(T)(1 + O(1/n)), with rho(T) = det(nI - S)/|det(inI - S)| the blow-up gain of THM-4755. Consequences: (a) the ternary towers have H(S_k) 2^(n-1)/n! -> kappa_S = prod_i ((3^(2i+1)+1)/(3^(2i+1)-1))^(3^i) = 2.791509888...; (b) the rotational tournaments have H(R_n) 2^(n-1)/n! -> cosh(1)/cos(1) = 2.855957892... (Wormald's Conjecture 5); (c) Paley tournaments have H(P_q) 2^(q-1)/q! -> e (HYP-9274 (A)); (d) H(R_q)/H(P_q) -> cosh(1)/(e cos(1)) = 1.050648... (sharpens HYP-9362 (c)); (e) v_3(H(S_k)) = (3^k - 1)/2 for every k != 2."
status: >
  OPEN, pre-audit (mac-mini-2026-10-22-sierp). The principle is not claimed for the ordering of classes at fixed n.
  Evidence (FINITE-EXACT data plus polynomial extrapolation in 1/n):
  - (a) exact H(S_k) for k <= 6 (729 vertices). Richardson on k = 3..6 gives 2.79150957576, against kappa_S = 2.79150988820
    (difference 3.1e-7). The deviation 1 - kappa(S_k)/rho(S_k) is 0.0159, 0.0055, 0.0019 at n = 81, 243, 729, i.e. about 1.35/n.
  - The principle across every circulant class: kappa/rho lies in [0.954, 0.982] at n = 25 (208 classes) and in
    [0.925, 0.964] at n = 13, so both edges tend to 1 roughly like 1/n. But the rank correlation between kappa and rho
    across classes is near 0 for n >= 15: rho fixes the first-order size, not the order (spectral_H_predictor.out).
  - (b) exact H(R_n) for odd n <= 31. Polynomial fits in 1/n give 2.85360 (3 points, n = 27..31) and 2.85715 (4 points,
    n = 25..31), bracketing cosh(1)/cos(1) = 2.85596. This contradicts the rotational half of HYP-9028, which predicts e.
  - (c) exact H(P_q), q = 19, 23, 31: the 3-point fit gives 2.718332, against e = 2.718282.
  - (d) the 3-point fit (q = 19, 23, 31) of H(R_q)/H(P_q) gives 1.04953, against 1.050648.
  - (e) holds for k = 1, 3, 4, 5, 6; at k = 2 the valuation is one higher (3159 = 3^5 * 13).
  Wormald (arXiv:2610.16245) proves the principle for his random blow-ups. Regular blocks multiply the gain by e^(1/d),
  which the lexicographic law reproduces as the within-block spectral factor. HYP-9364 is the deterministic, self-similar
  analogue.
session: mac-mini-2026-10-22-sierp
source: 01-canon/theorems/THM-4755-blow-up-gains-skew-spectrum-lexicographic-law-ternary-tower.md
scripts:
  - 04-computation/experiments/sierpinski_ocf_20261022/spectral_H_predictor.py (+ .out)
  - 04-computation/experiments/sierpinski_ocf_20261022/ternary_fast.py (+ ternary_fast_k6.out)
related:
  - 05-knowledge/hypotheses/HYP-9362-one-flip-conjecture-paley-peak-and-the-7-21-to-8-24-dictionary.md
  - 05-knowledge/hypotheses/HYP-9357-paley-hamiltonian-path-maximality-and-the-kantor-galois-exceptional-primes-3-7-11.md
---

# HYP-9364 — the self-blow-up principle

**The idea.** Wormald computes the expected Hamiltonian count of a random blow-up of a fixed regular tournament `T*`. Its gain over a random tournament is `ρ(T*) = det(dI − S)/|det(idI − S)|` (THM-4755). The principle says that a symmetric or self-similar regular tournament behaves like its own blow-up: its normalised count `κ(T) = H(T)·2^(n−1)/n!` is `ρ(T)` up to `O(1/n)`.

**Why it is plausible.**
* For the ternary tower the two sides are computed independently (exact counts against the level product) and agree to 7 digits.
* For blow-ups with regular blocks it is Wormald's Theorem 1.

**What it predicts.** The flat spectrum (Paley) gives `e`, and the most concentrated spectrum (rotational) gives `cosh(1)/cos(1)`.
* This settles in principle the conflict recorded in HYP-9362 (c). The ratio `H(R_q)/H(P_q)` tends to `cosh(1)/(e cos 1) ≈ 1.0506`, not to 1.
* The rotational half of HYP-9028 (limit `e`) should then be withdrawn, and HYP-9274 (A) (Paley limit `e`) kept.

**What it does not predict.**
* The ordering of circulant classes at a fixed `n`. That is set at order `1/n` by the disjoint-packing terms of the OCF (THM-4751 (7)).
* Which class maximises `H` at small `n`. The flat class wins for `n ≤ 11`, while `ρ` is always largest at `R_n`.

**Tests.**
* `H(R_33)`, `H(R_35)` by inclusion–exclusion (feasible) to tighten (b).
* `H(P_43)` is out of reach. `H(P_47)` could be attacked through the Paley automorphism group.
* The tower at `k = 7` needs a memory-lean DP.
* A proof of (a) via the exact recursion: the renormalisation `α ↦ (2 + α)/3` and the level product are its fixed-point form.

## Addendum (same session): the second-order law

The residual `y = n(1 − κ/ρ)` is almost a linear function of the normalised fourth spectral moment `M4 = Σ_k σ_k⁴/n⁴` (`principle_residual.out`).
* `corr(y, M4)` = 0.982, 0.983, 0.958, 0.987, 0.992 for `n = 17, 19, 21, 23, 25` (all circulant classes).
* The per-`n` fits `y = a + b·M4` have `a ≈ 0.24–0.35` and `b ≈ 4.3–5.3`.

**(f) Conjecture: second-order law.** For regular `T` on `n` vertices,

    H(T) 2^(n−1)/n!  =  ρ(T) · (1 − (a + b·M4(T))/n + o(1/n)),

with universal constants `a`, `b` (fitted: `a ≈ 0.35`, `b ≈ 5.2` near `n = 25`).
* On circulants it brings the rank correlation with `κ` from about 0 (with `ρ` alone) to 0.72–0.88.
* Its relative error is at most `2.6·10^(−3)` at `n = 23, 25`, and about `10^(−2)` at the composite `n = 15, 21`.
* It finds the maximiser every time.

For circulants, `M4` is an affine function of the 4-cycle count `c_4` (equivalently of the additive energy). So the second-order term is the energy correction of THM-4751 (6), now with a sign and a size. A fourth-cumulant (Edgeworth) correction to the local limit theorem is the natural source of an `M4` term. No derivation of `a`, `b` is claimed.

**Wormald's Conjecture 1** (regular bases, `d ≤ 13`) is verified exhaustively: see THM-4755, addendum.
