---
id: HYP-9364
title: "Self-blow-up principle: for a regular tournament T on n vertices, H(T) 2^(n-1)/n! = rho(T)(1 + O(1/n)), with rho(T) = det(nI - S)/|det(inI - S)| the blow-up gain of THM-4755. Consequences: (a) the ternary towers have H(S_k) 2^(n-1)/n! -> kappa_S = prod_i ((3^(2i+1)+1)/(3^(2i+1)-1))^(3^i) = 2.791509888...; (b) the rotational tournaments have H(R_n) 2^(n-1)/n! -> cosh(1)/cos(1) = 2.855957892... (Wormald's conjecture following his Conjecture 5, path form; the value is Bhattacharjee's); (c) Paley tournaments have H(P_q) 2^(q-1)/q! -> e (HYP-9274 (A)); (d) H(R_q)/H(P_q) -> cosh(1)/(e cos(1)) = 1.050648... (sharpens HYP-9362 (c)); (e) v_3(H(S_k)) = (3^k - 1)/2 for every k != 2."
status: >
  OPEN (mac-mini-2026-10-22-sierp; data and fits audited 2026-10-22, sierpinski_audit_20261022.md). The principle is not claimed for the ordering of classes at fixed n.
  UPDATE (same session): the principle is the one-loop term of the loop expansion THM-4757 (exact identity; one-loop
  sum = log rho PROVED); the second-order law (f) is superseded by the parameter-free two-loop law HYP-9365.
  Evidence (FINITE-EXACT data plus polynomial extrapolation in 1/n):
  - (a) exact H(S_k) for k <= 6 (729 vertices). Richardson on k = 3..6 gives 2.79150957576, against kappa_S = 2.79150988820
    (difference 3.1e-7). The deviation 1 - kappa(S_k)/rho(S_k) is 0.0159, 0.0055, 0.0018 at n = 81, 243, 729, i.e. about 1.35/n.
  - The principle across every circulant class: kappa/rho lies in [0.954, 0.982] at n = 25 (208 multiplier classes) and in
    [0.925, 0.964] at n = 13, so both edges tend to 1 roughly like 1/n. But the rank correlation between kappa and rho
    across classes is near 0 for n >= 15: rho fixes the first-order size, not the order (spectral_H_predictor.out).
  - (b) exact H(R_n) for odd n <= 33. Polynomial fits in 1/n give 2.85360 (3 points, n = 27..31) and 2.85715 (4 points,
    n = 25..31), but mixed fits are not stable: 5 and 6 points give 2.8707 and 2.9613, and with H(R_33) the 3-, 4- and 5-point
    fits give 2.8539, 2.8550, 2.8484 (audit). The instability is a small alternating term in n mod 4 (HYP-9365 (e)); fitted
    within each residue class mod 4 the fits are stable: 2.855930 (n = 1 mod 4, 7 points) and 2.85611-2.85631 (n = 3 mod 4),
    against cosh(1)/cos(1) = 2.855958. The data exclude e: n^2 (kappa(R_(n+2)) - kappa(R_n))/2 is 5.56 at n = 31 and still
    rising, while reaching e needs only 0.050 more. This contradicts HYP-9028.
  - (c) exact H(P_q), q = 19, 23, 31: the 3-point fit gives 2.718332, against e = 2.718282. With the Paley tournament on F_27
    (P27.out) the fits give 2.72566 (q = 19..31) and 2.72050 (q = 23, 27, 31): the 5e-5 agreement depends on the points chosen.
    The two-loop law (HYP-9365 (a)) gives the rate: kappa(P_q) = e(1 - 3/(2q) + O(q^-2)).
  - (d) the 3-point fit (q = 19, 23, 31) of H(R_q)/H(P_q) gives 1.04953, against 1.050648 (with q = 27: 1.04755 and 1.04894).
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
* For the ternary tower the two sides are computed independently (exact counts against the level product) and agree to about 6 digits (Richardson through 3, 4, 5 levels: −5.2e-6, −3.1e-7, +2.3e-7).
* For blow-ups with regular blocks it is Wormald's Theorem 1.

**What it predicts.** The flat spectrum (Paley) gives `e`, and the most concentrated spectrum (rotational) gives `cosh(1)/cos(1)`.
* If true, it resolves the conflict recorded in HYP-9362 (c). The ratio `H(R_q)/H(P_q)` tends to `cosh(1)/(e cos 1) ≈ 1.0506`, not to 1.
* The rotational half of HYP-9028 (limit `e`) should then be withdrawn, and HYP-9274 (A) (Paley limit `e`) kept. HYP-9028's universality clause ("every sequence of regular tournaments has exc → e") is already refuted, independently of this hypothesis: by Wormald's Theorems 1, 2 and 5, the regular blow-ups `R(C_3, n)` have `E[exc] → 2e^(1/3) = 2.7912`, so regular tournaments with `exc ≥ 2.79 − o(1)` exist for infinitely many `n`; also `exc(S_6) = 2.7826 > e`. This hypothesis predicts HYP-9028's control family `C_3[R_k, R_k, R_k]` → `2(cosh(1/3)/cos(1/3))^3 ≈ 2.7917`.

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

**(f) on all regular tournaments** (`principle_residual_regular.out`).
* At `n = 11` (1223 classes), `corr(y, M4) = 0.956` and the fit is `y = 0.118 + 5.82·M4`. The rank correlation with `κ` moves from −0.743 (`ρ` alone) to +0.849, and the corrected predictor picks the true maximiser `P_11`, which `ρ` alone misses.
* At `n = 9` (15 classes): from −0.721 to +0.796.

**Why Paley wins only at 3, 7, 11 (a quantitative reading).**
* The flat spectrum has `M4 ≈ 1/(2n)`.
* The rotational spectrum has `M4 → Σ_(j odd) (2/(πj))^4 = 1/6`.
* `R_n` gains `log(ρ(R_n)/ρ(flat))`, which tends to `log(cosh 1/(e cos 1)) ≈ 0.049`. It pays about `b·(1/6)/n`.
* With `b ≈ 4.3–5.8`, the two balance near `n ≈ 13`. That is where `R_13` first wins, by 0.22%.

This is the spectral form of the packing mechanism of THM-4751 (7). It is an empirical law (FINITE-EXACT fits), not a derivation.

**A proof route (exact identity).** Off the diagonal, `A_(uv) = (1 + S_(uv))/2`. So

    κ(T) = H(T) 2^(n−1)/n! = E_σ [ Π_(i=1)^(n−1) (1 + S_(σ_i σ_(i+1))) ],

the expectation over a uniformly random ordering `σ`.
* Expanding the product gives a sum over sets of steps, i.e. disjoint runs, each run a path weighted by its product of `S`.
* For a regular `T`, `S·1 = 0`, so isolated runs vanish to leading order. What survives comes from the distinctness constraints between runs: a cluster expansion in closed skew walks, the traces `tr S^(2m)`.
* The conjecture is that this cluster sum resums to `log ρ = −Σ tr(S^(4m+2))/((2m+1) n^(4m+2))` at leading order, with an `M4 = tr(S⁴)/(2n⁴)` correction at order `1/n`. That would match the second-order law.

## Addendum 2 (same session): the principle is the one-loop term

THM-4757 turns the proof route above into theorems.
* **Exact identity.** `κ(T) = Σ_L N_L(S)/(n)_L` over directed linear forests.
* **One loop.** The excess-0 clusters of `log Φ_T` sum exactly to `log ρ(T)`.
* **Two loops.** The excess-1 clusters (figure-eight, theta, dumbbell) have a closed matrix formula.

So the principle `κ ≈ ρ` is the one-loop approximation of a loop expansion with parameter `1/n`. Statement (f) is superseded by HYP-9365, the two-loop law with no fitted constants:
* corr 0.9984 with the exact residual over all circulant classes at `n = 25`;
* corr 0.978 over all regular 11-tournaments.

**(b) sharpened.** Exact `H(R_33) = 5394851076280218488551446699` (`R33.out`, 13.5 CPU-min).
* Fitted by residue mod 4, `n ≡ 1` (7 points) gives `κ(R_n) → 2.855930` and `n ≡ 3` gives 2.85611–2.85631, against `cosh 1/cos 1 = 2.855958` (`rot_limit_fit.out`).
* The mixed fit is unstable because of a small alternating term of non-perturbative origin (HYP-9365 (e)).
* The two-loop law gives the approach `κ(R_n)/ρ(R_n) = 1 − y/n + …` with `y → tan 1 tanh 1 + ½ tan² 1 − 1 = 1.398872`.

**(c) sharpened.** HYP-9365 (a) gives `κ(P_q) = e(1 − 3/(2q) + O(q^(−2)))` (`paley_y.out`).

**All regular `n = 13`** (`principle_regular_13.out`, 1495297 classes): `κ/ρ ∈ [0.9245, 0.9687]`, mean 0.9564, and the rank correlation of `κ` with `ρ` is −0.785.
