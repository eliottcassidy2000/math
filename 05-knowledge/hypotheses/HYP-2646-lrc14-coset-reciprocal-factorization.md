---
id: HYP-2646
title: LRC(14) the support-6 kernel factorizes EXACTLY as a finite mod-7 character coefficient times a pure reciprocal — K(n)=D7(n mod 7)/prod(n_j), with D7(-c)=conj(D7(c)); this makes "the signed/coset quotient is the ruler" (HYP-2640) literally precise and explains why the inf-norm box truncation plateaus
status: CONFIRMED (exact factorization verified to 1e-19 over 2000 random support-6 n; antipodal conjugate-symmetry verified to 3e-16); the precise signed-quotient form; does NOT close the cap (corroborates the far-element route HYP-2644)
source: kind-pasteur-2026-06-19-S12
depends_on:
  - THM-538   # support-6 floor (K=0 below support 6)
  - HYP-2606  # the signed relation-lattice Fourier identity
  - HYP-2632  # the finite chi_7 / additive-Fourier kernel (this is the per-coordinate atom)
related:
  - HYP-2640  # signed/coset quotient is the ruler (THIS is the precise form)
  - HYP-2644  # far-element plateau recursion (the route this feeds)
  - HYP-2645  # Poisson dual / box truncation plateaus (this EXPLAINS the plateau)
  - HYP-2633  # residue-lift equidistribution (the reciprocal sum S_c is the lift)
  - MISTAKE-078
  - OPEN-Q-108
---

# HYP-2646 — The exact coset/reciprocal factorization of the support-6 kernel

## The factorization (CONFIRMED, exact)

For a support-6 relation `n` (six coordinates nonzero and ≢0 mod 7; THM-538 floor), the
seven-sector signed Fourier kernel factorizes **exactly** as

> **`K(n) = D7(n mod 7) / ∏_j n_j`**,

where the **finite mod-7 character coefficient** is
```
D7(c) = Σ_{T⊆{1..6}} (−1)^|T| ∏_j h_T(c_j),   c_j ∈ F_7^*,
h_T(r) = −A(r)·Σ_{j∈T} ζ_7^{−rj},   A(r) = (1−ζ_7^{−r})/(2πi),   ζ_7 = e^{−2πi/7}.
```
This is the per-coordinate split of THM-538's `ĉ_T(m) = ŝ_T(m)`: writing
`ĉ_T(m) = h_T(m mod 7)/m` (magnitude `|sin(πm/7)|/(π|m|)`, phase a function of `m mod 7`), the
product over the six coordinates pulls out `1/∏ n_j` cleanly. (For support s>6 the same holds with
the product over all s effective coordinates; support <6 gives `K=0` by THM-538.)

**Verification** (`04-computation/lrc14_coset_quotient_decay_kps_s12.py`):
`max|K_direct − K_factored| = 1.6·10⁻¹⁹` over 2000 random support-6 `n`.

## The consequences (this is the *precise* signed-quotient ruler asked for by the endgame brief)

1. **The correction is a finite-coefficient–weighted reciprocal lattice sum:**
   ```
   corr(E) = meas(S7(E)) − M7(k) = Σ_{c ∈ (F_7^*)^6} D7(c) · S_c(E),
   S_c(E) = Σ_{n∈Λ(E), n≡c mod 7} 1/∏_j n_j   (+ the support-7,8 analogues).
   ```
   The **finite** outer object `D7` carries the algebraic cancellation; the **analytic**
   convergence lives entirely in the reciprocal lattice sums `S_c`. This is exactly the
   two-layer "expose the channel, then scalarize" structure of HYP-2632/2633/2640.

2. **The correction is REAL because `D7(−c) = conj(D7(c))` (verified 3·10⁻¹⁶).** The lattice is
   `n ↔ −n` symmetric and `∏(−n_j)=∏ n_j` (six coordinates, even), so antipodal pairs give
   `K(n)+K(−n) = 2·Re D7(c)/∏ n_j`. **The ruler is `Re D7(c)`**, with `|Re D7| ≤ 0.1431`. The
   imaginary part of every coset cancels under the antipodal pairing — this is the precise
   mechanism behind HYP-2632's `−108U+54U` signed cancellation.

3. **`D7` is NONZERO on ALL 46656 cosets** `(F_7^*)^6` (max `|Re|=0.1431`, max `|Im|=0.627`). There
   is **no coset-level vanishing** — the cancellation is not a sparsity of the finite kernel; it is
   the alternating sign of `Re D7` across cosets *plus* the conditional convergence of `S_c`.

## Why the inf-norm box truncation PLATEAUS (explains HYP-2645 / MISTAKE-078 sharply)

Because `K(n)=D7(c)/∏ n_j` and `∏ 1/|n_j|` over a rank-(k−2) lattice is harmonically
divergent, the series `Σ K(n)` is only **conditionally convergent**: the value depends on the
summation ORDER. The symmetric inf-norm box `[−L,L]^{k−1}` is the WRONG order. Measured
(`lrc14_signed_shell_decay`, `lrc14_coset_quotient_decay`):

| E | target `p0−M7` | box `|n|∞≤7` recovers | box `|n|∞≤9` recovers |
|---|---:|---:|---:|
| AP `[0..7]` | 0.30273 | 0.0249 (8.2%) | 0.0198 (6.5%) — **non-monotone** |
| odd-AP `[0,1,3,…,13]` | 0.21348 | 0.0228 (10.7%) | 0.0176 (8.3%) |

The per-shell SIGNED sum does NOT decay (AP: −0.001,+0.007,−0.017,+0.015,+0.015,+0.006); the per-shell
ABSOLUTE sum GROWS (0.014→0.44→1.32→2.97→3.19→1.82). Supports 6,7,8 each contribute the same order
and alternate sign (AP L=3: s6=−0.0111, s7=+0.0066, s8=−0.0009). **Conclusion: the box truncation is a
mathematically illegitimate way to evaluate this conditionally convergent sum** — which is exactly why
the correct route (HYP-2644) does NOT sum the lattice at all but uses the finite x-cell evaluation
(HYP-2645) and the single-frequency far-element plateau.

## Tie-in to the far-element direction (the live route)

The reciprocal sum `S_c(E)` over relations *involving the far coordinate `w`* is, by the linear
constraint `n_w w = −Σ_{j≠w} n_j e_j`, a sum with `1/n_w = −w / (Σ_{j≠w} n_j e_j)` — the `1/∏ n_j`
weight is *not* small, but such relations are SPARSE (they force a large companion coordinate). The
exact decorrelation defect `Δ_w = p0(E) − Plat(E')` (HYP-2644) is therefore a **single-frequency
(frequency `w`) signed object**, and the factorization shows its sign is governed by `Re D7`. Measured
uniform rate (`lrc14_weyl_through_w_coset_kps_s12b`, `lrc14_weyl_uniform_C_kps_s12c`): over consec,
near-AP, GAP, odd-value, and random cores at k=8,9,
```
sup_w  w·|Δ_w|  ≤  C  with  C ≈ 1.95 (k=8), 1.45 (k=9),
```
consistent with HYP-2644's measured `~1.25`. Since every core has `Plat(E') ≤ 0.362 < cap − 0.13`, the
far direction closes once `w > C/margin ≈ 15`, and `w ≤ 15` is bounded-spread (finite check).

## Honest status

CONFIRMED as an exact identity and the precise form of the "signed/coset quotient." It does NOT by
itself bound `corr` (the reciprocal sums `S_c` are still conditionally convergent, the same
obstruction). Its value: (i) it turns the vague "signed quotient is the ruler" into a literal
finite-coefficient formula `corr = Σ_c Re D7(c)·(2 Re S_c⁺)`; (ii) it explains, mechanistically, why
the lattice-box route plateaus and why the far-element x-cell route (HYP-2644/2645) is the correct one;
(iii) it identifies `Re D7` (|·|≤0.1431, antipodal-real) as the exact algebraic weight for any future
summation-by-parts on `S_c`.

Files: `04-computation/lrc14_{signed_shell_decay,coset_quotient_decay,conditional_truncation,weyl_through_w_coset_kps_s12b,weyl_uniform_C_kps_s12c}.py`;
outputs in `05-knowledge/results/`.
→ THM-538, HYP-2606, HYP-2632, HYP-2640, HYP-2644, HYP-2645, MISTAKE-078, OPEN-Q-108.
