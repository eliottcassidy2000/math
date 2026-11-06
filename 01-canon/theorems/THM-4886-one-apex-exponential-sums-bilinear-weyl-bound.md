---
id: THM-4886
title: "A bilinear Weyl bound for one-apex exponential sums. (A) If |2 theta - r/q| <= 1/q^2, then |sum over words x of length n beginning with 1 of e(theta H(x))| <= 2^(n-1) 2^(-L/2) (1 + 2(2^L/q + 1) sqrt((3 + 2^(3-a) q)(3 + 2^(2-l) q)))^(1/2) for every split n = a + L + l; hence 6.5 · 2^(n-1) q^(-1/2) for q <= 2^((n-3)/3), 12 · 2^(n-1) (q/2^n)^(1/4) for q <= 2^(n-6), and a power saving 2^(-eps n/4) on the minor arcs. (B) The same for every one-apex readout y_h H - y_e e + y_R (1 + A). (C) Fourier uniformity: the sum is o(2^n) for every theta outside (1/2)Z. Positive density of the spectrum would follow from linear energy E_n = O(2^n) (observed 3.25 |X_n|); density one and Chan-Pak OPEN"
status: >
  PROVED: (A), (B), (C); new, using Session 193's exact concatenation law (modular_sync/PROOF.md eq. (2)). (C) is new
  for irrational theta; for rational theta it is Session 193's Theorem 3 transferred to all words.
  FINITE-EXACT / numerical, pre-registered as HYP-9479 (commit 185c2b3efa) and CONFIRMED: square-root cancellation at
  five irrational theta to n = 28; E_n / |X_n| = 3.25 at n = 26, 28; Q / diagonal in [0.997, 1.069] at unseen splits.
  NOT CLAIMED: linear energy, positive density, density one, Chan-Pak (all OPEN). (A) is the minor-arc half of a circle
  method; the energy and pair-correlation half is not supplied.
  AUDIT: pending.
session: opus-2026-11-06-S83
source: 05-knowledge/results/one_apex_exponential_sums_bilinear_weyl_bound_20261106.md
related:
  - THM-4885
  - THM-4881
  - THM-4875
  - HYP-9479
  - HYP-9474
  - HYP-9478
  - 05-knowledge/results/session193_triangulated_proof_20261105/modular_sync/PROOF.md
  - 05-knowledge/results/feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md
  - 05-knowledge/results/seam_gadget_thickness_spectrum_20261104.md
verification:
  - 04-computation/experiments/one_apex_fourier_20261106/concat_check.py
  - 04-computation/experiments/one_apex_fourier_20261106/bound_check.py
  - 04-computation/experiments/one_apex_fourier_20261106/core_check.py
  - 04-computation/experiments/one_apex_fourier_20261106/const_check.py
  - 04-computation/experiments/one_apex_fourier_20261106/insample.py
  - 04-computation/experiments/one_apex_fourier_20261106/outsample.py
  - 04-computation/experiments/one_apex_fourier_20261106/outsample.out
---

# THM-4886 - a bilinear Weyl bound and Fourier uniformity for one-apex counts

**Notation** (Session 193).
- X_n is the set of binary words of length n beginning with 1, and H(x) = g(A(x)) is the one-apex Hamiltonian path count.
- A(w) is the binary value of w, and C(w) = sum_j (1 - w_j) 2^(j-1).
- p(w) is the last bit, and eps(w) = [w begins in 0].
- S_n(theta) = sum_(x in X_n) e(theta H(x)).

**Concatenation law** (Session 193 eq. (2)). For |v| = L:
H(uvw) = H(u) + H(v) + H(w) - 2 + A(u) C(v) + (2^L A(u) + A(v)) C(w) + p(u) eps(v) + p(v) eps(w).

## (A) The bilinear Weyl bound

For n = a + L + l (a >= 2, L, l >= 1) and r/q in lowest terms with |2 theta - r/q| <= 1/q^2:

    |S_n(theta)| <= 2^(n-1) 2^(-L/2) (1 + 2 (2^L/q + 1) sqrt((3 + 2^(3-a) q)(3 + 2^(2-l) q)))^(1/2).

Corollaries:
- **(A1)** |S_n| <= 6.5 · 2^(n-1) q^(-1/2) for 2 <= q <= 2^((n-3)/3).
- **(A2)** |S_n| <= 12 · 2^(n-1) (q/2^n)^(1/4) for 2^((n-6)/3) <= q <= 2^(n-6).
- **(A3)** If no r/q with q <= 2^(eps n) has |2 theta - r/q| <= 1/(q 2^((1-eps) n)) (0 < eps <= 1/3), then
  |S_n| <= 12 · 2^(n-1) 2^(-eps n/4).

**Proof idea.**
- Fix the boundary bits of u and w. Then A(u) and C(w) run over full arithmetic progressions, and the middle block v
  enters through the bilinear terms A(u) C(v) + A(v) C(w).
- Cauchy-Schwarz over (u, w), then geometric sums over the progressions.
- A counting lemma: in a block of q consecutive d, at most 2qt + 3 have ||beta d|| <= t. This bounds the off-diagonal
  pairs of middle words, because A and C are bijections onto [0, 2^L).

## (B) Readouts

The sum of e(theta (y_h H(x) - y_e [x ends in 0] + y_R (1 + A(x)))) over X_n satisfies (A), with 2 theta replaced by
2 y_h theta. These are the one-apex x bank families of THM-4885.

## (C) Fourier uniformity

S_n(theta) = o(2^n) for every real theta outside (1/2)Z, and S_n(theta + 1/2) = -S_n(theta).
- For irrational theta: (A1) along the convergents of 2 theta.
- For rational theta: Session 193 Theorem 3 (mixing of H mod M on strong words), transferred through x = x' 1^t.

## Toward density one (not claimed)

- If the energy E_n = #{(x, x') in X_n^2 : H(x) = H(x')} is O(2^n), Cauchy-Schwarz gives positive lower density of the
  spectrum. The observed value is 3.25 |X_n| at n <= 28.
- Density one needs, in addition, pair correlations of the bank families.
- (A)-(C) supply the minor-arc side only.

**Proofs.** Source note, sections 1-3.
