---
id: HYP-9782
title: "Generic-line large sieve: after the explicit lattice-ray families, the open core of the row-collision lemma (HYP-9752) is a large-sieve estimate for digital Riesz products prod_k (1 + e(h c_k/d)) summed across signatures, where p(x) - p(y) = c_0 + sum_k c_k z_k is affine in the shared bits with explicit coefficients"
status: >
  OPEN (CONJECTURE), klein-2026-11-15-llt.
  PROVED: the affine lemma for p(x) - p(y) at fixed signature, and the exact Riesz-product formula for the collision
  count.
  NUMERICAL: generic cores of 2U + p hit their lattice lines at 1.00-1.15 times the random-divisibility prediction
  (spans 10-13, common suffixes 0-4, prefix scales 2^3-2^5); d | rhs at rate 0.00301 against mean 1/d 0.00276.
  NOT CLAIMED: HYP-9752; positive density; Chan-Pak.
session: klein-2026-11-15-llt
date: 2026-11-15
source: 05-knowledge/results/chanpak_generic_lines_llt_20261115/GENERIC_LINES.md
scripts:
  - 05-knowledge/results/chanpak_generic_lines_llt_20261115/affine_lemma.py, equi_test2.py, residue_bias.py
related:
  - HYP-9752 (row collision lemma), THM-5232 (reverse-complement families), chanpak_cone_collision_placement_20261115/TWISTED_LATTICES.md
---

# HYP-9782 - the generic-line large sieve

**Setting.** A collision of G(u) = cu + kp(u) + φ(u mod 2) has:
- a two-ended core (x, y) of span a, with δ = x − y and d = rev_a x − rev_a y;
- a common suffix of length j;
- a prefix w.

Group the cores by signed difference σ ∈ {−1, 0, 1}^a, and let Z be the set of shared positions.

**Affine lemma (PROVED).** p(x) − p(y) = c₀(σ) + Σ_(k∈Z) c_k(σ) z_k, where

    c_k = −Σ_(j>k) σ_j 2^(j−k−1) − Σ_(i<k) σ_i 2^(k−i−1) − σ_(k±1).

So in the j = 0 sector the number of cores of signature σ on a line is

    #{z : d | A_σ + Σ c_k z_k, box} = (1/|d|) Σ_(h mod d) e(h A_σ/d) R_σ(h/d) B_σ(h),
    with R_σ(θ) = Π_(k∈Z) (1 + e(θc_k)).

**Conjecture.** Summed over generic signatures of span a (removing THM-5232's palindromic families and the
small-height slopes), the h ≢ 0 terms are o(Σ_σ 2^|Z(σ)| β_σ/|d_σ|), uniformly in the prefix scale. The analogue
holds for j ≥ 1, with the r-kernel.

**Why cancellation must come from across signatures.** Weighted by multiplicity, a signature has |Z| ≈ a/2 while
|d| ≈ 2^a. So its main term is about 2^(−a/2), and no single σ has an asymptotic.

**Why it should hold.** The geometric mean of |1 + e(φ)| is 1 (the Mahler measure of 1 + z), so off-resonance Riesz
products are e^(o(a)). Resonance needs the forward and reversed partial differences to align, which is the lattice-ray
structure already removed.

**Consequence.** With THM-5232 and the concurrent short-span bounds, it gives HYP-9752 in the form
N(l, w) ≪ 2^l/(1 + w). Hence the carrier rows have linear energy, and the Hamiltonian-path spectrum has positive
lower density.
