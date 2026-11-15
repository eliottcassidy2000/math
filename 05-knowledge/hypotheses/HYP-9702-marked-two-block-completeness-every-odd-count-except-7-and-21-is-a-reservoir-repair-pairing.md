---
id: HYP-9702
title: "Marked two-block completeness: every odd N other than 7 and 21 is H(T) for a tournament T made of an ordinal sum L => B of two blocks and one apex vertex with arbitrary arcs to them; equivalently, the bilinear pairing N = a_B (n_L - g) + g n_B + s_L t_B + d_L f_B of a marked reservoir's right-profile with a marked repair's left-profile covers every odd N >= 23. Blocks of order <= 7 already cover [23, 1,542,343), and orders (8,7) cover [23, 8,667,349)"
status: >
  OPEN (CONJECTURE), mac-mini-2026-11-15-slides, wave 26 of chanpak_waves_20261113.
  FINITE-EXACT (THM-5182 (E)):
  - block orders (<=6, <=6): every odd N in [23, 67,967);
  - block orders (<=7, <=7): every odd N in [23, 1,542,343);
  - block orders (<=8, <=7): every odd N in [23, 8,667,349);
  - 1, 3, 5, 9, ..., 19 are also attained, so every odd N < 8,667,349 other than 7 and 21 is a two-block count;
  - 160 sampled witnesses re-counted by direct DP.
  NUMERICAL: the coverage threshold grows by factors 22.7 (orders 6 -> 7) and 5.6 (7 -> (8,7)).
  Fixed-interface palettes have mean size 1.05-1.09, so the coverage is not a repair-window effect.
  Pre-registered P81 REFUTED (registered threshold >= 3*10^7 at (8,7)).
session: mac-mini-2026-11-15-slides
date: 2026-11-15
source: 05-knowledge/results/chanpak_waves_20261113/README.md#wave-26-fans-slides-and-reservoirs--unit-extensions-never-branch-twice-a-mod-3-law-for-passive-channels-and-marked-two-block-completeness
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/palette.c (profiles), twoblock.c (coverage),
    twowit.c and verify_two.py (witnesses), twomult.c (multiplicities)
related:
  - THM-5182 (pair expansion, mod-3 law), HYP-9692 (band fullness; the cyclic join is the case L all marked,
    B all unmarked), THM-5172 (congruences), chanpak_window_descent_20261115 (reservoir and repair language),
    THM-4094 (census to 2^32)
---

# HYP-9702 - marked two-block completeness

**The family.** Take two marked blocks L (the reservoir) and B (the repair). Form the ordinal sum L ⇒ B, and add
one apex that beats exactly the marked vertices. By THM-5182 (A) the count is

    N = a_B (n_L − g) + g n_B + s_L t_B + d_L f_B,

with:
- (g, n_L, s_L, d_L): L's ordinary count, apex count, and total and signed marked starts;
- (a_B, n_B, t_B, f_B): B's ordinary count, apex count, and total and signed unmarked ends.

So N is a bilinear pairing of the reservoir's right-profile with the repair's left-profile. This is the owner's
architecture without any slide. The reservoir's ordinary count g scales the repair's apex count. Its marked
endpoints (s_L, d_L) set the residue class and the absolute position.

**Conjecture.** Every odd N except 7 and 21 is a two-block count.

**Evidence (FINITE-EXACT, THM-5182 (E)).**

| block orders | distinct profiles per side | odd N covered |
|---|---:|---|
| (≤ 6, ≤ 6) | 2,831 | [23, 67,967) |
| (≤ 7, ≤ 7) | 52,313 | [23, 1,542,343) |
| (≤ 8, ≤ 7) | 1,693,156 × 52,313 | [23, 8,667,349) |

The values 1, 3, 5, 9, ..., 19 are attained too. Converse duality makes the (8,7) and (7,8) spectra equal.

For comparison, wave 25's cyclic join is the subfamily with L all marked and B all unmarked. With cores of order
≤ 10 it covers only (36,833, 2,398,737]. The marks fill everything below 36,833 and reach 8.7·10^6 with smaller
blocks.

**What carries the coverage.**
- *Not repair windows.* A last block's palette is its fiber over (a, t, f), and these fibers have mean size
  1.05–1.09 at orders 5–8 (longest run of consecutive odd apex counts: 4). Within a fixed interface there is
  essentially no room to slide.
- *The variety of interfaces.* For a fixed tournament the marks enter (t, f) linearly, as subset sums of the
  per-vertex end weights, and n quadratically, through the cut function of THM-5182 (D). So each tournament
  contributes up to 2^m profiles. The 4D cloud of left-profiles is much denser than the 2D (H, F₂) cloud of
  HYP-9692.
- *Local densities.* On odd N in [10^6, 1.4·10^6] the (7,7) multiplicity averages 215.9 at N ≡ 0 (mod 3) and 171
  at N ≡ ±1. Primes average 171.1, exactly their classes. Residues mod 5 and 7 move it by less than 2%. The only
  local bias is at 3, consistent with the mod-3 laws of THM-5172 (B) and THM-5182 (B).

**What a proof needs.** A fullness statement for the pairing: for every large odd N, some reservoir profile ρ whose
hyperplane ⟨ρ, ·⟩ = N meets an actual repair profile. The natural form is a local limit theorem for the
left-profile (a, n, t, f) of marked tournaments, in a slab transverse to ρ. The linear dependence of (t, f) on the
marks suggests Fourier analysis over the mark cube, with the quadratic cut function as the only non-linear
coordinate. HYP-9692 is the special case without marks. The conjecture is stronger than needed for Chan–Pak below
2^32 (THM-4094) and equivalent to it above that, within this family.

**Predictions.**
- *P84.* Blocks of order (8, 8) cover every odd N in [23, X) with X ≥ 5·10^7.
- *P85.* The multiplicity ratio between N ≡ 0 and N ≡ ±1 (mod 3) tends to a limit in [1.2, 1.35]. All other
  moduli are asymptotically unbiased.
