---
id: HYP-9256
title: "Finite Beckman-Quarles for the unit-circle Cayley graphs: for every prime power q not in {4, 8, 9}, every automorphism of G_q = Cay(F_(q^2), mu_(q+1)) is an affine semilinear map z -> a sigma(z) + b (a in mu_(q+1), sigma in Gal(F_(q^2)/F_p), b in F_(q^2)), so |Aut G_q| = q^2 (q+1) 2e for q = p^e"
status: "OPEN (2026-10-12, mac-mini-2026-10-12-logic4). FINITE-EXACT for every prime power q <= 256 (THM-4624 for q <= 81; audit B of 2026-10-12 for 83 <= q <= 256, all 38 prime powers, including the split-norm cases q = 1 mod 4). The exceptions q <= 256 are the basis/frame cases whose frame group exceeds Aff(q), namely q in {4, 8, 9} (Clebsch, H(3,4), ternary five-frame); the basis/frame cases are five (q in {2, 3, 4, 8, 9}; for q = 2, 3 the frame group equals Aff(q)), and THM-4624(3) shows no other q has a basis/frame structure, so this mechanism is exhausted. What is open is whether non-affine automorphisms can arise without it, for q > 256."
session: mac-mini-2026-10-12-logic4
related:
  - THM-4624 (the census and the exceptional identifications)
  - THM-4558 (kappa(q) = chi(G_q) in the residue-colouring bound)
  - M. Kiermaier, S. Kurz, Discrete Math. 309 (2009): Aut = AffAut for the integral-distance graph on F_q^2 when q = 3 mod 4 or q is a prime other than 5
---

# HYP-9256 — affine rigidity of the finite unit circles

## Statement

For every prime power `q = p^e` with `q ∉ {4, 8, 9}`, `Aut(G_q) = Aff(q)`, the group of affine
semilinear maps, of order `q^2 (q+1) · 2e`.

## Why it is plausible

* The census holds for all `70` prime powers `q ≤ 256` (THM-4624 and audit B), including all even
  `q ∉ {4, 8}`, all `q` with `p ≥ 5`, and all split-norm `q ≡ 1 (mod 4)`.
* The only mechanism found for extra symmetry — the unit circle being a basis or frame over a subfield —
  is exhausted by `q ∈ {2, 3, 4, 8, 9}` (elementary number theory, THM-4624(3)).
* Kiermaier–Kurz proved affine rigidity for the integral-distance graph (edges at square distances),
  a coarser relation containing `G_q` when `q ≡ 3 (mod 4)`, by the theory of directions determined by
  point sets and Paley cliques; their exception `q = 5` is not an exception here (`|Aut G_5| = 300`).

## Suggested proof route

Show that an automorphism of `G_q` maps the norm circles `{N(z - z_0) = n}` to norm circles (common
neighbourhood counting: two unit circles meet in at most two points, and the pattern of common
neighbours of two vertices determines the norm of their difference for large `q`), then that it
preserves lines, then apply the fundamental theorem of affine geometry.

## Decisive tests

The census to `q ≤ 256` is done (audit B, `audit_B_gq_census_*.out`; the split cases `q ≡ 1 (mod 4)`
show nothing special). Next: `256 < q ≤ 1024` by sparse nauty (up to `2^20` vertices), or a proof along
the common-neighbourhood route.
