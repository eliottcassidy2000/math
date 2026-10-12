---
id: HYP-9256
title: "Finite Beckman-Quarles for the unit-circle Cayley graphs: for every prime power q not in {4, 8, 9}, every automorphism of G_q = Cay(F_(q^2), mu_(q+1)) is an affine semilinear map z -> a sigma(z) + b (a in mu_(q+1), sigma in Gal(F_(q^2)/F_p), b in F_(q^2)), so |Aut G_q| = q^2 (q+1) 2e for q = p^e"
status: "OPEN (2026-10-12, mac-mini-2026-10-12-logic4). FINITE-EXACT for every prime power q <= 81 (THM-4624, nauty). The three exceptions are PROVED to be exactly the cases where mu_(q+1) is a basis or a frame over a subfield (Clebsch, H(3,4), ternary five-frame), and the subfield-line criterion admits no further exception for any q (THM-4624(3)); what is open is whether non-affine automorphisms can arise without a subfield-line structure."
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

* The census holds for all `32` prime powers `q ≤ 81` (THM-4624), including all even `q ∉ {4, 8}` and
  all `q` with `p ≥ 5`.
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

Extend the nauty census to `q = 83, ..., 256` (feasible: `q^2 ≤ 65536` vertices), and separately test the
Paley-type cases `q ≡ 1 (mod 4)` where the norm form is split over `F_q`.
