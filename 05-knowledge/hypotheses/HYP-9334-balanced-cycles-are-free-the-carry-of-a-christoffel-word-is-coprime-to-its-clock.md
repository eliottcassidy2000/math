---
id: HYP-9334
title: "Balanced cycles are free: for every Christoffel (balanced) necklace w of shape (q, p), gcd(2^q - 3^p, B_r) = 1 for every rotation r of w (B = the carry of the affine map); hence a balanced parity word is an integer cycle of 3x + d only when (2^q - 3^p) | d, there are no sporadic balanced cycles, the only balanced cycles of 3x+1 are the Gersonides ones (shapes (1,0), (1,1), (2,1), (3,2)), and the cardioid bulbs of the Mandelbrot set beyond the cusp, the 1/2-bulb and the 2/3-bulb carry no integer 3x+1 cycle"
status: >
  OPEN. FINITE-EXACT: gcd = 1 for all 686,732 rotations of all Christoffel necklaces with q <= 150 (balanced_free_150.out);
  by contrast random necklaces of the same shapes have gcd > 1 in 73 of 1194 trials (6%). The consequence for 3x+d is
  FINITE-EXACT in the census of THM-4715 (3) (every balanced cycle there is free). No proof: the carry
  B_w = sum_k 3^(p-k) 2^(ceil(kq/p) - 1) of the lower Christoffel word has no obvious structure modulo the primes of
  2^q - 3^p. What would refute: one Christoffel rotation with a common prime factor with its clock; what would prove it: an
  identity for B_w modulo 2^q - 3^p (e.g. a closed form of the free bulb cycle's elements).
session: mac-mini-2026-10-19-farey
source: 05-knowledge/results/farey_mandelbrot_20261019_collatz_angles_and_cross_conjecture_dictionary.md
related:
  - THM-4484 (free and sporadic cycles), THM-4714 (balanced = cardioid bulb), THM-4715 (free cycles realise every angle orbit)
  - THM-4633 (budget: large-minimum cycles have convergent shapes), THM-4591 (the five integer cycles)
  - 04-computation/experiments/farey_mandelbrot_20261019/balanced_free.py
---

# HYP-9334 — balanced cycles are free

## Statement

Let `w` be the Christoffel necklace of slope `p/q` (`gcd(p,q) = 1`, `q > p ≥ 1`), `D = 2^q − 3^p`, and for each rotation `r` of `w` let `B_r` be the carry of `F_r(x) = (3^p x + B_r)/2^q`. Then `gcd(D, B_r) = 1`.

## Consequences (PROVED from the statement)

* A balanced word of shape `(q, p)` is an integer cycle of `3x + d` iff `D | d` (the cycle element is `d B_r / D`).
* No sporadic balanced cycles exist for any `d`; the balanced cycles of `3x+1` are exactly those with `|D| = 1`, i.e. Gersonides' shapes `(1,0), (1,1), (2,1), (3,2)` (THM-4484), which are the cusp, the cusp, the `1/2`-bulb and the `2/3`-bulb of THM-4714.
* In Mandelbrot language: a positive-integer cycle of `3x+1` other than `{1, 2}` would have an unbalanced parity word, so its external angles land off the main cardioid.

## Why it might be true

* Modulo a prime `ℓ | D` the map `T_D` is linear (`n ↦ 3n/2` or `n/2`), so `ℓ | B_r` for one rotation iff for all; the statement says the balanced free cycle is never a multiple of a cycle of `3x + D/ℓ`.
* The balanced cycle is the one with the most evenly spread carries; its elements are the "rotation cycle" of Halbeisen–Hungerbühler, extremal for the minimum element. No divisibility mechanism is known to the session.

## Tests

* `balanced_free.py 150`: all rotations of all Christoffel words with `q ≤ 150`, exact big-integer gcds.
