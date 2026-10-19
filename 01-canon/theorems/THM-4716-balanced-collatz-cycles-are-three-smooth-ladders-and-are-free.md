---
id: THM-4716
title: "Balanced Collatz cycles are 3-smooth ladders and are free: for coprime (q, p), the q rotations of the Christoffel necklace of slope p/q are the lower mechanical words with intercepts t/q, consecutive intercepts differ by exactly one adjacent transposition 01 -> 10 (the one moves one step left), so the carries B_t (elements of the free cycle of 3x + (2^q - 3^p)) satisfy B_t - B_(t+1) = 2^i 3^j with i the new position of the moved one and j the number of ones to its right; hence the cycle elements sorted are a ladder whose consecutive differences are single 3-smooth numbers, and gcd(2^q - 3^p, B_t) = 1 for every t (a prime of the clock dividing one element divides all, hence a 3-smooth difference): no sporadic balanced cycle exists for any 3x + d, the balanced cycles of 3x+1 are exactly Gersonides' shapes (1,0), (1,1), (2,1), (3,2), and (THM-4714) the cardioid bulbs of the Mandelbrot set other than the cusp, the 1/2-bulb and the 2/3-bulb carry no integer 3x+1 cycle"
status: >
  PROVED (elementary: the intercept-shift lemma for mechanical words, the carry formula, and the linear dynamics modulo a
  prime of the clock). FINITE-EXACT confirmation: smooth_ladder_proof_check.py verifies the swap, the difference 2^i 3^j, the
  monotonicity and the coprimality for all Christoffel words with q <= 60 (44,230 checks); balanced_free.py the coprimality
  to q = 150 (686,732 rotations); THM-4715's census (120 balanced cycles, all free). Promotes HYP-9334 to a theorem. The
  3-smooth-ladder description is new to the repository; the session did not find it in Halbeisen-Hungerbuehler's rotation
  cycles, but no priority is asserted. Pre-audit.
session: mac-mini-2026-10-19-farey
source: 05-knowledge/results/farey_mandelbrot_20261019_collatz_angles_and_cross_conjecture_dictionary.md
scripts:
  - 04-computation/experiments/farey_mandelbrot_20261019/smooth_ladder_proof_check.py (+ .out), smooth_differences.py (+ .out), balanced_free.py (+ .out, balanced_free_150.out)
related:
  - THM-4484 (free and sporadic cycles), THM-4714 (balanced = cardioid bulb), THM-4715 (free cycles realise every angle orbit), HYP-9334 (now PROVED)
  - THM-4633 (large-minimum cycles have convergent shapes), THM-4591 (the five integer cycles of 3x+1)
  - Halbeisen-Hungerbuehler, Acta Arith. 78 (1997) (rotation cycles extremal for the minimum element)
---

# THM-4716 — balanced Collatz cycles are 3-smooth ladders and are free

## Setting

* `T_d(n) = n/2` (n even), `(3n + d)/2` (n odd). For a parity word `w = (w_0, …, w_(q−1))` with `p` ones, the affine map is `F_w(x) = (3^p x + B_w)/2^q` with carry `B_w = Σ_(i: w_i = 1) 3^(j_i) 2^i`, where `j_i` is the number of ones at positions `> i`. For `D = 2^q − 3^p` the free cycle of `T_D` with word `w` has elements `B_r`, `r` running over the rotations of `w` (THM-4484; `B_r = D·(fixed point of F_r)`).
* **Mechanical words.** For coprime `(q, p)` and intercept `t/q`, `0 ≤ t < q`, the lower mechanical word is `s^(t)_k = ⌊((k+1)p + t)/q⌋ − ⌊(kp + t)/q⌋`, `k = 0, …, q−1`. These `q` words are exactly the `q` rotations of the Christoffel (balanced) necklace of slope `p/q`.

## Statements

1. **Intercept-shift lemma (PROVED).** For `0 ≤ t < q − 1`, the words `s^(t)` and `s^(t+1)` differ in exactly two adjacent positions `k_0, k_0 + 1`, where `(s^(t)_(k_0), s^(t)_(k_0+1)) = (0, 1)` and `(s^(t+1)_(k_0), s^(t+1)_(k_0+1)) = (1, 0)`: one "1" moves one step to the left.

2. **3-smooth ladder (PROVED).** With `B_t = B_(s^(t))`, `B_t − B_(t+1) = 2^(k_0) 3^(j)`, where `j` is the number of ones to the right of position `k_0 + 1` in `s^(t)`. Hence `B_0 > B_1 > … > B_(q−1)`, the sorted elements of the balanced free cycle of `3x + (2^q − 3^p)` are `B_(q−1) < … < B_0`, and each consecutive difference is a single 3-smooth number `2^i 3^j`.

3. **Coprimality (PROVED).** `gcd(2^q − 3^p, B_t) = 1` for every `t`.

4. **Consequences (PROVED).**
   * (a) A balanced parity word of shape `(q, p)` is an integer cycle of `3x + d` iff `(2^q − 3^p) | d`: balanced cycles are free, and no sporadic cycle of any `3x + d` is balanced.
   * (b) The balanced integer cycles of `3x + 1` are exactly those with `|2^q − 3^p| = 1`, i.e. (Gersonides, THM-4484) the shapes `(1,0), (1,1), (2,1), (3,2)`: the cycles `{0}, {−1}, {1, 2}, {−5, −7, −10}`.
   * (c) By THM-4714 (1d): the root-angle orbits of the cardioid `p/q`-bulbs of the Mandelbrot set are integer `3x+1` cycles only for the cusp, the `1/2`-bulb and the `2/3`-bulb. A positive nontrivial `3x+1` cycle would have an unbalanced parity word, landing off the main cardioid.
   * (d) For every coprime `(q, p)` the free balanced cycle of `3x + (2^q − 3^p)` has `max − min = Σ_(t<q−1) 2^(i_t) 3^(j_t)`; for the convergent slopes `(q,p) = (3,2), (5,3), (8,5), (11,7), (13,8), …` the `q − 1` differences are the consecutive 3-smooth numbers of a window (e.g. `48, 54, 64, 72, 81, 96, 108` for `3x + 13`); in general they are the 3-smooth numbers `2^(i_t) 3^(j_t)` of the moved ones (FINITE-EXACT census in `smooth_differences.out`, 489 cycles with `q ≤ 40`, all differences 3-smooth; 14 windows consecutive).

## Proofs

**1.** Put `f_t(k) = ⌊(kp + t)/q⌋`. Increasing `t` by one increases `f_t(k)` by exactly one for those `k` with `(kp + t) mod q = q − 1`, and leaves the other values unchanged. Since `gcd(p, q) = 1`, the residues `kp + t mod q`, `k = 0, …, q`, take each value of `Z/q` once for `k = 0..q−1`, and `k = q` repeats `k = 0`; the value `q − 1` is attained at exactly one `k_1 ∈ {0, …, q−1}` (and again at `k_1 + q`). As `s^(t)_k = f_t(k+1) − f_t(k)`, the only letters that change are `s_(k_1 − 1)` (which gains 1) and `s_(k_1)` (which loses 1), indices mod `q`; since letters are 0 or 1, the change is `(0, 1) ↦ (1, 0)` at `(k_1 − 1, k_1)`. For `0 ≤ t < q − 1` the index `k_1 = q − 1 − t` (mod `q`) satisfies `1 ≤ k_1 ≤ q − 1`, so the two positions `k_0 = k_1 − 1` and `k_1` are adjacent without wrap-around. ∎

**2.** Only the one at position `k_0 + 1` moves, to position `k_0`; the ones to its right are unchanged, so its exponent `j` of 3 is unchanged, and the carry changes by `3^j(2^(k_0) − 2^(k_0+1)) = −2^(k_0) 3^j`. Strict monotonicity follows, and the sorted order of the elements is the intercept order. ∎

**3.** Let `ℓ` be a prime dividing `D = 2^q − 3^p`; `ℓ ≠ 2` (`D` is odd) and `ℓ ≠ 3` (`D ≡ 2^q mod 3`). Modulo `ℓ` the map `T_D` is linear: `n ↦ 3n·2^(−1)` or `n·2^(−1)`. So if `ℓ` divides one element of the cycle it divides every element, in particular `B_t` and `B_(t+1)` for some `t < q − 1` (there are `q ≥ 2` elements), hence their difference `2^(k_0) 3^j`, which is impossible. ∎

**4.** (a) The cycle element is `d B_t / D`, an integer for all `t` iff `D | d B_t`, iff `D | d` by 3. (b) `D | 1` iff `|D| = 1`, the Gersonides shapes of THM-4484. (c) THM-4714 (1d) identifies balanced words with cardioid-bulb root orbits; by (b) only `|D| = 1` shapes are `3x+1` cycles. (d) is 2 and the census. ∎

## Remark

The same argument applies verbatim to `qx + d` maps with the carry `Σ q^(j_i) 2^i` and clocks `2^K − q^X` (THM-4484's setting): balanced cycles of `y ↦ y/2, (qy + d)/2` are `{2, q}`-smooth ladders and are free.
