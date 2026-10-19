---
id: THM-4715
title: "Free cycles realise every external-angle orbit: for coprime (q, p) with clock D = 2^q - 3^p, the integer cycles of 3x + D of shape (q, p) are in bijection with the necklaces of length q with p ones, i.e. with the doubling orbits of period-q angles having p ones in the repetend (THM-4484 free case); exactly one of them, the Christoffel necklace, is the root-angle orbit of the cardioid p/q-bulb, the others land at period-q roots off the cardioid; 3x+13 realises all seven period-8 orbits with five ones (bulb cycle min 319), 3x+5 both period-5 orbits with three ones (bulb 23, other 19); the census of 3x+d, d odd <= 49, |n| <= 20000, has 219 cycles, 120 balanced, every balanced one free; the sporadic cycles sit at Stern-Brocot nodes of log_3 2 when their minimum is large (the -17 cycle at the semiconvergent 7/11, D = -139; the 27-cycles of 3x+5 and 3x+15 at 17/27, D = 5077565)"
status: >
  PROVED: statements 1-2 (THM-4484's free case plus THM-4714 (1d)). FINITE-EXACT: statement 3 (the census; the
  Stern-Brocot placements). CITED: THM-4633 (Legendre criterion: a cycle whose minimum is large has K/X a convergent
  of log2 3) for the "when the minimum is large" clause. Pre-audit.
session: mac-mini-2026-10-19-farey
source: 05-knowledge/results/farey_mandelbrot_20261019_collatz_angles_and_cross_conjecture_dictionary.md
scripts:
  - 04-computation/experiments/farey_mandelbrot_20261019/collatz_angles.py (+ .out), balanced_free.py (+ .out, balanced_free_150.out)
related:
  - THM-4484 (free iff (2^p - q^a) | d), THM-4714, THM-4633, THM-4591, HYP-9334
---

# THM-4715 — free cycles realise every external-angle orbit

## Statements

1. **Bijection (PROVED).** Let `q > p ≥ 1` be coprime and `D = 2^q − 3^p`. For every parity necklace `w` of length `q` with `p` ones, the affine map `F_w(x) = (3^p x + B_w)/2^q` has the integer fixed point `B_w` (the carry), and `{B_r : r a rotation of w}` is an integer cycle of `T_D(n) = n/2, (3n + D)/2` with parity word `w` (THM-4484, free case `D | d` with `d = D`). Every integer cycle of `T_D` of shape `(q, p)` arises this way. Hence the cycles of shape `(q,p)` of `3x + D` correspond bijectively to the doubling orbits of the angles `m/(2^q − 1)` with `p` ones in the repetend.

2. **The bulb cycle (PROVED).** Exactly one of these cycles has a balanced word, the Christoffel necklace of slope `p/q`; by THM-4714 (1d) its angle orbit is the root-angle orbit of the cardioid `p/q`-bulb. The others land at period-`q` roots off the cardioid. Examples: `3x + 13` (`(q,p) = (8,5)`) has seven cycles, one per necklace, with minima `211, 227, 251, 259, 283, 287, 319`; the bulb cycle is `319 → 485 → 734 → 367 → 557 → 842 → 421 → 638` (word `11011010`), the `5/8`-bulb's root angles being `181/255, 182/255`. `3x + 5` (`(5,3)`) has the two cycles `23` (balanced, the `3/5`-bulb, root angles `21/31, 22/31`) and `19`.

3. **Census and Stern–Brocot placement (FINITE-EXACT).** For odd `d ≤ 49` and `|n| ≤ 20000` the map `T_d` has 219 cycles, 120 of them balanced; every balanced cycle is free (`D | d`), in line with HYP-9334. The sporadic cycles with large minima sit at Stern–Brocot nodes of `log_3 2 = [0; 1, 1, 1, 2, 2, 3, 1, 5, 2, 23, …]`: the `−17` cycle of `3x+1` has shape `(11, 7)`, the semiconvergent `7/11` (`D = −139`); the cycles of `3x+5` with minima `187, 347` and of `3x+15` with minima `561, 1041` have shape `(27, 17)`, the semiconvergent `17/27` (`D = 5077565`); the free cycles of `3x+13` sit at the convergent `5/8`. The bulb cycles along the convergents `1/2, 2/3, 5/8, 12/19, 41/65, 53/84, …` are the free cycles of `3x + 1, 3x − 1, 3x + 13, 3x − 7153, 3x + 420491770248316829, …`; for none of them beyond `|D| = 1` does `D` divide the carry (so none is a `3x+1` cycle).

## Proofs

**1.** THM-4484: for `d = D` every parity word of shape `(q, p)` is integral, the element at the start of the word being `d·B_w/D = B_w`. Conversely a cycle of shape `(q, p)` of `T_D` has an odd denominator-free fixed point `B_w d/D = B_w`. Distinct necklaces give distinct cycles because the parity vector determines the 2-adic point (Lagarias). The angle of the element `B_r` is `0.\overline{r}`, so the cycle is the doubling orbit of its angle. ∎

**2.** A necklace of length `q` with `p` ones, `gcd(p,q) = 1`, is balanced iff it is the Christoffel necklace (the unique balanced necklace of its slope); THM-4714 (1d). The listed elements and angles are computed in `collatz_angles.py` (`bulb_angles(5, 8)` gives `181/255, 182/255`). ∎

**3.** Direct search (`cycles_3xd`) with cycle detection and bound `|n| ≤ 50·20000` on the orbit; balancedness by the cyclic factor test. The convergents and semiconvergents are computed from the continued fraction of `log_3 2`; `D | B_w` is tested exactly for the Christoffel words of the first eight convergents. ∎
