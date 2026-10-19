---
id: THM-4714
title: "Collatz parity numbers are Mandelbrot external angles and Minkowski-Farey addresses: (1) the parity number rho(x) of THM-4634 read as an angle makes every eventually periodic Terras orbit a periodic external angle of the Mandelbrot set, landing (Douady-Hubbard) at the root of a hyperbolic component whose period is the Terras period; the trivial cycle {1,2} is the root of the 1/2-bulb (angles 2/3, 1/3), the cycle {-5,-7,-10} is the root of the 2/3-bulb (6/7, 5/7) together with the airplane root c = -7/4 (3/7), the fixed points 0 and -1 are the cusp c = 1/4 (angle 0), the limit anchor -1/3 is the tip c = -2 (angle 1/2), the anchor orbits c_a = 1/(2^a-3) carry the root angles 1/(2^a-1), 2/(2^a-1) of the 1/a-bulbs accumulating at the cusp, and the -17 cycle is a period-11 orbit landing off the cardioid; a cycle word is balanced iff its angle orbit is a cardioid-bulb root orbit (the Christoffel necklace of slope X/K); (2) the Farey address alpha(x) = ?^-1(rho(x)) (Minkowski's question mark) conjugates the Terras map to the Farey-type map x/(1-x), 2 - 1/x, whose continued-fraction partial quotients are the parity run lengths; the two parabolic fixed points 0 and 1 are the even and odd anchors, the Syracuse map is the Gauss jump transformation, the trivial cycle is the golden 2-cycle {1/phi, 1/phi^2}, and Collatz says that every positive integer's address is a quadratic irrational in Q(sqrt 5) whose Farey orbit ends on that golden pair; (3) the addresses of the integers are distributed by the singular Minkowski measure"
status: >
  PROVED (elementary, given THM-4634 and the cited classical theorems): statements 1 (the dictionary of angles; the
  balanced/bulb equivalence relies on the CITED Douady-Hubbard landing theorem and the Sturmian description of bulb-root
  angles), 2 (the conjugacy via ?, verified exactly on every Farey fraction of denominator <= 60 and proved from Panti's
  tent conjugacy and the symmetry ?(1-x) = 1 - ?(x)), 3 (Terras equidistribution + continuity of ?). CITED: Douady-Hubbard
  (periodic angles land at roots of hyperbolic components; the p/q-limb root angles), Devaney, AMM 106 (1999) (the Farey
  tree of bulbs and the Fibonacci path), Minkowski (1904), Denjoy (1938), Salem (1943) (the ? function and its singular
  measure), Panti, Monatsh. Math. (2008) and Bonanno-Isola, Colloq. Math. 116 (2009) (? conjugates the Farey map to the
  tent map; Stern-Brocot trees and dynamics). FINITE-EXACT: the adjacent-pair rule for bulb-root angles checked against
  the classical table (nine bulbs) and its uniqueness for all p/q with q <= 15; the five integer cycles of 3x+1. Not
  claimed: any new statement about the Mandelbrot set itself, or any progress on Collatz beyond the reformulation.
  Pre-audit.
session: mac-mini-2026-10-19-farey
source: 05-knowledge/results/farey_mandelbrot_20261019_collatz_angles_and_cross_conjecture_dictionary.md
scripts:
  - 04-computation/experiments/farey_mandelbrot_20261019/collatz_angles.py (+ .out; 1671 exact checks)
related:
  - THM-4634 (the parity-number dictionary rho: Terras map semi-conjugate to doubling), THM-4633 (cycles as rational approximations of log2 3)
  - THM-4484 (free and sporadic cycles; Gersonides), THM-4591 (the five integer cycles to period 301,993), THM-4600 (the anchors c_a)
  - THM-4515 (fair splits: concatenation of blocks is the mediant of their slopes), THM-4556 (the Mersenne run as the odd anchor)
  - THM-4715 (free cycles of 3x+(2^q-3^p) as the period-q angle orbits), HYP-9334 (balanced cycles are free), HYP-9335
  - 05-knowledge/results/collatz_parity_cf_boundary_20261014a.md (codex: the parity-bit continued fraction beta_n and the Gauss map; the Farey address differs: its partial quotients are the run lengths)
  - 05-knowledge/results/parabolic_scaling_session_20261019.md (codex: pi-clocks at the cusp c = 1/4; the airplane c = -7/4 residue 51/49)
---

# THM-4714 — Collatz parity numbers are Mandelbrot external angles and Minkowski–Farey addresses

## Setting

* `T(x) = x/2` (x even), `(3x+1)/2` (x odd) on `Z_2`; THM-4634's parity number `ρ(x) = Σ_i b_i 2^(−i−1) ∈ [0,1]`, `b_i = T^i(x) mod 2`, satisfies `ρ(T x) = 2ρ(x) mod 1`. For eventually periodic orbits `ρ(x)` is rational; for a purely periodic orbit of Terras period `K` with parity word `w`, `ρ(x) = 0.\overline{w} = m/(2^K − 1)`.
* **External angles.** The Mandelbrot set `M` has external rays `R_θ`, `θ ∈ R/Z`; `θ` rational with odd denominator is periodic under doubling. **CITED (Douady–Hubbard):** a periodic angle of exact period `K` lands at the root (or co-root) of a hyperbolic component of period `K`; each root of a satellite `p/q`-bulb of the main cardioid is the landing point of exactly two angles, whose binary expansions are the two rotations of the Christoffel (balanced, Sturmian) necklace of slope `p/q` that are numerically adjacent, i.e. differ by `1/(2^q − 1)`; the angle `0` lands at the cusp `c = 1/4`, and the angle `1/2` at the tip `c = −2`.
* **Minkowski's question mark.** `?([0; a_1, a_2, …]) = 2 Σ_k (−1)^(k+1) 2^(−(a_1+…+a_k))`; `?` is an increasing homeomorphism of `[0,1]`, `?((a+c)/(b+d)) = (?(a/b) + ?(c/d))/2` on Farey neighbours, `?(1−x) = 1 − ?(x)`, `?` maps quadratic irrationals onto non-dyadic rationals (CITED: Minkowski, Denjoy, Salem).
* **Farey-type map.** `F̃(x) = x/(1−x)` on `[0, 1/2]`, `F̃(x) = 2 − 1/x` on `[1/2, 1]`. Its fixed points `0` and `1` are parabolic (`F̃'(0) = F̃'(1) = 1`).

## Statements

1. **The angle dictionary (PROVED from THM-4634 + CITED landing theorem).** Read `ρ(x)` as an angle `θ(x) ∈ R/Z`. Then `θ(Tx) = 2θ(x)`, and:
   * (a) `x` has an eventually periodic Terras orbit iff `θ(x)` is rational with odd denominator; a purely periodic orbit of period `K` gives a periodic angle of period dividing `K`, hence (Douady–Hubbard) the landing point of `R_θ` is the root of a hyperbolic component of that period.
   * (b) The five integer cycles of `3x+1` (THM-4591): `{0}` and `{−1}` have angles `0` and `1 ≡ 0` — the cusp `c = 1/4`; `{1, 2}` has angles `2/3, 1/3` — the root `c = −3/4` of the `1/2`-bulb; `{−5, −7, −10}` has angles `6/7, 5/7, 3/7` — `6/7, 5/7` are the root angles of the `2/3`-bulb and `3/7` lands at the airplane root `c = −7/4`; the `−17` cycle has eleven angles with denominator `2047 = 2^11 − 1` (`ρ(−17) = 1976/2047`), landing at period-11 roots off the main cardioid, because its word `11110111000` is not balanced.
   * (c) **Anchors.** The anchor `c_a = 1/(2^a − 3)` (THM-4600) has the periodic word `1 0^(a−1)`; its orbit carries the angles `{2^(a−1), 1, 2, …, 2^(a−2)}/(2^a − 1)`, among which `1/(2^a−1), 2/(2^a−1)` are the root angles of the `1/a`-bulb. As `a → ∞` the bulbs accumulate at the cusp, `c_a → −1/3` 2-adically, and `θ(−1/3) = 1/2`: the limit anchor is the tip `c = −2`. The whole real spine `[−2, 1/4]` of `M` is thus spanned by the two degenerate Collatz points `−1/3` (tip) and `0 ≡ −1` (cusp).
   * (d) **Balanced words are bulbs.** A periodic word `w` of shape `(K, X)` with `gcd(K, X) = 1` is balanced iff its angle orbit is the root-angle orbit of the cardioid `X/K`-bulb. Hence an integer cycle of any `3x + d` is a cardioid-bulb cycle iff its parity word is a Christoffel necklace; THM-4715 and HYP-9334 classify these.

2. **The Farey address (PROVED).** Put `α(x) = ?^(−1)(ρ(x)) ∈ [0, 1]`. Then

       α(T x) = F̃(α(x)),

   and the continued-fraction partial quotients of `α(x)` are the run lengths of the parity vector of `x` (`a_1 − 1` leading even steps, then `a_2` odd steps, then `a_3` even steps, …). Consequently:
   * (a) the even anchor `0` (all-even word) and the odd anchor `−1` (all-odd word; the Mersenne run) are the parabolic fixed points `α = 0` and `α = 1` of `F̃`; a run of `r` odd steps is a passage near the cusp `α = 1` at distance `≈ 1/r`;
   * (b) the Syracuse map `U(n) = (3n+1)/2^(v)` is the jump transformation of `F̃` over the cusp at 0, i.e. the Gauss map `x ↦ 1/x − ⌊1/x⌋` on the odd-run coordinates; the valuations of `U` are every other partial quotient;
   * (c) the trivial cycle `{1, 2}` has `ρ = {2/3, 1/3}` and `α = {1/φ, 1/φ²}` (`φ` the golden ratio), the unique 2-cycle of `F̃` in `(0,1)`;
   * (d) `α(x)` is a quadratic irrational iff the Terras orbit of `x` is eventually periodic (Lagrange + THM-4634); the Collatz conjecture says that for every positive integer `n` the address `α(n)` lies in `Q(√5)` with continued fraction ending in `\overline{1}`.

3. **Minkowski distribution (PROVED).** For `n` uniform in `{1, …, N}`, `α(n)` converges in distribution, as `N → ∞`, to the singular Minkowski measure `d?`: the empirical distribution function of `α(n)` tends to `?(α)`. (Numerically, for `N = 20000` and 40 parity bits, the Kolmogorov distance to `?` is `0.0001`, against `0.14` for the uniform law.)

## Proofs

**1.** `θ(Tx) = 2θ(x)` is THM-4634. (a): rationality with odd denominator is eventual periodicity of the parity vector; purely periodic orbits give `0.\overline{w}`. The landing statement is the cited theorem. (b): the parity numbers are computed exactly in `collatz_angles.py` (`ρ(1) = 2/3`, `ρ(2) = 1/3`, `ρ(−5) = 6/7`, `ρ(−7) = 5/7`, `ρ(−10) = 3/7`, `ρ(0) = 0`, `ρ(−1) = 1`, `ρ(−17) = 1976/2047`); the identifications of `1/3, 2/3` (root of the `1/2`-bulb), `5/7, 6/7` (root of the `2/3`-bulb), `3/7, 4/7` (airplane) and `0` (cusp) are the classical table, which the adjacent-pair rule reproduces for all nine listed bulbs. (c): `ρ(c_a) = 0.\overline{10^(a−1)} = 2^(a−1)/(2^a−1)` and `T` doubles the angle; `ρ(−1/3) = 0.1000… = 1/2` since `T(−1/3) = 0`. (d): the root angles of the `p/q`-bulb are rotations of the Christoffel necklace (cited), and the orbit of a periodic angle is the set of rotations of its repetend. ∎

**2.** By the cited theorem (Panti; Bonanno–Isola) `?` conjugates the Farey map `F(x) = x/(1−x)`, `(1−x)/x` to the tent map `t ↦ min(2t, 2 − 2t)`. On `[1/2, 1]`, `F̃ = 1 − F` composed appropriately: `F̃(x) = 2 − 1/x = 1 − (1−x)/x`, so `?(F̃(x)) = 1 − ?((1−x)/x) = 1 − (2 − 2?(x)) = 2?(x) − 1 = 2?(x) mod 1`, using `?(1−y) = 1 − ?(y)`; on `[0, 1/2]`, `?(F̃(x)) = ?(F(x)) = 2?(x)`. Hence `? ∘ F̃ = D ∘ ?` with `D` the doubling map, and `α ∘ T = F̃ ∘ α` follows from `ρ ∘ T = D ∘ ρ`. The identity is also verified exactly for every Farey fraction with denominator `≤ 60` in the script. The run-length statement is the definition of `?`: the binary expansion of `?([0; a_1, a_2, …])` is `0.0^(a_1−1) 1^(a_2) 0^(a_3) …`. (a) `F̃'(0) = 1 = F̃'(1)`; a run of `r` ones is `ρ ≈ 1 − 2^(−r)`, and `?^(−1)(1 − 2^(−r)) = [0; 1, r, …] ≈ 1 − 1/r`. (b) The jump transformation of `F̃` over `[0, 1/2]` is the Gauss map (classical; the Farey map's induced map on the cusp). (c) `?(1/φ) = ?([0; 1, 1, 1, …]) = 0.101010… = 2/3` and `?(1/φ²) = 1/3`; `F̃(1/φ) = 2 − φ = 1/φ²` and `F̃(1/φ²) = 1/φ`. (d) Lagrange's theorem and THM-4634. ∎

**3.** Terras: the parity vectors of length `L` of `n ≤ N` are asymptotically uniform over `{0,1}^L` (each residue class mod `2^L` once), so `ρ(n)` is asymptotically uniform on `[0,1]`; `?^(−1)` is continuous, so the push-forward of Lebesgue measure under `?^(−1)` is the law of `α`, i.e. the measure with distribution function `?`. ∎

## Remarks

* **Mediants.** Concatenating two parity blocks adds their `(K, X)` vectors, so the slope of a concatenation is the Farey mediant of the slopes (THM-4515's fair splits are equal-slope decompositions). The bulbs of `M` are ordered by the same mediant law (Devaney). Collatz cycle shapes approximate `log_3 2 = [0; 1, 1, 1, 2, 2, 3, 1, 5, 2, 23, …]`, whose first convergents `1/2, 2/3, 5/8` coincide with the Fibonacci path of bulbs to the golden rotation number `1/φ² = [0; 2, 1, 1, …]` — hence the golden coincidences recorded across the repository (THM-4590, THM-4656): `log_3 2 = 0.6309` and `1/φ = 0.6180` share their first three convergents.
* **Codex's continued fraction** `β_n = [0; 1 + b_0, 1 + b_1, …]` (parity bits as digits 1, 2; one Collatz step = one Gauss shift) is a different faithful code; the Farey address replaces digits by run lengths and the Gauss map by the Farey map, which is why the anchors become parabolic points and the Syracuse map the jump transformation.
