# Farey fractions, the mediant and the Mandelbrot set across five conjectures: Collatz parity numbers as external angles, Minkowski–Farey addresses, balanced cycles as cardioid bulbs, and what the mediant does not do for primes

Session mac-mini-2026-10-19-farey (2026-10-19). Owner prompt: continue Collatz, Hadwiger–Nelson, the planar Jacobian
conjecture, Goldbach and twin primes; merge the repository's Farey-fraction, mediant and Mandelbrot themes and any others.
**Status: THM-4714 PROVED (elementary on top of THM-4634 and cited classical theorems; pre-audit); THM-4715 PROVED +
FINITE-EXACT; THM-4716 PROVED (balanced cycles are 3-smooth ladders and are free; HYP-9334 promoted); HYP-9335 OPEN dictionary
with one row REFUTED. None of the five
conjectures is advanced beyond reformulation; every row is typed.**

## 0. What the repository already had, and what codex did today

* Farey: 180 canon files, almost all in the lonely-runner program (THM-998: the deep set is a Farey major-arc dissection;
  THM-501/503: the LRC singular series; THM-536: Sturmian partial-sum walks). Stern–Brocot: THM-4057/4059 (codex, 2026-08-24:
  rational-edge tournament gauge, Berggren/Calkin–Wilf). Christoffel words: THM-778/808/812/813. Collatz continued
  fractions: codex's `collatz_parity_cf_boundary_20261014a` (parity bits as partial quotients 1, 2; one step = one Gauss
  shift; exponent 2) and THM-4633/4634 (convergents of `log_2 3`; the parity number `ρ`). Minkowski's `?`: no mention.
* Mandelbrot: only codex's parabolic-scaling session of today (`π/√ε` escape clocks at the cusp `c = 1/4`, guarded
  integrals for twin-prime and Goldbach counts, the airplane residue `51/49` at `c = −7/4`, cotangent lifts for the
  Jacobian side). Goldbach and twin primes: old S5xx-era scripts and the twin-centre weld theorems (THM-2413..2447); thin.
* **Board.** Anchor: Collatz through the Farey address and external angles (THM-4714, THM-4715, HYP-9334). Niche: the
  cross-conjecture dictionary with a refuted prime row (HYP-9335). Wildcard: the golden 2-cycle and the cusp/tip spine.

## 1. The dictionary (THM-4714)

Read THM-4634's parity number `ρ(x)` as an external angle of the Mandelbrot set. Because `ρ ∘ T = 2ρ mod 1`, a
Collatz cycle is a periodic angle, and Douady–Hubbard place it at the root of a hyperbolic component of the same period.
Exact values (`collatz_angles.py`, 1671 checks):

| Collatz object | parity word | angle(s) | Mandelbrot landing point |
|---|---|---|---|
| fixed point 0 | `0` | 0 | cusp `c = 1/4` |
| fixed point −1 (the Mersenne run) | `1` | 1 ≡ 0 | cusp `c = 1/4` |
| trivial cycle {1, 2} | `10` | 2/3, 1/3 | root `c = −3/4` of the 1/2-bulb |
| {−5, −7, −10} | `110` | 6/7, 5/7; 3/7 | root of the 2/3-bulb; airplane root `c = −7/4` |
| −17 cycle (period 11, 7 odd) | `11110111000` | 1976/2047 and rotations | period-11 roots off the cardioid (unbalanced) |
| anchor `c_a = 1/(2^a − 3)` | `1 0^(a−1)` | `{2^(a−1), 1, 2, …, 2^(a−2)}/(2^a − 1)` | the 1/a-bulb (root angles `1/(2^a−1), 2/(2^a−1)`) |
| limit anchor −1/3 | `1 0 0 0 …` | 1/2 | tip `c = −2` |

So the real spine `[−2, 1/4]` of the Mandelbrot set runs from the Collatz point `−1/3` (whose image is 0) to the two
anchors `0, −1`, and the trivial cycle is the first bulb on it. A cycle word is balanced iff its angle orbit is a
cardioid-bulb root orbit; the bulb-root angles are the two numerically adjacent rotations of the Christoffel necklace,
a rule checked against the classical table for nine bulbs and unique for all `q ≤ 15`.

**Minkowski.** `α(x) = ?^(−1)(ρ(x))` turns the parity vector into a continued fraction of run lengths and the Terras map
into the Farey-type map `x/(1−x), 2 − 1/x` (exact conjugacy on every Farey fraction of denominator ≤ 60; general proof
from Panti's tent conjugacy and `?(1−x) = 1 − ?(x)`). The two anchors become the two parabolic fixed points (a run of
`r` odd steps is a cusp passage at distance `1/r`, the Collatz counterpart of codex's `π/√ε` cusp clocks), the Syracuse
map becomes the Gauss jump transformation, the trivial cycle becomes the golden pair `{1/φ, 1/φ²}`, and the integers'
addresses follow the singular Minkowski measure (Kolmogorov distance 0.0001 at `N = 20000`). Collatz reads: every
positive integer's address is a quadratic irrational of `Q(√5)` whose Farey orbit ends on the golden pair.

## 2. Free cycles are all the angle orbits; balanced cycles are free (THM-4715, HYP-9334)

For `D = 2^q − 3^p` the free cycles of `3x + D` of shape `(q, p)` are exactly the necklaces of length `q` with `p` ones
(THM-4484), i.e. all doubling orbits of period-`q` angles with `p` ones; the Christoffel necklace is the cardioid
`p/q`-bulb, the others are off-cardioid components. `3x + 13` realises all seven period-8 orbits with five ones (bulb cycle
`319 → 485 → 734 → 367 → 557 → 842 → 421 → 638`, root angles `181/255, 182/255`); `3x + 5` both period-5 orbits with
three ones. Census (`d ≤ 49` odd, `|n| ≤ 20000`): 219 cycles, 120 balanced, every balanced one free. Sporadic cycles with
large minima sit at Stern–Brocot nodes of `log_3 2`: the `−17` cycle at the semiconvergent `7/11` (`D = −139`), the
27-cycles of `3x + 5` and `3x + 15` at `17/27` (`D = 5077565`).

**THM-4716 (HYP-9334 proved).** The `q` rotations of the Christoffel necklace are the lower mechanical words with
intercepts `t/q`; consecutive intercepts differ by exactly one adjacent transposition `01 → 10`, so consecutive carries
differ by exactly `2^i 3^j` (`i` = new position of the moved one, `j` = ones to its right). Hence the sorted elements of
the balanced free cycle of `3x + (2^q − 3^p)` form a **3-smooth ladder** (`3x + 13`: `319, 367, 421, 485, 557, 638, 734,
842`, differences `48, 54, 64, 72, 81, 96, 108`, the consecutive 3-smooth numbers of that window; in general the
differences are the 3-smooth numbers of the moved ones), and since a prime of the clock dividing one element divides
all (the dynamics is linear modulo it), it would divide a 3-smooth difference: impossible. So `gcd(2^q − 3^p, B_r) = 1`
for every rotation (checked to `q = 150`, 686,732 rotations, while random necklaces fail in 6% of trials), no sporadic
cycle of any `3x + d` is balanced, the balanced `3x+1` cycles are Gersonides' four, and a positive nontrivial `3x+1`
cycle would land off the main cardioid.

## 3. The other four conjectures (HYP-9335)

* **Goldbach and twin primes.** The circle method is a Farey dissection and the singular series is a sum over Farey
  fractions of Ramanujan sums. The natural mediant guess — truncate the series by Stern–Brocot depth — is REFUTED as a
  convergence statement: for the twin-prime series the depth-16 partial sum is 0.885 against `2C_2 = 1.3203`, while
  denominators `≤ 30` already give 1.313; for `N = 210` depth 16 gives 4.205 against 4.225 (`singular_series_farey.out`).
  The cancellation of the circle method lives inside each denominator (the Ramanujan sum over all `a/q`), and the mediant
  tree cuts across it. The Farey structure relevant to primes is the conductor filtration `q`, not the mediant order.
  Codex's guarded parabolic integrals (same day) are the Mandelbrot-side formulation; they do not estimate the prime
  predicates either.
* **Hadwiger–Nelson.** The rational unit vectors are Farey-parametrised (THM-4057, THM-3333), the Farey graph itself is
  3-coloured by the parity class of `(a, b)` (the three cusps of `Γ(2)`), and the rational unit-distance graph is bipartite
  (Woodall 1973): the Farey skeleton contributes nothing to `χ(R²) ≥ 5`; the six-colour field work (THM-4700/4701) lives
  on irrational directions. No new HN statement.
* **Planar Jacobian.** `Aut(C²) = A *_(A∩E) E` (Jung–van der Kulk) has a Bass–Serre tree, the analogue of the Farey tree
  (which is the tree of `PSL_2(Z) = Z_2 * Z_3`, dual to the Farey tessellation); HYP-9070's "Euclidean depth" of the
  degree pair is the Stern–Brocot depth of the degree ratio. Codex's parabolic cotangent germs and the Collatz odd anchor
  are the same cusp (THM-4714 (1b)). No new JC statement.

## 4. Themes pulled in

Gersonides and the Catalan carries (THM-4484, THM-4640: the `|D| = 1` shapes are the cusp and the first two bulbs);
fair splits as mediants (THM-4515); the golden Wieferich and golden Collatz resonance (THM-4656, THM-4590: the first
three convergents of `log_3 2` are the Fibonacci bulbs `1/2, 2/3, 5/8`); the anchors `c_a` (THM-4600) as the `1/a`-bulb
sequence; the Mersenne run (THM-4556) as the cusp passage; codex's parabolic clocks; the LRC Farey dissection (THM-998)
as the prototype of row (b).

## 5. Reproduction

    python3 collatz_angles.py 49 20000        # 1671 exact checks; cycle census; Stern-Brocot path; Minkowski distribution
    python3 balanced_free.py 150              # HYP-9334: 686,732 rotations, all coprime
    python3 singular_series_farey.py          # the refuted mediant truncation

All in `04-computation/experiments/farey_mandelbrot_20261019/`.

Sources used for the cited facts: Devaney, *The Mandelbrot set, the Farey tree, and the Fibonacci sequence*, Amer. Math.
Monthly 106 (1999) 289–302; Douady–Hubbard (Orsay notes) and Milnor, *Periodic orbits, external rays and the Mandelbrot
set* (Astérisque 261, 2000); Panti, *Multidimensional continued fractions and a Minkowski function*, Monatsh. Math. 154
(2008) (arXiv:0705.0584: `?` is the unique homeomorphism conjugating the Farey map to the tent map); Bonanno–Isola,
*Orderings of the rationals and dynamical systems*, Colloq. Math. 116 (2009); Woodall, *Distances realized by sets
covering the plane*, J. Combin. Theory A 14 (1973).
