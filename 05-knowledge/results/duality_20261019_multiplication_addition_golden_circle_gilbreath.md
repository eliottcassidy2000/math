# The duality between multiplication and addition: Sturmian circles of the Collatz map, the golden circle and its constant, the Gilbreath depth law, and the dictionary

Session mac-mini-2026-10-19-duality (2026-10-19). Owner prompt: merge the repository's work on the duality between
multiplication and addition — the golden ratio, exponentiation, modular arithmetic, dynamical sequences, Gilbreath's
conjecture — and intertwine them in an investigation of structure. **Status: THM-4717 PROVED (pre-audit); HYP-9336 OPEN
(numerical law, new data to p < 2·10^8); HYP-9337 dictionary (typed rows). No conjecture is settled.**

## 0. Inheritance

* **Golden:** THM-4528 (parity words of `n/2, 3n+1` are golden-mean words; base-φ reading; golden β-map), THM-4714 (the
  Farey address of the trivial cycle is `1/φ`; the Fibonacci bulbs `1/2, 2/3, 5/8` are the first convergents of `log_3 2`),
  THM-4590 (golden resonance), THM-4656 (Wall–Sun–Sun = golden Wieferich), the Pisano theorems THM-486/491.
* **Gilbreath:** the S9/S10 opus sessions of 2026-09-26 (the 0/2 sea is the XOR automaton, fronts move left at speed 1,
  the first all-0/2 row is row 65 for `p < 2·10^5`; THM-4511 wall theorem; HYP-9163 extinction law). This session only
  adds the depth law at larger scale (HYP-9336).
* **Exponentiation / carries:** THM-4484 (free cycles iff the clock divides `d`), THM-4633 (integer cycles sit at the
  convergents of `log_3 2`), THM-4716 (yesterday: balanced cycles are 3-smooth ladders and free).
* **Modular:** THM-4641 (Fermat/Wilson quotients, Lerch), THM-4651 (Giuga at `p²`), THM-4556 (Mersenne words periodic in
  the exponent on the 2-adic order of 3).
* **Board.** Anchor: the ladder formula and the Sturmian circles (THM-4717). Niche: the Gilbreath depth law (HYP-9336).
  Wildcard: the dictionary (HYP-9337), including the critical-slope reading of why integer cycles are hard.

## 1. The explicit ladder and the Sturmian circles (THM-4717)

THM-4716 proved that consecutive intercepts of a Christoffel word differ by one adjacent swap, so the sorted elements of the
balanced free cycle of `3x + (2^q − 3^p)` differ by single 3-smooth numbers. Computing the moved position gives the closed
form: the differences are exactly `d_i = 2^i 3^(p−1−⌊(i+1)p/q⌋)`, `i = 0..q−2` (each power of two once; the 3-exponents a
reversed Beatty sequence), verified for all 1101 coprime shapes with `q ≤ 60`. Summing the carry formula and letting the
shape run through the convergents of an irrational slope `a` gives the limit circles:

* `a < log_3 2`: `x_a(θ) = ½ Σ_(r≥0) 3^r 2^(−⌊r/a + θ⌋)` (positive);
* `a > log_3 2`: `x_a(θ) = −Σ_(m≥1) 2^(⌈m/a − θ⌉ − 1) 3^(−m)` (negative);

with `θ ∈ R/(1/a)Z`, on which the Sturmian-coded Collatz step is the rotation `θ ↦ θ + 1` (rotation number `a`): an exact
telescoping identity, checked at `√2 − 1`, `1/φ`, `π − 2.5`, `1/√2` along their convergents and at random points. The
circle is the image of a strictly monotone pure-jump function — a measure-zero Cantor set, the Collatz counterpart of a
Denjoy exceptional minimal set. At the critical slope `a = log_3 2`, which is the slope of every integer `3x+1` cycle
(THM-4633), both geometric ratios equal 1 and the balanced rational cycles' minima blow up (`5, 24.5, −219, 867, −6143`).

**The golden circle.** For `a = 1/φ`: `c_min = Σ 3^r 2^(−⌈(r+1)φ⌉) = 10.370127141747049161823234577141…`,
`c_max = ½ Σ 3^r 2^(−⌊rφ⌋) = 31.610381425241147485469703731423…`, `c_max/c_min = 3.04821541656776…`; the Fibonacci-shape
cycles `(F_(n+1), F_n)` converge alternately (`|x_min − c_min| = 2.5·10^(−11)` at `(610, 377)`). It exists because
`log_2 3 = 1.585 < φ = 1.618`: the golden slope is subcritical by `0.0129`, the term ratio is `3/2^φ = 0.9773`. Two new
integer sequences fall out (neither in OEIS on 2026-10-19): the Fibonacci clocks `2^(F_(n+1)) − 3^(F_n) = 1, −1, 5, 13,
1631, 502829, 6719515981, 19351615319297399, …` and the golden Collatz minima `1, 5, 23, 319, 14501, 5805215, 68498816767,
204377029899901229, …` (the minimal elements of the balanced free cycles of `3x + clock`).

In THM-4528's coding (`1 → 10`) the golden Terras word becomes the Fibonacci word of slope `1/φ²`, so the golden circle is the
real avatar of THM-4528's golden-mean words: the same object seen through the base-φ reading `Θ` (2-adic/β-shift side) and
through the real series `x_a(θ)` (real side).

## 2. Gilbreath depth (HYP-9336)

For the first `N` primes, the depth `k(N)` is the first row that is `1` followed by `0/2` only on its first `N − k` entries;
the leading 1 then persists (XOR sea). New data (`gilbreath_depth.py`, 11,078,937 primes below `2·10^8`, 239 rows, 8 s):
`k = 35, 65, 97, 147, 175, 239` for `N = 10^3, …, 10^7, 1.1·10^7`, against the maximal prime gaps `G = 34, 72, 114, 154,
222, 248` below `p_N`: `k/G ∈ [0.79, 1.03]`. Odlyzko's 635 rows for `p < 10^13` (cited from memory) against `G = 674`
(T5K maximal-gap table: 674 first after 7,177,162,611,713) gives `0.94`. An isolated large gap dies in three rows, so the
long-lived defects are clusters; no mechanism is proved. The S9 note's `k = 65` for `p < 2·10^5` is reproduced.

## 3. The dictionary (HYP-9337)

Golden self-duality (`φ² = φ + 1`; Zeckendorf; the golden Farey pair; the golden circle), exponentiation (carry = additive
defect of the multiplicative word; clock = multiplicative defect of the additive shape; cycles iff one divides the other;
balanced cycles have multiplicative differences), modular arithmetic (Fermat quotient as a homomorphism from products to
sums; Wilson = Σ Fermat quotients; Pisano period = order of φ; Wall–Sun–Sun = golden Wieferich; Mersenne word period = 2-adic
order of 3), Gilbreath (the additive derivative of the multiplicative atoms; its sea is Pascal mod 2, the carry-free shadow
of `3^E = (1+2)^E`), Gersonides/Catalan (the only additive coincidences of `2^a, 3^b`), and the mediant (yesterday's HYP-9335,
with the refuted Stern–Brocot truncation of the singular series). Two conjectures: every nontrivial integer `3x+1` cycle is
unbalanced (half proved by THM-4716) and no Sturmian circle exists at the critical slope.

## 4. Reproduction

    python3 ladder_formula.py 60; python3 golden_limit.py; python3 sturmian_circles.py; python3 critical_slope.py
    python3 golden_circle.py          # 45-digit constants, conjugacy at 300 points, monotone/pure-jump checks
    python3 gilbreath_depth.py 200000000

All in `04-computation/experiments/duality_20261019/`.
