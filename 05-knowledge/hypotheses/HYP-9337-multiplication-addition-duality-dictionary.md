---
id: HYP-9337
title: "The multiplication/addition duality dictionary: (a) the golden ratio is the self-dual point (phi^2 = phi + 1; Zeckendorf; THM-4528's base-phi words; THM-4714's golden Farey pair; THM-4717's golden circle, which exists because log_2 3 < phi); (b) exponentiation carries sums to products and the Collatz carry B_w is the additive defect of the multiplicative word 3^p/2^q, the clock 2^q - 3^p the multiplicative defect of the additive shape, cycles existing iff the one divides the other (THM-4484), and balanced cycles having multiplicative (3-smooth) differences (THM-4716/4717); (c) modular arithmetic: the Fermat quotient is a homomorphism from products to sums, the Wilson quotient their sum (Lerch, THM-4641), the Pisano period is the multiplicative order of phi (THM-486/491, Wall-Sun-Sun = golden Wieferich THM-4656), the Mersenne word period is the 2-adic order of 3 (THM-4556); (d) Gilbreath's triangle is the additive derivative of the multiplicative atoms and its sea is Pascal mod 2, the carry-free shadow of 3^E = (1+2)^E; (e) Gersonides/Catalan are the only additive coincidences of the multiplicative atoms 2, 3 (THM-4484, THM-4640); conjectures: the critical slope log_3 2 is the unique irrational slope without a Sturmian Collatz circle, and every integer 3x+1 cycle word is unbalanced"
status: >
  OPEN as a dictionary; every row typed. (a), (b), (c), (e): KNOWN or PROVED in the cited theorems. (d): PROVED elementary
  (the XOR sea), the Gilbreath depth law is HYP-9336. The two conjectures at the end are consequences of Collatz (THM-4716 proves
  the balanced half: a balanced word is a 3x+1 cycle only for |2^q - 3^p| = 1) and of THM-4717 (the critical slope is where the
  series diverge; whether some other mechanism produces a circle there is open).
session: mac-mini-2026-10-19-duality
source: 05-knowledge/results/duality_20261019_multiplication_addition_golden_circle_gilbreath.md
related:
  - THM-4717, THM-4716, THM-4714, THM-4528, THM-4484, THM-4640, THM-4641, THM-4651, THM-4556, THM-4486 and THM-4491 (Pisano), THM-4656, THM-4511, HYP-9336
---

# HYP-9337 — the multiplication/addition duality dictionary

| theme | additive side | multiplicative side | bridge (repo object) | status |
|---|---|---|---|---|
| golden ratio | `φ + 1`, Zeckendorf sums, Beatty/Wythoff `⌊nφ⌋` | `φ²`, `φ^n = F_n φ + F_(n−1)` | THM-4528 (base-φ words), THM-4714 (`?^(−1)(2/3) = 1/φ`), THM-4717 (golden circle: `c_min = Σ 3^r 2^(−⌈(r+1)φ⌉)`) | PROVED |
| exponentiation | shape `(q, p)`, carry `B_w = Σ 3^j 2^i` | `3^p/2^q`, clock `2^q − 3^p` | THM-4484 (`D | d·B`), THM-4716/4717 (differences of balanced cycles are `2^i 3^j`) | PROVED |
| the critical slope | density of odd steps `a` | `3^a` vs `2` | THM-4633 (integer cycles have `a = log_3 2`), THM-4717 (circles exist iff `a ≠ log_3 2`) | PROVED / conjecture |
| modular arithmetic | Fermat quotient `q_p(ab) = q_p(a) + q_p(b)`, Wilson quotient `= Σ q_p` (Lerch) | `(Z/p²)^×` | THM-4641, THM-4651 | KNOWN |
| recurrences mod `p` | Pisano period, Fibonacci entry point | order of `φ` in `F_(p²)^×`, Wall–Sun–Sun | THM-486/491, THM-4656 | KNOWN |
| 2-adic words | exponent `E` of `2^E − 1` | `3^(E−1) mod 2^K` | THM-4556 (word period `2^(K−2)`) | PROVED |
| Gilbreath | iterated differences of primes | the primes as multiplicative atoms | S9 note / THM-4511 (XOR sea = Pascal mod 2 = carry-free `(1+2)^E`), HYP-9336 (depth ≈ max gap) | PROVED / OPEN |
| coincidences | `2^a − 3^b = ±1` | `2^a, 3^b` | Gersonides (THM-4484), Catalan carry (THM-4640) | KNOWN |
| mediant | `(a + c)/(b + d)` | Ford circles `1/(2q²)`, Ramanujan sums | THM-4714/4715, HYP-9335 (b) | PROVED / REFUTED row |

## Conjectures

1. **Unbalanced cycles.** Every integer cycle of `3x+1` other than the four Gersonides cycles has an unbalanced parity word (consequence of Collatz; THM-4716 proves that a balanced word is an integer cycle only for `|2^q − 3^p| = 1`).
2. **No circle at the critical slope.** The balanced rational cycles of slopes converging to `log_3 2` have unbounded minima (THM-4717's data: `5, 24.5, −219, 867, −6143`), and no Sturmian circle of the Collatz map exists at `a = log_3 2`.
3. **Gilbreath depth** is `(0.8–1.05)` times the maximal prime gap (HYP-9336).
