---
id: THM-4632
title: "The Collatz parity-number dictionary: reading the parity vector Q(x) of a 2-adic integer x as the binary expansion of a real number rho(x) in [0,1] conjugates the Terras map T(x) = x/2, (3x+1)/2 to the doubling map; rho(x) is rational (irrationality exponent 1) iff the T-orbit of x is eventually periodic, in which case x is rational with the explicit value x_0 = sum_i 3^(k-1-i) 2^(d_i)/(2^P - 3^k) (Lagarias; the converse is the periodicity conjecture); the trivial cycle {1, 2} is the doubling 2-cycle {2/3, 1/3}, the unique doubling cycle of denominator 3, and the Collatz conjecture says exactly that every rho(n), n a positive integer, has lowest-terms denominator 3 times a power of 2, while a nontrivial cycle of period P would be a rho-value with odd denominator part dividing 2^P - 1 and not 3, and a divergent trajectory an irrational rho-value; the irrationality exponent mu(rho(x)) stratifies Z_2: the stratum mu = 1 is the countable set of eventually periodic points, Haar-almost every x has mu = 2, and for every tau >= 2 the set {mu(rho(x)) >= tau} has Hausdorff dimension 2/tau (Jarnik-Besicovitch transported by the 2-adic isometry Q), so the integers m >= 2 label the strata of 'constructed continued fractions' with partial quotients a_(n+1) of size q_n^(m-2), and the Collatz conjecture places all of N in the stratum mu = 1"
status: >
  PROVED, mostly by assembling KNOWN results: Terras (1976) and Lagarias (1985) for the parity-vector
  bijection modulo 2^k, the homeomorphism/isometry Q and the rationality of eventually periodic points;
  the cycle formula and the rational cycles from binary words (Lagarias 1985; Halbeisen-Hungerbuehler
  1997); Jarnik (1931)-Besicovitch (1934) for the dimension 2/tau and Khinchin for mu = 2 almost
  everywhere; the standard dimension preservation of the binary coding map. NEW here: the packaging
  as a dictionary 'integers as irrationality exponents' (statements 4-6), including the exact
  denominator form of the Collatz conjecture (3 times a power of two) and the transported strata.
  FINITE-EXACT checks (parity_number_checks.py): the Terras bijection for k <= 14; lowest-terms
  denominators 3 2^j for all n <= 2000; all 247 period words of length <= 7 realised as rational cycles
  with integer ones {-1}, {1, 2}, {-10, -7, -5} (and {0} for the empty word). Inspired by the owner's
  proposal (2026-10-14) to view the positive integers as irrationality measures with the accelerated
  cycle 1 -> 2 -> 1 as the only abnormality. Not independently audited yet.
session: mac-mini-2026-10-14-irrmeasure
source: 05-knowledge/results/collatz_irrationality_measures_20261014.md
scripts:
  - 04-computation/experiments/irrmeasure_20261014/parity_number_checks.py (+ .out)
related:
  - THM-4633 (cycles as rational approximations of log2 3; the admissible-minimum budget)
  - THM-4591 (the five integer T-cycles to period 301,993), THM-4623 (order dictionary, Koopman torsion), THM-4523 (rigidity)
  - HYP-9260 / HYP-9261 (opus S29/S30: Moebius dictionary of parity words, Farey census of log2 3)
  - HYP-9262 (the owner's mu = 2 hypothesis for logarithms and odd zeta values)
  - J. C. Lagarias, The 3x+1 problem and its generalizations, Amer. Math. Monthly 92 (1985) 3-23 (parity vectors, Theorem B, the periodicity conjecture)
  - R. Terras, A stopping time problem on the positive integers, Acta Arith. 30 (1976) 241-252
  - D. J. Bernstein, A non-iterative 2-adic statement of the 3N+1 conjecture, Proc. AMS 121 (1994) 405-408; D. J. Bernstein, J. C. Lagarias, The 3x+1 conjugacy map, Canad. J. Math. 48 (1996) 1154-1169
  - L. Halbeisen, N. Hungerbuehler, Optimal bounds for the length of rational Collatz cycles, Acta Arith. 78 (1997) 227-239
  - V. Jarnik, Diophantische Approximationen und Hausdorffsches Mass, Mat. Sb. 36 (1929); A. S. Besicovitch, Sets of fractional dimensions IV, J. London Math. Soc. 9 (1934)
---

# THM-4632 — the parity-number dictionary

## Setting

`T(x) = x/2` for even `x` and `(3x + 1)/2` for odd `x`, on the 2-adic integers `Z_2` (the "accelerated" or
Terras map; `1 → 2 → 1`). The *parity vector* of `x` is `Q(x) = (v_0, v_1, ...)`, `v_i = T^i(x) mod 2`,
identified with the 2-adic integer `Σ v_i 2^i`. The *parity number* is the real

    ρ(x) = Σ_(i ≥ 0) v_i 2^(-i-1) ∈ [0, 1],

the parity vector read as a binary expansion. `μ(ξ)` is the irrationality exponent of a real `ξ`
(`μ = 1` for rationals, `μ ≥ 2` for irrationals by Dirichlet, `μ(ξ) = 2 + limsup_n log a_(n+1)/log q_n` in
terms of the continued fraction). `D_3 := {a/(3·2^j) in lowest terms : j ≥ 0, 0 < a < 3·2^j}` is the set of
reals in `(0, 1)` whose binary expansion ends in `(10)^∞` or `(01)^∞`.

## Statement

1. **Isometry (KNOWN).** For every `k`, `x mod 2^k ↦ (v_0, ..., v_(k-1))` is a bijection `Z/2^k → {0,1}^k`
   (Terras). Hence `Q : Z_2 → Z_2` is a 2-adic isometry (Lagarias: a homeomorphism), and `ρ = r ∘ Q` with
   `r` the binary reading is 1-Lipschitz from the 2-adic metric to `R`, onto `[0, 1]`, injective off the
   countable set of `x` whose orbit reaches the fixed point `0` or `-1` (eventually constant parity).
2. **Conjugacy.** `ρ(T x) = {2 ρ(x)}` (fractional part): `ρ` semi-conjugates `T` to the doubling map;
   the parity vector is the itinerary.
3. **Rationality (KNOWN: Lagarias).** For `x ∈ Z_2` the following are equivalent: `ρ(x) ∈ Q`; the parity
   vector is eventually periodic; the `T`-orbit of `x` is eventually periodic. They imply `x ∈ Q ∩ Z_2`
   (a rational with odd denominator): if the orbit of `x` enters, at `x_0 = T^L x`, a cycle of period `P`
   whose period word has `k ≥ 1` ones at positions `d_0 < ... < d_(k-1)` (counted from `x_0`), then

       x_0 = ( Σ_(i<k) 3^(k-1-i) 2^(d_i) ) / (2^P - 3^k),

   and every binary word of length `P` with `2^P ≠ 3^k` is the period word of exactly one `T`-cycle in
   `Q ∩ Z_2`. The converse implication "`x ∈ Q ∩ Z_2 ⇒ ρ(x) ∈ Q`" is Lagarias's periodicity conjecture
   (OPEN; it contains the Collatz conjecture).
4. **The trivial cycle and the Collatz conjecture.** `ρ(1) = 2/3`, `ρ(2) = 1/3`: the cycle `{1, 2}` of `T`
   is the 2-cycle `{2/3, 1/3}` of doubling, the unique doubling cycle with denominator 3. For `n ∈ N`:

       the T-orbit of n reaches {1, 2}  ⟺  ρ(n) ∈ D_3  ⟺  the doubling orbit of ρ(n) hits {1/3, 2/3}.

   Hence **the Collatz conjecture is the statement that the lowest-terms denominator of ρ(n) is
   `3·2^j` for every positive integer `n`.** A nontrivial cycle of `T` on `N` of period `P` would be an
   `n` with `ρ(n) = a/(2^P - 1)`, odd denominator part dividing `2^P - 1` and different from `3`; a
   divergent trajectory would be an `n` with `ρ(n)` irrational.
5. **Strata.** Put `S_τ = {x ∈ Z_2 : μ(ρ(x)) ≥ τ}`. Then `S_1 ∖ S_2 = {x : ρ(x) ∈ Q}` is the countable set
   of eventually periodic points (contained in `Q ∩ Z_2`); Haar-almost every `x` has `μ(ρ(x)) = 2`; and
   for every `τ ≥ 2`, `dim_H S_τ = 2/τ` (2-adic metric), with `S_∞` (Liouville parity numbers) of
   dimension `0` but dense. In particular each integer `m ≥ 2` labels the stratum `{μ(ρ(x)) = m}`, of
   dimension `2/m`, carried by the parity numbers whose continued fractions have partial quotients
   `a_(n+1) ≍ q_n^(m-2)` along a subsequence — the "constructed continued fractions" of the proposal —
   and the Collatz conjecture says `N ⊆ S_1 ∖ S_2` (with the denominator form of 4), while Lagarias's
   periodicity conjecture says `Q ∩ Z_2 ⊆ S_1 ∖ S_2`.
6. **Small cycles (FINITE-EXACT).** Among the `247` period words of length `≤ 7` the integer cycles are
   `{-1}` (word `1`), `{1, 2}` (word `10`), `{-10, -7, -5}` (word `110`), and `{0}` (empty word); all others are
   rational non-integral cycles such as `1/5 → 4/5 → 2/5 → 1/5` (word `100`, `ρ = 1/7`).

## Proofs

**(1).** The first `k` parities of `x` depend only on `x mod 2^k` (each step is affine with the parity
decided by the current residue), and the map `Z/2^k → {0,1}^k` is injective because `T` restricted to a
residue class mod `2^k` with prescribed parity is injective onto a class mod `2^(k-1)` (both branches
are bijections `2Z_2 → Z_2` and `1 + 2Z_2 → Z_2`); counting gives bijectivity (checked for `k ≤ 14`).
So `|Q(x) - Q(y)|_2 = |x - y|_2`. The binary reading satisfies `|r(u) - r(v)| ≤ 2^(-k)` when `u, v`
agree in the first `k` digits, and identifies exactly the pairs `u = w0111..., v = w1000...`, i.e. the
parity vectors eventually constant `1` (orbit reaching `-1`, since `T(-1) = -1` and an eventually
constant-1 parity vector means the tail orbit is a fixed point with parity 1, which is `-1`) or
eventually constant `0` (tail orbit the fixed point `0`); those `x` form a countable set.

**(2)** is the definition: the parity vector of `Tx` is the shifted parity vector of `x`.

**(3).** `ρ(x) ∈ Q` iff the binary expansion is eventually periodic iff `Q(x)` is eventually periodic.
If the orbit is eventually periodic so is the parity vector; conversely an eventually periodic parity
vector with preperiod `L` and period `P` gives `Q(T^(L+P) x) = Q(T^L x)`, hence `T^(L+P) x = T^L x` by
injectivity of `Q`. For the formula, iterate `2^(e) y_(i+1) = 3 y_i + 1` along the period word (with
`2^P` total halvings and `k` odd steps) to get `2^P x_0 = 3^k x_0 + Σ_i 3^(k-1-i) 2^(d_i)`; the
denominator `2^P - 3^k` is odd, so `x_0 ∈ Q ∩ Z_2`, and `x = T^(-L)`-preimage data are rational as well
(the branches of `T` are rational affine maps). Conversely, given a word, the displayed `x_0` has the
word as its period word (direct substitution; all 247 words of length `≤ 7` checked exactly).

**(4).** `Q(1) = (1, 0, 1, 0, ...)` and `Q(2) = (0, 1, 0, 1, ...)`, so `ρ(1) = 2/3`, `ρ(2) = 1/3`; a
doubling cycle `{a/(2^P - 1)}` has denominator `3` only if the period word is `10` or `01` repeated.
The orbit of `n` reaches `{1, 2}` iff from some step on the parities alternate, iff
`ρ(n) = (a + 2/3)/2^j` or `(a + 1/3)/2^j = (3a + 2)/(3·2^j)` or `(3a + 1)/(3·2^j)`, which are in lowest
terms (the numerator is prime to 3; dividing out 2's only lowers `j`). A cycle of period `P` has
`ρ = a/(2^P - 1)` with `a` the period word read in binary; divergence means a non-eventually-periodic
parity vector, i.e. an irrational `ρ(n)`.

**(5).** `S_1 ∖ S_2 = ρ^(-1)(Q)` is countable by (3) (eventually periodic points are determined by a
word and a preperiod). For `τ ≥ 2`, `{ξ ∈ [0,1] : μ(ξ) ≥ τ}` has Lebesgue measure zero for `τ > 2` and
Hausdorff dimension `2/τ` (Jarník–Besicovitch), and `μ = 2` almost everywhere (Khinchin). The binary
coding map preserves Hausdorff dimension between `{0,1}^N` (metric `2^(-first disagreement)`) and
`[0, 1]` (dyadic intervals of length `2^(-k)` are the cylinders; net measures), and `Q` is an isometry,
so the same holds for `S_τ ⊆ Z_2`; Haar measure on `Z_2` corresponds to the fair-coin measure on parity
vectors (Terras) and to Lebesgue measure on `[0,1]`. The description of the stratum `μ = m` through
partial quotients is the formula `μ = 2 + limsup log a_(n+1)/log q_n`; the Liouville set is a dense
`G_δ` of dimension zero, and `S_τ` is a `G_δ` (a limsup of open conditions), so all strata are low Borel.

**(6)** is the exact enumeration in `parity_number_checks.py` (the `x_0` of (3) for every word of length
`≤ 7`, with the period word recovered by iterating `T` on rationals with odd denominators).  ∎

## Reading

* **The owner's dictionary made exact.** "Integers as irrationality measures": the positive integers sit
  inside `Z_2`, whose points are stratified by the irrationality exponent of their parity number;
  "`1 → 2 → 1` is the only abnormality in the connected tree" is the statement that the only odd
  denominator a parity number of a positive integer ever has is `3`, i.e. `N ⊆ ρ^(-1)(D_3)`; "the integers
  above 2 correspond to constructed continued fractions" is the Jarník stratification, where the level
  `m` is realised by parity vectors with partial quotients of size `q_n^(m-2)`, a set of 2-adic integers
  of dimension `2/m` that contains no eventually periodic point and, conjecturally, no positive integer.
* **Dirac and Fisher information.** A rational parity number is a deterministic, eventually periodic
  itinerary (zero entropy, "known exactly"); `μ = 2` is the law of large numbers of the Gauss map
  (`log a_(n+1)/log q_n → 0` almost surely); the strata `μ = m > 2` are large-deviation sets. For a
  divergent `n` the parity number is irrational with binary digit density `liminf (odd steps)/(steps)
  ≥ log 2/log 3` (from `n_N = n (3^k/2^N) Π(1 + 1/(3n_i)) ≥ 1`), a non-normal number whose digit law is at
  Kullback–Leibler distance at least `1 - H(log_3 2) ≈ 0.049` bits per step from the fair coin: divergence
  costs information at every step. No finite irrationality-exponent bound for such `ρ(n)` follows from
  this (good rational approximations need not come from the binary expansion).
* **What the dictionary does not do.** It is a conjugacy, not a new invariant: every statement about
  `ρ(n)` is a statement about the parity vector. Its value is to name the three failure modes of the
  conjecture by Diophantine type (denominator `3·2^j`; odd denominator not `3`; irrational) and to make
  THM-4633's "budget" reading of cycles available.
