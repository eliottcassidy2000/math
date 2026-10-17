---
id: THM-4666
title: "Class numbers are ranking excesses: the natural order of the Kronecker oriented graph of a negative fundamental discriminant D has forward-minus-backward excess (2/w_D)|D|h(D), the area under the Legendre walk; for the Paley tournament P_p a single-vertex insertion changes the excess by -2S(t), so the natural order is insertion-locally optimal exactly for Kalmynin's primes L+ (all partial sums of (./p) nonnegative), which by Fekete's identity have no real zero of L(s, chi_(-p)); h(-p) = 1 (the least excess) exactly for p = 7, 11, 19, 43, 67, 163, i.e. Euler's lucky numbers"
status: >
  PROVED (elementary): statements 1, 3, 5, 6 (the tournament readings) from Dirichlet's class number formula and quadratic
  reciprocity. CITED: Dirichlet's half-window formula (statement 2), Heegner-Baker-Stark and Rabinowitsch (statement 4),
  Fekete-Polya (statement 5, proof included), Baker-Montgomery and Kalmynin (density of L+), Alon (order of the spectral
  bound). FINITE-EXACT: statement 7. No novelty claim for the identities themselves (they rewrite Dirichlet's formula);
  the insertion/walk dictionary and the reading of L+ as local optimality of a ranking were not found in a quick search.
session: opus-2026-10-17-S40
source: 05-knowledge/results/binary_information_class_groups_largest_anomalies_20261017.md
scripts:
  - 04-computation/experiments/paley_class_number_ranking_20261017.py (+ .out, .json)
  - 04-computation/experiments/paley_legendre_walk_positivity_20261017.py (+ .out)
  - 04-computation/experiments/paley_mas_exact_20261017.py (+ p31 .out, .json)
related:
  - HYP-9287 (pre-registered: optimality of the natural order is governed by L(1, chi_(-p)))
  - THM-871 (rotational tournaments at the Fermat rungs; all rotational tournaments are self-converse), THM-4644 (Fermat = rigid circulant)
  - prime463_conductor_clock_20261017e.md (clock_duality_oct15: h(-463) = 7 and the half-sum 7 at p = 463 as separate predicates)
  - HYP-9286 (the Paley/Legendre local system)
---

# THM-4666 — Class numbers as ranking excesses of Kronecker tournaments

## Setting

`D < 0` a fundamental discriminant, `χ_D` its Kronecker character (odd: `χ_D(-1) = -1`), `h(D)` the class number,
`w_D = 6, 4, 2` for `D = -3, -4, D < -4`. The **Kronecker oriented graph** `Γ_D` has vertex set `Z/|D|` and an arc
`i → j` whenever `χ_D(j - i) = +1`; each pair with `gcd(j - i, D) = 1` carries exactly one arc. For `D = -p`,
`p ≡ 3 (mod 4)` prime, `χ_D(n) = (n/p)` and `Γ_D` is the Paley tournament `P_p`. The **natural order** lists the
vertices `0 < 1 < ... < |D| - 1`; its excess is (forward arcs) − (backward arcs). `S_D(t) = Σ_(1 <= j <= t) χ_D(j)` is
the Legendre walk.

## Statements

1. **(Excess = class number, PROVED.)** The natural order of `Γ_D` has excess
   `Σ_(d=1)^(|D|-1) (|D| - d) χ_D(d) = Σ_(t=1)^(|D|-1) S_D(t) = (2/w_D) |D| h(D)`.
   Every cyclic shift has the same excess, and the multiplier order `i ↦ a i (mod |D|)` (`a` prime to `D`) has excess
   `χ_D(a) (2/w_D) |D| h(D)`. For `P_p`, `p > 3`: the excess is `p h(-p)`.
2. **(Half window, CITED: Dirichlet.)** For `p ≡ 3 (mod 4)`, `p > 3`: vertex 0 of `P_p` has out-degree minus in-degree
   `S((p-1)/2) = (2 - (2/p)) h(-p)` on the half circle `{1, ..., (p-1)/2}`.
3. **(Arcs at a vertex are split primes, PROVED.)** For every prime `q ≠ p`, `0 → (q mod p)` in `P_p` iff `q` splits in
   `Q(sqrt(-p))` (reciprocity: `(-p/q) = (q/p)` for odd `q`; `q = 2` splits iff `p ≡ 7 (mod 8)` iff `(2/p) = 1`). Hence
   `h(-p) = (sqrt(p)/π) Σ_(n >= 1) ε(n)/n` with `ε(n) = +1` if `0 → n (mod p)`, `-1` if `n → 0`, `0` if `p | n`: the class
   number is the harmonically weighted out-flow of one vertex.
4. **(The least excess, CITED + reading.)** The natural order of `P_p` (`p > 3`) has the least possible excess `p`
   exactly when `h(-p) = 1`, i.e. for `p = 7, 11, 19, 43, 67, 163` (Heegner-Baker-Stark); equivalently (Rabinowitsch
   1913) for `p = 4q - 1` with `q = 2, 3, 5, 11, 17, 41`, Euler's lucky numbers (`n^2 + n + q` prime for
   `0 <= n <= q - 2`). `P_163` is the largest Paley tournament whose natural order is this poor; the statement is
   unconditional.
5. **(Local optimality = positive walk; no real zero, PROVED; Fekete-Polya CITED.)** In the natural order of `P_p`,
   moving one vertex forward or backward over `t` consecutive vertices changes the excess by exactly `-2 S(t)`, for every
   vertex and every admissible `t` (and `S(p - 1) = 0`, `S(p - 1 - t) = S(t)`). Hence the natural order admits no
   improving single-vertex insertion iff `S(t) >= 0` for all `t`, i.e. iff `p` lies in Kalmynin's set `L+`
   (arXiv:2105.06910). For such `p`, Fekete's identity `f_p(x) = Σ_t (t/p) x^t = (1 - x) Σ_t S(t) x^t` makes `f_p`
   positive on `(0, 1)`, and `Γ(s) L(s, χ_(-p)) = ∫_0^∞ u^(s-1) f_p(e^(-u)) / (1 - e^(-pu)) du > 0` for real `s > 0`:
   **a real (Siegel) zero of `L(s, χ_(-p))` forces an improving single-vertex move of the natural order of `P_p`.**
   `L+` has relative density 0 among the primes (Baker-Montgomery), with counting function
   `x / (log x (log log x)^(c + o(1)))`, `c ≈ 0.0368` (Kalmynin; CITED).
6. **(Ranking bounds, PROVED.)** `MAS(P_p) >= (p(p-1)/2 + p h(-p))/2`. Every linear order of `P_p` has excess at most
   `(sqrt(p)/2) Σ_(k=1)^p |cot((2k-1)π/(2p))| = (p^(3/2)/π)(log p + γ + log(4/π) + o(1))`: the excess is
   `(1/2)<S, P T_p P^T>` with `S` the skew adjacency matrix (`||S||_op = sqrt(p)`) and `T_p` the sign matrix of a linear
   order, whose eigenvalues are `±i cot((2k-1)π/(2p))`; trace duality gives the bound (order of magnitude as in Alon,
   "Ranking tournaments", with an explicit constant). The natural order reaches the fraction
   `L(1, χ_(-p)) / (log p + γ + log(4/π) + o(1))` of it.
7. **(FINITE-EXACT.)** Statement 1 at all 9006 primes `p ≡ 3 (mod 4)` below `2·10^5` and all 6079 fundamental
   `D >= -20000` (class numbers from reduced forms); statement 2 at the same primes; the area form at 9005 primes; the
   `h = 1` list below `2·10^5` is exactly `3, 7, 11, 19, 43, 67, 163`. `L+ ∩ [7, 2·10^5]` has 2937 of 9005 primes, and the
   proportion in `L+` rises with `L(1, χ)`: 0.000, 0.060, 0.451, 0.674, 0.879, 0.996 on `L < 1`, `[1, 1.5)`, `[1.5, 2)`,
   `[2, 2.5)`, `[2.5, 3)`, `L >= 3`. Exact `MAS(P_p)` (subset DP with vertex 0 first, brute force at 11; at 31 two DPs
   with the second vertex a residue or a non-residue): 2, 14, 35, 107, 161, 285 for `p = 3, 7, 11, 19, 23, 31`; the
   natural order (2, 14, 33, 95, 161, 279) is optimal at 3, 7, 23 only; at 43 the certified 543 (Chindelevitch-
   Harutyunyan, arXiv:2607.26690, Thm 7.1) exceeds the natural 473.

## Proof

1. Each pair at distance `d = j - i` (`1 <= d <= |D| - 1`) occurs `|D| - d` times and contributes `χ_D(d)`; since
   `Σ_d χ_D(d) = 0`, the total is `-Σ_d d χ_D(d) = (2/w_D)|D|h(D)` by `h(D) = -(w_D/(2|D|)) Σ_(a=1)^(|D|) a χ_D(a)`. Summation
   by parts gives `Σ_t S_D(t)`. A cyclic shift preserves all differences; a multiplier multiplies every `χ_D(j - i)` by
   `χ_D(a)`.
3. Reciprocity as displayed; `χ_(-p)(n) = (n/p)`; the series for `L(1, χ_(-p))` and `h(-p) = sqrt(p) L(1, χ_(-p))/π`.
5. In the natural order the vertex `v = i` and the vertices it jumps over differ by `1, ..., t`. Jumping forward, each
   pair changes from `(v/i+j)` order `χ(j)` to `χ(-j) = -χ(j)`; backward likewise. The change is `-2 Σ_(j <= t) χ(j)`.
   Vertex 0 moving forward realises every `t <= p - 1`. Fekete: Abel summation with `S(p - 1) = 0`, and
   `Σ_n χ(n) x^n = f_p(x)/(1 - x^p)` in the Mellin integral for `Γ(s)L(s, χ)`.
6. `MAS` bound from statement 1. For the upper bound, trace duality `|<A, B>| <= ||A||_op ||B||_(S_1)`; the skew matrix of
   `P_p` is circulant with eigenvalues `0` and `±i sqrt(p)`; the eigenvalues of `T_n` are verified in the script and
   classical; the asymptotic of the cotangent sum is checked numerically to 10 digits at `n = 10^5`.

## Remarks

* **Extreme values.** Joshi (J. Number Theory 2, 1970) proved `limsup L(1, χ)/log log q >= e^γ` over real characters of
  prime modulus. If this holds along `q ≡ 3 (mod 4)` (his method prescribes residue classes; the restriction has not
  been checked here), statement 6 gives Paley tournaments whose best ranking has excess at least
  `(e^γ/π + o(1)) p^(3/2) log log p`, unlike random tournaments (`Θ(n^(3/2))`, Spencer; de la Vega). Typed CONDITIONAL.
* **What is not claimed.** No statement about MAS beyond the computed and cited values; HYP-9287 records the observed
  dependence of global optimality on `L(1, χ_(-p))`.
