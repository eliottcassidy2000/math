---
id: THM-4666
title: "Class numbers are ranking excesses: the natural order of the Kronecker oriented graph of a negative fundamental discriminant D has forward-minus-backward excess (2/w_D)|D|h(D), the area under the Legendre walk; for the Paley tournament P_p a single-vertex insertion changes the excess by -2S(t), so the natural order is insertion-locally optimal exactly for Kalmynin's primes L+ (all partial sums of (./p) nonnegative), which by Fekete's identity have no real zero of L(s, chi_(-p)); h(-p) = 1 (the least excess) exactly for p = 7, 11, 19, 43, 67, 163, i.e. Euler's lucky numbers; the inverse order has excess 3h(-p) and every pair-block order (1 - 2(2/p))h(-p)"
status: >
  PROVED (elementary): statements 1, 1b (discrete-log, inverse and pair-block orders; the residue sub-tournament), 1c
  (the Farey sweep over all rotation orders and the Farey class number identity), 3, 5, 5b, 6, from Dirichlet's class number formula, quadratic reciprocity and Jacobsthal sums. CITED: Dirichlet's half-window
  formula (statement 2), Heegner-Baker-Stark and Rabinowitsch (statement 4), Fekete-Polya (statement 5, proof included),
  Kalmynin (density 0 of L+: his Thm 1, via Baker-Montgomery; upper bound: Thm 2), Alon-Spencer section 9.1 and Alon,
  Ranking tournaments, Cor. 3.2 (order of the spectral bound). FINITE-EXACT: statement 7 and the identity checks of 1b.
  NUMERICAL: the Fekete-positivity census in 5b. Independently audited 2026-10-17
  (04-computation/experiments/s40_audit_20261017/REPORT.md): every statement and every FINITE-EXACT number reproduced
  with independent code; citation fixes applied (MISTAKE-626). No novelty claim for the identities of statements 1, 2
  (they rewrite Dirichlet's formula). Repository prior art: statement 2 as a Paley half-row imbalance
  (ogg_triangular_chowla_20260927.md section 2; special_primes_curling_genus_zero_langlands_20261015.md section 4), the
  character form of statement 1 (collatz_mod6_20260917_floor_sums_odd_functions.md S5), both identities at p = 463
  (prime463_conductor_clock_20261017e.md P4), MFAS(P_7) = 7 attained by the natural order
  (tournament_heptagon_flipdepth_kps.out), Euler's lucky numbers and class number one (HYP-2226, HYP-3730). Not found in
  the repository or in a quick literature search: the insertion/walk dictionary, the reading of L+ and of Fekete
  positivity as local optimality of a ranking, and the inverse-order and pair-block identities.
session: opus-2026-10-17-S40
source: 05-knowledge/results/binary_information_class_groups_largest_anomalies_20261017.md
scripts:
  - 04-computation/experiments/paley_class_number_ranking_20261017.py (+ .out, .json)
  - 04-computation/experiments/paley_legendre_walk_positivity_20261017.py (+ .out)
  - 04-computation/experiments/paley_mas_exact_20261017.py (+ p31 .out, .json)
  - 04-computation/experiments/paley_arithmetic_orders_20261017.py (+ .out, .json)
  - 04-computation/experiments/paley_mobius_orders_20261017.py (+ .out, .json)
  - 04-computation/experiments/paley_fekete_geometric_moves_20261017.py (+ .out, .json)
  - 04-computation/experiments/paley_rotation_farey_sweep_20261017.py (+ 480_20000 .out, .json)
related:
  - HYP-9287 (REFUTED in its threshold form by the audit: rotation orders beat the natural order at high-L primes)
  - THM-871 (rotational tournaments at the Fermat rungs; clock_duality's scope repair: all rotational tournaments are self-converse), THM-4644 (Fermat = rigid circulant)
  - prime463_conductor_clock_20261017e.md (clock_duality_oct15 P4: sum a(a/463) = -463*7 and sum_(a<=231)(a/463) = 7, i.e. statements 1-2 at p = 463; it keeps the mu14 all-odd Paley property separate from h(-463) = 7)
  - HYP-9286 (the Paley/Legendre local system); keller_marked_cubic_fourier_20261017d.md (the marked-coefficient mechanism used in 1b)
---

# THM-4666 — Class numbers as ranking excesses of Kronecker tournaments

## Setting

`D < 0` a fundamental discriminant, `χ_D` its Kronecker character (odd: `χ_D(-1) = -1`), `h(D)` the class number,
`w_D = 6, 4, 2` for `D = -3, -4, D < -4`. The **Kronecker oriented graph** `Γ_D` has vertex set `Z/|D|` and an arc
`i → j` whenever `χ_D(j - i) = +1`; each pair with `gcd(j - i, D) = 1` carries exactly one arc. For `D = -p`,
`p ≡ 3 (mod 4)` prime, `χ_D(n) = (n/p)` and `Γ_D` is the Paley tournament `P_p`. The **natural order** lists the
vertices `0 < 1 < ... < |D| - 1`; its excess is (forward arcs) − (backward arcs). `S_D(t) = Σ_(1 <= j <= t) χ_D(j)` is
the Legendre walk; `<x>` is the least nonnegative residue.

## Statements

1. **(Excess = class number, PROVED.)** The natural order of `Γ_D` has excess
   `Σ_(d=1)^(|D|-1) (|D| - d) χ_D(d) = Σ_(t=1)^(|D|-1) S_D(t) = (2/w_D) |D| h(D)`.
   Every cyclic shift has the same excess, and the multiplier order `i ↦ a i (mod |D|)` (`a` prime to `D`) has excess
   `χ_D(a) (2/w_D) |D| h(D)`. For `P_p`, `p > 3`: the excess is `p h(-p)`.

1b. **(Arithmetic orders carry class numbers, PROVED.)** For `P_p`, `p > 3`, with vertex 0 first:
   * the **discrete-log order** `g^0, g^1, ..., g^(p-2)` has excess exactly `1`, for every primitive root `g`;
   * the **inverse order** `1^(-1), 2^(-1), ..., (p-1)^(-1)` has excess `3 h(-p)`; equivalently
     `T3 := Σ_(1 <= a < b <= p-1) χ(ab(b - a)) = -3 h(-p)`;
   * every **pair-block order** (each pair `{x, -x}` in consecutive positions, its smaller representative in
     `[1, (p-1)/2]` first, the pairs in any order), for instance the order by `x^2 mod p`, has excess
     `(1 - 2(2/p)) h(-p)`: `3h` for `p ≡ 3 (mod 8)`, `-h` for `p ≡ 7 (mod 8)`;
   * the natural orders of the **residue and non-residue sub-tournaments** both have excess `(p - 3) h(-p)/4`.
   The natural order (`p h`) is the best affine order. Among the Moebius orders (rank `x` by `<φ(x)>`,
   `φ ∈ PGL_2(F_p)`, the pole first) the affine ones give `±p h(-p)` and the non-affine ones only small odd excesses
   (FINITE: exhaustive for `p <= 23`, extremes `±3, ±3, ±5, ±9 = ±3h` at `p = 7, 11, 19, 23`; samples of 300 at
   `31 <= p <= 71`, extremes `±3h` when `h >= 3` and `±7, ±9` at the class-number-one primes 43, 67).
1c. **(All rotation orders at once: the Farey sweep, PROVED; added after the audit.)** For `α ∈ (0, 1)` outside the Farey
   sequence `F_(p-1)`, let `R_α` list `x ∈ {0, ..., p-1}` by `frac(α x)` (vertex 0 first; the keys are distinct). Crossing a
   Farey fraction `a/b` changes the excess of `R_α` by
   `Δ_b = 2 χ(b) [S(q) + (b - s) A(q) + s A(q + 1)]`, `p - 1 = q b + s`, `0 <= s < b`, `A(k) = Σ_(t<k) S(t)`,
   independently of `a`. So the excess on the interval after the `n`-th Farey fraction is `p h(-p)` plus the first `n`
   jumps, and the best rotation order is the maximal prefix sum, an `O(p^2)` computation. Since `R_α` tends to the
   reversed natural order as `α -> 1`, this gives the **Farey class number identity**
   `Σ_(b=2)^(p-1) φ(b) Δ_b = -2 p h(-p)`. The natural order is the best rotation order iff every prefix sum of the jump
   sequence is `<= 0`; call this set of primes `R+`, the rotation analogue of `L+`. Since
   `frac((1 - α)x) = 1 - frac(αx)`, the excess is antisymmetric, `E(1 - α) = -E(α)`. In particular the order
   `(0, p-2, p-4, ..., 1, p-1, p-3, ..., 2)` (just below `α = 1/2`) has excess
   `-Δ_2/2 = -χ(2)[S((p-1)/2) + 2A((p-1)/2)]`. The first `(p-1)/2` jumps (`b = p - k`, `q = 1`) are `-2kχ(k)`, so `p ∈ R+`
   requires the weighted walk positivity `Σ_(k<=K) kχ(k) >= 0` for `K <= (p-1)/2`. This is necessary, but it fails to
   decide at 118 of the 245 scanned primes.
2. **(Half window, CITED: Dirichlet.)** For `p ≡ 3 (mod 4)`, `p > 3`: vertex 0 of `P_p` has out-degree minus in-degree
   `S((p-1)/2) = (2 - (2/p)) h(-p)` on the half circle `{1, ..., (p-1)/2}`.
3. **(Arcs at a vertex are split primes, PROVED.)** For every prime `q ≠ p`, `0 → (q mod p)` in `P_p` iff `q` splits in
   `Q(sqrt(-p))` (reciprocity: `(-p/q) = (q/p)` for odd `q`; `q = 2` splits iff `p ≡ 7 (mod 8)` iff `(2/p) = 1`). Hence,
   for `p > 3`, `h(-p) = (sqrt(p)/π) Σ_(n >= 1) ε(n)/n` with `ε(n) = +1` if `0 → n (mod p)`, `-1` if `n → 0`, `0` if
   `p | n`: the class number is the harmonically weighted out-flow of one vertex.
4. **(The least excess, CITED + reading.)** The natural order of `P_p` (`p > 3`) has the least possible excess `p`
   exactly when `h(-p) = 1`, i.e. for `p = 7, 11, 19, 43, 67, 163` (Heegner-Baker-Stark); equivalently (Rabinowitsch
   1913) for `p = 4q - 1` with `q = 2, 3, 5, 11, 17, 41`, Euler's lucky numbers (`n^2 + n + q` prime for
   `0 <= n <= q - 2`). `P_163` is the largest Paley tournament whose natural order is this poor; the statement is
   unconditional.
5. **(Local optimality = positive walk; no real zero, PROVED; Fekete-Polya CITED.)** In the natural order of `P_p`,
   moving one vertex forward or backward over `t` consecutive vertices changes the excess by exactly `-2 S(t)`, for every
   vertex and every admissible `t` (and `S(p - 1) = 0`, `S(p - 1 - t) = S(t)`). Hence the natural order admits no
   improving single-vertex insertion iff `S(t) >= 0` for all `t`, i.e. iff `p` lies in Kalmynin's set `L+`
   (Mathematika 69 (2023); arXiv:2105.06910). For such `p`, Fekete's identity `f_p(x) = Σ_t (t/p) x^t = (1 - x) Σ_t
   S(t) x^t` makes `f_p` positive on `(0, 1)`, and `Γ(s) L(s, χ_(-p)) = ∫_0^∞ u^(s-1) f_p(e^(-u)) / (1 - e^(-pu)) du > 0`
   for real `s > 0`: **a real (Siegel) zero of `L(s, χ_(-p))` forces an improving single-vertex move of the natural order
   of `P_p`.** `L+` has relative density 0 among the primes (Kalmynin, Thm 1, credited there to Baker-Montgomery, whose
   theorem is stated for almost all fundamental discriminants), and `#(L+ ∩ [1, x]) << π(x)/(log log x)^c`,
   `c ≈ 0.0368` (Kalmynin, Thm 2; CITED). This is an upper bound only; it is not known whether `L+` is infinite.

5b. **(Randomised moves = Fekete positivity, PROVED; census NUMERICAL.)** Moving vertex 0 of the natural order forward by
   a random distance with `P(T = t) = (1 - x) x^(t-1)` (`t <= p - 1`; with the remaining probability `x^(p-1)` it goes
   to the end, a change of 0) changes the excess by `-(2/x) f_p(x)` in expectation. So `f_p > 0` on `(0, 1)`, Fekete's
   sufficient condition for no real zero, says exactly that every geometrically randomised insertion strictly lowers the
   expected excess; `L+` is the deterministic version. Numerically 173 of the 217 primes `p ≡ 3 (mod 4)` in
   `[7, 3000]` are Fekete-positive against 97 in `L+` (every `L+` prime is Fekete-positive), with slowly decreasing
   shares (0.84 to 0.76 and 0.53 to 0.40 from `[7, 500)` to `[2000, 3000)`).
6. **(Ranking bounds, PROVED.)** `MAS(P_p) >= (p(p-1)/2 + p h(-p))/2`. Every linear order of `P_p` has excess at most
   `(sqrt(p)/2) Σ_(k=1)^p |cot((2k-1)π/(2p))| = (p^(3/2)/π)(log p + γ + log(4/π) + o(1))`: the excess is
   `(1/2)<S, P T_p P^T>` with `S` the skew adjacency matrix (`||S||_op = sqrt(p)`) and `T_p` the sign matrix of a linear
   order, whose eigenvalues are `±i cot((2k-1)π/(2p))`; trace duality gives the bound. The order `p^(3/2) log p` is
   classical (Alon-Spencer section 9.1; Alon, "Ranking tournaments", Cor. 3.2 gives `|excess| <= p^(3/2) ceil(log2 p)`);
   the trace-norm bound lowers the leading constant from `1/log 2` to `1/π`. The natural order reaches the fraction
   `L(1, χ_(-p)) / (log p + γ + log(4/π) + o(1))` of it.
7. **(FINITE-EXACT.)** Statement 1 at all 9006 primes `p ≡ 3 (mod 4)` below `2·10^5` and all 6079 fundamental
   `D >= -20000` (class numbers from reduced forms); statement 2 at the same primes; the area form at 9005 primes; the
   `h = 1` list below `2·10^5` is exactly `3, 7, 11, 19, 43, 67, 163`. `L+ ∩ [7, 2·10^5]` has 2937 of 9005 primes, and the
   proportion in `L+` rises with `L(1, χ)`: 0.000, 0.060, 0.451, 0.674, 0.879, 0.996 on `L < 1`, `[1, 1.5)`, `[1.5, 2)`,
   `[2, 2.5)`, `[2.5, 3)`, `L >= 3`. Exact `MAS(P_p)` (subset DP with vertex 0 first, brute force at 11; at 31 two DPs
   with the second vertex a residue or a non-residue): 2, 14, 35, 107, 161, 285 for `p = 3, 7, 11, 19, 23, 31`; the
   natural order (2, 14, 33, 95, 161, 279) is optimal at 3, 7, 23 only. At 43 the natural order (473) is not optimal: a
   witness of excess 183 (forward 543) exists, and 543 is the certified maximum (Chindelevitch-Harutyunyan,
   arXiv:2607.26690, Thm 7.1, an unrefereed preprint). The 1b identities hold at all 217 primes `p ≡ 3 (mod 4)` in
   `[7, 3000]`, and the marked-coefficient lemma below at every `t` for the ten primes `p <= 71`.

## Proof

1. Each pair at distance `d = j - i` (`1 <= d <= |D| - 1`) occurs `|D| - d` times and contributes `χ_D(d)`; since
   `Σ_d χ_D(d) = 0`, the total is `-Σ_d d χ_D(d) = (2/w_D)|D|h(D)` by `h(D) = -(w_D/(2|D|)) Σ_(a=1)^(|D|) a χ_D(a)`. Summation
   by parts gives `Σ_t S_D(t)`. A cyclic shift preserves all differences; a multiplier multiplies every `χ_D(j - i)` by
   `χ_D(a)`.

1b. *Discrete log.* Pairs with 0 contribute `Σ_k χ(g^k) = 0`. For positions `i < j` the sign is
   `χ(g^j - g^i) = (-1)^i χ(g^(j-i) - 1)`; summing `(-1)^i` over the admissible `i` leaves exactly the odd gaps
   `d = j - i`, so the excess is `Σ_(d odd) χ(g^d - 1) = Σ_(n non-residue) χ(n - 1)
   = (1/2)[Σ_(n ≠ 0) χ(n - 1) - Σ_n χ(n(n - 1))] = (1/2)(1 + 1) = 1` (Jacobsthal).
   *Inverse order.* The excess is `Σ_(a<b) χ(b^(-1) - a^(-1)) = -T3`, since `χ(b^(-1) - a^(-1)) = -χ(a)χ(b)χ(b - a)`.
   Put `b = <a t>`: `T3 = Σ_(t ≠ 0, 1) χ(t(t - 1)) N(t)` with `N(t) = Σ_a χ(a)[<a t> > a]` (**marked coefficient**). The
   order condition is linear in least residues, `[<a t> > a] = 1 - (<a(t - 1)> - <a t> + a)/p`, and
   `Σ_a χ(a)<a c> = χ(c) Σ_m m χ(m) = -χ(c) p h(-p)`; hence `N(t) = h(-p)(1 + χ(t - 1) - χ(t))`. With
   `Σ_t χ(t(t - 1)) = -1`, `Σ_(t ≠ 0, 1) χ(t) = -1` and `Σ_(t ≠ 0, 1) χ(t - 1) = 1` this gives
   `T3 = h(-p)(-1 - 1 - 1) = -3h(-p)`. (Sum out the free variable and keep the coefficient that can vanish: the
   marked-cubic mechanism of keller_marked_cubic_fourier_20261017d.md in a tournament setting.)
   *Pair blocks.* Between two pairs `{x, -x}`, `{y, -y}` the four cross terms are
   `χ(y - x) + χ(x - y) + χ(x + y) + χ(-x - y) = 0`, so only the within-pair terms `χ(-2s) = -(2/p)χ(s)` remain; their sum
   is `-(2/p) S((p-1)/2) = -(2/p)(2 - (2/p)) h(-p)` by statement 2.
   *Sub-tournaments.* With `T0 = Σ_(a<b) χ(b - a) = p h`, `T1 = Σ_(a<b) χ(a)χ(b - a) = Σ_a χ(a)S(a) = (p - 1)/2`
   (telescoping `S(a)^2 - S(a)S(a - 1)`) and `T2 = Σ_(a<b) χ(b)χ(b - a) = -T1` (reflection `(a, b) -> (p - b, p - a)`), the
   residue part is `(T0 + T1 + T2 + T3)/4 = (p h - 3h)/4` and the non-residue part `(T0 - T1 - T2 + T3)/4` the same.
1c. At `a/b` the keys coincide exactly on residue classes mod `b`. For `α = a/b ∓ ε` a class with nonzero key `r/b` is a
   contiguous block, ordered by decreasing `x` before and increasing `x` after. Its internal pairs have differences
   `-(j - i)b`, so it contributes `-χ(b)A(k)` before and `+χ(b)A(k)` after. The class of 0 (without the vertex 0) moves
   from the end, keys `1 - εx`, to just after 0, keys `εx`, and is reversed. Its cross terms with the other vertices
   change by `2 Σ_(x ∈ B) χ(x) = 2χ(b)S(k_0)`, because `Σ_(r ∈ Z/p) χ(x - r) = 0` and the double sum over the block is
   antisymmetric; the pairs with 0 are unchanged. There are `s` classes of size `q + 1`, `b - 1 - s` of size `q`, and
   `k_0 = q`. As `α -> 1^-` the order is `0, p-1, ..., 1`, of excess `-p h(-p)`, and each denominator `b >= 2` occurs
   `φ(b)` times. Checks: direct excess on every Farey interval for all `p ≡ 3 (mod 4)`, `p <= 47`; the identity at all
   338 primes `p ≡ 3 (mod 4)` up to 5000; the `α = 1/2` corollary at 154 primes below 2000. The sweep reproduces every audited rotation certificate (167: 1861, 263: 4079, 383: 7519, 439: 9313,
   503: 11657, 599: 15507, 647: 17661, 887: 28235, 983: 34225).
3. Reciprocity as displayed; `χ_(-p)(n) = (n/p)`; the series for `L(1, χ_(-p))` and `h(-p) = sqrt(p) L(1, χ_(-p))/π`
   (`p > 3`; at `p = 3` the unit factor gives 1/3).
5. In the natural order the vertex `v = i` and the vertices it jumps over differ by `1, ..., t`. Jumping forward, each
   pair changes from `χ(j)` to `χ(-j) = -χ(j)`; backward likewise. The change is `-2 Σ_(j <= t) χ(j)`. Vertex 0 moving
   forward realises every `t <= p - 1`. Fekete: Abel summation with `S(p - 1) = 0`, and `Σ_n χ(n) x^n = f_p(x)/(1 - x^p)`
   in the Mellin integral for `Γ(s)L(s, χ)`.
5b. `E[Δ_x] = Σ_(t=1)^(p-1) (1 - x) x^(t-1) (-2 S(t)) = -(2/x)(1 - x) Σ_t S(t) x^t = -(2/x) f_p(x)` (exact rational
   check in the script).
6. `MAS` bound from statement 1. For the upper bound, trace duality `|<A, B>| <= ||A||_op ||B||_(S_1)`; the skew matrix of
   `P_p` is circulant with eigenvalues `0` and `±i sqrt(p)`; the eigenvalues of `T_n` are verified in the script and
   classical; the asymptotic of the cotangent sum is checked numerically to 10 digits at `n = 10^5`.

## Remarks

* **Extreme values.** Joshi (J. Number Theory 2, 1970) proved `limsup L(1, χ)/log log q >= e^γ` over real characters of
  prime modulus; his Thm 2 (as stated in Brink-Moree-Osburn, Thm 3.1) gives the same bound along `q ≡ 7 (mod 8)` (local
  factor 1). So statement 6 gives Paley tournaments whose best ranking has excess at least
  `(e^γ/π + o(1)) p^(3/2) log log p`, unlike random tournaments (`Θ(n^(3/2))`, Spencer; de la Vega). Typed CITED
  (secondary source).
* **The natural order is usually not the best ranking.** Rotation (Beatty) orders `x -> {αx}`, which contain the affine
  orders, beat it at many primes of large `L(1, χ_(-p))` (167, 263, 383, 439, 503, 599, 647, 887, 911, 983; HYP-9287). It
  is exactly optimal at 3, 7, 23 and unbeaten by any known order at 47, 71, 191, 239, 311, 359, 431, 479 (audit). By the
  Farey sweep (1c), `R+ ∩ (479, 20000]` (among the 245 primes with `L >= 2.4`) is the 46 primes 719, 839, 1151, 1319, 1511, 1559, 2351, 2999, 3191, 3671, 3911, 4079, 4679, 4751, 4871, 5591, 5711, 5879, 6551, 6599, 8111, 8951, 9239, 9719, 10391, 10559, 10631, 10799, 11399, 11471, 12239, 12479, 12911, 13751, 14159, 14759, 14951, 15791, 16631, 17471, 18119, 18191, 18311, 18959, 19031, 19319, all with
  `L >= 3.58`. It is never optimal for `h(-p) = 1`, `p > 7` (11, 19 by exact MAS;
  43, 67, 163 by witnesses).
