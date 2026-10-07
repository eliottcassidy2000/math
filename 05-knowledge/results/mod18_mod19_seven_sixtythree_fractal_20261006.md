# Mod 18, mod 19, 7 and 63 on one 3-adic clock tower: the backward sieve keeps a positive proportion (D62 answered, limit in [0.288, 0.299]) and mirrors the forward glide law under 2 <-> 3; the rewrite residuals and the odd Mersenne numbers are backward-minimal; the Sierpiński tower's 63 = 8·7 + 7 and its self-similarity

**Session:** mac-mini-2026-10-06-mod1819 (worktree `math-wt-chessboard-20261006`), 2026-10-06.

**Owner's prompt (verbatim, attached to the line "I had called all 63 of P₇'s minima Hall obstructions, but 7 are not."):**
"you should merge past work with numbers mod 18 and mod 19 in with past ideas regarding 7 and 63 and fractal recursion systems. keep focusing on progressing next steps toward a collatz proof"

**Canon:** [THM-4554](../../01-canon/theorems/THM-4554-backward-sieve-keeps-a-positive-proportion-moran-duality.md), [THM-4555](../../01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md) (section 6b). **Hypothesis updated:** [HYP-9162](../hypotheses/HYP-9162-sierpinski-tournament-tower-frobenius-21.md) (a reduction and new structure; still OPEN).

**Scripts:**

* [backward sieve](../../04-computation/experiments/mod1819_20261006_backward_sieve.py);
* [clock tower, Mersenne, residual seeds, Sierpiński tower](../../04-computation/experiments/mod1819_20261006_clock_tower.py);
* [uniform switches = collisions at −1](../../04-computation/experiments/mod1819_20261006_trailing_ones_switch.py) (section 6b).

All three `.out` files sit beside their scripts and end in `ALL CHECKS PASSED`.

**Inputs merged** (read in full for this note, with four survey passes):

* the mod-19 cluster of 2026-10-04: `mod19_recursive_observers`, `mod19_route_lifts`, `checked_switch_phase19`, `inverse_ray_torus_clock`, `collatz_partitioned_completion`, `difference_families`, `triplet_crt_223_233`;
* the mod-6/mod-18 cluster: `collatz_mod6_20260917_*`, `collatz_blueprint_20260921_*`, `thirtysix_signed`, `procgen_cdiff`;
* the fractal-recursion cluster: `collatz_procgen_20260924_inverse_tree_mod192`, `procgen_family27`, THM-4504, HYP-9162;
* the 7/Paley thread: note 19, THM-4553;
* the opus S16 and S17 notes of today.

**Status.**

* **PROVED:**
  * THM-4554 (i), (ii), (iv), (v): backward-sieve positivity, the Moran duality, Mersenne sampling;
  * the clock-tower identities of section 1;
  * the explicit arc rule of the Sierpiński tower, its self-similarity (`q(8x, 8y) = q(x, y)`), and `F_21 <= Aut(T_k)`;
  * the reduction of HYP-9162 to one lemma (section 5).
* **FINITE-EXACT:**
  * the bracket `0.28820 <= s_∞ <= 0.29912`;
  * the Mersenne census;
  * all 239 residual seeds backward-minimal;
  * `|Aut(T_k)| = 21` (`k <= 8`);
  * the lemma of HYP-9162's reduction (`k <= 9`).
* **Numerical fit:** the census law `m_d ~ C(d) d^(-3/2) 3^-(1-h)d` (only the upper bound `m_d <= (3/2) 3^-(1-h)d` is proved).
* **NUMEROLOGY:** typed where it occurs.
* **Collatz:** OPEN.
* **Audit:** independently audited 2026-10-06 (section 9); the corrections are in MISTAKE-569.

## 0. Answers in brief

1. **One clock.** The appearances of 18, 19, 7 and 63 surveyed here (section 1's table) are levels of a single tower. Not every occurrence is structural: `checked_switch_phase19` itself flags its count of 19 even debt heights as non-structural. The level-`n` clock is `L_n = ord_(3^n)(2) = 2·3^(n-1)` (2 is a primitive root mod every `3^n`):

   | level | clock | `2^(L_n) - 1` |
   |---|---|---|
   | 2 | 6 | `63 = 9·7`, adds the prime 7 |
   | 3 | 18 | `2^18 − 1 = 27·7·19·73`, adds 19 and 73 |
   | 4 | 54 | adds 87211 and 262657 |

   **Mod 18 is the level-3 clock, and mod 19 lives on it:** `ord_19(2) = ord_27(2) = 18`, so `2^a mod 27 -> 2^a mod 19` is a group isomorphism `C_18 -> C_18`. The exponent `a mod 18` drives the 27-adic digit and the 19-adic phase simultaneously. This is why the mod-19 observers, the `a ≡ 13 (mod 18)` translation branch and the refinement of `a mod 6` to `a mod 18` keep appearing together.

   **63 = 2^6 − 1 is the level-2 Mersenne number**, giving the row braid, `R^3(n) − n = 21(3n+1)`, the hexagon operator `63λ^6 + λ^3 − 1` (THM-4520, THM-4521) and Tao's mod-9 law over 63. **7 is its new prime**, and the clock of the trivial cycle (code `QR_7`).
2. **Fractal recursion: the two sieves share one Moran function (THM-4554).**
   * The backward (3-adic) predecessor tree of a unit has Chernoff generating function `(3/2)ρ(θ)^d` with `ρ(θ) = 3^(θ−1)/(2^θ − 1)`. Its zero set is that of the inverse tree's Moran function `g(s) = 2^-s + (1/3)(3/2)^s` (THM-4504).
   * Its minimum is `3^-(1-h)`, with `h = H(log_3 2) = 0.949956` the forward exceptional dimension.
     * The backward first-passage mass is at most `(3/2) 3^-(1-h)d` (proved).
     * Numerically it decays like `d^(-3/2) 3^-(1-h)d`, where the forward glide census decays like `2^-(1-h)` per step (THM-4495).
   * The two are mirror families on opposite sides of the critical line `2^K = 3^d`. Glide words keep their density of ones above `log_3 2` on every prefix; backward first-passage words keep it below until the last step. The 2-adic weight `2^-K` and the 3-adic weight `3^-d` coincide on the line, and the shared exponent is the Cramér rate there.
3. **Progress toward the proof program: D62 is answered.**
   * Note 18's minimal counterexample is either divisible by 3 or survives a backward sieve (no smaller predecessor). That sieve keeps a positive proportion: `s_∞ >= 0.22425` from the first moment alone, and `0.28820 <= s_∞ <= 0.29912` with exact depth-16 data.
   * So no argument using only the 3-adic backward sieve can exclude a minimal counterexample. An exclusion must also use the forward side (a null set of dimension `h`) or integrality and size.
4. **The residual obligations are backward-minimal.**
   * All 239 residual seeds of the rewrite compiler (`n <= 10000`) have no smaller ancestor: 153 are multiples of 3, and the other 86 were checked to depth 40 (and 60). The base rate among comparable sources is `245/834 = 29.4%`.
   * The odd-exponent Mersenne numbers `2^a − 1` sample the 3-adic classes uniformly. So at each depth `k` the backward-minimal ones (no smaller ancestor through words of length `<= k`) have density exactly `2s(k)`, which tends to `2s_∞ ≈ 0.59`.
   * Their exceptions are governed by the higher clocks: the class `13 mod 54` (word `(2,2,1,1)`) and so on.
   * Hence every residual obligation below `10^4` needs a **forward** join. The residuals are a proper subset of the backward-minimal sources: 153 of the 416 comparable multiples of 3 and 86 of the 245 comparable others. For the multiples of 3, which have no ancestors at all, only forward joins can exist.
5. **7 and 63 in the Sierpiński tower.**
   * The tower of HYP-9162 has the explicit arc rule `q(x,y) = Σ_(l : x_l = 1) (1 + y_l + [x, y differ below bit l])` over `F_2`.
   * `F_21 = Aut(P_7)`, acting on the low three bits, is a group of automorphisms of every `T_k` (PROVED).
     * It has `2^(k−3)` orbits of size 7 and `2^(k−3) − 1` fixed points, and the fixed points induce `T_(k−3)` (PROVED).
     * `|Aut(T_k)| = 21` for `k <= 8` (nauty).
     * At 63 vertices this is `63 = 8·7 + 7`: eight heptagon orbits around a fixed copy of the Paley heptagon `P_7`.
   * HYP-9162 reduces to one lemma: no vertex of the second copy has a doubly regular out-neighbourhood. The lemma is verified for `k <= 9`.
   * The resemblance of `63 = 56 + 7` to `P_7`'s blocking census (56 Hall obstructions + 7 exotic stars, the line you attached) is NUMEROLOGY: one comes from `2^6 − 1 = 7·2^3 + (2^3 − 1)`, the other from `Aut(P_7)`-orbit counts `21 + 21 + 14 + 7`.

## 1. The clock tower (PROVED identities; merge table)

| level `n` | modulus | clock `L_n` | `2^(L_n) − 1` | where it lives in the repo |
|---|---|---|---|---|
| 1 | 3 | 2 | 3 | parity gate; `n ≡ 2 (mod 3)` has the smaller predecessor `(2n−1)/3` |
| 2 | 9 | 6 | `63 = 9·7` | row braid (rows 2, 5, 8 mod 9); `R(n) = 4n+1`, `R^3(n) − n = 21(3n+1)`; hexagon operator `63λ^6+λ^3−1` (derived in THM-4520, stated in THM-4521); Tao's law `(8,16,11,4,2,22)/63` on the units `1,2,4,5,7,8 mod 9` (`procgen_tcpc` §5); trunk dead end `a_3 = 21 = 63/3`; `4 mod 9` rule `(8n−5)/9`; reset switch `63 ⇒ 31` |
| 3 | 27 | **18** | `262143 = 27·7·19·73` | the mod-18 exponent clock; **19**: `ord_19(2) = 18`, `19 = 3^3 − 2^3`, branches `C_a(x) = 2^-a(3x+1) mod 19` depend on `a mod 18`, `a ≡ 13 (mod 18)` is the order-19 translation; the `−17` cycle's code denominator `2^18 − 1` (unaccelerated map; `2^11 − 1` for the shortcut map); trunk dead end `a_9 = 87381 = (2^18−1)/3 = 9·7·19·73` |
| 4 | 81 | 54 | adds 87211, 262657 | the boundary-compiler family `(2^(30+54h) − 73)/81`; Mersenne exception class `13 mod 54` (section 4) |
| 5 | 243 | 162 | adds 163, 2593, … | `ord_163(2) = 162` |

**Facts (PROVED; `clock_tower` A–B):**

* `v_3(2^(L_n) − 1) = n` (lifting the exponent).
* `2^a mod 27 -> 2^a mod 19` is a group isomorphism `(Z/27)^* -> F_19^*`, and the joint residue `2^a mod 513` has period exactly 18.
* The primes with 2-clock dividing 18 are exactly 3, 7, 19, 73.
* The trunk `a_j = (4^j − 1)/3 = R^(j−1)(1)` is the sibling ladder of 1 (the inverse tree's recursion `E D^2 = S E`, `S = R`). Its member `a_j` is a dead end (a multiple of 3, with no odd predecessor) iff `3 | j`, and `v_3(a_(3^k)) = k`. So the clock tower's Mersenne numbers divided by 3 are exactly the ladder's dead ends at the tower levels: `21`, `87381`, `(2^54 − 1)/3`, …

**NUMEROLOGY (recorded, not used).**

* The primes `2·3^k + 1` with 2 a primitive root, `19` and `163` (and `3`), join the tower at levels `k + 1`. They and `7` are Heegner numbers, the class-number-one discriminants whose `Q(sqrt(−d))` S17's note studies (`j = −15^3` at `d = 7`).
* The coincidence stops at 163. The next such prime is `86093443 = 2·3^16 + 1`, which is not a Heegner number (there is also a probable prime at `k = 320`).
* `11`, `43` and `67` are Heegner numbers that are not of this form. No mechanism is known.

## 2. Fractal recursion: the backward sieve and its Moran duality (THM-4554)

The inverse tree (`collatz_procgen_20260924_inverse_tree_mod192`, Prop. 11) is a doubling chain carrying the backward trees of one sibling ladder `S^j(p_0)`, `S(p) = 4p+1`.

* **Its forward shadow:**
  * the Moran function `g(s) = 2^-s + (1/3)(3/2)^s` (roots 1, 2; THM-4504);
  * the exceptional set `Bad` of dimension `h = H(log_3 2)` (THM-4504);
  * the glide census `W_k = Θ(2^(hk) k^(−3/2))` (THM-4495; THM-4504 gives `log_2 W_m / m -> h`).
* **Its backward shadow (new):** inverse words from a 3-adic unit, each legal on exactly one class mod `3^d` (weight `(1/2)3^(1−d)`), smaller iff `2^K < 3^d`.
  * Generating function `(3/2)ρ(θ)^d`, with `ρ(θ) = 3^(θ−1)/(2^θ − 1)` and `ρ − 1 = 2^θ(g − 1)/(2^θ − 1)`.
  * `min ρ = 3^-(1-h)` at `θ* = log_2(log_2 3/(log_2 3 − 1))`.
  * The first-passage masses are at most `(3/2) 3^-(1-h)d` (proved). Numerically they follow `C(d) d^(−3/2) 3^-(1-h)d` with `C(d)` in `[0.22, 0.62]`: the fitted slope agrees with `−(1−h) ln 3` to `10^-5`, and a free fit gives the exponent `−1.52`. This numerically matches the forward glide law with `2 -> 3`.

**Why the exponents agree.**

* A forward parity word of length `K` with `d` odd steps has 2-adic weight `2^-K`. A backward word (`d` inverse odd steps, `K` doublings) has 3-adic weight `3^-d`.
* The two families are not the same paths. Glide words keep their density of ones above `log_3 2` on every prefix; backward first-passage words keep it below until the last step. They are mirror families on opposite sides of the critical line `2^K = 3^d`.
* On the line the two weights are equal, and the number of paths along it is `2^(K H(d/K)) = 3^(h d)`. So both families decay at the Cramér rate of the line: `2^-(1-h)` per forward step, `3^-(1-h)` per backward step.

This makes precise the mod-6 synthesis's "3-adic mirror (reverse the arrows and swap 2 and 3)".

## 3. D62 answered: the backward sieve keeps a positive proportion (THM-4554)

Note 18 (Proposition 3) showed that a minimal counterexample `n_0` is divisible by 3 or has no smaller predecessor, and it observed survival fractions `0.500, …, 0.305`. Direction D62 asked whether the limit is positive.

* **First moment (PROVED).** `P(smaller predecessor) <= Σ_d (1/2)3^(1−d) F_d = 0.7757408148`. Here `F_d` counts first-passage compositions (`F_1..F_14 = 1, 1, 0, 2, 0, 7, 30, 0, 113, 0, 525, 2652, 0, 11433`). The sum is computed in exact rationals to `d = 1500`.
  * The dropped window is bounded by `Σ_d N_(d−1)·(1/2)3^(1−d)·R(θ)·2^(−θ e_d)/(1−2^(−θ)) = 4.1·10^-30`. Here `N_(d−1)` counts the in-window prefixes, `e_d` is the excess at the window edge, `θ = 3/2`, `R = ρ/(1−ρ)` and `ρ(3/2) = 0.94729`.
  * The depth tail is at most `(3/2)ρ^1501/(1−ρ) = 1.4·10^-34`.
  * So `s_∞ >= 0.22425`. The audit's untruncated recomputation gives the sum `0.77574081479760749602…`.
* **Exact depths (FINITE-EXACT).** `s(r)`, `r = 1..16`: `1/2, 1/3, 1/3, 17/54, 17/54, 151/486, 445/1458, 445/1458, 3977/13122, 3977/13122, 35669/118098, 106405/354294, 106405/354294, 955147/3188646, 955147/3188646, 4292002/14348907` (`= 0.299117`).
  * These come from the vectorised recursion. A brute force agrees for `r <= 8`, and the audit's independent class marking reproduces all sixteen rationals.
  * The repo's own table (`collatz_connectivity_from_rigidity_20261001.out`, section D) uses standard-map depth `D = K + d`. Every first-passage word of length `d >= 2` has `K = ⌊d log_2 3⌋`, so its value at `D` is `s(max{d : d + ⌊d log_2 3⌋ <= D})`, which reproduces all ten entries.
* **Bracket.** `0.28820 <= s_∞ <= 0.29912` (`s(16)` minus the tail is `0.2882095614`).
* **Meaning.**
  * The last-dip note's set `D` (targets with a smaller ancestor) is a set of integers. Its 3-adic analogue has Haar measure `1 − s_∞ ∈ [0.70088, 0.71180]` among units.
  * Integers are consistent: `992/3332 = 29.77%` of the odd `n` in `[5, 10^4]` prime to 3 are backward-minimal (`29.79%` with `n = 1`). The last-dip note observes `0.703004` with a smaller ancestor among units below `2^28`.
  * Equality of the integer and 3-adic densities is not proved.

**Typed consequence for the proof program.** The two-sided trap of note 18 is lopsided:

* the forward 2-adic survivors form a null set of dimension `h`;
* the backward 3-adic survivors have measure about 0.29;
* by CRT the two sieves are independent at every finite level and couple only through size.

So no argument using only the 3-adic backward sieve can exclude a minimal counterexample; an argument must also use the forward side or size, which is what the repo's "two walls" diagnosis says. The theorem is sheet-blind (the `3n−1` sheet has the same backward sieve under `x -> −x`), as it must be by note 18's Theorem B.

## 4. Where the obligations live: residual seeds and Mersenne numbers (FINITE-EXACT + PROVED periodicity)

**The rewrite compiler's residual seeds.** In `checked_switch_phase19_20261004` the adaptive compiler on odd `n <= 10000` leaves 239 seeds. All are first-reset-2, and their open target is the "debt boundary" rewrites.

* **All 239 are backward-minimal.** 153 are multiples of 3 (no ancestors at all; 64%), and the other 86 have no smaller ancestor (exact integers, depth 40 and 60). All 239 are `≡ 3 mod 4` with first reset 2.
* Among comparable sources (odd, prime to 3, `≡ 3 mod 4`, first reset 2) only `245/834 = 29.4%` are backward-minimal.
* So below `10^4` the compiler's ancestor and suffix machinery handles every source with a smaller ancestor. **The residual task is a finite statement:** forward joins for the 239 seeds.
* The seeds are a proper subset of the backward-minimal sources: 153 of the 416 comparable multiples of 3, and 86 of the 245 comparable others. Backward-minimality is necessary for a residual, not sufficient.

**Mersenne numbers.** `2^a − 1` is the canonical hostile family (`T^j(2^k − 1) = 3^j 2^(k−j) − 1`):

* **even `a`:** a multiple of 3; reset switch `2^a − 1 ⇒ 2^(a−1) − 1` (e.g. `63 ⇒ 31`);
* **`a ≡ 5 (mod 6)`:** `≡ 4 (mod 9)`, with smaller ancestor `(8n−5)/9` at depth 2;
* **`a ≡ 1, 3 (mod 6)`:** backward-minimal for every odd `3 <= a <= 119` except `13, 67, 69, 103`. Here `13` and `67 = 13 + 54` share the word `(2,2,1,1)` (class `13 mod 54`, `54 = ord_81(2)`); `103` uses a depth-6 word, and `69` a depth-14 word.

In total 35 of the 59 odd `3 <= a <= 119` are backward-minimal (`0.593`).

**PROVED (THM-4554(v)).** Odd `a mod 2·3^(k−1)` maps bijectively onto the classes `≡ 1 (mod 3)` of `2^a − 1 mod 3^k`. So for every depth `k` the exponents with a smaller ancestor via a word of length `<= k` are periodic with density exactly `1 − 2s(k)`, tending to `1 − 2s_∞`. The census `0.593` sits inside `2s_∞ ∈ [0.5764, 0.5983]`. A density at infinite depth would need an exchange of limits and is not claimed.

`7 = 2^3 − 1` is the first backward-minimal Mersenne number with odd `a >= 3` (`3 = 2^2 − 1`, a multiple of 3, also qualifies). Its debt state 13 then joins the root.

## 5. 7 and 63 in the Sierpiński tower (HYP-9162: reduction and structure)

The tower `T_(k+1) = T_k + {0'} + T_k'` (doubly regular, order `2^k − 1`, `T_3 = P_7`).

* **Arc rule (PROVED, checked `k <= 9`).** `H_(2^k)(x,y) = (−1)^(q(x,y))`, where `q(x,y) = Σ_(l : x_l = 1) (1 + y_l + [x ≢ y mod 2^l])` over `F_2` (bit 0 lowest). Equivalently, `x -> y` iff `#{l : x_l = 1, y_l = 0} + #{l > v_2(x ⊕ y) : x_l = 1}` is even.
  * *Proof.* `H` is skew off the diagonal: `H_2 + H_2^T = 2I`, and `H_2n + H_2n^T = [[H + H^T, 0], [0, H^T + H]]`.
  * Induct on the top bit, with `x'`, `y'` the lower bits. If `x`'s top bit is 0, the entry is `H_n(x', y')`, and `q` gains no term.
  * If it is 1, the entry is `−H_n(y', x')` (for `y`'s top bit 0) or `+H_n(y', x')` (for 1). By skewness `H_n(y', x') = (−1)^[x' ≠ y'] H_n(x', y')`, which is exactly the new term `1 + y_top + [x' ≠ y']`. ∎
* **Self-similarity and the Frobenius action (PROVED).**
  * Shifting bits gives `q(8x, 8y) = q(x, y)`, so for every `k` the multiples of 8 induce exactly `T_(k−3)`: the tower is a heptagon bundle over a smaller copy of itself.
  * `q(x,y)` is `q_3(x mod 8, y mod 8)` plus terms that depend only on the high bits and on `[x ≢ y (mod 8)]`. So every automorphism of `T_3 = P_7`, acting on the low three bits and fixing the multiples of 8, is an automorphism of every `T_k`.
  * Hence `F_21 <= Aut(T_k)`, with `2^(k−3)` orbits `{8m+1, …, 8m+7}` of size 7 and `2^(k−3) − 1` fixed points (H-index `≡ 0 mod 8`).
* **Equality (FINITE-EXACT `k <= 8`, nauty).** `|Aut(T_k)| = 21`. At `k = 6`: **63 = 8·7 + 7**, eight 7-orbits around a fixed Paley heptagon `T_3 = P_7`.
* **The distinguished set (FINITE-EXACT `k <= 9`).** Let `D` be the set of vertices whose out-neighbourhood is doubly regular.
  * For `k <= 9`, `D` is exactly the base heptagon plus the apex chain `0', 0'', …`.
  * The apexes form a transitive tournament dominating the base, and the top apex `0'` dominates all of them.
* **Reduction (PROVED).** HYP-9162 follows from the lemma: *no vertex of the second copy `T_k'` lies in `D`*. The argument does not use the finite characterization above:
  * `N^+(0') = T_k` is doubly regular, so `0' ∈ D`;
  * the lemma gives `D ⊆ T_k ∪ {0'}`;
  * since `0' -> T_k`, `0'` is the unique source of the `Aut`-invariant set `D`, so every automorphism fixes `0'`;
  * HYP-9162's audit proved that the stabiliser of `0'` is the diagonal extensions of `Aut(T_k)`, so `Aut(T_(k+1)) ≅ Aut(T_k)`;
  * by induction from `Aut(P_7) = F_21`, `Aut(T_k) = F_21` for all `k >= 3`.
* **The lemma's failing pairs.** Computationally, in the out-neighbourhood of `i'` every pair with the wrong common-out-neighbour count has the form `(x, y')`, `x ∈ N^+(i)`, `y ∈ N^-(i)`. The count is `|N^+(x) ∩ N^+(y) ∩ N^+(i)| + |N^+(x) ∩ N^-(y) ∩ N^-(i)|` against `λ = 2^(k−2) − 1`. Pairs of other types provably have count `λ`, and the failing type's average is exactly `λ` (PROVED, double counting). So the lemma is a variance statement about triple intersections. The arc rule above makes these intersections computable level by level; this is the proposed route to a proof.
* **Typed.** `F_21` is the Borel of `PSL(2,7)`, and its point stabiliser `C_3` is the Collatz Frobenius on the trivial cycle (THM-4553). The tower is a tournament object; its tie to Collatz is that dictionary, not a map.

## 6. What this changes for the Collatz program

**Settled:**

* D62: positive, with bracket.
* The backward/forward exponent duality.
* The residual obligations of the rewrite program below `10^4` are backward-minimal.

**Not changed:**

* Collatz is OPEN.
* None of the repo's barriers is crossed: these results are local and sheet-blind, as Theorem B requires.

**Next steps (ranked):**

1. **Forward joins for the backward-minimal reset-2 sources**, starting with the 153 multiples of 3. For those, only forward joins exist, and their orbits begin `3m -> oddpart(9m+1)`. Two questions:
   * Do the 153 split into finitely many word families with uniform joins, as the reset switch does for `(1^r, a >= 3)`?
   * Is there a "3-multiple switch"?

   *Attempted in section 6b (THM-4555).* Uniform switches are exactly collisions at −1, so no uniform switch is special to multiples of 3. Sporadic collisions resolve 16 of the 239 seeds (10 multiples of 3). The rest need run-length-dependent rewrites, and the debt states remain the live route.
2. **Prove HYP-9162 via the lemma of section 5**, using the arc rule.
3. **Sharpen D62's bracket** by a first moment conditioned on survival to depth 16. The audit measured the survival-conditioned tail at about 21–24% of the unconditioned mass at `d = 18, 19`, so the lower end should rise well above `0.2882`. Also identify the oscillating constant `C(d)` of the census law in terms of `{d log_2 3}`.
4. **Inherited and unchanged:** the last-dip sweep, H1 (HYP-9176) and the Mahler-side attack on HYP-9134/HYP-9127.

## 6b. Step 1 attempted: uniform switches are collisions at −1 (THM-4555)

The first next step asked for forward joins for the backward-minimal reset-2 sources, and for a "3-multiple switch". The attempt produced a classification of every switch that is uniform in the run length, the kind the reset switch is.

* **Endpoint progression (PROVED).** Write `n = 2^(r+1)t − 1` with word `1^r u`, `u` reduced (first letter `>= 2`). The endpoints `z = U^j(n)` are exactly `z ≡ f_u(−1) (mod 3^(r+|u|))`, where `f_c(x) = (3x+1)/2^c`. The reason: the word map is affine with slope `3^j/2^A`, and `f_1` fixes `−1`, so `n + 1 = 2^A (z − f_u(−1))/3^j`.
* **Uniform families are collisions (PROVED).** Two families `1^r u` and `1^(r+s) u'` share endpoints for infinitely many `r` iff `f_u(−1) = f_u'(−1)`. Colliding reduced words have equal totals, because the denominator of `f_u(−1)` is exactly `2^(A_u − 1)`.
* **The switch (PROVED).** If `u ~ u'` and `|u'| = |u| + D`, `D >= 1`, then every `n` with word `1^r u`, `r >= D`, merges at equal length `j = r + |u|` with
  * `m = (n+1)/2^D − 1`, i.e. `n` with `D` trailing binary ones deleted.
  * The reset switch `m = (n−1)/2` is the root collision `(a) ~ (2, a−2)`.
* **Reset 2 (PROVED).**
  * For a reset-2 source the root collision `(2, c, ·) ~ (c+2, ·)` points upward: it is the reset switch read backwards.
  * Any longer partner fires only on `n ≡ −1 (mod 3^(k−j))`, where `(2n−1)/3 < n` already. Shifted by 3, the root collision fires exactly on `n ≡ 8 (mod 9)`: 1/9 of reset-2 sources, none of them backward-minimal, hence none of the 239 seeds.
* **No 3-multiple switch (PROVED).** Within any family the multiples of 3 form one class of the endpoint progression mod `3^(j+1)`. The switch applies to them like any other source, so D-α has a negative answer in the uniform sense.
* **Sporadic collisions (FINITE-EXACT).**
  * `(8, c) ~ (4, 1, 1, c+2)`, both `125/2^(7+c)` (`125 = 2^7 − 3 = 2^5 + 3·2^4 + 9·2^3 − 27`). Composed with the root collision, `(2, 6, c) ~ (4, 1, 1, c+2)`, so reset-2 sources `1^r (2, 6, c)` switch to `(n−1)/2` at depth `r + 3`.
  * Among reduced words of length `<= 6`, letters `<= 8`, 774 of the 28014 collision values are not explained by the root collision.
* **Census (FINITE-EXACT).**
  * Reset-2 sources `n < 2·10^4`: 377 of 2500 (15.1%) have a uniform trailing-ones switch. Every equal-length merge of `n` with some `(n+1)/2^D − 1` is a collision.
  * The 239 residual seeds: 16, including 10 of the 153 multiples of 3.
  * Odd Mersenne numbers: 37 of the 60 odd `a` in `[3, 121]` switch to `2^(a−D) − 1`, always with `D` odd. The partner is an even Mersenne number, which the reset switch then shortens.
* **What it means.**
  * The reset-2 "debt" is the absence of a downward root collision. Uniform relief must come from sporadic collisions: Diophantine coincidences `−3^(p−1) + Σ 3^(p−1−i) 2^(A_i)` with two representations. These cover about 15% of reset-2 sources in the census.
  * The remaining rewrites must depend on the run length (the debt states of `checked_switch_phase19` (4)–(5)) or use partners other than trailing-ones deletions.
  * Mirror: the forward words of `2^a − 1` run on the 2-adic clock of 3 (`a mod 2^(K−2)`), while their backward-minimality runs on the 3-adic clock of 2 (`a mod 2·3^(k−1)`, THM-4554(v)).

## 7. Verdicts

| claim | status |
|---|---|
| clock tower identities, `(Z/27)^* ≅ F_19^*` by `2 -> 2`, primes with 2-clock dividing 18 are `3, 7, 19, 73`, trunk dead ends at `3 | j`, `v_3(a_(3^k)) = k` | PROVED |
| backward sieve: `s_∞ >= 0.22425` (first moment) | PROVED (THM-4554(ii)) |
| `0.28820 <= s_∞ <= 0.29912` | FINITE-EXACT + PROVED tail |
| Moran duality: `ρ − 1 = 2^θ(g−1)/(2^θ−1)`, `min ρ = 3^-(1-h)`, `m_d <= (3/2)3^-(1-h)d`; census `~ C(d) d^(−3/2) 3^-(1-h)d` | PROVED; census a numerical fit |
| Mersenne sampling: density `1 − 2s(k)` at each depth `k` | PROVED; census `35/59` (odd `3 <= a <= 119`) FINITE-EXACT |
| all 239 compiler residual seeds backward-minimal (a proper subset of the backward-minimal sources) | FINITE-EXACT |
| Sierpiński arc rule; fixed set `≅ T_(k−3)`; `F_21 <= Aut(T_k)`; `|Aut| = 21`, 63 = 8·7 + 7 | PROVED, except `|Aut| = 21`: FINITE-EXACT (`k <= 8`) |
| HYP-9162 ⇐ the second-copy lemma; lemma verified `k <= 9` | PROVED reduction; FINITE-EXACT lemma |
| uniform switches ⟺ collisions at −1; switch `m = (n+1)/2^D − 1`; reset switch = root collision; reset 2 has no downward root collision; no uniform 3-multiple switch | PROVED (THM-4555) |
| sporadic collisions; 377/2500 reset-2 sources, 16/239 seeds, 37/60 odd Mersenne exponents have uniform switches | FINITE-EXACT (THM-4555(vi)) |
| Heegner primes `7, 19, 163` on the tower (stops at 163); `63 = 56 + 7` in two places | NUMEROLOGY |
| Collatz | OPEN |

## 8. Directions

* **D-α.** A "3-multiple switch": a uniform common-future rule for reset-2 multiples of 3 (the 153 seeds). *Answered negatively in the uniform sense by THM-4555(iv): uniform switches are collisions at −1, and multiples of 3 are not special within a family.*
* **D-ε (new).** The sporadic collisions of the rational tree of −1: is their set described by finitely many parametric families? The numerators `N = −3^(p−1) + Σ 3^(p−1−i) 2^(A_i)` resemble the sums in the Collatz cycle equation; a collision is one `N` with two such representations. Their density decides how much of the reset-2 debt is uniformly payable.
* **D-β.** HYP-9162 via the triple-intersection variance on the arc rule.
* **D-γ.** The census constant `C(d)`: is it a continuous function of `{d log_2 3}` (a Sturmian modulation)?
* **D-δ.** The level primes of the tower (`7`, `19 & 73`, `87211 & 262657`, `163 & …`): is there a uniform role for the "Artin primes" `p = 2·3^k + 1` (2 primitive) in the mod-`p` observers, beyond `p = 19`?

## 9. Audit record

**Independent adversarial audit (2026-10-06, subagent, own code only).** The audit used Python, a GMP C search, PARI and dreadnaut, and did not run the session's scripts. **Verdict: all the mathematics checks out.** It required corrections to wording, citations and four bounds rounded the wrong way, all applied above and logged as MISTAKE-569.

* **Legality (PASS).** 8953 words with `d <= 6` follow the literal rule; each word is legal on exactly one class, `x_0 ≡ B_w 2^-K (mod 3^d)`.
* **First moment (values PASS).**
  * An untruncated dynamic program gives `0.77574081479760749602…`; the sums to `d <= 1500` and `d <= 3000` agree to about `10^-40`.
  * The prose tail bound was false as written: the measured dropped mass is about `0.152·2^-W`. It was replaced by the script's actual (valid) bound and the clean depth-tail bound `(3/2)ρ^1501/(1−ρ)`. Also `ρ(3/2) = 0.94729`.
* **Exact depths (values PASS).**
  * Independent class marking mod `3^16` gives the sixteen rationals of section 3, and a per-class search agrees for `r <= 10`.
  * Four bounds had been rounded inward; they are now `0.28820`, `0.71180`, `0.5983` and `0.22425`.
  * The repo's table was relabelled as standard-map depth (the script now reproduces all ten entries).
* **Moran duality (PASS).**
  * The identity is exact, the generating function checks to `10^-25`, `θ* = 1.4380327` and `min ρ = 0.9465045768 = 3^-(1-h)`.
  * The census fit is off by `9.4·10^-6` in slope (`3.0·10^-6` to `d = 3000`), with free-fit exponents `−1.526` and `−1.516` and `C(d)` in `[0.22, 0.62]`.
  * Corrections: cite THM-4495; type the census law as a fit; replace "same lattice paths" by mirror families.
* **Mersenne sampling (PASS).**
  * The density is `1 − 2s(k)` exactly for `k <= 12` (`10/27` at `k = 4`, `70742/177147` at `k = 12`).
  * `a ≡ 5 (mod 6)` gives `2^a − 1 ≡ 4 (mod 9)`.
  * The carry threshold is finite. The session's own follow-up bounds it by `max B_w/(2^K − 3^d) = 381.3` for words of length `<= 12`, and a direct search finds no odd unit `n < 3000` with a smaller supercritical ancestor at all.
* **Censuses (PASS).**
  * The exact-integer GMP search matches BFS on all odd `n < 400`.
  * Mersenne exceptions are exactly `a ≡ 5 (mod 6)`, plus 13 and 67 (word `(2,2,1,1)`), 103 (depth 6) and 69 (depth 14).
  * The 239 seeds are 153 + 86, all `≡ 3 mod 4` with first reset 2, and backward-minimal at depths 40 and 60.
  * The base rate is `245/834`, and `992/3332` excludes `n = 1`.
* **Clock tower (PASS, minor corrections).**
  * Every identity checks in PARI. The level-5 primes are `163, 2593, 71119, 135433, 97685839, 272010961`.
  * Tao's law is `(0, 8, 16, 0, 11, 4, 0, 2, 22)/63`, and the boundary family `(2^(30+54h) − 73)/81` satisfies `U^4 = 1` for `h <= 2`.
  * Corrections: cite THM-4520 as well as THM-4521; the `−17` cycle's `2^18 − 1` is for the unaccelerated map; the Heegner coincidence stops at 163.
* **Sierpiński tower (computations PASS).**
  * The arc rule matches for `k <= 10`, and `|Aut| = 21` for `k <= 8` with the stated orbits.
  * `D` is the heptagon plus the apex chain through `T_10`, and `Stab(0') = 21` with diagonal generators for `T_4..T_9`.
  * Failing pairs are only of type `(x, y')`, with exact average `λ`.
  * Corrections: write out the arc-rule proof (induction on the top bit using skewness); type the self-similarity PROVED (`q(8x, 8y) = q(x, y)`); and make the reduction independent of the `k <= 9` characterization of `D`. All three are done in section 5 and the HYP-9162 update.
* **Overclaims (corrected).** These were "every repo appearance", "a proof cannot come from the backward side", the "interlock only through size" sentence (deleted), "the residual task is exactly", "density `2s_∞`", the setting of note 18's Proposition 3 (or divisible by 3), and the last-dip measure sentence.
* **Confirmed:** the survival-conditioned tail at `d = 18, 19` is about 21–24% of the unconditioned mass (section 6, step 3).

## 10. Reproduction

```bash
python3 04-computation/experiments/mod1819_20261006_backward_sieve.py 16   # ~20 s, ~2 GB
python3 04-computation/experiments/mod1819_20261006_clock_tower.py         # ~3 min, needs dreadnaut
```
