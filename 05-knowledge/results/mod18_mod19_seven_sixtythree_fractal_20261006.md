# Mod 18, mod 19, 7 and 63 on one 3-adic clock tower: the backward sieve keeps a positive proportion (D62 answered, limit in [0.288, 0.299]) and mirrors the forward glide law under 2 <-> 3; the rewrite residuals and the odd Mersenne numbers are backward-minimal; the Sierpiński tower's 63 = 8·7 + 7 and its self-similarity

**Session:** mac-mini-2026-10-06-mod1819 (worktree `math-wt-chessboard-20261006`), 2026-10-06.

**Owner's prompt (verbatim, attached to the line "I had called all 63 of P₇'s minima Hall obstructions, but 7 are not."):**
"you should merge past work with numbers mod 18 and mod 19 in with past ideas regarding 7 and 63 and fractal recursion systems. keep focusing on progressing next steps toward a collatz proof"

**Canon:** [THM-4554](../../01-canon/theorems/THM-4554-backward-sieve-keeps-a-positive-proportion-moran-duality.md). **Hypothesis updated:** [HYP-9162](../hypotheses/HYP-9162-sierpinski-tournament-tower-frobenius-21.md) (a reduction and new structure; still OPEN).

**Scripts:**

* [backward sieve](../../04-computation/experiments/mod1819_20261006_backward_sieve.py);
* [clock tower, Mersenne, residual seeds, Sierpiński tower](../../04-computation/experiments/mod1819_20261006_clock_tower.py).

Both `.out` files sit beside their scripts and end in `ALL CHECKS PASSED`.

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
  * the explicit arc rule of the Sierpiński tower;
  * the reduction of HYP-9162 to one lemma (section 5).
* **FINITE-EXACT:**
  * the bracket `0.28821 <= s_∞ <= 0.29912`;
  * the census law;
  * the Mersenne census;
  * all 239 residual seeds backward-minimal;
  * the tower's orbit structure and self-similarity (`k <= 8`);
  * the lemma of HYP-9162's reduction (`k <= 9`).
* **NUMEROLOGY:** typed where it occurs.
* **Collatz:** OPEN.

## 0. Answers in brief

1. **One clock.** Every repo appearance of 18, 19, 7 and 63 in the Collatz work is a level of a single tower. The level-`n` clock is `L_n = ord_(3^n)(2) = 2·3^(n-1)` (2 is a primitive root mod every `3^n`):

   | level | clock | `2^(L_n) - 1` |
   |---|---|---|
   | 2 | 6 | `63 = 9·7`, adds the prime 7 |
   | 3 | 18 | `2^18 − 1 = 27·7·19·73`, adds 19 and 73 |
   | 4 | 54 | adds 87211 and 262657 |

   **Mod 18 is the level-3 clock, and mod 19 lives on it:** `ord_19(2) = ord_27(2) = 18`, so `2^a mod 27 -> 2^a mod 19` is a group isomorphism `C_18 -> C_18`. The exponent `a mod 18` drives the 27-adic digit and the 19-adic phase simultaneously. This is why the mod-19 observers, the `a ≡ 13 (mod 18)` translation branch and the refinement of `a mod 6` to `a mod 18` keep appearing together.

   **63 is the level-2 clock**, giving the row braid, `R^3(n) − n = 21(3n+1)`, the hexagon operator `63λ^6 + λ^3 − 1` and Tao's mod-9 law over 63. **7 is its new prime**, and the clock of the trivial cycle (code `QR_7`).
2. **Fractal recursion: the two sieves share one Moran function (THM-4554).**
   * The backward (3-adic) predecessor tree of a unit has Chernoff generating function `(3/2)ρ(θ)^d` with `ρ(θ) = 3^(θ−1)/(2^θ − 1)`. Its zero set is that of the inverse tree's Moran function `g(s) = 2^-s + (1/3)(3/2)^s` (THM-4504).
   * Its minimum is `3^-(1-h)`, with `h = H(log_3 2) = 0.949956` the forward exceptional dimension. The backward census decays like `3^-(1-h)` per step where the forward glide census decays like `2^-(1-h)`.
   * Both count the same lattice paths, weighted `2^-K` and `3^-d`, and these weights coincide on the critical line `2^K = 3^d`.
3. **Progress toward the proof program: D62 is answered.**
   * Note 18's minimal counterexample must survive a backward sieve (no smaller predecessor). That sieve keeps a positive proportion: `s_∞ >= 0.2243` from the first moment alone, and `0.28821 <= s_∞ <= 0.29912` with exact depth-16 data.
   * So the backward side cannot exclude a minimal counterexample. All exclusion power sits on the forward side (a null set of dimension `h`) and in integrality.
4. **The residual obligations are backward-minimal.**
   * All 239 residual seeds of the rewrite compiler (`n <= 10000`) have no smaller ancestor: 153 are multiples of 3, and the other 86 were checked to depth 40. The base rate among comparable sources is 29%.
   * The odd-exponent Mersenne numbers `2^a − 1` (7, 127, 511, …) are backward-minimal with density `2s_∞ ≈ 0.59`, because they sample the 3-adic classes uniformly.
   * Their exceptions are governed by the higher clocks: `a ≡ 13 (mod 54)` (word `(2,2,1,1)`) and so on.
   * Hence the rewrite program's remaining task is to supply **forward** joins for a positive-density set. For the multiples of 3, which have no ancestors at all, only forward joins can exist.
5. **7 and 63 in the Sierpiński tower.**
   * The tower of HYP-9162 has the explicit arc rule `q(x,y) = Σ_(l : x_l = 1) (1 + y_l + [x, y differ below bit l])` over `F_2`.
   * `Aut(T_k)` has `2^(k−3)` orbits of size 7 and `2^(k−3) − 1` fixed points, and the fixed points induce `T_(k−3)`. At 63 vertices this is `63 = 8·7 + 7`: eight heptagon orbits around a fixed copy of the Paley heptagon `P_7`.
   * HYP-9162 reduces to one lemma: no vertex of the second copy has a doubly regular out-neighbourhood. The lemma is verified for `k <= 9`.
   * The resemblance of `63 = 56 + 7` to `P_7`'s blocking census (56 Hall obstructions + 7 exotic stars, the line you attached) is NUMEROLOGY: one comes from `2^6 − 1 = 7·2^3 + (2^3 − 1)`, the other from `Aut(P_7)`-orbit counts `21 + 21 + 14 + 7`.

## 1. The clock tower (PROVED identities; merge table)

| level `n` | modulus | clock `L_n` | `2^(L_n) − 1` | where it lives in the repo |
|---|---|---|---|---|
| 1 | 3 | 2 | 3 | parity gate; `n ≡ 2 (mod 3)` has the smaller predecessor `(2n−1)/3` |
| 2 | 9 | 6 | `63 = 9·7` | row braid (rows 2, 5, 8 mod 9); `R(n) = 4n+1`, `R^3(n) − n = 21(3n+1)`; hexagon operator `63λ^6+λ^3−1` (THM-4521); Tao's law `(8,16,11,4,2,22)/63`; trunk dead end `a_3 = 21 = 63/3`; `4 mod 9` rule `(8n−5)/9`; reset switch `63 ⇒ 31` |
| 3 | 27 | **18** | `262143 = 27·7·19·73` | the mod-18 exponent clock; **19**: `ord_19(2) = 18`, `19 = 3^3 − 2^3`, branches `C_a(x) = 2^-a(3x+1) mod 19` depend on `a mod 18`, `a ≡ 13 (mod 18)` is the order-19 translation; the `−17` cycle's code denominator `2^18 − 1`; trunk dead end `a_9 = 87381 = (2^18−1)/3 = 9·7·19·73` |
| 4 | 81 | 54 | adds 87211, 262657 | the boundary-compiler family `(2^(30+54h) − 73)/81`; Mersenne exception class `a ≡ 13 (mod 54)` (section 4) |
| 5 | 243 | 162 | adds 163, 2593, … | `ord_163(2) = 162` |

**Facts (PROVED; `clock_tower` A–B):**

* `v_3(2^(L_n) − 1) = n` (lifting the exponent).
* `2^a mod 27 -> 2^a mod 19` is a group isomorphism `(Z/27)^* -> F_19^*`, and the joint residue `2^a mod 513` has period exactly 18.
* The primes with 2-clock dividing 18 are exactly 3, 7, 19, 73.
* The trunk `a_j = (4^j − 1)/3 = R^(j−1)(1)` is the sibling ladder of 1 (the inverse tree's recursion `E D^2 = S E`, `S = R`). Its member `a_j` is a dead end (a multiple of 3, with no odd predecessor) iff `3 | j`, and `v_3(a_(3^k)) = k`. So the clock tower's Mersenne numbers divided by 3 are exactly the ladder's dead ends at the tower levels: `21`, `87381`, `(2^54 − 1)/3`, …

**NUMEROLOGY (recorded, not used).**

* The primes `2·3^k + 1` with 2 a primitive root, `19` and `163` (and `3`), join the tower at levels `k + 1`. They and `7` are Heegner numbers, the class-number-one discriminants whose `Q(sqrt(−d))` S17's note studies (`j = −15^3` at `d = 7`).
* `11`, `43` and `67` are Heegner numbers that are not of this form. No mechanism is known.

## 2. Fractal recursion: the backward sieve and its Moran duality (THM-4554)

The inverse tree (`collatz_procgen_20260924_inverse_tree_mod192`, Prop. 11) is a doubling chain carrying the backward trees of one sibling ladder `S^j(p_0)`, `S(p) = 4p+1`.

* **Its forward shadow (THM-4504):**
  * the Moran function `g(s) = 2^-s + (1/3)(3/2)^s` (roots 1, 2);
  * the exceptional set `Bad` of dimension `h = H(log_3 2)`;
  * the glide census `W_k = Θ(2^(hk) k^(−3/2))`.
* **Its backward shadow (new):** inverse words from a 3-adic unit, each legal on exactly one class mod `3^d` (weight `(1/2)3^(1−d)`), smaller iff `2^K < 3^d`.
  * Generating function `(3/2)ρ(θ)^d`, with `ρ(θ) = 3^(θ−1)/(2^θ − 1)` and `ρ − 1 = 2^θ(g − 1)/(2^θ − 1)`.
  * `min ρ = 3^-(1-h)` at `θ* = log_2(log_2 3/(log_2 3 − 1))`.
  * The first-passage masses follow `C(d) d^(−3/2) 3^-(1-h)d`; the fitted slope agrees with `−(1−h) ln 3` to `10^-5`. This is exactly the forward glide law with `2 -> 3`.

**Why the exponents agree.** A forward parity word of length `K` with `d` odd steps has 2-adic weight `2^-K`. The same path read backward (`d` inverse odd steps, `K` doublings) has 3-adic weight `3^-d`. On the critical line `2^K = 3^d` the weights are equal, and the path count there is `2^(K H(d/K)) = 3^(h d)`. This is the exact form of the mod-6 synthesis's "3-adic mirror (reverse the arrows and swap 2 and 3)". It is also the quantitative content of note 18's remark that the two sieves "interlock only through size".

## 3. D62 answered: the backward sieve keeps a positive proportion (THM-4554)

Note 18 (Proposition 3) showed that a minimal counterexample `n_0` has no smaller predecessor, and it observed survival fractions `0.500, …, 0.305`. Direction D62 asked whether the limit is positive.

* **First moment (PROVED).** `P(smaller predecessor) <= Σ_d (1/2)3^(1−d) F_d = 0.7757408148`. Here `F_d` counts first-passage compositions (`F_1..F_14 = 1, 1, 0, 2, 0, 7, 30, 0, 113, 0, 525, 2652, 0, 11433`). The sum is computed in exact rationals to `d = 1500`, with Chernoff-bounded tails below `10^-29`. So `s_∞ >= 0.22426`.
* **Exact depths (FINITE-EXACT).** `s(r)`, `r = 1..16`: `0.5, 0.3333, 0.3333, 0.31481, 0.31481, 0.31070, 0.30521, 0.30521, 0.30308, 0.30308, 0.30203, 0.30033, 0.30033, 0.29955, 0.29955, 0.29912`. These come from the vectorised recursion; an independent brute force agrees for `r <= 8`, and so does the repo's T-depth table.
* **Bracket.** `0.28821 <= s_∞ <= 0.29912`.
* **Meaning.** Among 3-adic units, the last-dip note's set `D` (targets with a smaller ancestor) has measure `1 − s_∞ ∈ [0.70088, 0.71179]`. Integers agree: 29.77% of odd `n <= 10^4` prime to 3 are backward-minimal.

**Typed consequence for the proof program.** The two-sided trap of note 18 is lopsided:

* the forward 2-adic survivors form a null set of dimension `h`;
* the backward 3-adic survivors have measure about 0.29;
* by CRT the two sieves are independent at every finite level and couple only through size.

So a proof cannot come from the backward side, and any proof must exploit forward size, which is what the repo's "two walls" diagnosis says. The theorem is sheet-blind (the `3n−1` sheet has the same backward sieve under `x -> −x`), as it must be by note 18's Theorem B.

## 4. Where the obligations live: residual seeds and Mersenne numbers (FINITE-EXACT + PROVED periodicity)

**The rewrite compiler's residual seeds.** In `checked_switch_phase19_20261004` the adaptive compiler on odd `n <= 10000` leaves 239 seeds. All are first-reset-2, and their open target is the "debt boundary" rewrites.

* **All 239 are backward-minimal.** 153 are multiples of 3 (no ancestors at all; 64%), and the other 86 have no smaller ancestor (exact integers, depth 40).
* Among comparable sources (odd, prime to 3, `≡ 3 mod 4`, first reset 2) only 29.4% are backward-minimal.
* So the compiler's ancestor and suffix machinery already handles every source with a smaller ancestor. **The residual task is exactly: forward joins for backward-minimal reset-2 sources.**

**Mersenne numbers.** `2^a − 1` is the canonical hostile family (`T^j(2^k − 1) = 3^j 2^(k−j) − 1`):

* **even `a`:** a multiple of 3; reset switch `2^a − 1 ⇒ 2^(a−1) − 1` (e.g. `63 ⇒ 31`);
* **`a ≡ 5 (mod 6)`:** `≡ 4 (mod 9)`, with smaller ancestor `(8n−5)/9` at depth 2;
* **`a ≡ 1, 3 (mod 6)`:** backward-minimal for every `a <= 120` except `13, 67, 69, 103`. Here `13` and `67 = 13 + 54` share the word `(2,2,1,1)` (class `13 mod 54 = ord_81(2)`); `103` uses a depth-6 word, and `69` a depth-14 word.

In total 35 of the 59 odd `a <= 120` are backward-minimal (`0.593`).

**PROVED (THM-4554(v)).** Odd `a mod 2·3^(k−1)` maps bijectively onto the classes `≡ 1 (mod 3)` of `2^a − 1 mod 3^k`. So for every depth `k` the exponents with a smaller ancestor via a word of length `<= k` are periodic with density exactly `1 − 2s(k)`. The census `0.593` sits inside `2s_∞ ∈ [0.5764, 0.5982]`.

7 = `2^3 − 1` is the first backward-minimal Mersenne number. Its debt state 13 then joins the root.

## 5. 7 and 63 in the Sierpiński tower (HYP-9162: reduction and structure)

The tower `T_(k+1) = T_k + {0'} + T_k'` (doubly regular, order `2^k − 1`, `T_3 = P_7`).

* **Arc rule (PROVED, checked `k <= 9`).** `H_(2^k)(x,y) = (−1)^(q(x,y))`, where `q(x,y) = Σ_(l : x_l = 1) (1 + y_l + [x ≢ y mod 2^l])` over `F_2` (bit 0 lowest). Equivalently, `x -> y` iff `#{l : x_l = 1, y_l = 0} + #{l > v_2(x ⊕ y) : x_l = 1}` is even.
* **Orbits and self-similarity (FINITE-EXACT `k <= 8`).**
  * `|Aut(T_k)| = 21`, with `2^(k−3)` orbits of size 7 (the low three bits nonzero) and `2^(k−3) − 1` fixed points (H-index `≡ 0 mod 8`).
  * The fixed points induce exactly `T_(k−3)`: the tower is a heptagon bundle over a smaller copy of itself.
  * At `k = 6`: **63 = 8·7 + 7**, eight 7-orbits around a fixed Paley heptagon `T_3 = P_7`.
* **The distinguished set (FINITE-EXACT `k <= 9`).**
  * The vertices whose out-neighbourhood is doubly regular are exactly the base heptagon plus the apex chain `0', 0'', …`.
  * The apexes form a transitive tournament dominating the base, and the top apex `0'` dominates all of them.
  * Since this set is `Aut`-invariant, `0'` (its unique source) is fixed by every automorphism.
* **Reduction (PROVED).** HYP-9162 follows from the lemma: *no vertex of the second copy `T_k'` has a doubly regular out-neighbourhood*. The argument:
  * the lemma puts the distinguished set inside `T_k ∪ {0'}`, and `0'` dominates `T_k`;
  * so every automorphism fixes `0'`;
  * HYP-9162's audit proved that the stabiliser of `0'` is the diagonal extensions of `Aut(T_k)`;
  * so `Aut(T_(k+1)) ≅ Aut(T_k)`, and by induction from `Aut(P_7) = F_21`, `Aut(T_k) = F_21` for all `k >= 3`.
* **The lemma's failing pairs.** Computationally, in the out-neighbourhood of `i'` every pair with the wrong common-out-neighbour count has the form `(x, y')`, `x ∈ N^+(i)`, `y ∈ N^-(i)`. The count is `|N^+(x) ∩ N^+(y) ∩ N^+(i)| + |N^+(x) ∩ N^-(y) ∩ N^-(i)|` against `λ = 2^(k−2) − 1`. Its average over all such pairs is exactly `λ` (PROVED), so the lemma is a variance statement about triple intersections. The arc rule above makes these intersections computable level by level; this is the proposed route to a proof.
* **Typed.** `F_21` is the Borel of `PSL(2,7)`, and its point stabiliser `C_3` is the Collatz Frobenius on the trivial cycle (THM-4553). The tower is a tournament object; its tie to Collatz is that dictionary, not a map.

## 6. What this changes for the Collatz program

**Settled:**

* D62: positive, with bracket.
* The backward/forward exponent duality.
* The residual obligations of the rewrite program are backward-minimal.

**Not changed:**

* Collatz is OPEN.
* None of the repo's barriers is crossed: these results are local and sheet-blind, as Theorem B requires.

**Next steps (ranked):**

1. **Forward joins for the backward-minimal reset-2 sources**, starting with the 153 multiples of 3. For those, only forward joins exist, and their orbits begin `3m -> oddpart(9m+1)`. Two questions:
   * Do the 153 split into finitely many word families with uniform joins, as the reset switch does for `(1^r, a >= 3)`?
   * Is there a "3-multiple switch"?
2. **Prove HYP-9162 via the lemma of section 5**, using the arc rule.
3. **Sharpen D62's bracket** by a first moment conditioned on survival to depth 16 (the restricted tail is much smaller than `0.0109`). Also identify the oscillating constant `C(d)` of the census law in terms of `{d log_2 3}`.
4. **Inherited and unchanged:** the last-dip sweep, H1 (HYP-9176) and the Mahler-side attack on HYP-9134/HYP-9127.

## 7. Verdicts

| claim | status |
|---|---|
| clock tower identities, `(Z/27)^* ≅ F_19^*` by `2 -> 2`, primes with 2-clock dividing 18 are `3, 7, 19, 73`, trunk dead ends at `3 | j`, `v_3(a_(3^k)) = k` | PROVED |
| backward sieve: `s_∞ >= 0.22425` (first moment) | PROVED (THM-4554(ii)) |
| `0.28821 <= s_∞ <= 0.29912` | FINITE-EXACT + PROVED tail |
| Moran duality: `ρ − 1 = 2^θ(g−1)/(2^θ−1)`, `min ρ = 3^-(1-h)`; census `~ d^(−3/2) 3^-(1-h)d` | PROVED; census FINITE-EXACT |
| Mersenne sampling: density `1 − 2s(k)` at each depth `k` | PROVED; census `35/59` FINITE-EXACT |
| all 239 compiler residual seeds backward-minimal | FINITE-EXACT |
| Sierpiński arc rule; orbit structure; fixed set `≅ T_(k−3)`; 63 = 8·7 + 7 | arc rule PROVED; structure FINITE-EXACT (`k <= 8`) |
| HYP-9162 ⇐ the second-copy lemma; lemma verified `k <= 9` | PROVED reduction; FINITE-EXACT lemma |
| Heegner primes `7, 19, 163` on the tower; `63 = 56 + 7` in two places | NUMEROLOGY |
| Collatz | OPEN |

## 8. Directions

* **D-α.** A "3-multiple switch": a uniform common-future rule for reset-2 multiples of 3 (the 153 seeds).
* **D-β.** HYP-9162 via the triple-intersection variance on the arc rule.
* **D-γ.** The census constant `C(d)`: is it a continuous function of `{d log_2 3}` (a Sturmian modulation)?
* **D-δ.** The level primes of the tower (`7`, `19 & 73`, `87211 & 262657`, `163 & …`): is there a uniform role for the "Artin primes" `p = 2·3^k + 1` (2 primitive) in the mod-`p` observers, beyond `p = 19`?

## 9. Audit record

*(independent audit pending; to be filled)*

## 10. Reproduction

```bash
python3 04-computation/experiments/mod1819_20261006_backward_sieve.py 16   # ~20 s, ~2 GB
python3 04-computation/experiments/mod1819_20261006_clock_tower.py         # ~70 s, needs dreadnaut
```
