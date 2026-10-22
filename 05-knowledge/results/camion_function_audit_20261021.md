# Audit: THM-4750 (the Camion function) and HYP-9361

**Independent adversarial audit, 2026-10-21.** Pre-audit packet of mac-mini-2026-10-21-fanopsp. The auditor is an independent Claude subagent. It did not edit any canon, hypothesis or results file and ran no git write command. Its scripts and outputs are in `04-computation/experiments/camion_wieferich_20261021/audit/`. C code was compiled into a scratch directory outside the repository.

**Files audited.**

* THM-4750 = `01-canon/theorems/THM-4750-the-camion-function-...-and-the-camion-limit.md`
* HYP-9361 = `05-knowledge/hypotheses/HYP-9361-camion-function-open-laws-....md`
* HYP-9358 = `05-knowledge/hypotheses/HYP-9358-camion-pseudoprimes-open-laws-....md`: the "Continuation (THM-4750 (5))" bullets and the 10^4 statistics
* RES = `05-knowledge/results/tournament_pseudoprimes_camion_fermat_20261021.md`, §13
* the scripts and outputs in `04-computation/experiments/camion_wieferich_20261021/`

References of the form `THM-4750:194` give a file and line number.

## 0. Verdict

**The core mathematics is sound.** I re-derived every proof step of statements 1, 2, 3, 4(a)–(d), 5, 6(a) and 7, and found no error in any of these arguments. That includes:

* the Leibniz/Mahler estimate;
* the nonunit-node vanishing;
* the `F′` filtration at Wieferich primes;
* the grade-1 computation behind the curvature formula;
* the existence and uniqueness of the Camion function and its restricted-power-series bound;
* the Lipschitz limit;
* Moon parity;
* the Vandermonde period and the threshold for hypergraphs.

**Independent code reproduces every FINITE-EXACT number I tested.** The audit code uses different recurrences and different algorithms from the session's:

* the Camion–Wieferich list `p ≤ 10^5` = {3, 13, 22409}, by a quasi-linear FLINT search;
* `v_22409(D(22409)) = 2` and `h(22409) = 2`, by a new 3-term evaluation of `Λ_p mod p³`;
* `n*` to 140 bits by Kitamasa evaluation, and to 90 bits by a closed form through the unit roots of `T`;
* the full 4906-class root table, and the r = 3 table.

**Two MAJOR findings.** Both concern statements, not the core proofs.

* **M1. The split/inert dichotomy of statement 4 is false in the double-root case, inside its stated range `p ≥ 5`.**
  * At `p = 5` the triangle class (`m = 3`) has a ramified pair of second-order zeros, and so does class 2 at `p = 7` (`m = 8`). Neither is split nor inert.
  * So the title clause "the triangle class m = 3 is always split" is false, and (4e) omits the case `χ(Δ) = 0`. Finding 11 has the details.
* **M2. The rational-reconstruction evidence for `n*` is false as written.**
  * The session's own output shows `765490658/1598160473`, with numerator and denominator below `2^31`, matching `n*` mod `2^64`.
  * The intended conclusion is true. This audit certifies it more strongly: that fraction fails mod `2^80`, and there is no rational of height `≤ 2^69.5`. Finding 22 has the details.

**MINOR findings.** Thirteen findings in §1 carry a MINOR tag (5, 9, 12, 13, 15, 18, 21, 25–29, 31); several bundle more than one issue. They include:

* a wrong range of validity for the (4e) valuation formula when `n` is even;
* missing hypotheses in 5(e) and in the general-base "exactly";
* a false sub-claim in 6(b) (the `mod 2^(j+3)` statement fails at `n = 1, 4`);
* a wrong expectation value (5.5 should be about 4.0);
* an unproved mechanism presented as an explanation (the r = 3 forced root);
* provenance gaps (the 17158 checks have no output file; one histogram is cited to the wrong file);
* one overclaim in RES ("sharp at every tested non-Wieferich prime").

**Upgrades the session can take (OK+).**

* `ε^a D(p^a) ≡ Λ_p` holds mod `p^(a+1)`, not just mod `p^a`.
* `Λ_p mod p³` needs only `s` mod `p²` up to `3p`, plus `s(p) mod p³`.
* The extra divisibility at `p = 3, 5, 7` equals `v_p((2k−1)!!)` exactly.
* The HEURISTIC closed form of HYP-9361(e) for `n*` is correct to 90 bits, and the 2-adic statement 6(c) is provable along the lines sketched in finding 23.

**Literature.** I found no prior statement of the higher Kummer congruences, of the curvature formula, of the 2-adic statements, or of the strong-class parity law. The nearest prior art is listed in §3.

**Recommended labels after the fixes.**

* PROVED: statements 1, 2, 3, 4(a)–(d), 5, 6(a), 7, with the M1 rewording of the title; and 4(e) with the corrected trichotomy and range.
* FINITE-EXACT: 6(b) and 6(c), with the corrected sub-claim and corrected evidence. Upgradable.
* HYP-9361 stays OPEN, with the (e) evidence line rewritten.

## 1. Findings

### Statement 1: half-speed bound (checklist 1)

**1. OK. The proof is correct** (THM-4750:86–100). Each step:

* *Linearity.* `S` is a finite `R`-combination of generators, and coefficients in `R` do not lower valuations.
* *Nonunit nodes.* `v ≥ j + n − d ≥ n − j ≥ (2k−1) − (k−1) = k`, using `d ≤ 2j` and `j ≤ k − 1`. `β = 0` is included because `n > d`.
* *Unit nodes.* `β^(p−1) = ε(1 + py)` since `b^(p−1) ≡ 1` in `R_0`. The identity `ε^(−l) g(m + l(p−1)) = p^j β^(m−d) P(l) Q(l)` is right.
* *Leibniz rule.* The discrete Leibniz rule `Δ^i(PQ)(0) = Σ_a C(i,a) Δ^a P(0) Δ^(i−a) Q(a)` is correct (write `Δ = Δ_P E_Q + Δ_Q`). With `Δ^c Q = (py)^c Q` and `Δ^a P(0) = 0` for `a > d`, it gives `v ≥ j + max(0, i − d) ≥ max(j, i − j) ≥ ⌈i/2⌉`.
* *Circulants.* `2^((m−1)/2)(pz)^i` has `v ≥ i`.
* *Index shift.* There is no issue: the class indices are `n ≥ m ≥ 2k − 1 ≥ 1`, where `S_n = s(n)`, and `S_0 = 0` is never used.
* *Integer congruence.* `Z/p^k → R` is injective, so the congruence holds for integers.

The higher Kummer congruence is `ε^(2k−1) a_(2k−1)(m)`, so its sign convention is right.

**2. OK. Independent exact check** (`audit/b1_half_speed_exact.py`, `.out`).

* *Method.* Exact big integers. `s` comes from the Hurwitz inverse `W(n) = −Σ_(k≥1) C(n,k) T(k) W(n−k)`, not Wright's recurrence. This is cross-checked against Wright for `n ≤ 60` and against brute force for `n ≤ 5`.
* *Bound.* 28 pairs `(p,k)`: `p = 3` (`k = 2, 3, 4, 6`), `p = 5` (`k ≤ 5`), `p = 7, 11, 13` (`k ≤ 4`), `p = 17, 19, 23` (`k ≤ 3`), `p = 29..43` (`k = 2`). 34400 coefficient checks at `m ≥ 2k − 1`, **0 violations**.
* *Higher Kummer congruence.* 0 failures at `m ≥ 2k − 1`, for both `D` and `s`.
* *Failures below the threshold* are as in `higher_kummer.out`. For example `p = 5, k = 5` fails at `m = 1, 2, 3, 4, 6` (`m = 5` passes), and `p = 7, k = 2` and `p = 13, k = 2` fail only at `m = 1`.
* *Strongest exact form.* `v(a_i(m)) ≥ min(⌊(m+1)/2⌋, ⌈i/2⌉)` (all admissible `k` at once) holds for every `m ≥ 1`.
* *Sharpness.* The minimum over classes equals `min(k, ⌈i/2⌉)` exactly for `p ≥ 7`, and for `p = 5` with `k ≤ 3`.

### Statement 2: full speed at Wieferich primes (checklist 2)

**3. OK. `F′` is a ring, `T ∈ F′`, and the inversion carries over** (THM-4750:113–121).

* *Ring.* `(p^J δ_d ∗ E_β)(p^J′ δ_d′ ∗ E_β′) = p^(J+J′) C(d+d′, d) δ_(d+d′) ∗ E_(β+β′)` with `d + d′ ≤ J + J′`.
* *`T ∈ F′`.* At a Wieferich prime `pu = p²u′`, so the grade-`j` Newton term is `p^(2j)` times degree `≤ 2j`: take `J = 2j`.
* *Inversion.* `A` is the same invertible group-algebra element, and `U` has `J ≥ 2 > 0`, so the finite geometric series stays in `F′`.
* *Unit nodes.* `J + max(0, i − d) ≥ i` because `d ≤ J`. For unit `β` the formula `C(n,d) β^(n−d)` holds for every `n ≥ 0`, so no condition on `m` arises.
* *Nonunit nodes.* `v ≥ J + n − d ≥ n`, so the threshold `m ≥ k` is right.

**4. OK. Independent test, with the threshold `m ≥ k` sharp** (`audit/b2_full_speed.c`, Hurwitz-inverse recurrence, `.out`). The session tested only `m ≥ 2K − 1`.

| prime and base | `K` | result |
|---|---|---|
| `p = 1093` | 3, 4, 5 | full speed `min(K,i)` holds for every base `m ≥ K`; failures at exactly `m = 1..K−1` |
| `p = 3511` | 3, 4 | same |
| base 3 at `p = 11` | 4, 5 | full speed from `m ≥ K` |
| base 7 at `p = 5` | 4, 5 | full speed from `m ≥ K` |
| controls `p = 13` (bases 2, 3), `p = 7` (base 5), `p = 1091` | | half speed only, as claimed |

**5. MINOR. The general-base claim "full speed occurs exactly at the base-`a` Wieferich primes" (THM-4750:123) needs its converse, which is not argued.**

* The missing step is that at a non-Wieferich prime some class has `v(a_2) = 1`, which needs some `X_a(m−2) ≢ 0`.
* A two-line fix: `X_a(1) = a³ − 2a` and `X_a(2) = a⁶ − 4a³ − 2a² + 6a` (check: 4 and 36 at `a = 2`, and 621 at `a = 3`).
* If `p ∤ a` and `X_a(1) ≡ 0`, then `a² ≡ 2`, so `X_a(2) ≡ 4 − 2a`. That vanishes only if `a ≡ 2`, which gives `4 ≡ 2`, impossible for odd `p`. So class 3 or class 4 has `v(a_2) = 1`.

### Statement 3: Fermat-quotient curvature (checklist 3)

**6. OK. The derivation is correct** (THM-4750:159–168).

* `T ≡ t_0(1 + pu·C(n,2))` mod `p²` holds exactly for all `n`.
* `C(n,2) t_0(n) = (δ_2 ∗ Δ²t_0)(n)`.
* `T^(−1) ≡ A^(−1) − U_1 A^(−2)`, so the grade-1 part of `S` is `+pu δ_2 ∗ X̃`. The sign is right.
* `X̃ ≡ X` mod `p`, and the quasi-periodicity of `X` mod `p` needs index `≥ 1`, hence `m ≥ 3`.
* `Δ_l² C(m + l(p−1), 2) = (p−1)² ≡ 1`.
* Grade 0 and the circulants contribute `(py)² ≡ 0`, and nonunit nodes vanish mod `p²` for `n ≥ 3`.
* `u ≡ −q_2`, from `2^(p−1) = (1+pu)^(p−1)`.

**7. OK. Independent check** (`audit/b3_curvature.py`, `.out`).

* `X` is computed through the other identity, `X = S″ + 2S′²T`, from `s`. The session used `Δ²T ∗ (e − S)²`.
* The formula holds in all 13806 classes of all primes `5 ≤ p ≤ 400`, with 0 failures.
* The direct quadratic law mod `p²` (`l ≤ p + 1`, all classes, `p ≤ 61`) passes 20928 checks with 0 failures.
* At `p = 1091` the unique class with `p² | a_2` is `m = 603`, where `X(601) ≡ 0`, as claimed (THM-4750:143). At `p = 1097` there are two such classes, `m = 558` and `843`, both with `X ≡ 0`.

**8. OK. The combinatorial meaning of `X` is proved and verified.**

* *Proof sketch (THM-4750:64).* With `{u,v}` made two-way, the digraph is semicomplete, so its condensation is transitive and `u, v` share a strong block. Hence `T(m)/2 = (T ∗ Ỹ ∗ T)(m − 2)` with `Ỹ = X/2`, and each strong digraph comes from two tournaments.
* *Checks.* My own brute force (`m ≤ 6`) gives 2, 4, 36, 680, 25152. The exact `X(0..6)` includes `X(6) = 244023264`. The base-derivative identity `C(n,2)X(n−2) = 2 ds_a(n)/da` at `a = 2` was checked exactly for `n ≤ 7` with polynomial arithmetic in `a`.

**9. MINOR. The 17158 direct checks have no recorded output** (THM-4750:14, 172).

* All four committed `second_order_*.out` files say `direct-law checks=0`.
* Rerunning the session's `camion_second_order.c 3 59 60` (compiled in scratch) gives `direct-law checks=17158 failures=0`, so the claim is true; the output should be committed.
* The rerun of `camion_second_order 3 4000` is byte-identical to the committed `second_order_3_4000.out`.

### Statement 4: second-order zeros (checklist 4)

**10. OK. Parts (a)–(d) are correct, and so are most of (e)** (THM-4750:182–198).

* *Existence and uniqueness.* Compatibility of the level-`k` polynomial functions on the dense sets `{l ≥ L_k}` gives existence; continuity and density give uniqueness.
* *Bound at the origin.* "Re-expanding at the origin keeps the bound" is right: with `C(l − L, i) = Σ_j C(−L, i−j) C(l, j)`, the coefficient of `C(l,j)` is a `Z_p`-combination of `a_i(m′)` with `i ≥ j`, so `v ≥ min(k, ⌈j/2⌉)` for every `k`.
* *Power series.* `v(b_j) ≥ min_(i≥j)(⌈i/2⌉ − v(i!)) ≥ (j+1)/4` for `p ≥ 5`, so `F_m` is a restricted power series.
* *Reduction.* `⌈i/2⌉ − v(i!) ≥ 2` for `i ≥ 3` (checked case by case for `i ≤ 6`, and from `(i+1)/4` for `i ≥ 7`), so the reduction of `F_m/p` is `g_m`.
* *Exact data.* `audit/b10_valuation_law.py` uses exact integers for `n ≤ 800`. For every zero class with simple roots at `p = 7, 11, 13`, one consistent `ℓ ∈ Z_p` reproduces every exact valuation `v_p(D(m + l(p−1)))` on its disc (the available data fix `ℓ` mod `p` to mod `p³`, depending on the class), and the valuation is exactly 1 off the root discs.

**11. MAJOR (M1). The "split or inert" dichotomy is false when `χ(Δ) = 0`, and "the triangle class m = 3 is always split" is false, both inside the stated range `p ≥ 5`.**

*Where.* THM-4750:3 (title of (4): "two p-adic zeros (split in Z_p or an inert conjugate pair) … the triangle class m = 3 is always split") and THM-4750:193 ((4e) lists only `χ(Δ) = ±1`).

*Counterexamples* (`audit/b11_double_root_ramified.py`, `.out`; exact integers). The session's own `second_order_3_4000.out` tags both classes as double roots: `p=5 … 3:1D` and `p=7 … 2:1D`.

* `p = 5`, triangle class (`m = 3`). `g_3 = 4l + 3C(l,2)` has a double root at `l = 0`, and `q_2 X(1) ≢ 0`. Then `v_5(D(3 + 4l)) = 2` for all 39 values `l ∈ {5, 10, …, 195}`, which cover every residue `≡ 0 (mod 5)` mod 25.
* `p = 7`, class 2 (`m = 8`). `g` has a double root at `l = 6`. Then `v_7(D(8 + 6l)) = 2` for all 19 values `l ≡ 6 (mod 7)`, which cover every residue mod 49.

*Why the zeros are ramified in both cases.*

* By (4e) itself, `F_m/p = P·U`, so `v(D) = 1 + v(P(l))` at these valid points.
* Split zeros would give `v(D) ≥ 4` somewhere on the disc, and an unramified pair would give `v(D) ≥ 3` on the whole disc.
* A constant valuation 2 means `P(r + X) = X² + bX + c` with `v(c) = 1` and `v(b) ≥ 1`. That polynomial is Eisenstein, so the two zeros are a **ramified** conjugate pair at distance `p^(−1/2)` from `r`.
* At `p = 5` this also shows that `F_3(0) ≠ 0`. The exact zero `D(3) = 0` lies below the level-3 threshold, so it does not force a rational zero.

*Correction.* For `χ(Δ) = 0` the two zeros are congruent mod `p` and may be split, unramified or ramified; the valuations of `D` on the disc decide which. The title should read "the triangle class is split whenever `p` is not Wieferich and the root `l = 0` is simple (at `p = 5` it is ramified)". The body (4d) already has the correct qualifier, "unless the root is double".

**12. MINOR. The range of the (4e) valuation formula is wrong by one for even `n`** (THM-4750:194: "right side below `(m + l(p−1) + 1)/2`").

* `F_m ≡ A_m` mod `p^k` requires `n ≥ 2k − 1`, so the agreement precision is `⌊(n+1)/2⌋ = ⌈n/2⌉`.
* For even `n`, "below `(n+1)/2`" allows `1 + v(l − ℓ) = n/2`, where only `v(D) ≥ n/2` follows.
* Correct form: for `l ≡ ℓ (mod p)` (replacing the vague "near ℓ"), `v_p(D(n)) = 1 + v_p(l − ℓ)` whenever `1 + v_p(l − ℓ) < ⌈n/2⌉`. No such boundary case occurs in the b10 data.

**13. OK, with one MINOR. The statistics re-derive exactly from the outputs** (`audit/b4_parse_second_order.py`, `.out`; my own parser of the four `second_order_*.out` files).

* *Coverage.* 1228 primes, none missing in [3, 10^4]. That is 1226 primes `≥ 7` and 4906 zero classes.

| class | 0 roots | 1 root | 2 roots |
|---|---|---|---|
| 1 | 633 | 2 | 591 |
| 2 | 612 | 5 | 609 |
| 3 | 0 | 2 | 1224 |
| sporadic | 618 | 5 | 605 |

* *Special classes.* The Wieferich lines are as stated. The 13 linear classes are the 8 Wieferich ones plus (23: 16), (23: 2), (97: 24), (389: 2) and (3989: 2275). The single double root is (7: 2).
* *Odd `λ_C`.* `p = 7, 97, 389, 3989`.
* *Class 1 against class 2.* Counts 320 / 312 / 292 / 297, `χ² = 0.136`.
* *Split rate by `ε`.* 0.483 at `ε = +1` and 0.501 at `ε = −1` (`z ≈ 1.05`, not significant).
* **MINOR (THM-4750:212).** "5 with `X(m−2) ≡ 0`, where about 5.5 are expected" is wrong. Class 3 can never be linear, because `X(1) = 4`. The 1/p model gives **4.02** expected linear classes. For double roots it gives 5.47 (all classes), and `P(≤ 1) = 0.027` matches the stated "about 3%".

### Statement 5: Camion limit and tower height (checklist 5)

**14. OK, with an upgrade available. Statements (a)–(d) are correct** (THM-4750:218–230).

* *Lipschitz bound.* `|C(x,i) − C(y,i)| ≤ p^⌊log_p i⌋ |x − y|`, by Vandermonde and `|C(h,j)| ≤ |h|/|j|`.
* *Signs and parametrisation.* `l_a = (p^a − 1)/(p−1) ≡ a (mod 2)`, so `ε^(l_a) = ε^a` and `A_1(l_a) = ε^a D(p^a)`. Relative to base `m = p`, `F_p(z) = ε F_1(z + 1)` and `z_∞ = −p/(p−1) = p + p² + … ≡ p (mod p²)`. So (d) and "`h ≥ 3` iff `ℓ ≡ p (mod p²)`" (for a simple root) are right.
* *Height.* (c) follows as stated.
* **Upgrade.** Since `⌈i/2⌉ ≥ 1 + ⌊log_p i⌋` for all `i ≥ 1`, `F_1` is a contraction: `|F_1(x) − F_1(y)| ≤ p^(−1)|x − y|`. Because `p^a ≥ 2a + 1`, this gives **`ε^a D(p^a) ≡ Λ_p (mod p^(a+1))`**.
* *Support for the upgrade.* My values of `Λ_p mod p³` (finding 16) equal the session's direct `ε²D(p²) mod p³` (`tower_height_direct_checks.out`) for all of `p = 7, 11, 13, 17, 19, 23`.

**15. MINOR. Statement 5(e) needs `p ≥ 2K − 1` (THM-4750:224), and less data suffices.**

* For `p < 2K − 1` the base `m = p` is below the level-`K` threshold, so the stated index bound is wrong; for example `p = 3`, `K = 3` needs `n ≤ 13`, not `11`. The division by `i!` also needs `p > 2K − 2`. `tower_height.c` itself requires `p > 2K`.
* *Improvement.* `v(C(z,i)) ≥ 1` at `z = −p/(p−1)` for `1 ≤ i < p`, so for `p ≥ 2K − 1`, `Λ_p mod p^K` needs only:
  * `D(p) mod p^K`;
  * `D(p + j(p−1)) mod p^(K−1)` for `1 ≤ j ≤ 2K − 4`.
* For `K = 3` this is **`Λ_p ≡ ε[D(p) + a_1 z + a_2 C(z,2)] (mod p³)`**, with `a_1, a_2` needed only mod `p²`.

**16. OK. 22409 is Camion–Wieferich and `h(22409) = 2`, independently** (`audit/b5_cw_tower_direct.c`, `.out`).

* *Method.* `s mod p³` to `n = p` (u128) and `s mod p²` to `3p` (u64), both from the Hurwitz inverse. The two runs are cross-checked against each other on `n ≤ p`. Total run time 5 s.
* *Values.* `D(22409) ≡ 6885662909072 (mod p³)`, so `v_22409(D(22409)) = 2`. `F_p(z) ≡ 5060801545918`, so **`v(Λ_22409) = 2` and `h(22409) = 2`**.
* *Simple root.* `c_1 − c_2/2 = 18775 ≢ 0`, so the root `l = 0` is simple.
* *Other primes.* `h(13) = 2`, and `h = 1` at 5, 7, 11, 17, 19, 23, 29, 1093 and 3511. The Camion quotients for `p ≤ 29` equal `cw_fast_validation_5_10000.out`.

**17. OK. The `cw_fast.c` identity holds, and the whole search reproduces with an independent algorithm.**

* *Derivation of the identity.*
  * `s(p) = −p![x^p] 1/T` and `T = T_<p + t_p x^p + …`, so `[x^p] 1/T = u_p − t_p`.
  * Hence `s(p) = 2^C(p,2) − p·(p−1)!·u_p ≡ 2^C(p,2) + p·u_p (mod p²)`, by Wilson.
  * `2^C(p,2) = (ε + p e′)^p ≡ ε (mod p²)`.
  * So `D(p)/p ≡ u_p − e′`, with `u_p` the reduction of a `p`-integral number. That is exactly what `cw_fast.c:29–30` computes.
* *Search.* `audit/b6_cw_search_flint.py` uses FLINT's Newton power-series inversion, an `O(p log p)` algorithm per prime, with coefficients built through discrete logarithms. It covered all 9590 primes `5 ≤ p ≤ 10^5` in 168 s and found **{13, 22409}**.
* *Validation of the search code.*
  * against the direct mod-`p²` recurrence for `p ≤ 400`;
  * against a plain-loop coefficient builder for `p ≤ 3000`;
  * against the session's printed `q_C(5..37)`.
* The session's ranges [5, 10^4], [10001, 63000], [63001, 79400], [79401, 90850] and [90851, 100000] have no gaps.
* `cw_fast.c` has no overflow: the u64 accumulator is `< p³` for `p < 2.6·10^6`.

**18. MINOR. "Two independent codes" overstates the search** (THM-4750:16).

* Only `cw_fast.c` searched (10^4, 10^5]. `tower_height.c` confirmed the hit 22409, not the absence of others.
* The claim is now true by finding 17. The wording should name `audit/b6`.
* The same applies to HYP-9358:59, "end of the feasible range for `O(p²)` methods". That is literally true, but a quasi-linear method makes `10^6` routine.

### Statement 6: the prime 2 (checklist 6)

**19. OK. 6(a) is PROVED correctly** (THM-4750:245–251).

* `f(n) = Σ 1/|Aut|` by orbit–stabilizer, and Moon's odd-order theorem makes each term a 2-adic unit.
* `T ≡ Φ_3 (mod 2)` coefficientwise in `Z_(2)[[x]]`, and `1/Φ_3 ≡ (1+x)/(1+x³)`.
* The OEIS data agree with the law up to `n = 19`, not just 18: A051337 odd iff `n ≢ 2 (mod 3)` for `n = 1..19`, and A000568 even for `n = 3..19`.
* OEIS has no parity comment on either sequence.

**20. OK. 6(b) and 6(c) reproduce** (`audit/b7_two_adic.py`, `.out`). These are exact big-integer values of `s(n)` for `n ≤ 1500` (Hurwitz inverse; 387 s):

* the parity law holds for `n ≤ 1500`;
* `v_2(f(n)) = 2 + v_2(n − 3250)` holds for all 499 values `n ≡ 2 (mod 3)`, `5 ≤ n ≤ 1500`;
* the exact data alone determine `n* ≡ 3250 mod 2^12`, as in `two_adic_zero_1500.out`;
* `f(n + 3·2^j) ≡ f(n) (mod 2^(j+2))` fails only at `n ≤ 3` for `j ≤ 7`, and at `n ≤ 4` for `j = 8` (tested to `n ≈ 1500 − 3·2^j`).

**21. MINOR. "It holds mod `2^(j+3)` on the class `n ≡ 1 (mod 3)`" is false as stated** (THM-4750:253).

* It fails at `n = 1` for every `j ≥ 4`, and at `n = 4` for `j ≥ 7`.
* The session's evidence was two values of `n` per `j`, both near `10^5`.
* Corrected statement: it holds for all `n ≡ 1 (mod 3)` with `7 ≤ n ≤ 1500 − 3·2^j` and `j ≤ 8`.

**22. MAJOR (M2). The rational-reconstruction evidence contradicts the session's own output.**

*Where.*

* THM-4750:261: "Rational reconstruction at `2^64` finds no fraction with numerator and denominator below about `2^31`."
* HYP-9361:13: "no rational with numerator and denominator below about `2^31` matches `n*` mod `2^64`."
* THM-4750:19 (status): "not a rational with numerator and denominator below about `2^31`."

*Contradiction.* `two_adic_zero_fast.out` prints `rational reconstruction: 765490658 / 1598160473  check: True`. That fraction has numerator about `2^29.5` and denominator about `2^30.6`, both below `2^31` and inside the uniqueness box `√(2^63)`.

*Audit result* (`audit/b7`). A Kitamasa evaluation, a different algorithm from the session's companion matrix, at 160-bit precision gives:

    n* ≡ 1090019148193376184792522869885760632192178 (mod 2^140).

* This agrees with the session's 64 bits and with `3250 mod 4096`.
* The law holds at 300 random huge `n`.
* The candidate `765490658/1598160473` **fails** mod `2^80`.
* The `2^128` reconstruction (height about `2^64`) also fails mod `2^140`.
* At `2^140` no valid reconstruction exists. So **`n*` is not a rational with numerator and denominator below `2^69.5`**.

The status clause is therefore true, but only because of this audit's check. HYP-9361:13 and THM-4750:261 must be rewritten.

**23. OK+. The HEURISTIC log formula is correct to 90 bits, and 6(b)/(c) are provable** (`audit/b9_nstar_log_formula.py`, `.out`).

*Numerical check.* Arithmetic was done in `Z/2^150[ζ_3]`.

* `ρ ≡ ζ (mod 2)` is a root of `T` mod `2^150`, found by Newton iteration (`T′(ζ) ≡ 1` is a unit).
* `f(n) = Tr(cρ^(−n))` mod `2^100` for all tested `n ∈ [1500, 2200]`, with `c = 1/(ρT′(ρ))`.
* The closed form of HYP-9361:81, `n* = log(−cζ/σ(cζ)) / log(w/σ(w))` with `ρ = ζw`, agrees with `n*` in **all 90 available bits**. The valuations are `v_2(log(w/σw)) = 2` and `v_2(n*) = 1`.

*Proof sketch for the upgrade.*

* `T` is entire on `C_2`.
* By its Newton polygon, its only zeros on `|x| = 1` are `ρ` and `σρ`, and all others have `|x| ≥ 2^(3/2)`.
* So `f(n) = Tr(cρ^(−n)) + r_n` with `v_2(r_n) ≥ (3/2)n − C`.
* On `n ≡ 2 (mod 3)`, `Tr(cρ^(−n)) = cζ w^(−n)(1 − v^(n − n*))` with `v = w/σw`, and `v_2(1 − v^t) = 2 + v_2(t)`.
* This proves that the 2-adic interpolation of `f` on the class is `n ↦ Tr(cζ w^(−n))`, that its unique zero is the closed form, and that `v_2(f(n)) = 2 + v_2(n − n*)` whenever `2 + v_2(n − n*)` is below the remainder's decay rate.
* The same argument gives the quasi-period in 6(b).

The HYP-9361(e) reading can be relabelled FINITE-EXACT now, and PROVED once the remainder constant is written out.

### Statement 7: r-uniform hypergraphs (checklist 7)

**24. OK. Statement 7 is proved.**

* *Vandermonde period.* `t = ord_p(2)` divides `P = p − 1`, and `gcd(t, i!) = 1` gives `t | C(P, i)`.
* *Fourier expansion.* The period holds from `n = 0`, so the Fourier expansion over `μ_(p−1)` is valid, with multiplier 1.
* *Filtration.* `d ≤ rj` at grade `j`.
* *Unit nodes.* `v ≥ max(j, i − (r−1)j) ≥ ⌈i/r⌉`.
* *Nonunit nodes.* They vanish for `n ≥ r(k−1) + 1`, and the threshold is `r(k−1) + 2` after the shift `Δb = ΔT ∗ T^(−1)`.
* *Class 1.* `b_3(p) ≡ 2^((p−1)(p−2)/6) − 1 (mod p)`, because rotations act freely except on unions of the `(p−1)(p−2)/6` free orbits of 3-sets, and any nonempty union is connected. The right side is `0` because `t | (p−1)(p−2)/6` when `gcd(t, 6) = 1`.

**25. OK, with one MINOR. Independent replication** (`audit/b8_hyper3.c`, `.out`; computed through `Δb = ΔT ∗ T^(−1)`).

* *Reproduced exactly.* 60 primes, 60444 classes, 58 not genuinely cubic. The table matches (class 1: 22/25/1/12; class 2: 0/31/0/29; sporadic: 21/21/0/14), and the class-1 formula holds with 0 failures.
* *Thresholds are sufficient, not sharp* (FINITE-EXACT). At `k = 2` the cubic law fails exactly at the bases `m ∈ {1, 3}` at every one of the 60 primes, against the proved threshold 5. At `k = 3` it fails exactly at `m ≤ 6` (15 primes `≤ 400`), against the proved 8, and speed `1/3` is attained.
* **MINOR (THM-4750:18, 291; HYP-9361:63).**
  * The class-2 "forced root" is real: the cubic vanishes at `l ≡ −1` (the point `n = 2`) in 60/60 cases.
  * But its mechanism is that base `m = 2` already obeys the cubic law (60/60). That is outside the proved threshold `m ≥ 5`, so "exact zero `b_3(2) = 0`: forced root" is a FINITE-EXACT observation, not a consequence of the theorem.
  * The tournament analogue confirms that a threshold matters. `D(1) = D(2) = 0` lie below the `k = 2` threshold and force nothing.

### Status labels, wording and provenance (checklist 8)

**26. MINOR. HYP-9361:9 says "split in all 1224 non-Wieferich cases (as THM-4750 (4d) forces)".**

* (4d) forces a split only when the root `l = 0` is simple; finding 11 shows it can be ramified.
* The data contain no class-3 double root for `7 ≤ p ≤ 10^4`.

**27. MINOR. HYP-9358:12–14 cites the wrong file, and the dispersion is understated.**

* The joint sporadic histogram `423, 473, 250, 67, 11, 2` is cited to `camion_irregular_10000.out`. That file contains only the separate odd and even histograms. The joint one is derivable from `second_order_*.out`, and `audit/b4` reproduces it exactly: mean 1.0016, variance 0.8744.
* "Slightly underdispersed" understates the case. Against the binomial model `N_p ~ Bin(p−4, 1/p)` of HYP-9358(b), the dispersion statistic is `z = −2.47` (two-sided about 1.4%).
* The quoted `χ² ≈ 10.7 on 5 df` depends on the binning, and should use `5 − 1 − 1 = 3` df once the mean is fitted. Bins 0, 1, 2, 3, `≥ 4` give `χ² = 10.86`.
* This is mild evidence against the exact Poisson form of HYP-9358(b), and should be recorded as such.

**28. MINOR. RES:358 says the half-speed bound "is sharp at every tested non-Wieferich prime (`camion_speed.out`)".**

* That is false at `p = 3, 5, 7`, where the same output shows extra divisibility. THM-4750 §2 states the exception correctly.
* *OK+.* The extra divisibility is exactly `v_p((2k−1)!!)` at `i ∈ {2k−1, 2k}` wherever the cap allows:
  * `p = 3`: +1, +1, +1, +3 at `k = 2, 3, 4, 5` (`v_3(9!!) = 3`);
  * `p = 5`: +1 at `k = 3, 4, 5`;
  * `p = 7`: +1 at `k = 4`.
* This sharpens the double-factorial mechanism of HYP-9358(e2) quantitatively.

**29. MINOR. Wording.**

* "`p`-adic analytic" (THM-4750:49, 299) is proved only for `p ≥ 5`, as a restricted power series. For `p = 3` the bound gives continuity, and local analyticity on discs of radius `1/3` by Amice's criterion.
* "agree to `n = 18`" (THM-4750:252) should read 19.

### Script review (checklist 9)

**30. OK. No arithmetic bugs that affect any reported number.** Overflow margins checked:

* `camion_second_order.c`: unreduced u128 accumulation `< 3·10^4 · 10^16`.
* `camion_speed.c`: reduction every 64 terms, `Q ≤ 1093^4`.
* `tower_height.c`: every 16 terms.
* `cw_fast.c`, `camion_irregular.c`: u64 `< p³`.
* `hyper3_second_order.c`: u128.

Other points checked:

* The `ε` handling (`A[l]` negated for odd `l` when `ε = −1`) is correct in every C file.
* Class ranges cover each residue class once:
  * `m ∈ [3, p+1]`, with class 1 at `m = p` and class 2 at `m = p + 1`;
  * `m ∈ [2K−1, 2K+p−3]`;
  * `m ∈ [6, p+4]` for r = 3.
* The discriminant `(c_1 − c_2/2)² − 2c_0c_2` is right.
* `tower_height.c` evaluates `C(z,i)` correctly for `i < p`.

**31. MINOR. Cosmetic and provenance issues in the scripts.**

* `two_adic_zero_fast.out` is not byte-identical to a rerun of `two_adic_zero_fast.py`: the committed file repeats the line `shift 3*2^8` twice, and the script prints it once.
* `camion_speed_base.c:1` still carries the header `camion_speed.c`.
* `tower_height.c:34` contains dead code (`Dpa`).
* `two_adic_face.py` tests the refuted guess, which is documented as such.

## 2. Independent checks run

| id | script (in `audit/`) | method | result |
|---|---|---|---|
| b1 | `b1_half_speed_exact.py` | exact integers, Hurwitz inverse; 28 pairs `(p,k)` | statement 1 and the higher Kummer congruence: 0 violations in 34400 checks at `m ≥ 2k−1`; failure sets below; sharpness for `p ≥ 7` |
| b2 | `b2_full_speed.c` | mod `p^K`, Hurwitz inverse | full speed from `m ≥ K` at 1093 (`K ≤ 5`) and 3511 (`K ≤ 4`); base 3 at 11, base 7 at 5; threshold sharp; non-Wieferich controls |
| b3 | `b3_curvature.py` | numpy; `X = S″ + 2S′²T`; brute force | curvature formula in 13806 classes (`p ≤ 400`), 0 failures; direct law 20928 checks, 0 failures; `X` meaning (`m ≤ 6`); base-derivative identity; 1091/1097 exceptions |
| b4 | `b4_parse_second_order.py` | parser of the session outputs | root table, linear/double lists, `λ_C`, contingency, `ε` split rates, sporadic histogram and dispersion (`z = −2.47`); expected linear 4.02 |
| b5 | `b5_cw_tower_direct.c` | `s` mod `p³` (u128) and mod `p²` (u64); 3-term `Λ_p` | `v(D(22409)) = 2`, `h(22409) = 2`, `h(13) = 2`, `h = 1` at the other tested primes |
| b6 | `b6_cw_search_flint.py` | FLINT Newton series inversion, `O(p log p)` per prime | CW primes in [5, 10^5] = {13, 22409}; identity validated against direct recurrences |
| b7 | `b7_two_adic.py` | exact `s(n)`, `n ≤ 1500`; Kitamasa at 160 bits | parity law; `v_2(f) = 2 + v_2(n − 3250)` (499 values); quasi-period; `n*` mod `2^140`; rational candidate refuted; no rational of height `≤ 2^69.5` |
| b8 | `b8_hyper3.c` | rooted identity `Δb = ΔT ∗ T^(−1)` | r = 3 table reproduced; observed thresholds `{1,3}` (`k = 2`) and `≤ 6` (`k = 3`); class-2 root at `n = 2`; speed 1/3 attained |
| b9 | `b9_nstar_log_formula.py` | unit root of `T` in `Z_2[ζ_3]`, 2-adic logarithms | `f = Tr(cρ^(−n))` at large `n`; closed form equals `n*` to 90 bits |
| b10 | `b10_valuation_law.py` | exact integers, `n ≤ 800` | (4e) valuation law with a consistent `ℓ` at `p = 7, 11, 13`; no even-`n` boundary case met |
| b11 | `b11_double_root_ramified.py` | exact integers | ramified second-order zeros at `p = 5` (triangle) and `p = 7` (class 2): counterexamples M1 |

**Reruns of session code.** All were compiled and run in the scratch directory.

* `camion_second_order 3 4000` is byte-identical to the committed output.
* `camion_second_order 3 59 60` gives 17158 checks with 0 failures, matching the claim, but it was never committed.
* `two_adic_zero_fast.py` reproduces its output apart from the duplicated line.

## 3. Literature

No source found states the higher Kummer (Mahler half-speed) congruences for labelled counts, the Fermat-quotient curvature formula, the Wieferich collapse, the Camion limit, the strong-class parity law, or the 2-adic zero. All of them appear to be new. The parity law is a two-line consequence of Moon's theorem and may be folklore. The nearest prior art:

* **A. P. Mani and R. J. Stones, "The number of labeled connected graphs modulo prime powers", SIAM J. Discrete Math. 30(2) (2016) 1046–1057** (doi:10.1137/15M1024615).
  * As quoted by Kotek–Makowsky (Prop. 1.11), they prove `C_n(b) ≡ b^(φ(p^k)/2) C_(n−φ(p^k))(b) (mod p^k)` for `n ≥ p^k` and `b ≢ 1 (mod p)`.
  * This is the first-order quasi-period with multiplier `(b/p)`, for connected graphs. It contains no higher differences, curvature or interpolation.
  * Same authors: **"Congruences for weighted number of labeled forests", Integers 16 (2016) #A17**, which has first-order recurrence congruences relating `f_(n+p^k)` to `f_n`.
* **T. Kotek and J. A. Makowsky, "On the Tutte and matching polynomials for complete graphs", Fundamenta Informaticae (2022)**, arXiv:2112.06581.
  * They prove ultimate periodicity of `T(K_n; a, b)` mod `μ` via Specker–Blatter, and state the Mani–Stones conjecture.
* **E. Specker and C. Blatter (1981/1984)**, and **E. Fischer and J. A. Makowsky, "Extensions and limits of the Specker–Blatter theorem", J. Symbolic Logic (2024)**, arXiv:2206.12135. These give MC-finiteness (periodicity mod every `m`) for MSOL-definable classes, a qualitative statement only.
* **J. W. Moon, "Tournaments with a given automorphism group", Canad. J. Math. 16 (1964) 485–489**, and *Topics on Tournaments* (1968): the automorphism group of a tournament has odd order.
* **G. Royle, C. Praeger, S. Glasby, S. Freedman and A. Devillers, "Tournaments and even graphs are equinumerous", J. Algebraic Combin. (2023)**, arXiv:2204.01947: Cauchy–Frobenius counting over odd-order permutations. They do not state parity.
* **A. O'Desky and D. H. Richman, "Derangements and the p-adic incomplete gamma function", Trans. AMS 376 (2023)**, arXiv:2012.04615: a `p`-adic continuity criterion for combinatorial sequences (derangement-type), not graph counts. This is the closest "p-adic interpolation of a counting sequence" prior art found.
* **J. Archer, I. Gessel, C. Graves and X. Liang, "Counting acyclic and strong digraphs by descents"**, arXiv:1909.01550: generating functions for strong tournaments and digraphs, not congruences.
* **Background.** Kummer congruences and Kubota–Leopoldt interpolation (classical); Mahler's theorem and Amice's analyticity criterion; the Gessel umbral-calculus literature on Kummer congruences for EGF coefficients. None of it covers `2^C(n,2)`-type counts.

**Sources consulted:**

* [Mani–Stones, SIAM J. Discrete Math. (2016)](https://epubs.siam.org/doi/10.1137/15M1024615)
* [Mani–Stones forests, Integers 16 (2016) A17](https://emis.univie.ac.at/journals/INTEGERS/papers/q17/q17.pdf)
* [Kotek–Makowsky, arXiv:2112.06581](https://arxiv.org/abs/2112.06581)
* [Kotek–Makowsky, Fundamenta Informaticae PDF](https://fi.episciences.org/8839/pdf)
* [Fischer–Makowsky, arXiv:2206.12135](https://ar5iv.arxiv.org/html/2206.12135)
* [Moon 1964, Canad. J. Math.](https://core-cms.cambridgecore.org/core/services/aop-cambridge-core/content/view/CF0B648313DD8BE07D726932193D9044/S0008414X00035598a.pdf/tournaments-with-a-given-automorphism-group.pdf)
* [Moon, Topics on Tournaments (Project Gutenberg)](https://www.gutenberg.org/files/42833/42833-pdf.pdf)
* [Royle et al., arXiv:2204.01947](https://arxiv.org/pdf/2204.01947)
* [O'Desky–Richman, arXiv:2012.04615](https://arxiv.org/abs/2012.04615)
* [Archer–Gessel–Graves–Liang, arXiv:1909.01550](https://arxiv.org/pdf/1909.01550)
* OEIS [A051337](https://oeis.org/A051337), [A000568](https://oeis.org/A000568) and [A054946](https://oeis.org/A054946), queried through the JSON API: no parity or congruence comments.

## 4. Required corrections (summary for the session)

1. **M1.**
   * Title of (4): replace "(split in Z_p or an inert conjugate pair) … the triangle class m = 3 is always split" by the trichotomy, with `χ(Δ) = 0` covering split, unramified or ramified pairs congruent mod `p`, and give the qualified triangle statement.
   * Add the counterexamples `p = 5` (`m = 3`) and `p = 7` (`m = 8`), both ramified (`audit/b11`).
2. **M2.** Rewrite THM-4750:261 and HYP-9361:13 with the audit result: the `2^64` reconstruction `765490658/1598160473` fails mod `2^80`, and there is no rational of height `≤ 2^69.5` (`audit/b7`).
3. **Finding 12.** In (4e), state the range as `1 + v_p(l − ℓ) < ⌈n/2⌉` with `l ≡ ℓ (mod p)`.
4. **Finding 15.** In 5(e), add `p ≥ 2K − 1`. Optionally add the cheaper formula and the sharpened congruence mod `p^(a+1)` (finding 14).
5. **Finding 5.** In the general-base paragraph, add the converse argument through `X_a(1)` and `X_a(2)`.
6. **Finding 13.** Expected linear classes: about 4.0, not 5.5.
7. **Finding 21.** In 6(b), correct the `mod 2^(j+3)` sub-claim (exceptions at `n = 1, 4`).
8. **Finding 25.** In §6b, label the class-2 forced-root mechanism FINITE-EXACT and record the observed thresholds.
9. **Findings 9, 18, 27, 31 (provenance).**
   * Commit the `VMAX = 60` output.
   * Cite `audit/b6` for the second independent search code.
   * Cite `second_order_*.out` / `audit/b4` for the joint histogram.
   * Regenerate `two_adic_zero_fast.out`.
10. **Findings 26, 28, 29 (wording).** HYP-9361:9 ("forces"), RES:358 ("sharp at every tested prime"), "analytic" for `p = 3`, and "n = 18".
