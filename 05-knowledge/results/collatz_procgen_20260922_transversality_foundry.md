# Transversality foundry: the divergence target is Lagarias's Periodicity Conjecture, a typed 2-versus-3 catalogue, and a proof that no rational has a Sturmian parity vector

**Status.**
* **AUDIT (2026-09-23, session orchestrator, independent code path):** the proofs of Theorem R, Lemma J, Lemma H and Theorem S were re-derived line by line, and no gap was found. [`..._transversality_audit.py`](../../04-computation/experiments/collatz_procgen_20260922_transversality_audit.py) recomputes `Phi_T(s) mod 2^N` from the closed form `-r sum_l 2^(d_l) q^(-l)` rather than the affine recursion, and sums every approximant `Phi(u v^inf)` exactly as a geometric series. On 9 rows at `N = 20000` it confirms, at every convergent, the parity prefix, Lemma J (period, gaps, first defect), both lambda lower bounds and the exact isometry `v_2 = lambda`. Gains track `q_(n+1)/mu - (mu-1) q_n` up to `O(log)`. The rows cover `3x+1` at critical (three intercepts), supercritical and subcritical slopes, `3x-1`, and `5x+1` at `mu = 1.435` and `1.642`. A `5x+1` control with `mu(mu-1) = 2.29 > phi` shows negative gains at small partial quotients, as the theorem's scope predicts. Output: [`..._transversality_audit.out`](collatz_procgen_20260922_transversality_audit.out).
* **PROVED** (this note; elementary; audited as above):
  * Theorem R, a 2-adic repetition criterion for the irrationality of Bernstein numbers. It holds for every affine 2-adic shift map.
  * Theorem S: no rational with odd denominator, positive or negative, has an eventually Sturmian parity vector under any `3x+r` map, for every irrational slope and every intercept. The same holds for `5x+1` at slopes below `0.804`, and for Mahler's map `a -> ceil(3a/2)`, so no Z-number has an eventually Sturmian carry word.
  * Proposition B: the Periodicity Conjecture holds on words of bounded critical discrepancy. It is assembled from the in-house discrepancy theorem; its positive-tail step uses that note's section 4 extension to `3n+b`, `b > 0`, which is written there as a short sketch.
  * The exact count `N_j = 2^(j-1)`.
  * The kappa formula `v_3(2^K - r) = [K = e mod 2](1 + v_3(K - kappa(r)))`.
* **FINITE-EXACT:**
  * Sturmian certificates at `N = 20000` bits for 58 words, plus one deep row at `N = 10^5` bits. For every word, no rational of height `<= 2^9999` (deep row: `2^49999`) equals its Bernstein number.
  * Lemma J checked on 580 cases; Theorem R's inequality and the exact isometry `v_2(x - P/Q) = lambda` checked on 17,909 pairs (genuine rational `x`, approximant `P/Q`) for `3x+1`, `5x+1` and Mahler's map, with no violation.
  * Erdős scan for `K <= 2*10^5`, base-5 analogue, ternary and binary digit profiles.
  * Loop-clock landings up to `s = 2*10^7`.
* **CITED:** the 35-item catalogue (`[P]` means a primary source was read; `[R]` means a named secondary source was used) and the formulations of the Periodicity Conjecture.
* **OPEN:**
  * Collatz divergence (T1) and Collatz cycles.
  * E-SCC Q1, and Q2 (T2).
  * Erdős's ternary problem (T3) and Mahler's Z-numbers (T4).
  * The Periodicity Conjecture (PC).
  * The existence of a divergent `5x+1` orbit.
  * Furstenberg's `x2 x3` measure conjecture.
  * Every hypothesis candidate of section 6.
* **REFUTED** (generated statements, each with an elementary counterexample): B1 and B2 (3-adic-window Erdős), D1, E1 (hostile-free loop-clock landings), H0 (lambda-uniform Erdős), M3 (PC for Mahler's map), and the eventually-periodic rows.

Session `collatz-procgen-20260922`, transversality lane (machine `mac-mini`, 2026-09-22).

* **Generator:** [`collatz_procgen_20260922_transversality.py`](../../04-computation/experiments/collatz_procgen_20260922_transversality.py).
* **Helpers:** `collatz_procgen_20260922_transversality_{catalogue,engines,sturmian}.py`.
* **Output:** [`collatz_procgen_20260922_transversality.out`](collatz_procgen_20260922_transversality.out).
* **Literature:** three sub-lanes, whose notes are in `scratch/procgen_transversality/lit_{A,B,C}/`. Downloaded texts are git-ignored. No bot-check or paywall was bypassed.

## 0. Answers in brief

1. **The target.** The synthesis's "p-adic irrationality of the Bernstein series" is exactly **Lagarias's Periodicity Conjecture**, `Q_inf(Q_2) = Q_2` (Lagarias 1985 section 2.8; `Phi(Q cap Z_2) = Q cap Z_2` in Bernstein–Lagarias 1996).
   * The synthesis's own sentence, "such a series is never a positive integer", is the weaker, positive-integer, `k=1` instance. That is Bernstein's 1994 non-iterative form of "no divergent positive orbit" (T1).
2. **Catalogue** (section 2, 35 items).
   * Every counting, density or dimension theorem is DEFECT-blind.
   * Every every-element 2-versus-3 theorem excludes only a **zero-entropy** exceptional class: few nonzero digits, bounded p-adic closeness, runs at the ends of an expansion, or balanced words.
   * The only proved every-orbit constraints on **positive-entropy** classes of parity words use growth or capacity, not transversality:
     * Monks–Yazinski: subcritical density;
     * the in-house theorem: bounded critical discrepancy.
   * The untouched core of T1, PC and T4 is the class of supercritical, positive-entropy, non-repetitive words.
3. **Generator** (section 3). 49 statements over 8 objects, 5 properties and 3 quantifiers. The minimal unblocked statements are:
   * T1: the HARD-class form (G[HARD]).
   * T2: the steering statement E3.
   * T3: the Erdős statement itself, A1. No 3-adic-window statement can serve (B1 and B2 are FALSE).
   * T4: H4, M1, or H3 ("spread > 1/2", which the drift control `5/2` shows to be DRIFT-sensitive).
4. **Sturmian test** (section 4).
   * The literature knew only slopes below `log_3 2`.
   * **Theorem S settles every slope and every intercept for all rationals**, by an elementary 2-adic Liouville argument.
   * Computations at `2*10^4` and `10^5` bits agree, with repetition gains growing like `q_(n+1)/mu`.
5. **Cross-problem matrix** (section 5).
   * One statement type, "a Bernstein number is not a positive integer (not rational) on a word class", serves T1, PC and T4. Theorem S proves all three on the Sturmian class at once.
   * The `2^K`-in-`Z_3` family serves T2 and T3, but with incompatible windows.
   * The `5x+1` question needs the opposite transfer.
   * Furstenberg's conjecture serves none of these problems.

## 1. The target, named (task 1)

**Setting.** `T(x) = x/2` (x even) or `(3x+1)/2` (x odd), on `Z_2`. The parity vector is `v(x) = (T^i(x) mod 2)_(i>=0)`.

* **Lagarias 1985 [P], (2.33) and Theorem L:** `Q_inf(x) = sum v_i(x) 2^i` is a continuous, measure-preserving bijection of `Z_2`.
* **Bernstein–Lagarias 1996 [P], App. A, Cor. A.3:** such solenoidal bijections are exactly the 2-adic isometries.
* **Bernstein 1994 [P]:** the inverse `Phi = Q_inf^(-1)` has the non-iterative form
  `Phi(sum_l 2^(d_l)) = -sum_(l>=1) 3^(-l) 2^(d_l)` for `0 <= d_1 < d_2 < ...` (also B–L (1.6)).
  His Theorems 1–3 show that the 3x+1 conjecture is equivalent to `Z^+ subset Phi((1/3)Z)`.

**The Periodicity Conjecture, exact statement.**
* Lagarias 1985, section 2.8 [P], names it but does not number it: `Q_inf(Q_2) = Q_2`. Here `Q_2` is the set of rationals with odd denominator, equivalently the 2-adic integers with eventually periodic expansion.
* Bernstein 1994, Cor. 1 [P], proved the easy half, `Phi(Q_2) subset Q_2`: every eventually periodic (EP) word is the parity vector of a rational, given by the cycle formula. He wrote that this is half of the conjecture.
* What remains is: **for every non-EP word `w`, the Bernstein number `Phi(w) = -sum_l 2^(d_l) 3^(-l)` is 2-adically irrational.** This is word for word the synthesis's target.

**Equivalent forms (CITED).**
* **B–L 1996 [P]:** PC holds iff no `3x+k` map with `k = +-1 (mod 6)` has a divergent trajectory on `Z` (either sign).
  * They derive this from Lagarias 1990, Cor. 2.1b, which we did not read ([R]: bot-blocked).
  * PC also holds iff `Phi^k(Q cap Z_2) = Q cap Z_2`, for any fixed `k >= 1`.
* **Monks–Yazinski 2004 [P], Thm 2.1:** PC holds iff the Autoconjugacy Conjecture `Omega(Q_odd) subset Q_odd` holds, iff no rational 2-adic integer has a divergent `T`-orbit.
  * `Omega` complements the parity vector.
  * Also, `Omega(Z^+) subset Q_odd` iff there is no divergent positive orbit.
* **Rozier 2019 [P]:** Conj. 1 and 5 restate PC.
* **Akin 2004 [R]** (via the zbMATH review and Lagarias's bibliography II) poses the analogue for `tau_a(x) = x/2, (ax+1)/2`:
  * it holds for `a = +-1` and fails for every non-integer odd rational `a`;
  * Akin gives a heuristic that it holds for `a = +-3` and fails for odd `|a| >= 5`.
  * That heuristic is the foundry's DRIFT barrier.

**Partial results before this note.**
1. **The easy half** (Bernstein 1994).
2. **The density constraint.** A rational with a non-EP parity vector has `liminf a_s/s >= log_3 2`, where `a_s` counts the 1s among the first `s` letters.
   * Stated without proof for integers in Lagarias 1985 (2.31).
   * Proved for rationals in Monks–Yazinski 2004, Thm 2.7(b) [P].
   * We call it Lemma L below and re-prove it.
3. **No arithmetic class of rationals is known to satisfy PC.** The case of `Z` is the Divergent Trajectories Conjecture.
4. **Sturmian words.**
   * López–Stoll 2009 [P] (Integers 9, A13) give a 2-adic series for `Phi(1c_alpha)` and a generalized continued fraction for `-1/Phi(1c_alpha)`. They state that it is unknown whether any aperiodic word has an eventually periodic `Phi`.
   * López–Stoll 2021 (arXiv 2101.12747, unrefereed) claim `liminf = ln 2/ln 3` for divergent rationals. The literature sub-lane found an apparent gap: the argument uses the real sum `Phi_R` where the 2-adic value is needed.
   * zbMATH searches for "Sturmian & Collatz" and "Beatty & Collatz" returned no hits.
   * So before this note the Sturmian case was settled **only for slopes below `log_3 2`**, as a consequence of the density constraint.
5. **Balanced cycles.** Knight 2026 [P] (Discrete Math. 349, Art. 114812; HAL preprint read): no positive-integer cycle other than `{1,2}` has a Christoffel (balanced periodic) parity word.
   * The negative cycle `{-5,-7,-10}`, with word `110`, is the exception. It comes from the Catalan unit gap `2^3 - 3^2 = -1`.
   * So this is a sheet-aware, zero-entropy cycle result: the periodic twin of Theorem S.
6. **In-house.** [collatz_guards_20260921_discrepancy](collatz_guards_20260921_discrepancy.md) proves that no positive-integer orbit has bounded critical discrepancy. This covers the critical-slope Sturmian halving words (section 4.4).

**Verdict.**
* "p-adic irrationality of the Bernstein series" **is** PC.
* "Such a series is never a positive integer" is **T1** in Bernstein's form: the `k = 1`, positive-integer instance of the B–L equivalence. PC implies T1; the converse is not known.

The synthesis (section 5, item 5) should cite PC (Lagarias 1985 section 2.8; B–L 1996; Monks–Yazinski 2004) and name T1 as its positive-integer instance. This note does not edit the synthesis.

## 2. Catalogue of PROVED 2-versus-3 transversality theorems (task 2)

**Typing rule.** The barrier names follow the [foundry](collatz_procgen_20260922_foundry.md) and the [atlas](collatz_procgen_20260922_barrier_atlas.md).

| barrier | B (blind) | O (overcomes) |
|---|---|---|
| SHEET | the conclusion also holds for the negated or `3n-1` analogue | it fails there |
| DRIFT | the theorem holds for every multiplicatively independent pair, e.g. `(2,5)` or `5/2`, where the Collatz drift `3 < 2^2` (resp. Mahler's threshold `p < q^2`) has the other sign | it uses `3 < 2^2` (resp. `p < q^2`) |
| DEFECT | counting, density, measure and dimension statements | every-element statements |
| INTEGRAL | the statement tolerates extra rational or periodic exceptions | it isolates them |
| UNIFORM | the proof runs unchanged over an explicit infinite family | per-instance: an explicit finite computation, or bounded complexity |

Kind: COUNT (a counting/density statement) or EVERY (every element or orbit). The `.out` (S2) prints all 35 items with full statements; the table below keeps the ones that were asked for, plus the key additions.

| item | exact statement (short) | kind | SH | DR | DE | IN | UN | access |
|---|---|---|---|---|---|---|---|---|
| Furstenberg 1967 | `p,q` mult. independent: a closed `xp, xq`-invariant `X subset R/Z` is finite (rational) or everything; `{p^m q^n alpha}` is dense for irrational `alpha`. The measure form is OPEN | EVERY | B | B | O | O | B | [R] via Rudolph 1990, Lindenstrauss 2005 Thm 1.1, zbMATH |
| Senge–Straus 1973 | `{n : s_a(n) <= c, s_b(n) <= c}` is finite for all `c` iff `log a/log b` is irrational (ineffective) | EVERY | B | B | O | - | B | [R] via Stewart 1980, Spiegelhofer, Dimitrov–Howe |
| Stewart 1980 | `L_a(n) + L_b(n) > loglog n/(logloglog n + C) - 1` for `n > 25` (effective); `3^n` has more than `log n/(loglog n + C_0) - 1` nonzero bits, `n > 4` | EVERY | B | B | O | - | B | [P] Thms 1–2 |
| Mahler 1957 | rational non-integer `alpha > 1`: `dist(alpha^n, Z) < theta^n` only finitely often. Ineffective `dist((3/2)^n, Z) > (3/4)^n`; effective `0.5803^k` (Zudilin 2007) | EVERY | B | B | O | O | B | [R] via Corvaja–Zannier [P], Zudilin [P] |
| Mahler 1968 | each `[g, g+1)` has at most one Z-number, in `[g, g+1/2)`; at most `x^0.7` of them up to `x` (Flatto: exponent `log_2(3/2)`) | COUNT | - | **O** | B | - | O | [P] reprint |
| Flatto–Lagarias–Pollington 1995 | for every `xi > 0` and coprime `p > q > 1`, `limsup - liminf {xi(p/q)^n} >= 1/p` (so at least `1/3` for `3/2`); `-3/2`: at least `11/27` (Lu–Zheng 2026) | EVERY | B | B | O | - | B | [R] via Lu–Zheng [P], Farhi [P] |
| Lagarias 2009 | the Thm 1.1 count; Thm 1.2 exceptional `lambda`; Thm 1.4 count; Thm 1.5 dimensions (details below the table) | COUNT | - | B | B | B | B | [P] |
| Narkiewicz 1980 | `#{n <= X : 2^n omits 2 in base 3} <= 1.62 X^(log_3 2)`; HKR 2015: `1.3 a^(log_3 2)` | COUNT | - | B | B | B | B | [R] via Lagarias 2009; HKR [P] |
| Dupuy–Weirich 2016 | Cesàro average over `n` of the frequency of a digit among the `m` lowest base-`q` digits of `p^n` tends to `1/q` | COUNT | - | B | B | B | B | [R] via Li–Zhao 2026 |
| Hochman–Shmerkin 2012 | `dim(X + sY) = min(1, dim X + dim Y)` for closed `x2`-, `x3`-invariant `X, Y` and `s != 0`; the measure form Thm 1.3 | COUNT | B | B | B | B | B | [P] |
| Shmerkin 2019; Wu 2019 | `upper-box-dim(A cap g(B)) <= max(0, dim A + dim B - 1)` for `xp`-, `xq`-invariant closed `A, B` and every invertible affine `g` | COUNT | B | B | B | B | B | [P] |
| Tao 2022 | Prop 1.17: `abs(E e(-2 pi i xi Syrac(Z/3^n)/3^n)) << n^(-A)` for `3` not dividing `xi`, with `Syrac = sum_(i<=n) 3^(i-1) 2^(-(a_1+..+a_i))`, `a_i` iid `Geom(2)`; Thm 1.3 holds in log density | COUNT | B | B (Prop 1.17); O (Thm 1.3, via GGM's `q < p^(p/(p-1))`) | B | B | B | [P] |
| Yu; Bugeaud–Laurent | explicit bounds for `v_p(prod alpha_i^(b_i) - 1)`. Consequence: `v_3(2^n - r) <= C(r) log n` for rational `r` not a power of 2 | EVERY | B | B | O | O | B | [R] via Palojärvi–Seppälä, Pink–Ziegler |
| Saye 2022 | every `16 <= n <= 2*3^45 = 5.9e21`: `2^n` contains each ternary digit, so Erdős holds for `8 < n <= 2*3^45` | FINITE | - | - | O (range) | - | O | [P] |
| Dimitrov–Howe 2021/25 | the only powers of 3 with at most 22 binary ones are `3^0..3^25`; for `x not in {0,2,8}`, `2^x` has a ternary 2 or at least 26 ternary ones | EVERY | - | B | O | - | O | [P] |
| Abram–Lagarias; ABL 2017 | `dim E(Z_3) <= log_3 phi = 0.438` | COUNT | - | B | B | B | B | [P] |
| Ren–Roettger 2025 | the run of 0s after the leading ternary digit of `2^n` is `O(log n)` (Baker) | EVERY | - | B | O | - | B | [P] |
| Monks–Yazinski 2004 | Thm 2.1 (PC equivalences); Thm 2.7(b): divergent rationals have `liminf kappa_n/n >= ln2/ln3` | EVERY | B | B | O | O | B | [P] |
| Knight 2026 | no positive-integer Christoffel cycle except `{1,2}`; `{-5,-7,-10}` exists | EVERY | **O** | - | O | O | O | [P] preprint |
| Calegari 2005 | `zeta_2(3)`, `zeta_3(3)` and `L_2(2, chi_4)` are irrational (p-adic Beukers method, criterion Lemma 2.2) | EVERY | - | - | O | - | O | [P]; refereed (IMRN) |
| in-house discrepancy | no positive-integer orbit has bounded critical discrepancy (section 4.4) | EVERY | B | **O (proof)** | O | O | B | [repo] |
| Theorem S (this note) | no rational has an eventually Sturmian parity vector (section 4.3) | EVERY | B | B | O | O | B | [repo] |

**Lagarias 2009 in detail.**
* Thm 1.1: for each `lambda > 0`, the number of `n <= X` for which the ternary expansion of `floor(lambda 2^n)` omits 2 is at most `25 X^0.9725`.
* Thm 1.2: uncountably many `lambda` omit 2 along an infinite sparse sequence of `n`.
* Thm 1.4: for nonzero `lambda in Z_3`, the corresponding count is at most `2 X^(log_3 2)`.
* Thm 1.5: `dim E^(1) = log_3 2` and `(1/2) log_3 2 <= dim E^(2) <= 1/2`.

**The further items**, in the `.out` only:
* Rudolph–Johnson, Host / Hochman–Shmerkin 2015, Cassels–Schmidt, Ellison 1971, Rhin / Simons–de Weger, BLMV 2009, Spiegelhofer 2023, Bugeaud–Corvaja–Zannier 2003, Koksma, Akin 2004, López–Stoll 2009, ADQZ 2001 / BHZ 2006, Dubickas–Mossinghoff 2009.
* Ellison 1971 [P]: `|2^x - 3^y| > 2^x e^(-x/10)`, except for listed small `x`.
* ADQZ 2001 [R via BHZ, P]: every Sturmian word begins with infinitely many squares.
* Dubickas–Mossinghoff 2009 [R via Andrieu–Eliahou–Vivion]: no Z-number below `2^57`.

**What the pattern says.**
1. **COUNT theorems are DEFECT-blind by construction.** A dimension-0 intersection (Shmerkin–Wu), a log-density-one set (Tao) or a `X^(log_3 2)` count (Narkiewicz) never excludes a single integer. Lagarias's Thm 1.2 even shows that the lambda-uniform version of Erdős's statement is FALSE (generated statement H0).
2. **Every-element 2-versus-3 theorems reach only zero-entropy classes.** Examples:
   * Senge–Straus and Stewart: few nonzero digits;
   * Yu: bounded p-adic closeness to one rational;
   * Ren–Roettger and LTE: runs at the two ends of an expansion;
   * Dimitrov–Howe: at most 25 ternary ones;
   * Knight and Theorem S: balanced words.

   The targets need the avoidance of positive-dimensional Cantor sets:
   * dimension `0.95` for Collatz;
   * `log_3 2` for Erdős;
   * `log_2(3/2)` for Mahler.
3. **The only every-orbit results on positive-entropy classes of parity words use growth or capacity, not transversality.**
   * Monks–Yazinski 2.7(b) handles subcritical words, whose orbit would shrink.
   * The in-house theorem handles bounded critical discrepancy, whose orbit would grow at most linearly and so fill a positive density of integers.
   * Neither says anything about **supercritical** words, where the orbit grows exponentially.
4. **DRIFT.** Almost every item is uniform over multiplicatively independent pairs. The exceptions are:
   * Mahler 1968's count;
   * Tao's Theorem 1.3 (not Prop 1.17);
   * Akin's heuristic;
   * the proof of the in-house theorem.

   For T4 the drift control is **proved**: for `p/q = 5/2 > q^2`, Z-numbers exist (Tijdeman 1972, Flatto 1992; [R] via Andrieu–Eliahou–Vivion). For Collatz, the `5x+1` control is only conjectural.
5. **SHEET.** Only Knight's cycle theorem is sheet-sensitive, through the Catalan gaps `|2^k - 3^x| = 1`. That matches the synthesis: the sign enters through the unit-gap clocks.

## 3. The generator (task 3)

**Grammar.** A statement is `S = (object, property, quantifier)`, with a parameter.

* **Objects:**
  * `2^K` (ternary expansion, or `Z_3`);
  * `2^K mod 3^D` (the 3-adic window);
  * `2^(K0(s)) w` along the loop clocks `K0(s) = floor((s+1) log_2 3)`;
  * `3^a` (binary, or `Z_2`);
  * Bernstein numbers `Phi_T(w)` for `w` in a word class `W`:
    * EP: eventually periodic;
    * STURM: Sturmian;
    * SUB: subcritical;
    * BCD: bounded critical discrepancy;
    * REP: strong repetitions;
    * LIN: linear complexity;
    * HARD: the rest;
    * ALL: every non-EP word.
  * the Mahler-map numbers `Phi_M(c)`;
  * `lambda 2^n` and `xi (3/2)^n` (real systems).
* **Properties:** digit avoidance, number of nonzero digits, p-adic distance to a rational, run lengths, rationality/integrality.
* **Quantifiers:** all, almost all (density one, or Haar), infinitely many.

**Evaluation.** Each statement carries the following.
* **Status:** PROVED (catalogue item or a proof here), OPEN (with FINITE-EXACT support to a stated bound where computed), or FALSE (with a counterexample).
* **Barrier profile:** the same statement for the `3n-1` sheet and for `5n+1`, where that makes sense. For Mahler, the analogue of `5n+1` is `5/2`; for Erdős, it is base 5. Plus DEFECT (EVERY/COUNT), INTEGRAL and UNIFORM.
* **Targets:** `<=>`, `=>` or `part` for T1–T4, with PC as an extra column.

A statement is **blocked** for a target when one of the following holds:
* its analogue holds on a control where the target fails (DRIFT, SHEET);
* it is a counting statement (DEFECT);
* it only sees a periodic 3-adic window (INTEGRAL);
* it is uniform over a family that contains counterexamples (UNIFORM);
* it excludes only a zero-dimensional set (DIM).

The output has 49 statements: 25 PROVED, 13 OPEN, 11 FALSE.

**Selected rows** (the full matrix is `.out` S3.2).

| id | statement | status | T1 | T2 | T3 | T4 | PC | blocked by |
|---|---|---|---|---|---|---|---|---|
| A1 | every `K >= 9`: `2^K` has a ternary 2 | OPEN; FINITE-EXACT `K <= 2*10^5` (CITED to `2*3^45`) | . | . | <=> | . | . | – |
| A2 | almost all `K` | PROVED (Narkiewicz) | . | . | . | . | . | DEFECT |
| A4 | `K not in {0,2,8}`: a ternary 2 or at least 26 ternary 1s | PROVED (Dimitrov–Howe) | . | . | part | . | . | bounded complexity |
| B1[D] | every `K`: the lowest `D` ternary digits contain a 2 | **FALSE** (`K = 2*3^(D-1)`) | . | . | . | . | . | INTEGRAL (periodic in `K`) |
| B2[D] | almost all `K`: same | **FALSE** (exceptional density `(1/2)(2/3)^(D-1)`) | . | . | . | . | . | INTEGRAL |
| R1 | the 0-run after the leading ternary digit is `O(log K)` | PROVED (Ren–Roettger) | . | . | part | . | . | does not force a 2 |
| R2 | the run above the lowest ternary digit has length exactly `v_3(K)` (`K` even) or `1 + v_3(K)` (`K` odd) | PROVED (LTE); checked `K <= 2*10^4` | . | part | . | . | . | only the balls of `+-1` |
| D1 | `v_3(2^K - 1)` bounded | **FALSE** (`K = 2*3^j` gives `j+1`) | . | . | . | . | . | INTEGRAL |
| D2 | `v_3(2^K - r) <= C_r log K` | PROVED (Yu) | . | part | . | . | . | one terminal landing only |
| D3 | `v_3(2^K - r) <= (1+eps) log_3 K + C` | OPEN (with `eps = 0`: expected FALSE, by Borel–Cantelli) | . | part | . | . | . | as D2 |
| E1 | every `s`: `2^(K0(s)) w` avoids the depth-13 balls at 1 and 1/2 | **FALSE** (Weyl; hits found) | . | . | . | . | . | equidistribution |
| E3 | every hostile chain ends in a descending class | OPEN | . | <=> | . | . | . | – |
| G[STURM] | Sturmian: `Phi` irrational | **PROVED (Thm S)** | part | . | . | . | part | DRIFT (the `5x+1` analogue is proved too) |
| G[SUB] | subcritical: `Phi` irrational | PROVED (Monks–Yazinski) | part | . | . | . | part | DRIFT |
| G[BCD] | bounded critical discrepancy | PROVED (in-house; Prop B) | part | . | . | . | part | critical slope only |
| G[HARD] | HARD class: `Phi` not in `Z^+` | OPEN | **<=>** | . | . | . | . | – |
| G[ALL].irr | PC | OPEN | => | . | . | . | <=> | – |
| G[ALL].irr (Haar) | almost all words | PROVED (trivial) | . | . | . | . | . | DEFECT |
| M1 | `Phi_M(S)` contains no positive integer | OPEN | . | . | . | <=> | . | – |
| M2 | Sturmian carry words: `Phi_M` irrational | **PROVED (Thm S)** | . | . | . | part | . | zero entropy |
| M3 | PC for the Mahler map | **FALSE** (`Phi_M(parity word of 1) = 1`) | . | . | . | . | . | – |
| H0 | every `lambda > 0`: eventually a ternary 2 | **FALSE** (Lagarias Thm 1.2) | . | . | . | . | . | UNIFORM in `lambda` |
| H2 | every `xi`: spread at least `1/3` | PROVED (FLP) | . | . | . | part | . | needs `> 1/2` |
| H3 | every `xi`: spread `> 1/2` | OPEN | . | . | . | => | . | – |
| H4 | every `xi`: some `frac(xi (3/2)^n) >= 1/2` | OPEN | . | . | . | <=> | . | – |

**Minimal unblocked statements** (`.out` S3.3).

* **T1:** G[HARD]. It is OPEN and equivalent to T1, since SUB, BCD, STURM and REP are proved. Its `3n-1` analogue is OPEN and expected true, and its `5n+1` analogue is expected FALSE (Kontorovich–Lagarias), so it sees DRIFT. PC implies it.
* **T2:** E3 (steering). D2 bounds the depth of a terminal landing by `C log K`, and R2 gives the exact depth at the balls of `+-1`. E1 shows why no "for all K" landing statement can work.
* **T3:** A1 itself.
  * The 3-adic-window strengthenings (B1, B2) are FALSE, so every proof must use the top (archimedean) digits.
  * The lambda-uniform strengthening (H0) is FALSE, so the proof must use `lambda = 1`.
  * A4 and R1 are partial.
* **T4:** H4 or M1, or the stronger H3. All three have a **FALSE** `5/2`-analogue (Tijdeman), so each must see `p < q^2`. The partial results are FLP, Mahler's countability, Dubickas–Mossinghoff, and M2.

**New exact facts used by the generator** (PROVED, elementary, and checked by the engines).
1. **Low-digit classes.** Let `N_j` be the number of classes `K mod 2*3^(j-1)` for which the lowest `j` ternary digits of `2^K` avoid 2. Then `N_j = 2^(j-1)`, checked for `j <= 16`.
   *Proof.* `2^(2*3^(j-1)) = 1 + c 3^j (mod 3^(j+1))` with `3` not dividing `c`, because 2 is a primitive root mod 9. So the three lifts of a class change digit `j` by `t c (2^k mod 3)`, and digit `j` takes each value once. ∎
2. **The kappa formula.** For a 3-adic unit `r`, write `r = (-1)^e r_1` with `r_1 = 1 mod 3`, and put `kappa(r) = log_3(r_1)/log_3(-2)`, which lies in `Z_3`. Then
   `v_3(2^K - r) = 0` if `K != e (mod 2)`, and `v_3(2^K - r) = 1 + v_3(K - kappa(r))` if `K = e (mod 2)`.
   The formula was checked for seven `r` and `K <= 3000`.
   *Proof.* `2^K = (-1)^K (-2)^K`, and the 3-adic logarithm is an isometry from `1 + 3Z_3` onto `3Z_3`, with `v_3(log_3(-2)) = 1`. ∎

   So a hostile landing at `h` is exactly `K = kappa(h/w) (mod 3^(D-1))`. The Q2 endgame asks whether the chain's total halving count `K` avoids finitely many 3-adic balls: a "3-adic digits of an integer versus a Cantor-type target" question.
3. **Loop-clock landings.** `K0(s) = floor((s+1) log_2 3)` is equidistributed mod `2*3^12` (Weyl). So each depth-13 hostile ball is hit with density `1/(2*3^12)`; the hit counts for `s <= 2*10^7` are in the `.out`, and they refute E1.

## 4. The Sturmian test (task 4)

### 4.1 Computations (FINITE-EXACT; `.out` S4)

`Phi_T(w) mod 2^N` was computed exactly from the first `N` letters, with `N = 20000`.

The words were exact mechanical words for 14 slopes: the critical `log_3 2`; sub-, near- and supercritical quadratic slopes; `log_5 2` and quadratic slopes near it; and Mahler-safe slopes. Each word is generated from a 520-bit rational approximation of its slope, with a margin check that certifies every letter.

There are three rationality tests per word.

1. **Rational reconstruction.** Wang's algorithm at `N` bits returns the unique candidate of height `<= 2^9999`. Its parity vector is then followed exactly for `8N` letters and compared with `w`.
2. **Repetition certificates.** These come from Theorem R, with the options A and B of the proof of Theorem S. Each certificate is cross-checked through the isometry: `v_2(Phi(w) - P/Q) = lambda` exactly.
3. **Bit length.** The bit length of `Phi mod 2^N` is at least `N - 8 = 19992` in every row, so `Phi(w)` is not a positive integer below `2^19991`.

| map | slopes (regime) | intercepts | words | rational certificate | best gains (bits) |
|---|---|---|---|---|---|
| `3x+1` | `log_3 2` (crit), 4 subcritical, 4 supercritical incl. `0.6309297540` just above critical | `0`, `0` upper, `1/2`, random `k/2^60` | 36 | no rational of height `<= 2^9999` for every word | `2*10^4` to `1.2*10^5` |
| `3x-1` | `log_3 2`, `sqrt2/2` | `0` | 2 | same (`Phi_(3,-1) = -Phi`) | as `3x+1` |
| `5x+1` | `log_5 2` (crit), `[0;2,3,(9)]` and `[0;2,(3)]` (super), `sqrt2-1` (sub), `1/phi`, `sqrt2/2` | `0`, `0` upper, `1/2` | 18 | same | `3.5*10^4` to `9.8*10^4` |
| Mahler | `[0;3,(3)] = 0.303`, `sqrt2-1` | `0` | 2 | same | `3.5*10^4`, `6.1*10^4` |
| deep row | `3x+1`, `log_3 2`, `N = 10^5` | `0` | 1 | no rational of height `<= 2^49999` | `1.26*10^5` at `p = 50508` |

**Findings.**

* **Gains track the proof bound.** At each convergent denominator `q_n` the best gain is above, or within `O(log q_(n+1))` of, `q_(n+1)/mu - (mu-1) q_n`, as proved. For example, at the critical slope, `q = 1054` gives a gain of `25771` bits (predicted `24727`). The partial quotient 23 of `log_3 2` makes that repetition 24 periods long. The last listed period of a row can show a small or negative gain; that is a window artifact, because its long repetition runs past the `8N` letters examined.
* **A small-height candidate that agrees for a long stretch is not evidence of rationality.** A rational candidate can agree with `w` far beyond `N`. For `log_3 2` the candidate has height `2^1063` and agrees for `25780` letters. It is exactly the approximant `Phi(v^inf)` of a long repetition.

  So a finite rationality test must follow the candidate until it leaves `w`. An earlier version that checked only `N + 256` bits returned false "RATIONAL" verdicts for many words at `N = 4000`, including every critical-slope row; section 7 keeps this as a warning.
* **Finite windows cannot separate nearby slopes.** The supercritical slope `[0;1,1,1,2,2,3,1,5,2,23,2,(3)] = 0.6309297540` differs from `log_3 2` by `4*10^-10`. Its mechanical words agree with the critical ones on the whole window, and its rows coincide with the critical rows. Only Theorem S covers both.
* **Controls behave as predicted.**
  * EP words reconstruct exactly: `223/45`, and `27`.
  * Random words have best gains around `8` to `14` bits (logarithmic).
  * The random bounded-critical-discrepancy word has gain `13.6`: Theorem R does not reach it, and the in-house theorem does.
  * The in-house critical halving words gain `563` (they are Sturmian).
  * Thue–Morse gains grow like `p` along `p = 2^k`, but Thue–Morse is subcritical anyway.

### 4.2 Literature on Sturmian parity vectors

This is summarized in section 1, items 4–6, and rests on lit_A's searches (zbMATH Open, arXiv API, web).

* **Known before this note:**
  * the subcritical case, via Monks–Yazinski 2.7(b);
  * the periodic balanced case for positive integers (Knight 2026);
  * the critical-slope case for positive integers (in-house).
* **Open before this note:** rationals with Sturmian parity vectors of slope `>= log_3 2`. López–Stoll 2009 state the general aperiodic question as open; López–Stoll 2021 is unrefereed and has an apparent gap.
* **Not found:** any statement of Theorem S. The proof below is elementary, experts may know it, and no priority is claimed.

### 4.3 The theorems (PROVED here)

**Setting.** An *affine 2-adic shift map* is `T(x) = (m_e x + r_e)/2` for `x = e (mod 2)`, with `m_0, m_1` odd and `r_e = e (mod 2)`. Examples:
* `3x+r`: `m_0 = 1, r_0 = 0, m_1 = 3, r_1 = r` odd;
* `5x+1`;
* Mahler's `a -> ceil(3a/2)`: `m_0 = m_1 = 3, r_0 = 0, r_1 = 1` (THM-2228).

For a finite word `z`, `T^|z|(x) = (M_z x + R_z)/2^|z|` on the cylinder of `z`. Here `M_z = prod m_(z_i)`, `R_() = 0`, and `R_(ze) = m_e R_z + r_e 2^|z|`.

*Sturmian words* are the aperiodic balanced binary words, equivalently the lower and upper mechanical words of irrational slope (Morse–Hedlund 1940; Lothaire, *Algebraic Combinatorics on Words*, Thm 2.1.13; standard, [R], not re-read).

The cylinder of `z` is one class mod `2^|z|`, so the parity-vector map is a bijective isometry `Z_2 -> {0,1}^N` (sibling ladder Lemma 1; THM-2228 section 1). Write `Phi_T` for its inverse. Then `|Phi_T(w) - Phi_T(w')|_2 = 2^(-lambda)`, where `lambda` is the length of the common prefix of `w` and `w'`. `T` maps `Q cap Z_2` into itself.

**Periodic approximants.** For finite `u` and nonempty `v`, `Phi_T(u v^inf) = P/Q`, where:
* `Q = M_u D_v`;
* `P = 2^|u| R_v - R_u D_v`;
* `D_v = 2^|v| - M_v`, which is odd and nonzero.

Indeed, `R_v/D_v` is the fixed point of the `v`-branch, and `M_u x' + R_u = 2^|u| R_v/D_v`.

**Lemma G (2-adic gap).** If `a/b != P/Q` with `b` and `Q` odd, then `|a/b - P/Q|_2 >= 1/(|a||Q| + |P| b)`.
*Proof.* `aQ - Pb` is a nonzero integer, and an integer `N != 0` has `|N|_2 >= 1/|N|`. ∎

**Theorem R (repetition criterion).** Let `w` be a word that is not eventually periodic, and suppose `Phi_T(w) = a/b` with `b` odd. Then for every pair `(u, v)`:

`max(|a|, b) >= 2^g`, where `g(u,v) = lambda(w, u v^inf) - log_2(|P| + |Q|)`.

So if `g` is unbounded, `Phi_T(w)` is irrational.

*Proof.* `u v^inf != w`, so the two Bernstein numbers differ (`Phi` is injective). The isometry and Lemma G give `2^(-lambda) >= 1/(max(|a|, b)(|P| + |Q|))`. ∎

FINITE-EXACT control: on 17,909 pairs (genuine rational `x = a/b`, approximant `P/Q` from a prefix of its own parity word) for `3x+1`, `5x+1` and Mahler's map, `2^lambda <= |a||Q| + |P|b` and `v_2(x - P/Q) = lambda` held in every case.

Theorem R is Calegari's p-adic irrationality criterion (2005, Lemma 2.2) with Padé approximants replaced by **periodic-extension approximants**. That is why it reaches exactly the words with long, early repetitions.

**Lemma L** (Monks–Yazinski 2004, Thm 2.7(b); re-proved). For `T_(q,r)`, a nonzero rational with a non-EP parity vector has `liminf a_s/s >= log_q 2`.

*Proof.* Let the denominator be `b`. The orbit lies in `(1/b)Z`, and its values are distinct and nonzero. An odd step is `y -> (q/2) y (1 + r/(qy))`, and at most `2Mb + 1` orbit values satisfy `|y| <= M`. Hence

`1/b <= |T^s(x)| <= |x| C_M q^(a_s) (1 + |r|/(qM))^(a_s) 2^(-s)`.

Let `s -> infinity`, then `M -> infinity`. ∎

**Lemma J (periodic stretches of Sturmian words).** Let `s` be the lower (floor) or upper (ceiling) mechanical word of irrational slope `alpha` and intercept `rho`. Let `p_n/q_n` be a convergent of `alpha` and `delta = q_n alpha - p_n`. Let `J` be the set of `j` for which the floor (or ceiling) of `j alpha + rho + delta` differs from that of `j alpha + rho`.
1. If neither `i` nor `i+1` is in `J`, then `s(i + q_n) = s(i)`.
2. Distinct elements of `J` differ by at least `q_(n+1)`.
3. `min J < q_n + q_(n+1)`, and `J` is infinite.

*Proof.*
1. `s(i + q_n) - s(i) = chi(i+1) - chi(i)`, where `chi` is the difference in the definition of `J`.
2. `j` is in `J` iff `{j alpha + rho}` lies in a fixed half-open arc of length `||q_n alpha||`. Two such `j` would give `||(j - j') alpha|| < ||q_n alpha||`, which is impossible for `0 < |j - j'| < q_(n+1)` by best approximation.
3. The first `q_n + q_(n+1)` points of the orbit leave gaps of at most `||q_n alpha||` (three-distance theorem). `J` is infinite because the rotation is minimal. ∎

FINITE-EXACT check: 580 triples (slope, intercept, `n`) over 14 slopes and 5 intercepts, on words of length 60,000. Defect clusters have size at most 2, start before `q_n + q_(n+1)`, and are at least `q_(n+1) - 1` apart. There were no failures.

**Lemma H (heights of balanced approximants).** A factor `z` of a Sturmian word of slope `alpha` has `|a(z) - |z| alpha| < 1`, and its `l`-th 1 sits at a position `d_l <= l/alpha`. Hence, for `T_(q,r)`, with `mu = max(1, alpha log_2 q)`:
* `M_z <= q 2^(mu |z|)`;
* `|R_z| <= |r| 2^(1/alpha) a(z) 2^(mu |z|)`;
* `log_2(|P| + |Q|) <= mu (|u| + |v|) + log_2(|u| + |v| + 2) + C(alpha, q, r)`.

For Mahler's map, `M_z = 3^|z|` and `0 <= R_z < 3^|z|`, so `log_2(|P| + |Q|) <= log_2 3 (|u| + |v|) + 2`.

*Proof.* `R_z/r = sum_(l <= a) q^(a-l) 2^(d_l) <= q^a sum_(l <= a) theta^l`, with `theta = 2^(1/alpha)/q`. That sum is at most `a max(theta, theta^a)`. The other bounds follow from `|P| <= 2^|u| |R_v| + |R_u||D_v|` and `|D_v| <= max(2^|v|, M_v)`. ∎

**Theorem S.** Let `T` be an affine 2-adic shift map, and `s` a Sturmian word (any irrational slope, any intercept, lower or upper) for which Lemma H holds with `mu (mu - 1) < phi`, i.e. `mu < 1.8668`. Then `Phi_T(u s)` is irrational for every finite `u`. In particular:
1. **`3x+r`, every odd `r`, every slope.** Here `mu <= log_2 3 = 1.585`.
   * No rational with odd denominator has an eventually Sturmian parity vector, whether it is positive or negative, integer or not.
   * The Periodicity Conjecture holds on the class of eventually Sturmian words.
   * No divergent orbit of a rational has an eventually Sturmian parity vector.
2. **`5x+1`: slopes `alpha < 0.8040`.** This range contains the drift threshold `log_5 2 = 0.4307` and the typical density `1/2`. For `7x+1` the range is `alpha < 0.6650`.
3. **Mahler's map** (`mu = log_2 3`). The carry word of a Z-number `xi` is the parity word of `floor(xi (3/2)^n)`, and THM-2228 gives `Phi_M(c) = A`, a positive integer. So for a Z-number these parities are never eventually Sturmian.

*Proof.* Fix `n`, put `q = q_n`, and let `j_0 < j_1` be the two least elements of `J`. There are two options.
* **Option A:** `u = ()` and `v = s[0, q)`. Then `lambda_A >= j_0 - 1 + q`.
* **Option B:** `u = s[0, j_0]` and `v = s[j_0 + 1, j_0 + q]`. Then `lambda_B >= j_1 - 1 + q >= j_0 + q_(n+1) - 1 + q`.

By Lemma H:
* `g_A >= j_0 - (mu - 1) q - O(log q_(n+1))`;
* `g_B >= q_(n+1) - (mu - 1)(j_0 + q) - O(log q_(n+1))`.

Split at `j_0 = q_(n+1)/mu`. Either way,

`G_n := max(g_A, g_B) >= (q_n/mu)(q_(n+1)/q_n - mu(mu - 1)) - O(log q_(n+1))`.

For every irrational `alpha`, `limsup q_(n+1)/q_n >= phi`: the ratio is at least 2 whenever `a_(n+1) >= 2`, and otherwise it tends to `phi`. Since `mu(mu - 1) < phi`, `G_n >= c q_n` along infinitely many `n`, and Theorem R applies. For `u s`: `T^|u|` maps a rational `Phi_T(u s)` to the rational `Phi_T(s)`. ∎

**A second route (CITED input).** ADQZ 2001 [R via BHZ 2006, P] show that every Sturmian word begins with infinitely many squares `v^2`. With `u = ()` this gives `g >= (2 - mu)|v| - O(log |v|)`, so Theorem S holds for all `mu < 2`. For `5x+1` that means `alpha < 0.8614`.

### 4.4 What this covers, and the in-house discrepancy theorem (coordinator item 1)

**What the in-house theorem covers.** Its statement is in [collatz_guards_20260921_discrepancy](collatz_guards_20260921_discrepancy.md), sections 2a, 3 and 4.

* **Covered:**
  * positive-integer `3n+1` orbits whose halving word has bounded discrepancy `K_j - j log_2 3` (any width), also after dropping a finite prefix;
  * `3n+b` for every positive odd `b` (section 4, in prose);
  * as a corollary, critical-slope Sturmian halving words. In `T`-coding these are Sturmian parity vectors of slope `log_3 2`;
  * for `b < 0`, the result (D8) that every positive orbit has discrepancy tending to `-infinity`.
* **Not covered:**
  1. **Supercritical slopes.** The orbit grows exponentially, so the capacity count (distinct images below a linear bound) and the density-zero lemma both give nothing. Covered here for Sturmian words by Theorem S, all slopes; OPEN in general (candidate C2).
  2. **Critical slope with unbounded discrepancy**, for example excursions of size `sqrt j` or `log j`. This is OPEN. Theorem S covers only its zero-entropy Sturmian part.
  3. **Rationals.** These are covered by **Proposition B** below.
  4. **`5n+1`.** The density-zero lemma needs mean halving `2 > log_2 q`, which fails for `q = 5`. Only the width bound (D7), `B/A >= 12.5`, survives. Theorem S covers the Sturmian `5n+1` words.
  5. **Mahler's map.** Both branches expand, so there is no critical slope. Theorem S covers the Sturmian carry words.

**Proposition B** (PC on bounded critical discrepancy; PROVED, assembled). No rational with odd denominator has a non-EP parity vector with `sup_s |a_s - s log_3 2| < infinity`.

*Proof.* The orbit takes infinitely many distinct values in `(1/b)Z`, so it is unbounded (Monks–Yazinski Cor. 3.3). Some tail of it has constant sign:
* positive values stay positive;
* a value in `[-1/3, 0)` either stays in `[-1/3, 0)` forever, which would make the orbit bounded, or jumps into `[0, 1/2)`;
* so some tail lies entirely in `(0, infinity)` or entirely in `(-infinity, -1/3)`.

Write `x = c/b`.
* **Positive tail.** The numerators form a positive `3n+b` orbit with the same halving word. This contradicts the in-house section 2a as extended in section 4. That extension is stated there in prose, with the changed constants, and was not re-derived line by line here.
* **Negative tail.** The negated numerators form a positive `3n-b` orbit, and `3n - b > 0` holds along it because the values lie below `-1/3`. By (D8), `q_j -> 0`, so the discrepancy tends to `-infinity` and is unbounded. ∎

**Typing of the in-house theorem.**
* **SHEET: B.** On `3n-1` it holds in the stronger form (D8). The `3n-1` cycles have discrepancy that drifts linearly, because `2^K != 3^L`.
* **DRIFT: O for the proof, B/OPEN for the statement.** The density-zero stopping-time lemma uses `2 > log_2 3`. The `5n+1` statement is OPEN; heuristically it is also true, since bounded-discrepancy orbits are non-generic. Only (D7) transfers.
* **DEFECT: O.**
* **INTEGRAL: O.** It uses the discreteness of the integers.
* **UNIFORM: B.** It holds for all `3n+b`, `b > 0`.
* **Structural:** it is the **only** every-orbit result that excludes a positive-entropy class of words. Strips wider than 1 have positive entropy; strips narrower than 1 contain only Sturmian words.

So the two results are complementary:
* the in-house theorem: positive entropy, critical slope only, integers, by capacity;
* Theorem S: zero entropy, all slopes, all rationals, by transversality.

**The remaining class.** Supercritical, positive-entropy words without strong repetitions (HARD) are untouched by both. This is the precise missing piece of T1 and PC.

### 4.5 Padé and holonomy p-adic methods (coordinator item 2)

**Judgment: they cannot reach lacunary or Bernstein-type 2-adic series of general words.** They need the target to be a special value of one holonomic or modular function, which supplies approximants with controlled denominators. But a 0/1 word's series `sum w_n z^n` is D-finite only when `w` is eventually periodic (Pólya–Carlson; [R] via Bell–Miles–Ward 2014 and Bell–Chen–Nguyen–Zannier), which is the INTEGRAL case where `Phi(w)` is rational anyway, and PC quantifies over uncountably many words.

The refereed anchor is Calegari 2005 (IMRN 2005, no. 20, 1235–1249; arXiv text read, publication checked via zbMATH and Crossref), which proves `zeta_2(3)`, `zeta_3(3)` and `L_2(2, chi_4)` irrational. Its criterion (Lemma 2.2) is exactly the one Theorem R uses, with Padé approximants replaced by periodic-extension approximants. The unrefereed 22-value draft ([audit](../reference/p-adic-zeta-irrationality-source-audit-20260825.md)) is not used anywhere in this note.

### 4.6 Transcendence (hypothesis candidate C1, with a sketch)

The claim is that `Phi(s)` is transcendental for every Sturmian `s` under `3x+1`.

*Sketch.* Take the approximants of the proof above, with gain `>= eps (|u| + |v|)`. Consider the integer points `X = (3^(a(u)) 2^|v|, 3^(a(u)+a(v)), P)` and the places `{inf, 2, 3}`.
* At the 2-adic place, use the forms `X_1`, `X_2` and `x(X_1 - X_2) - X_3`. Elsewhere use the coordinates.
* The product over places is at most `2^(-g + O(1))`, which is at most `H(X)^(-eps')`.
* By Schlickewei's p-adic subspace theorem, infinitely many of the points lie on one rational plane `c_1 X_1 + c_2 X_2 + c_3 X_3 = 0`.
* If `c_3 != 0`, letting `|v| -> infinity` 2-adically forces `x = c_2/c_3`, which contradicts Theorem S.
* If `c_3 = 0`, then `2^|v|/3^(a(v))` is constant, which is impossible.

This mirrors Adamczewski–Bugeaud 2007 section 6 [P] and Ferenczi–Mauduit 1997 [R]. It is not audited, so it stays a candidate.

## 5. Cross-problem matrix (task 5)

| problem | status | single missing ingredient | proved partial ingredients | shared statement |
|---|---|---|---|---|
| Collatz divergence (T1) | OPEN | every-orbit 2-adic non-integrality of `Phi` on HARD (supercritical, positive entropy, no strong repetitions) | SUB (Monks–Yazinski), BCD (in-house), STURM and REP (this note) | Bernstein family |
| Collatz cycles | OPEN | integrality exclusion of cycle points `R_v/(2^len(v) - 3^(a(v)))` at the convergent clocks (carry anti-concentration plus Baker) | `m <= 91` (Hercher), circuits (Steiner), Christoffel cycles (Knight), Ellison, Rhin | outside the grammar (Theorem R is vacuous on EP words) |
| E-SCC Q1 | OPEN | a sheet-sensitive escape family for the hostile points `-1 - 2^i c/3^j` | dimension evidence favours 0; the A–L escape is sheet-blind | mirror of the D/F families |
| E-SCC Q2 (T2) | OPEN | chain control plus terminal steering E3 (3-adic balls at `kappa(h/w)`) | terminal depth `<= C log K` (D2), exact LTE depth at `+-1` (R2) | `2^K`-in-`Z_3` family (low window) |
| Erdős ternary (T3) | OPEN | an every-`n` statement about the **top** ternary digits of `2^n` (`1 not in E(Z_3)`) | Narkiewicz/HKR counts, `dim E(Z_3) <= log_3 phi`, Dimitrov–Howe, Saye, Ren–Roettger | `2^K`-in-`Z_3` family (full expansion) |
| Mahler `3/2` (T4) | OPEN | 2-adic non-integrality of `Phi_M` on the `beta`-shift `S` (entropy `log(3/2)`) | Mahler countability, FLP, Dubickas–Mossinghoff, M2 (this note) | Bernstein family (map `ceil(3a/2)`) |
| a divergent `5x+1` orbit exists | OPEN (expected) | the **reverse** transfer: one rational (e.g. 7) with a non-EP parity vector | none. Theorem S only says such an orbit is not Sturmian for `alpha < 0.804` | opposite direction |
| Periodicity Conjecture | OPEN | `Phi` irrational on HARD, for all `3x+k` | SUB, BCD (Prop B), STURM, REP | Bernstein family |
| Furstenberg `x2x3` measure | OPEN | measure rigidity at zero entropy | Rudolph–Johnson, Host, Hochman–Shmerkin, Shmerkin–Wu | none: even a proof is COUNT-type, DEFECT-blind for T1–T4 |

**Does one statement serve several problems?**

* **Yes, in form: the Bernstein family.** "`Phi_T(w)` is not a positive integer (or not rational) for `w` in a class" is T1 and PC for `3x+1`, and T4 for `ceil(3a/2)` (THM-2228). Theorem S proves all three on the Sturmian class with one argument.
* **The drift decides the direction.**
  * For `3x+1` the statement should hold (PC).
  * For `5x+1` its failure is exactly the existence of a divergent orbit.
  * For `5/2` its Mahler analogue is FALSE (Tijdeman).

  So a single drift-blind proof cannot serve T1 and T4 while their controls fail. Akin's heuristic draws the same line at `|a| = 3` versus `|a| >= 5`.
* **Only partly for T2 and T3.** Both concern `2^K` in `Z_3`, but T3 needs the top digits (every 3-adic-window version is FALSE), while T2 needs the low digits along dynamically chosen `K` (the kappa formula). The shared ingredient would be a 2-versus-3 statement about `K` versus a Cantor-type set in `Z_3`, which is not in the catalogue.
* **Furstenberg's conjecture serves none of them.** It is a measure statement.

## 6. Hypothesis candidates (for `INDEX.md` numbering by the integrator; no HYP files created)

* **C1 (transcendence).** `Phi_T(s)` is transcendental for every Sturmian `s`, under `3x+r`, `5x+1` (`alpha < 0.804`) and Mahler's map. Sketch in section 4.6.
* **C2 (supercritical strips).** No rational has a non-EP parity vector with `sup_s |a_s - alpha s| < infinity` for some `alpha in (log_3 2, 1)`.
  * This is the positive-entropy extension of Theorem S (strips of width at least 1) and the supercritical extension of Proposition B.
  * Neither capacity (the orbit grows exponentially) nor repetitions (positive entropy) reaches it.
  * **Recommended next target.** It is the smallest statement that needs a new mechanism.
* **C3 (quasi-Sturmian words and rotation codings).** PC for words of complexity `n + C`, and for codings of irrational rotations by finitely many intervals. The Lemma J argument should generalize. Probably PROVABLE.
* **C4 (3-adic Lang–Waldschmidt).** `v_3(2^K - r) <= (1 + eps) log_3 K + C(r, eps)` (D3). With `eps = 0` it is expected FALSE for generic `r`. The FINITE-EXACT excess up to `3*10^5` is in the `.out`.
* **C5 (Q2 steering, E3).** Every hostile see-saw chain ends in a descending class. Via the kappa formula this is a statement about the 3-adic digits of the total halving count.
* **C6 (the zero-entropy class is settled).** Together, Knight 2026 (cycles, positive integers) and Theorem S (divergence, all rationals) settle both halves of Collatz on balanced words. Check whether Knight's difference argument extends to all rational Christoffel cycles of `3x+k`. This is sheet-sensitive through `|2^k - 3^x| = 1`.

## 7. Reproduction

```bash
cd <worktree>
python3 04-computation/experiments/collatz_procgen_20260922_transversality.py \
    > 05-knowledge/results/collatz_procgen_20260922_transversality.out
# about 10 minutes, one process, peak memory well under 900 MB; --quick for a 1-minute smoke run.
# Timing goes to stderr; the output is deterministic (fixed seeds).
```

* **Sections of the `.out`:**
  * S1: targets and grammar.
  * S2: catalogue.
  * S3.0: engines.
  * S3.1–S3.3: statements, matrix, and minimal statements.
  * S4.1–S4.3: Sturmian rows, controls, and the deep row.
  * S5: cross-problem matrix.
* **Every certificate is exact.**
  * The mechanical words carry a 520-bit margin check.
  * Rational candidates are followed letter by letter until they leave `w`.
  * Each repetition certificate is cross-checked through the isometry, `v_2(Phi(w) - P/Q) = lambda`.
  * The low-digit counts are checked by brute force for `j <= 9`.
  * The kappa formula is checked against direct valuations.
* **Helper scripts:** `scratch/procgen_transversality/check_lemmaJ.py` (the Lemma J census; 580 checks, about 5 s) and `check_theoremR.py` (17,909 inequality and isometry checks on genuine rationals, under 1 s).
* **Warning.** Rational reconstruction checked only `N + 256` bits past the window gives false "RATIONAL" verdicts for Sturmian words with large partial quotients, because the approximant of a long repetition agrees far beyond `N`. The generator therefore follows the candidate until it leaves `w`.

## Sources

`[P]` means a primary text was read by the lane's sub-lanes; `[R]` means a named secondary source was used. Blocked sources were not bypassed.

* **Collatz and the Periodicity Conjecture.**
  * J. C. Lagarias, Amer. Math. Monthly 92 (1985) 3–23 [P, CECM reprint].
  * D. J. Bernstein, Proc. AMS 121 (1994) 405–408 [P].
  * Bernstein–Lagarias, Canad. J. Math. 48 (1996) 1154–1169 [P].
  * Monks–Yazinski, Discrete Math. 275 (2004) 219–236 [P, preprint].
  * E. Akin, Contemp. Math. 356 (2004) [R].
  * López–Stoll, Integers 9 (2009) A13 [P]; arXiv 2101.12747 [P, unrefereed].
  * O. Rozier, Integers 19 (2019) A8 [P].
  * B. Knight, Discrete Math. 349 (2026) 114812 [P, HAL preprint].
  * Lagarias 1990, Acta Arith. 56 [R].
  * Halbeisen–Hungerbühler, Acta Arith. 78 (1997) [P].
  * Tao, Forum Math. Pi 10 (2022) e12 [P].
* **Digits and Diophantine results.**
  * Senge–Straus 1973 [R].
  * Stewart 1980 [P].
  * Mahler 1957 [R]; Mahler 1968 [P].
  * Corvaja–Zannier 2004 [P]; Zudilin 2007 [P].
  * Flatto–Lagarias–Pollington 1995 [R]; Lu–Zheng arXiv 2603.16794 [P].
  * Dubickas–Mossinghoff 2009 [R]; Andrieu–Eliahou–Vivion arXiv 2510.11723 [P].
  * Narkiewicz 1980 [R]; Holdum–Klausen–Rasmussen, INTEGERS 15 (2015) A43 [P].
  * Lagarias, J. London Math. Soc. 79 (2009) [P].
  * Abram–Lagarias 2014 [P]; Abram–Bolshakov–Lagarias 2017 [P].
  * Dupuy–Weirich 2016 [R].
  * Saye 2022 [P].
  * Dimitrov–Howe 2025 [P].
  * Yu 2007 [R]; Bugeaud–Laurent 1996 [R].
  * Ellison 1971 [P].
  * Ren–Roettger arXiv 2511.03861 [P].
* **Dynamics and measure.**
  * Furstenberg 1967 [R].
  * Rudolph 1990 [P].
  * Host 1995 [R]; Hochman–Shmerkin 2012 [P] and 2015 [P].
  * Cassels 1959 [R]; Schmidt 1960 [P].
  * Shmerkin 2019 [P]; Wu 2019 [P].
  * BLMV 2009 [R/P].
* **Transcendence and irrationality.**
  * Calegari 2005 [P].
  * Pólya–Carlson [R].
  * Mahler 1929 [R]; Loxton–van der Poorten 1977 [P].
  * Ferenczi–Mauduit 1997 [R]; Adamczewski–Bugeaud 2007 [P].
  * Allouche–Davison–Queffélec–Zamboni 2001 [R]; Berthé–Holton–Zamboni 2006 [P].
* **Repository.**
  * [synthesis](collatz_procgen_20260922_synthesis.md), [barrier atlas](collatz_procgen_20260922_barrier_atlas.md), [foundry](collatz_procgen_20260922_foundry.md).
  * [sibling ladder](collatz_procgen_20260922_sibling_dimension_ladder.md), [order laws](collatz_procgen_20260922_order_laws.md), [loops and escapes](collatz_procgen_20260922_loops_and_escapes.md).
  * [in-house discrepancy](collatz_guards_20260921_discrepancy.md), THM-2228, THM-3848.
