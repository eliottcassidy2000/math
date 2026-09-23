# Nine outside items against the Periodicity-Conjecture frontier: an interplay map, the moonshine chain from 6 = 5+1 to 196884 = 196883+1, and the seven-hexagon packing

**Status.**
* **AUDIT (orchestrator, 2026-09-23).** Re-derived by hand, with no gap found:
  * **Corollary M.** `j(rho^2) = 256(1 - l + l^2)^3/(l^2(1-l)^2)` with `l = 16 rho psi(rho^2)^4/theta_3(rho)^4`. This is a formal q-series identity with integer coefficients, valid 2-adically. BDGP (the p-adic Mahler–Manin theorem) applies at the algebraic `rho^2` with `0 < |rho^2|_2 < 1`. Both swap numbers are affine in `theta_3(rho)` and `psi(rho^2)`.
  * **Proposition T′.** The coordinates `X1 = 3^(a_U) 2^p` and `X2 = 3^(a_U + a_V)` are S-units. The forms are independent at each place (`b, d != 0`). The product is `<= C p |Q| 2^(-eta lambda) = H^(1 - rho' + o(1))`. S-unit gcds are harmless, and all three exceptional-plane cases close, using `Phi_2(w) != Phi_R(w)` for non-EP words. It depends on the CITED Schlickewei p-adic Subspace Theorem.
* **PROVED** (hand proofs in this note; audited as above):
  * **Proposition T′** (§2.1). Assume the CITED rational p-adic Subspace Theorem. Let `w` be strictly supercritical, of bounded discrepancy, not eventually periodic (EP), with `Dio(w) > 1`. Then `Phi_2(w)` and `Phi_R(w)` are not both rational. The hard note's Proposition T needed `Dio > 2`.
  * **Corollary M** (§2.3). Assume the CITED p-adic Mahler–Manin theorem (Barré-Sirieix–Diaz–Gramain–Philibert 1996). Then the 2-adic Bernstein numbers of the square-swap word `Y` and of the pronic-swap word `Y_pr` are not both algebraic.
  * **No-go for HC1** (§2.1, elementary). Two approximant families have natural heights: the periodic approximants and the partial sums of `sum rho^(k^3)`. For both, the Subspace product over `S = {inf, 2, 3}` is `H^(1 - 1/eta + o(1))`, where `1 - 1/eta = 0.299 > 0`. So a three-place Subspace cover, such as the Lean one in Erdős 1062(ii), cannot settle the cube-swap word `Y3` (HYP-9127) through them.
* **FINITE-EXACT** (`.out`; independent code paths where stated):
  * The whole moonshine chain (§3): the outer automorphism of `S_6`; the hexacode, the ternary Golay code and the binary Golay code as extended quadratic-residue (QR) codes; `S(5,6,12)` and `S(5,8,24)`; `|M_12| = 95040`; `|M_24| = 244823040`; the hexad stabilizer `S_6` and the dodecad stabilizer `M_12`, each acting on its set and on the complement through an outer twist; the 196560 Leech minimal vectors and their inner-product distribution; `Theta_Leech/eta^24 - 24 = j - 744`; and every decomposition of 196884 and its successors.
  * The seven-hexagon packing at side `5/sqrt3`, in exact `Q(sqrt3)` arithmetic (§4).
  * Kruer–Kohlmeyer's finite formula for Erdős 1062(ii) against an exact ILP brute force for `n <= 70`, and their density `L = 0.672965651199489347471...` (§2.1).
  * The 2-adic identities for `Y`: the Jacobi triple product and the swap identities mod `2^3000`–`2^4000`, and `j(rho^2)` computed two ways mod `2^3940`.
  * The `N = 5` finite quantum dilogarithm, numerically to `1e-40`.
  * The `Out(S_6)` twist on all `2^15` six-vertex Rédei–Berge functions.
  * Statistics of the `{log n}` construction.
* **CITED**:
  * Radchenko–Wheeler, Karingula–Lovett and Ezeunala–Jiang [P, PDFs].
  * The accepted Lean results for Erdős 859 and 1062(ii) [P, expositions only; Lean sources not read].
  * BDGP 1996 [R, via Waldschmidt's Bourbaki exposé, read].
  * Mason 1985 [R]; McKay/FLM/Borcherds [R].
  * Friedman's Packing Center and Berthold et al. 2026 [P].
* **OPEN**:
  * PC, HC1 (`Y3`) and C2.
  * Whether `theta_3(rho)` is transcendental (this would need a p-adic Nesterenko theorem).
  * Whether `5/sqrt3` is optimal for 6 or 7 hexagons; in every source found it is a record, not a proven optimum.
  * Whether the log-table number is normal.
  * The source of the `{log n}` construction, which was **not located**.
* **REFUTED**:
  * "The `Out(S_6)` twist is a symmetry of the six-vertex Rédei–Berge family." It maps `U_T` into the family only when it fixes `U_T` (3 of 29 cycle-count triples).
  * "The natural-log reading of the `{log n}` construction has a column-1 limit law." Its frequencies oscillate from round to round.

Session `collatz-procgen-20260922` (wave 5), sources lane, 2026-09-23. It builds on the [synthesis](collatz_procgen_20260922_synthesis.md) (§2c, §4) and the [hard-class note](collatz_procgen_20260922_hard_class.md) (§1.5, §2.4).
* **Scripts:** [`procgen_sources_20260923_moonshine.py`](../../04-computation/experiments/procgen_sources_20260923_moonshine.py), [`_hexagons.py`](../../04-computation/experiments/procgen_sources_20260923_hexagons.py), [`_logdigits.py`](../../04-computation/experiments/procgen_sources_20260923_logdigits.py) and [`_interplay.py`](../../04-computation/experiments/procgen_sources_20260923_interplay.py). The runner is [`_run.py`](../../04-computation/experiments/procgen_sources_20260923_run.py).
* **Output:** [`procgen_sources_20260923.out`](procgen_sources_20260923.out). It is deterministic and identical under a different `PYTHONHASHSEED`.

## 0. Answers in brief

1. **Erdős 1062(ii) and the missing mechanism.** Their proof and Theorems S/D run on the same engine:
   * rationality forces infinitely many small, nonzero S-unit linear forms, each with a single non-S-unit integer coordinate `P`;
   * a Diophantine rigidity statement then forbids them.

   The two settings differ in one number, the **exchange rate** `eta`: the bits of height paid per bit of precision gained.
   * **Their series.** `L = (32 + sum a_k 4^-k + sum b_k 3^-k)/180`, with bounded coefficients, has `eta = 1`. So any repetition exponent above 1 suffices, and the hard part is nonvanishing.
   * **Supercritical parity words.** Here `eta = beta log2 3 > 1`, so the threshold is `Dio > eta` (Theorem D, and HARD §1.5).
   * **Rotation.** Their coefficient words `a_k`, `b_k` are codings of the rotations by `log_3 4 = 2 log_3 2` (twice the critical `3x+1` slope) and by its reciprocal `log_4 3` (FINITE-EXACT, `k < 3000`).
   * **Formalization.** Theorems S and D do not need the cover: they are elementary 2-adic Liouville arguments.
   * **What the cover would give.** A rational three-place cover would formalize the new Proposition T′.
   * **Y3.** It cannot reach `Y3`: the no-go above.
2. **Theorem Y against the Radchenko–Wheeler finite quantum dilogarithms.** The shared predicate is a *first-order q-difference equation*. The Padé polynomials `(qT;q)_n` of Theorem Y and Faddeev's `Phi_b` both come from one. Radchenko–Wheeler prove *algebraicity* at `|q| = 1`; Theorem Y proves *irrationality* at `|q|_2 < 1`. Cubes have no such equation. Label: ANALOGY. There is no transfer to HC1.
3. **Tate curve.**
   * `theta_3(rho)` is a theta-null of the Tate curve `E_(rho^2)` over `Q_2`.
   * `j(rho^2) = rho^-2 + 744 + 196884 rho^2 + ...` is a rational function of `theta_3(rho)` and `psi(rho^2)`, checked two ways mod `2^3940`.
   * The citation is verified: BDGP, Invent. Math. 124 (1996) 1–9, solved Manin's p-adic case together with Mahler's complex case. So `j(rho^2)` is transcendental, which gives Corollary M.
   * Transcendence of `theta_3(rho)` alone would need a p-adic Nesterenko theorem, not available in the sources read.
   * `sum rho^(k^3)` is not modular, so it has no Tate-curve handle.
4. **Komlós and Kadison–Singer.** They share the predicate "bounded prefix discrepancy" with the C2 strips (HYP-9123), with the opposite quantifier. Discrepancy theory proves that balanced signings *exist*. C2 asks that no rational *realize* a balanced non-EP word. By Terras, every finite strip word is realized by exactly one class mod `2^s`, so C2 is a pure 2-adic limit statement. Label: ANALOGY. There is no transfer.
5. **Erdős 859.** Ford's `delta = 1 - (1 + ln ln 2)/ln 2 = 0.0860713` equals `Q(1/ln 2)`, a Poisson relative entropy. The session's critical rate `D(log_3 2 || 1/2) = 0.0346882` nats is a Bernoulli relative entropy. The shape is shared; the numbers are unrelated (PSLQ finds no relation with coefficients `<= 10^4`). Label: NUMEROLOGY.
6. **Moonshine.** §3 gives a verified chain.
   * **Reading of the hint.** "S(5), S(6)" are `S_5` and `S_6`: `S_5 = PGL(2,5)` acts on the `5+1` points of `P^1(F_5)`, and the outer automorphism of `S_6` swaps the two classes of `S_5`.
   * **The rungs.** The hexad stabilizer in `M_12` is `S_6`, acting through that twist, and the dodecad stabilizer in `M_24` is `M_12`, acting through *its* twist. Then come the Golay code (`24 = 23+1`), Leech (196560), `V_Leech` (300+24+196560) and FLM ((300+98280)+98304). Finally `2^(1+24).Co_1` gives 1+299+98280+98304, and the Griess algebra gives 196884 = 1+196883.
   * **Side results.** The `M_24` Frame shape `1^2 11^2` is the repository's level-11 newform `eta(t)^2 eta(11t)^2`, a REAL link. The `24 = 23+1` rung is THM-481.
7. **Hexagons.** Seven unit hexagons do fit in side `5/sqrt3` (exact check, 12 contacts): it is the honeycomb flower, rotated 30 degrees against the container.
   * The record is Morandi's, from April 2015, listed on Friedman's page as "smallest known".
   * The 2026 global-optimization paper did not improve it, and its dual bounds for `n >= 6` stay at the area bound `sqrt7 = 2.6458`.
   * Optimality: OPEN.
8. **The `{log n}` construction.** The source was not found.
   * **Reading.** Read with `log_b` and row order, it is a logarithm table: round `r` appends the `r`-digit mantissa logs of the `r`-digit numbers.
   * **Statistics.** The checkpoints are exact rectangles. The statistics fit normality with bias about `C_b/r`, like Champernowne. Column 1 follows Benford's law exactly.
   * **Relation to the session.** It lives on the *top-digit* (archimedean) side. The session's endgames read *low* p-adic digits.
   * **2-versus-3 slice.** In base 3, the rows `n = 2^k` carry the top ternary digits of `2^k`. Their column 1 is a 3-arc coding of the rotation by `log_3 2`.

## 1. The items: exact statements and methods

**1.1 Radchenko–Wheeler, arXiv 2609.21892v1 (math.NT, 18 Sep 2026)** [P].
* **Theorem 1.** For hyperbolic `gamma` in `SL_2(Z)` (trace `N+2`, `c > 0`), the values `F_gamma^pm(u)` of Faddeev's *modular quantum dilogarithm* `Phi_(gamma,m,n)(z;tau) = (q^m e(z); e(tau))_inf / (q~^n e(z/(c tau+d)); e(gamma tau))_inf` at the fixed point `tau` of `gamma` are algebraic.
* **Stark units.** These values give the Stark–Shintani ray class invariants: `eps_c = |F^+_(gamma_r)(u)|^(-t)`, with `t` in `{1,2}`. So Stark units of real quadratic fields are algebraic.
* **Theorem 2.** The values satisfy an explicit, overdetermined polynomial system over `Z[zeta_2N, 1/N]`: `F(0) = eps^(+-1/2)`, a reflection identity, and discrete Fourier identities.
* **Definition 4.** A *finite quantum dilogarithm* on a metric group `G` is a solution of the pentagon relation with a nonzero defect `C` (Andersen–Kashaev have `C = 0`).
* **Theorem 8.** Its values are algebraic. The proof goes through Izumi near-group fusion categories and Ocneanu rigidity. For prime `N` (Theorem 9) there is a short argument by Tao's uncertainty principle.
* **Theorem 7 and Remark 6.** These prove the Appleby–Flammia–Kopp quadratic (twisted convolution) conjecture [AFK, Conj. 1.35], which is motivated by Zauner's SIC-POVM conjecture.
* **Byproduct.** An infinite family of irrational near-group fusion categories, for every `G = Z/n x Z/m`.
* **Finite objects.** The paper's "finite" objects are functions on finite groups. Its only finite Pochhammers are the standard ones in Appendix A.1.

**1.2 Karingula–Lovett, arXiv 2609.20979v1** [P]. Vectors `v_i` in `R^d` with `||v_i||_2 <= 1` admit signs with `||sum eps_i v_i||_inf <= 36`.
* **Method.** An elementary simplification of Guo–Fang–Lu (Sept 2026, `C = 3 sqrt(2 pi)`).
  * The central quantity is the *shift distance* `Delta(P,u) = d_TV(P, P+u)`.
  * **Lemma 1.4.** If `Delta(P, 6 v_i) <= 1/3` for all `i`, then some signs put `mu(P) + sum eps_i v_i` in `conv supp P`. The proof is an induction with one split per vector, keeping an extra bit coordinate.
  * **Lemma 1.5.** Rounding a continuous product density to a rational grid gives a mean-zero `P` on `[-6,6]^d` with `Delta(P, v_i) <= 1/3`.
* There is no polynomial-time algorithm, and the strong (prefix) Komlós conjecture stays open.

**1.3 Ezeunala–Jiang, arXiv 2609.17266v1** [P]. A deterministic polynomial-time algorithm gives signs with `||sum s_i H_i|| <= 13 ||sum H_i^2||^(1/2)` for rational Hermitian `H_i` of rank at most 1. It solves Weaver's `KS_2` with constant `13/2`.
* The walk runs on fractional signings `x` in `[-1,1]^N`.
* It is guided by a potential `R(x) + lambda sum psi_i(x)`: `R` is the value of an SDP with resolvent-weighted reservoirs, and `psi_i = (1 - x_i^2)^(1/3)`.
* The analysis uses Bansal's covariance lemma.
* Marcus–Spielman–Srivastava gives a better constant non-algorithmically.

**1.4 Erdős 859** (Kruer–Kohlmeyer exposition of an accepted 36,197-line Lean proof, conjectures.io, dated 21 Sep 2026) [P; Lean not read].
* **Statement.** Let `d_t` be the density of the `n` for which `t` is a sum of distinct divisors of `n`. Then `liminf (-log d_t)/log log t = delta` and `(log t)^delta d_t -> 0`, with `delta = 1 - (1 + log log 2)/log 2`. So `d_t ~ c_1 (log t)^(-c_2)` is impossible.
* **Upper bound.** A reduction to divisor intervals `D(y, 2y)` (Ford-type bounds, log-support measure `L_h(S)`), with a `(log log)^(-5/4)` saving.
* **Lower bound.** A Poisson tilt at the critical `lambda* = 1/log 2`, together with:
  * Hölder change of measure;
  * conditioning on adaptively chosen prime windows;
  * `H`-proper divisor cubes;
  * Fourier inversion;
  * practical-number lifting.

  The rate optimizes to `lim_(q -> 1) R(lambda*, q) = delta`.

**1.5 Erdős 1062(ii)** (Kruer–Kohlmeyer exposition of an accepted Lean proof whose final declaration sits at line 74,198, dated 22 Sep 2026) [P; Lean not read].
* **Statement.** Let `f(n)` be the size of the largest *fork-free* subset of `{1..n}`, meaning no member divides two others. Then `f(n)/n -> L`, and `L` is irrational.
* **Finite formula (4).** `f(n) = sum_((q,6)=1, q <= n) g(floor(n/q))`, where `g = R - B` is a 2-3-smooth rank count with a window correction.
* **Density.** `L = (1/3) sum Delta_j/(j+1) = (32 + sum a_k/4^k + sum b_k/3^k)/180`, where `a_k` and `b_k` are read from tables comparing `4^k` with `3^floor(log_3 4^k)` and `3^k` with `4^floor(log_4 3^k)`.
* **Irrationality.**
  * Rationality gives fixed-dimension (`m <= 404`) *return approximations*: `0 < |sum d_i u_i - P| <= e^(-T - s/200)` with `u_i` in `U_(2,3)`, height `<= e^(27 s)` and `u_i >= e^(-T)`.
  * These contradict the rigidity bound `|E| >= c t H^(-1/10800)` (`oneForm_integer_rigidity`), which is derived from the Lean theorem `rational_three_place_subspace_cover`: a finite rational-hyperplane cover for small forms at `inf`, 2 and 3.
  * Nonvanishing of the errors is proved separately.
  * **The exact scope of the cover theorem is UNVERIFIED here.**

**1.6 The pasted `{log n}` construction.** It interleaves the base-`b` digits `D(n, j)` of `{log n}` round by round. After round `r` the prefix is exactly the rectangle `2 <= n <= b^r`, `1 <= j <= r`, with `T_r = r(b^r - 1)` digits. Source not located; see §5.

**1.7 The Moonshine exercise and 1.8 the hexagon claim.** See §3 and §4.

**1.9 AMM 12592** (cross-reference only; a separate lane handles it). `long-mathematics/deterministic-von-neumann-fair-extractor`, README [P]:
* `T(L) <= C_* L - log_5 L + 6`, with `C_* = 1 + 2 log(phi)/log 5 = 1.5979874...`.
* The constant is THM-3027's `gamma*` plus 1.
* The manuscript is attributed to the same author as the p-adic zeta preprint audited in [p-adic-zeta audit](../reference/p-adic-zeta-irrationality-source-audit-20260825.md).

## 2. The interplay map

**Typing** follows the repository convention: map / preserved predicate / loss / sidecar / test. **Labels:** REAL (a proved or exactly checked link), ANALOGY (shared structure, no transfer) and NUMEROLOGY (shared numbers only).

| # | item | target | type | precise connection | label | test |
|---|---|---|---|---|---|---|
| 1 | Erdős 1062(ii) | Theorems S, D (PC) | map | `L` is a pair of rotation-coding series at `log_3 4 = 2 log_3 2` and `log_4 3`; Sturmian and rotation-coding Bernstein numbers are the 2-adic analogue | ANALOGY | T1 |
| 2 | Erdős 1062(ii) | foundry v4 missing mechanism; HARD §1.5 | preserved predicate + **loss** | preserved: rationality forces nonzero small S-unit forms with one integer non-S-unit `P`; lost: the exchange rate (`eta = 1` there, `eta = beta log2 3 > 1` here) | **REAL** | T1, §2.1 |
| 3 | Erdős 1062(ii) Lean cover | formalizing S, D | sidecar | S and D need only Theorem R plus word combinatorics; a rational three-place cover formalizes Proposition T′; it cannot reach `Y3` (no-go) | REAL (T′ conditional) | §2.1 |
| 4 | Erdős 1062(ii) | Mahler 3/2 (THM-2228) | map | Mahler's `||(3/2)^n||` bound uses Ridout's theorem, the one-form, three-place case of the same cover (Mahler 1957; citation UNVERIFIED, from memory) | REAL | – |
| 5 | Radchenko–Wheeler | Theorem Y, HC1 | preserved predicate + loss | preserved: first-order q-difference equations (Padé `(qT;q)_n`, `Phi_b`); lost at cubes; also algebraicity versus irrationality, and `|q| = 1` versus `|q|_2 < 1` | ANALOGY | T2 (Jacobi triple product), T3 |
| 6 | Radchenko–Wheeler | AMM 12592 | numerology | the `N = 5` finite quantum dilogarithm has `E(0) = phi^2`, and `C_* = 1 + log_5 E(0)` | NUMEROLOGY | T3 |
| 7 | moonshine / Tate curve | Theorem Y | map | `Phi_2(Y)` is affine in `theta_3(rho)`; `j(rho^2)` is rational in `theta_3(rho)` and `psi(rho^2)`; BDGP gives Corollary M | **REAL** | T2 |
| 8 | moonshine / Tate curve | HC1 | loss | `sum rho^(k^3)` is not modular: no Tate-curve or Mahler–Manin handle | – | – |
| 9 | Komlós; Kadison–Singer | C2 (HYP-9123) | preserved predicate + loss | bounded prefix discrepancy is shared; the quantifiers are opposite (a signing exists versus no arithmetic realization); Terras makes C2 a limit-only statement | ANALOGY | T5 |
| 10 | Erdős 859 | session LD rates | numerology | `delta = Q(1/ln 2)` (Poisson) against `D(log_3 2 || 1/2) = 0.0347`, `-ln 0.758751 = 0.2761` and `c* = 0.1144`; no relation found | NUMEROLOGY (shape: ANALOGY) | T4 |
| 11 | Erdős 859 method | HYP-9122, chain budget | sidecar | two possible tools: (a) the finite identity for conditioning on an adaptively chosen window; (b) log-support measures `L_h(S)` against the covering of ratios `2^K/3^s` by loops | ANALOGY | – |
| 12 | moonshine chain | level-11 thread | map | `M_24` Frame shape `1^2 11^2` (`x -> 2x` on `P^1(F_23)`) gives `eta(t)^2 eta(11t)^2`, the level-11 newform | **REAL** | M4 |
| 13 | moonshine chain | THM-481 (Paley and Golay) | map | the `24 = 23+1` rung is the extended QR code of the Paley tournament on `F_23` | **REAL** | M2 |
| 14 | moonshine (`Out(S_6)`) | W(n), Rédei–Berge | test | the twist `tau` on `Lambda^6` keeps `U_T` in the family only when `c3 = 2 d33`, and then fixes it | ANALOGY (plus a FINITE-EXACT refutation) | T6 |
| 15 | hexagons | LRC14 hexagonal witness (HYP-3706/3726) | numerology | the flower centres form a complete residue system of `Z[w]/(3+w) = Z/Phi_6(3)`, and `x3` rotates them; `3k(k-1)+1 = Phi_6(n)` is equivalent to `X^2 - 3Y^2 = -2` | NUMEROLOGY | T7 |
| 16 | `{log n}` | top versus low digits; Theorem S | map | base-3 rows `n = 2^k` carry the top ternary digits of `2^k`; column 1 is a 3-arc rotation coding by `log_3 2` | ANALOGY (the identity is REAL) | L5 |
| 17 | AMM 12592 | p-adic zeta audit | sidecar | same author, another p-adic irrationality program; its holonomy method needs modular or D-finite input, which the lacunary `sum z^(k^3)` lacks | ANALOGY | – |
| 18 | sixty clocks (THM-4035) | moonshine | numerology | `60 = |A_5|` (`PSL(2,5)`, and the index-11 `A_5` in `PSL(2,11)`, M3) against `lcm(1..6) = 60` | NUMEROLOGY | M3 |
| 19 | Pell/Beatty crossover (HYP-2456) | Erdős 1062(ii); log digits | analogy | rotation codings by finitely many arcs: `a_k` has 7 arcs, `b_k` 5, the L5 word 3 | ANALOGY | T1, L5 |

### 2.1 Erdős 1062's three-place cover against the missing mechanism

**The shared engine.**
* Theorem R/D. A rational `Phi_T(w)` and a periodic approximant `P/Q` agree to `lambda` bits (the isometry), and heights are controlled by Lemma H′.
* Erdős 1062. A rational `L` and the return approximations give `sum d_i u_i - P` (small, nonzero) with `u_i` in `U_(2,3)`.
* **The Subspace form.** In both, the coordinates other than `P` are S-units for `S = {inf, 2, 3}`, so their product over `S` is 1. The only cost is the non-S-unit coordinate `P` at the archimedean place.

**The exchange rate.** Write `eta` for the height bits paid per bit (or letter) of precision.
* **Bounded-coefficient series** such as `sum a_k 4^-k`: each letter buys `log 4` of precision and costs `log 4` of height, so `eta = 1`.
  * Liouville or Subspace then needs only `Dio > 1`, which every rotation coding has (Lemma J′: `Dio >= 1 + 1/k`).
  * Redundant digit sets admit non-EP words with rational sums, so the real difficulty is *nonvanishing*. The exposition says the errors are proved nonzero.
* **Supercritical parity words.** Each letter buys 1 bit of 2-adic precision and costs `eta = beta log2 3` bits of height (`M_z = 3^(#1s)`).
  * The Subspace product is `|P| 2^(-lambda) ~ H^(1 - lambda/(eta|UV|))`, which is HARD §1.5's `Dio > eta`.
  * The 3-adic place contributes nothing, because `P` has no forced 3-divisibility.

So the missing mechanism of foundry v4 is not a missing *place*: the preserved predicate is already exactly the S-unit/Subspace form. What is missing is approximants whose height beats `eta`. Theorem Y's Padé rationals have height about `2^(0.95 lambda)` against the natural `2^(1.43 lambda)`, and so evade the loss. Row 2 of the table: **REAL**, with the loss exactly located.

**T1 (FINITE-EXACT).**
* Formula (4) equals an exact ILP brute force of `f(n)` for every `n <= 70`.
* `L = 0.6729656511994893474712279198522096447528`, and `f(10^6)/10^6 - L = -1.65e-6`.
* `a_k = A(frac(k log_3 4))` for `k < 3000`, with arcs at `log_3` of `16/15, 10/9, 4/3, 3/2, 2, 8/3`. The only endpoint hit is `k = 1`, where `4/3` goes to the next row.
* `b_k = B(frac(k log_4 3))`, with arcs at `log_4` of `4/3, 8/5, 5/3, 2`.
* The factor complexities are linear.

**Can their Lean cover formalize Theorem S/D?** It is not needed.
* Theorem S and Theorem D (with the Sturmian options A/B) are elementary: Theorem R (the 2-adic isometry plus Lemma H) and combinatorics on words.
* Corollary D1 also needs Bugeaud–Kim's Sturmian repetition theorem. That is the heavy part to formalize, and the cover does not help with it.
* The cover can serve the *two-place* strengthening below. It is rational, so it cannot give the transcendence version C1, which needs algebraic coefficients.

**Proposition T′ (PROVED modulo the CITED rational p-adic Subspace Theorem; not audited).** Let `T = x/2, (3x+1)/2`. Let `w` be strictly supercritical (`eta > 1`), of bounded discrepancy, not EP, with `Dio(w) > 1`. Then `x := Phi_2(w)` and `y := Phi_R(w)` are not both rational.

*Proof.* Suppose `x = a/b` and `y = c/d`.
1. **Setup.** Take patterns `U V^e` with `lambda >= rho |UV|`, where `rho > 1`. Set `X = (X1, X2, X3) = (M_U 2^p, M_U M_V, P)`, with `p = |V|`. Then `Q = X1 - X2` and `r = P/Q = Phi(U V^inf)`.
2. **Error bounds.** `|x - r|_2 = 2^-lambda` (isometry), and `|y - r| <= K p 2^(-(eta - 1) lambda)`.
   * The tail of `w` past `lambda` is at most `C 2^lambda/3^(a_lambda)` (bounded discrepancy), and the tail of `U V^inf` is at most `2^lambda/3^(a_lambda)` times `p/(1 - 2^p/M_V)`. Both are `O(p 2^(-(eta - 1) lambda))`, since `3^(a_lambda) >= 2^(eta lambda - O(1))`.
   * That factor is uniformly bounded, because every witnessing `V` is supercritical. For long `V`, bounded discrepancy gives `a_V >= alpha p - 2C`. A short `V` repeated unboundedly (`e -> infinity`) must have `a_V = alpha p` exactly, since otherwise the discrepancy along `V^e` grows linearly. So `M_V = 2^(eta p) > 2^p`.
   * Only finitely many short `V` occur.
3. **The forms.**
   * At `inf`: `X1`, `X2`, `c(X1 - X2) - d X3 = d Q (y - r)`.
   * At 2: `X1`, `X2`, `a(X1 - X2) - b X3 = b Q (x - r)`; `Q` is odd.
   * At 3: `X1`, `X2`, `X3`.
4. **The product.** `X1` and `X2` are S-units, so the product is at most `|d| K p |Q| 2^(-(eta - 1) lambda) 2^(-lambda)`.
   * The gcd of `X` divides `3^(a_U)`, an S-unit, so passing to primitive vectors changes nothing.
   * With `H(X) = 2^(eta |UV| + O(log))` (by `M_(UV)` and Lemma H′), the product is `H^(1 - rho + o(1)) <= H^(-eps)`.
5. **The Subspace Theorem** puts all these `X` in finitely many rational planes `z . X = 0`. One plane contains infinitely many of them.
   * **Case `z3 != 0` and `p -> infinity`.** Here `r = -(z1 2^p + z2 M_V)/(z3 (2^p - M_V))` tends to `z2/z3` both in `R` (`2^p/M_V -> 0`) and in `Q_2` (`2^p -> 0`). So `x = y = z2/z3`. But `x != y` (hard note Proposition M(iv): `L = 0` would force an EP orbit).
   * **Case `z3 != 0` and `p` bounded.** `r` takes finitely many values and `r -> x` 2-adically, so `r = x`, contradicting `v_2(x - r) = lambda < infinity`.
   * **Case `z3 = 0`.** Then `z1 2^p = -z2 3^(a_V)` fixes `(p, a_V)`. Restrict to the plane with forms `{X1, c(1 + z1/z2) X1 - d X3}` at `inf`, and similarly at 2 and 3. The product is at most `|Q| 2^(-eta lambda)` with `H ~ 2^(eta |U|)` and `lambda > |UV| > |U|`. The Subspace Theorem in dimension 2 leaves finitely many lines `X3 = kappa X1`, so again finitely many `r`, and we are in the previous case. ∎

**Collatz reading (sidecar).** A divergent integer `n` whose parity vector `w` is a strip word with `1 < Dio(w) <= eta` must have irrational `L = n - Phi_R(w)`, the constant of Proposition M. Theorem D already excludes `Dio > eta`. Random strip words and `Y3` have `Dio = 1`, so this says nothing about them.

**No-go for HC1 (PROVED, elementary).** Let `x = sum_(k >= 1) rho^(k^3)`, `rho = 2^10/3^9` and `eta = log2(3^9)/10 = 1.4265`.
* **Partial sums.** `S_K = N_K/3^(9K^3)` has an S-unit denominator and a 3-adic unit numerator of size about `3^(9K^3)`. The three-place product for `(N_K, 3^(9K^3))` is `3^(9K^3) 2^(-10 (K+1)^3) = H^(1 - 1/eta + o(1))`.
* **Periodic approximants** (`Dio(Y3) = 1`) give the same `H^(1 - 1/eta + o(1))`.
* Since `1 - 1/eta = 0.299 > 0`, no Subspace cover over `{inf, 2, 3}` applied to natural-height approximants reaches `Y3`. This matches HYP-9127's "Liouville would need `log2(3^9)/10 < 1`".

### 2.2 Theorem Y's finite q-Pochhammers against Radchenko–Wheeler

**T2 (FINITE-EXACT).** The 2-adic Jacobi triple product holds mod `2^4000`:

`sum_(k in Z) rho^(k^2) = prod_(m >= 1) (1 - rho^(2m)) (1 + rho^(2m-1))^2`.

So the square-swap Bernstein number is an *infinite* q-Pochhammer value, with `q = rho^2` and `|q|_2 = 2^-20`. Theorem Y's Padé polynomials are the finite Pochhammers `R_n(T) = (qT;q)_n`, with `q = rho^-2`.
* In Radchenko–Wheeler, `Phi_b` is a ratio of infinite Pochhammers evaluated at `|q| = 1` (real quadratic `tau`), where it "should" be 1. The special values satisfy a pentagon system with defect, which forces algebraicity.
* **T3.** For `E_a` on `Z/5` (Gaussian `zeta_10^(3x(x+5))`), identity (29) holds with constant `C = -phi^2`, and Theorem 5(i)–(ii) hold (40 digits).

**Preserved predicate.** "The value is controlled by a first-order q-difference equation." Tschakaloff's `F(x) = 1 + rho x F(rho^2 x)`, Faddeev's pair of difference equations in `q` and `q~`, and the Jacobi triple product all live there.

**Losses.**
1. **Cubes.** The coefficient ratio `rho^(3k^2 + 3k + 1)` has an exponent quadratic in `k`, so it is not `R(q^k)` for any `q`. Squares are q-hypergeometric: `rho (rho^2)^k`.
2. **Direction.** Radchenko–Wheeler prove *algebraicity*; HC1 needs *irrationality*.
3. **Place.** Their `|q| = 1`; ours `|q|_2 < 1`. Their concluding remarks name the behaviour of solutions of (7) over non-archimedean fields as an open question.

Label: ANALOGY. Radchenko–Wheeler supply no approximants for `Y3`.

### 2.3 The Tate curve at `q = rho^2` and the verified Mahler–Manin citation

Take nome `rho`, so the Tate parameter is `q = rho^2 = 2^20/3^18`. Note that `sqrt(rho)` is not in `Q_2`, since `3` is not a 2-adic square. The curve `E_(rho)` of the owner's suggestion is 2-isogenous, and `j(rho)` and `j(rho^2)` satisfy the level-2 modular equation.

**T2 (FINITE-EXACT, mod `2^3940`).** `j(q) q = 1 + 744 q + 196884 q^2 + ...` evaluated at `q = rho^2` equals the classical formula `256 (1 - l + l^2)^3/(l^2 (1 - l)^2)`. Here:
* `l = theta_2^4/theta_3^4 = 16 rho psi(rho^2)^4/theta_3(rho)^4 = 2^14 U`, with `U` a 2-adic unit;
* `theta_3(rho) = sum_(k in Z) rho^(k^2)`;
* `psi(x) = sum_(n >= 0) x^(n(n+1)/2)`.

Also `v_2(j(rho^2)) = -20`, and the 196884 term has valuation 22.

**Swap identities (T2, mod `2^3000`, from the parity words themselves).**
* `Phi_2(Y) = -1 - 512/(3^9 - 2^10) - (256/3^9) sum_(k >= 1) rho^(k^2)`.
* `Phi_2(Y_pr)` is the same with `sum rho^(k(k+1)) = psi(rho^2) - 1`. Here `Y_pr` swaps the block `1^8 0 1` in at the pronic indices `k(k+1)`.

**The citation, verified.** K. Barré-Sirieix, G. Diaz, F. Gramain, G. Philibert, *Une preuve de la conjecture de Mahler–Manin*, Invent. Math. 124 (1996) 1–9, DOI 10.1007/s002220050044.
* Waldschmidt's Bourbaki exposé 824 (Astérisque 245, 1997), read: Théorème 2 says `J(alpha)` is transcendental for algebraic `alpha` with `0 < |alpha| < 1`. The introduction says the Saint-Étienne team "ont résolu en même temps le problème analogue p-adique, qui avait été posé par Manin".
* The same exposé states Nesterenko's `tr.deg >= 3` for **complex** `q` only (Théorème 4). The p-adic analogue given is the Chudnovsky–Bertrand `tr.deg >= 2` (Théorème 3).

**Corollary M (PROVED modulo BDGP).** `j(rho^2)` is transcendental, and `l` is a rational function of `theta_3(rho)` and `psi(rho^2)` over `Q`. So `theta_3(rho)` and `psi(rho^2)` are not both algebraic. Equivalently, **`Phi_2(Y)` and `Phi_2(Y_pr)` are not both algebraic.** Label: REAL. This is a transcendence-type statement that Lemma P (irrationality only) does not give.

**Losses.** Transcendence of `Phi_2(Y)` alone would need a p-adic Nesterenko theorem (OPEN). `sum rho^(k^3)` is not a modular form, so HC1 has no Tate-curve counterpart.

### 2.4 Komlós and Kadison–Singer against balanced parity words

The C2 strips (HYP-9123) are the words with bounded prefix discrepancy `|a_s - alpha s| <= C`. Komlós in its prefix ("strong") form asks for signings with bounded prefix sums. In dimension 1 that is trivial (greedy), and in general it is still open (Karingula–Lovett).

**T5 (FINITE-EXACT).**
* The Terras bijection (each parity word of length `s` is realized by exactly one class mod `2^s`; Terras 1976 [R]) was rechecked for `s <= 16`.
* Strip words of length 60 number 64 for slope `9/10` and width 1 (0.100 bits per letter, matching the hard note's C1), 46656 for width `3/2` (0.258), and `1.6e11` for slope `3/4` and width 2 (0.612).

So every finite balanced word is realized arithmetically. Any finite-horizon discrepancy statement is automatically satisfied, and C2 is purely the 2-adic limit, the coupled Z-number condition of Proposition M. In foundry terms, discrepancy-existence theorems are DEFECT-blind. The only REAL discrepancy input for Collatz remains the in-house capacity theorem at the critical slope.

Ezeunala–Jiang's potential-guided deterministic walk is an algorithmic cousin of the foundry's "sound certificate search". That is only a design idea for potential-guided chain certificates in E-SCC. Label: ANALOGY.

### 2.5 Erdős 859's `delta` against the session's large-deviation rates

**T4.**
* **Ford's constant.** `delta = 0.086071332055934206888 = Q(1/ln 2)`, with `Q(l) = l ln l - l + 1`. This is the Poisson relative entropy at the tilt where `2^omega` matches the log scale.
* **The session's rate.** `D(log_3 2 || 1/2) = ln 2 - H(ln 2/ln 3) = 0.0346882` nats `= 0.0500445` bits `= 1 - h(log_3 2)` (with `h = 0.9499555`). This is the Bernoulli relative entropy at the frequency where `3^a` matches `2^s`.
* **Ratios.** `delta/D = 2.4813` and `delta/c* = 0.7524`.
* **PSLQ.** No integer relation among `(delta, D, c*, 1)` with coefficients up to `10^4`. The sanity check recovers `delta ln 2 = ln 2 - 1 - ln ln 2`.

Label: NUMEROLOGY. The shape, "relative entropy at a log-2 threshold", is an ANALOGY.

Methodological sidecar (row 11). Their exact finite identity for conditioning on an *adaptively chosen* window is the tool that the E-SCC chain-budget heuristics (Q2, where escapes are chosen after observing digits) lack.

### 2.6 Past threads

* **Level-11 eigenform (REAL).**
  * **Frame shape (M4).** In `M_24 = <PSL(2,23), delta>`, the map `x -> 2x` has Frame shape `1^2 11^2` and `x -> x+1` has `1.23`.
  * **Newform check.** The eta product `eta(t)^2 eta(11t)^2` has `a_p = p + 1 - #X_0(11)(F_p)` for all primes `p < 60`, `p != 11`. Also `sum_r a_(2^r) t^r = 1/(1 + 2t + 2t^2)`, which is the input of the level-11 note's `1/15` law.
  * **Mason's correspondence** [R: Mason, Contemp. Math. 45 (1985) 223–244; Dummit–Kisilevsky–McKay 1985]. It sends 21 `M_24` Frame shapes to multiplicative eta products. So the repository's level-11 form is the `M_24` "order-11" form.
* **THM-481 (REAL).** The Paley-tournament extended QR code at `q = 23` is the binary Golay code of rung 2 in §3.
* **W(n) and Rédei–Berge, T6.**
  * Let `tau` be the `Out(S_6)` twist on `Lambda^6`: `p_(3,1,1,1) <-> p_(3,3)`, `p_(2,1^4) <-> p_(2^3)` and `p_6 <-> p_(3,2,1)`. It commutes with `omega` and preserves `H` (since `zeta(p_lambda) = 1`).
  * For `n = 6`, `U_T = p1^6 + 2 c3 p3 p1^3 + 2 c5 p5 p1 + 4 d33 p3^2` (Grinberg–Stanley via THM-467). The OCF `H = 1 + 2 c3 + 2 c5 + 4 d33` was checked on all 32768 labelled tournaments.
  * `tau(U_T)` is again a Rédei–Berge function **only** for the three `tau`-fixed triples `(c3, c5, d33) = (0,0,0), (2,0,1), (8,6,4)`, which cover 1040 labelled tournaments. So `tau` is not a symmetry of the family (REFUTED). This is the only `n` where such a twist exists.
* **Hexagonal threads (NUMEROLOGY, T7).**
  * The LRC14 covering-min witness is multiplication by `n` (a sixth root of unity) on `Z/Phi_6(n) = Z[w]/(n + w)` (HYP-3706, HYP-3726). At `n = 3` this is exactly the 7-flower: its centres `{0, +-1, +-w, +-w^2}` form a complete residue system mod `3 + w`, and `x3 = x(-w)` rotates the petals by 60 degrees.
  * The flower sizes `3k(k-1)+1` equal `Phi_6(n)` exactly on the Pell sequence `X^2 - 3Y^2 = -2`: 1, 7, 91, 1261, 17557, and so on. `183 = Phi_6(14)` is not among them.
  * There is no transfer to packing or to LRC.
* **p-adic zeta audit (ANALOGY or sidecar).** Long's Kubota–Leopoldt irrationality program is the nearest "2-adic irrationality of an explicit number" effort in the repository. Its arithmetic-holonomy engine needs a modular or holonomic source. `theta_3` has one (it is modular in `tau`), while the lacunary `sum z^(k^3)` does not: the same loss as rows 5 and 8.
* **Sixty clocks and Kakeya (NUMEROLOGY).** `60 = |A_5|` appears at the `p = 5` and `p = 11` rungs (M3). The `lcm(1..6) = 60` of THM-4035 is an unrelated 60.
* **Pell/Beatty crossover (ANALOGY).** All the words met here, `a_k` (7 arcs), `b_k` (5 arcs) and the L5 word (3 arcs), are codings of a single rotation by finitely many arcs, the class of HYP-2456's address words. Only the rotation-coding lemmas (Lemma J′) transfer.

## 3. The moonshine exercise: from 6 = 5+1 to 196884 = 196883+1

**Reading the hint.** "S(5)" and "S(6)" are the symmetric groups. Their relation is the transitive embedding `S_5 = PGL(2,5)` in `Sym(P^1(F_5)) = S_6`, which is the source of `6 = 5 + 1`. The "exceptional automorphism" is the outer automorphism of `S_6`, the only symmetric group that has one.

If "S(5), S(6)" were instead read as the Steiner systems `S(5,6,12)` and `S(5,8,24)`, the chain below from rung 1 on is unchanged, but then the exceptional automorphism belongs to `M_12`, not to `M_24`.

Every numerical statement below is FINITE-EXACT (`.out`, M1–M6). Group-theoretic facts that are not recomputed are marked [R].

0. **Rung 0: `6 = 5 + 1` (M1).**
   * `PGL(2,5)` acts on the 6 points of `P^1(F_5)` sharply 3-transitively, with no transposition. It is a transitive `S_5` inside `S_6`, not conjugate to a point stabilizer.
   * The coset action of `S_6` on `S_6/PGL(2,5)` is an automorphism `psi`. It swaps the classes `2.1^4 <-> 2^3`, `3.1^3 <-> 3^2` and `6 <-> 3.2.1`, so it is outer.
   * The natural permutation character is `1 + chi^(5,1)`, and the exotic one is `1 + chi^(2,2,2)`. We have `<pi, pi'> = 1`: **`6 = 1 + 5` in two inequivalent ways, exchanged by `Out(S_6)`.**
   * The same `p + 1` pattern gives the hexacode: the extended QR code of length `5 + 1` over `F_4`, `[6,3,4]`, with weights `1 + 45 y^4 + 18 y^6`.
1. **Rung 1: `12 = 11 + 1 = 6 + 6` (M2, M3).**
   * **Ternary Golay code.** The extended QR code of length `11 + 1` over `F_3`: `[12,6,6]`, with weights `1 + 264 y^6 + 440 y^9 + 24 y^12`. Its 132 weight-6 supports form `S(5,6,12)`.
   * **`M_12`.** The group `<PSL(2,11), (2 10)(3 4)(5 9)(6 7)>` has order `95040 = 12·11·10·9·8` and is sharply 5-transitive, hence `M_12` [Jordan, R]. It preserves an `S(5,6,12)`.
   * **The key link.** The setwise stabilizer of a hexad is `S_6` (order 720). It acts on the hexad and on the complementary hexad (itself a hexad) by two faithful representations whose cycle-type correspondence is *exactly* the table of `Out(S_6)` from rung 0. So `12 = 6 + 6` carries `S_6` in its two actions. The 720 elements of `M_12` exchanging the two halves complete a pair stabilizer of order 1440. It contains elements of order 8 and 10, so it is not `S_6 x 2`: it is the extension of `S_6` by its exceptional automorphism.
   * `M_12` is 2-transitive, so `12 = 1 + 11`.
   * **The `p = 11` analogue of the exceptional 5-point action of `PSL(2,5) < PGL(2,5) = S_5`.** `PSL(2,11)` contains Galois's exceptional `A_5` of index 11, so it acts on 11 points as well as on the 12 points of `P^1(F_11)`. This `A_5` is transitive on `P^1(F_11)`. The two actions extend to the two classes of `M_11` in `M_12`, which `Out(M_12)` exchanges [R].
2. **Rung 2: `24 = 23 + 1 = 12 + 12` (M2, M4, M4b).**
   * **Binary Golay code.** The extended QR code of length `23 + 1` over `F_2` (THM-481's Paley code), with weights `1 + 759 y^8 + 2576 y^12 + 759 y^16 + y^24`. Its 759 octads form `S(5,8,24)`: `42504 = 759·56`.
   * **`M_24`.** The group `<PSL(2,23), delta>` has order `244823040`. Here `delta` is `x -> x^3/9` on the squares and `x -> 9x^3` on the non-squares, fixing `0` and `inf`.
   * **Dodecad stabilizer.** `M_24` is transitive on the 2576 dodecads, so a dodecad stabilizer has order `95040`; the subgroup generated by random stabilizing elements has that order, so it is the full stabilizer, `M_12`. Its actions on the dodecad and on the complement differ in cycle type (`4^2 1^4 <-> 4^2 2^2` and `8.2.1^2 <-> 8.4`), so they are twisted by the outer automorphism of `M_12`. This is exactly how rung 0 sat inside rung 1.
   * **Direct route from 6 to 24** [R: Conway's MOG, SPLAG ch. 11]. The hexacode builds the Golay code on a `6 x 4` array.
3. **Rung 3: Leech (M5).** The Conway–Sloane conditions [R: SPLAG ch. 4 §11] give, in coordinates `/sqrt8`, minimal vectors of three shapes:
   * `(+-4, +-4, 0^22)`: 1104;
   * `(+-2^8, 0^16)` on octads with an even number of minus signs: `759·2^7 = 97152`;
   * `(-+3, +-1^23)` from Golay codewords: `24·2^12 = 98304`.

   The total is **196560**, and every vector satisfies the membership conditions. The inner products with a fixed minimal vector are distributed as `1, 4600, 47104, 93150, 47104, 4600, 1`. Finally `Theta_Leech = E_4^3 - 720 Delta = 1 + 196560 q^2 + 16773120 q^3 + ...`, and `179280 + 720·24 = 196560`.
4. **Rung 4: `V_Leech` (M6).** The lattice VOA has character `Theta_Leech/eta^24`, and the identity `Theta_Leech/eta^24 - 24 = j - 744` holds through `q^11`.
   * Its weight-2 space is `Sym^2(h)` (`C(25,2) = 300`), plus `h_(-2)` (24), plus `e^alpha` with `alpha^2 = 4` (196560): **`300 + 24 + 196560 = 196884`**.
   * The count `324 = 300 + 24` is `p_24(2)`, the `q^2` coefficient of `1/prod(1-q^n)^24`.
5. **Rung 5: FLM orbifold** [R: Frenkel–Lepowsky–Meurman 1988]. `V♮ = V_Leech^+ + (V_Leech^T)^+`. Its weight-2 space is `(300 + 98280) + 24·2^12 = 98580 + 98304 = 196884`.
   * `h_(-2)` is `theta`-odd, so it drops out; `98280 = 196560/2` counts the pairs `+-alpha`.
   * The twisted sector contributes `h(-1/2)` tensored with the `2^12`-dimensional module of `2^(1+24)`.
   * The weight-1 space vanishes, which gives exactly `j - 744`.
6. **Rung 6: the Monster** [R: Griess 1982; Conway 1985; FLM]. `Aut(V♮) = M`. The involution centralizer `2^(1+24).Co_1` acts on the weight-2 space as `1 + 299 + 98280 + 98304`, and `300 = 1 + 299`, since the invariant form is the trivial summand of `Sym^2` of the 24-dimensional `Co_0` module. The Griess algebra is **`196884 = 1 + 196883`** as an `M`-module, with `196883 = 299 + 98280 + 98304 = 47·59·71`.
7. **Moonshine** [R: McKay 1978; Thompson 1979; Conway–Norton 1979; Borcherds 1992].
   * `j - 744 = q^-1 + 196884 q + 21493760 q^2 + 864299970 q^3 + ...`
   * `21493760 = 1 + 196883 + 21296876`.
   * `864299970 = 2·1 + 2·196883 + 21296876 + 842609326`.

**The "+1" ladder.** Each `+1` is the trivial summand of a 2-transitive or form-preserving action:
* `6 = 5 + 1` for `S_6`;
* `12 = 11 + 1` for `M_12`;
* `24 = 23 + 1` for `M_24`;
* `300 = 299 + 1` for `Co_1`;
* `196884 = 196883 + 1` for `M`.

The outer automorphisms of rungs 0 and 1 are realized inside the next rung by complementary halves: a hexad and its complement, then a dodecad and its complement.

NUMEROLOGY, flagged and not used: `5 -> 11 -> 23 -> 47` is a Cunningham chain (`p -> 2p+1`), and `47` divides `196883`.

**Side links.** Rungs 1–2 connect REALLY to two repository threads: the level-11 newform through the `M_24` Frame shape `1^2 11^2`, and THM-481 (§2.6). The `q`-expansion of `j`, with its 196884, converges 2-adically at `q = rho^2`. That is the Tate-curve bridge of §2.3.

## 4. The seven-hexagon packing

**Coordinates.**
* **Tiles.** Unit hexagons with vertices at angles 30 + 60k degrees, i.e. `c + (+-sqrt3/2, +-1/2)` and `c + (0, +-1)`. The centres `c` are `0` and `sqrt3 e^(i pi k/3)`, i.e. `(+-sqrt3, 0)` and `(+-sqrt3/2, +-3/2)`: the honeycomb flower.
* **Container.** A regular hexagon of side `s = 5/sqrt3`, with vertices at angles `0, 60, ..., 300` degrees, rotated 30 degrees relative to the tiles. Its apothem is `5/2`, so it is the set `n_k . v <= 5/2` for the six unit normals at 30 + 60k degrees.

**Exact check (all decisions in `Q(sqrt3)`, `.out` H1–H4).**
* All 42 tile vertices satisfy the six inequalities, with **12 equalities**: each outer tile touches two sides. For example, `(3 sqrt3/2, 1/2)` gives `(3 sqrt3/2)(sqrt3/2) + 1/4 = 5/2`, and `(sqrt3/2, 5/2)` touches the top side.
* The tile interiors are pairwise disjoint (exact separating axes).
* Shrinking the container by any `t > 0` makes *this* configuration infeasible (tested at `10^-6` and `10^-30`, exactly).
* **Answer: yes, 7 non-overlapping unit hexagons fit in side `5/sqrt3 = 2.886751345948`.**

**Literature status.**
* **Friedman's Packing Center**, "Hexagons in Hexagons" (read 2026-09-23). The row "6-7" lists `s = 5/sqrt3 = 2.886+`, "Found by Maurizio Morandi in April 2015". The page describes its entries as the "smallest known" containers. The picture matches the configuration above.
* **Berthold–Kamp–Mexi–Pokutta–Polik**, *Out-of-the-box global optimization for packing problems*, arXiv 2605.04850 (May 2026) [P].
  * Their Table 2 improves the `(6,6,n)` records for `n = 11, 12, 14, 15, 16, 17, 23`, and not for 6 or 7.
  * For polygon packing they state: "proving optimality for one and two inner elements and several instances with n ≤ 5, while for n ≥ 6, the dual bound generally remains at R_min". For `n = 7` the area bound is `R_min = sqrt7 = 2.6458`, and the record exceeds it by 9.1%.
* No source found proves optimality, and no 2026 optimality paper was found. **Status: record, optimality OPEN.**

**Family.** For `N_k = 3k(k-1) + 1` hexagons, the rotated `k`-ring flower fits exactly in `s_k = (3k-1)/sqrt3` (FINITE-EXACT, `k = 2, 3, 4`: 12, 18 and 24 contacts).
* `k = 3` is Friedman's 18–19 record, `8/sqrt3`.
* `11/sqrt3 = 6.3508530` agrees with the value listed for `n = 34` (6.35085+) to the displayed digits.

**Repository.** There is no hexagon-packing thread. The hexagonal threads (LRC14's `zeta_6` witness) relate only numerologically (§2.6, T7).

## 5. The `{log n}` construction

**Source.** Not located.
* Six web searches (normal-number, Champernowne/Copeland–Erdős, Benford and "digits of log n" phrasings) found nothing.
* A scan of the conjectures.io catalogue (`/problems`, `/results`) found nothing either.
* **UNVERIFIED.** What follows is analysis, not a source claim.

**What the rectangle is (L1).** Round `r` adds column `r` for the old rows `n <= b^(r-1)` (`b^(r-1) - 1` digits), plus all `r` digits of the new rows `b^(r-1) < n <= b^r`. This totals `T_r = r(b^r - 1)`, verified for `b = 2` (`r <= 20`), `b = 3` (`r <= 12`) and `b = 10` (`r <= 6`).

**The natural reading is `log_b`, which makes it a logarithm table.** The bound `n <= b^r` is aligned with `floor(log_b n) < r`. The new rows of round `r` are the `r`-digit numbers, each carrying its `r`-digit mantissa logarithm.
* **Resolution (L2).** Consecutive values of `log_b n` differ by at least `1/((n+1) ln b)`. For `b = 2` (`ln 2 < 1`), `r` digits therefore separate every `n <= 2^r`: 100% distinct 20-bit prefixes among the `2^19` new rows. The figure is 99.4% for `b = 3` and 77.4% for `b = 10`. For `b = 2` this is "the largest log table that `r` digits resolve".
* **Column 1 (L4).** It follows Benford's law exactly: `P(d) = (b^((d+1)/b) - b^(d/b))/(b - 1)`, e.g. `0.4142/0.5858` for `b = 2`.
* **Frequencies (L3).** In row order the block frequencies approach uniform with bias about `C_b/r`: `r · dev = 0.34, 0.57, 1.41` for `b = 2, 3, 10`, constant across the checkpoints. That is the Champernowne-type rate `1/log T`, consistent with normality in base `b` (**OPEN**). The exact-rectangle checkpoints make such frequencies exactly computable.
* **Column order.** It is not normal: the 2-block bias is 0.76, 1.6 and 9.7 and does not decay. Consecutive `n` share long runs of equal digits.
* **Natural-log reading.** Its column-1 law **oscillates** between rounds (for `b = 2`: 0.336 at round 19, 0.391 at round 20). This is the classical failure of `{ln n}` to have a distribution. It is another reason to prefer `log_b`.

**Best reconstruction of the purpose.** An explicit Champernowne-type ("log-table") number, built so that every prefix `T_r` is a complete rectangle. It is either a test of normality from a non-u.d. source, or a computable real that encodes every mantissa log at explicit positions.

**Relation to the digit and transversality themes (L5).**
* The construction lives wholly on the *top-digit*, archimedean side: mantissas and Benford.
* **The 2-versus-3 slice.** In base 3 the row `n = 2^k` is `{k log_3 2}`, which determines the leading ternary digits of `2^k` (FINITE-EXACT: top digit `= floor(3^frac(k log_3 2))` for `k <= 400`). Erdős's ternary problem asks about *all* ternary digits of `2^k`.
* Column 1 along `n = 2^k` is a 3-arc coding of the rotation by `log_3 2`: complexity exactly `3m` for `m <= 55`, and repetition estimate 5.69 on 6000 letters, driven by the partial quotient 23 of `log_2 3`. The 2-arc coding of the same rotation is the critical `3x+1` Sturmian word of Theorem S.
* The session's endgames read *low* p-adic digits (synthesis §4, item 2). So the construction is on the other side of the top/low divide. Label: ANALOGY; the rotation-coding identity is REAL.

## 6. Hypothesis candidates (for `INDEX.md` numbering by the integrator; no HYP files created)

* **PS1 (Proposition T′; PROVED modulo the Subspace Theorem, not audited).** For strictly supercritical, bounded-discrepancy, non-EP `w` with `Dio(w) > 1`, `Phi_2(w)` and `Phi_R(w)` are not both rational. Formal route: a rational `{inf, 2, 3}` Subspace cover of the Erdős 1062(ii) type, whose exact statement is UNVERIFIED. The algebraic version (HC4) stays a sketch: the exceptional-subspace case only yields "`x`, `y` are two embeddings of one algebraic number".
* **PS2 (Corollary M; PROVED modulo BDGP).** `Phi_2(Y)` and `Phi_2(Y_pr)` are not both algebraic. **Conjecture PS2′ (OPEN):** `Phi_2(Y)` is transcendental. This is a p-adic Nesterenko statement at `q = rho^2`.
* **PS3 (no-go; PROVED, elementary).** Any approximation scheme with natural heights gives a three-place Subspace product `H^(1 - 1/eta + o(1))`, with `1 - 1/eta = 0.299` for `Y3`. HC1 requires anomalously small heights, as in Theorem Y.
* **PS4 (OPEN).** The base-`b` log-table number (`log_b`, row order) is normal in base `b`. Evidence: bias about `C_b/r` with constant `r · dev`.
* **PS5 (OPEN).** `5/sqrt3` is optimal for 6 and for 7 unit hexagons, and more generally `(3k-1)/sqrt3` is optimal for `3k(k-1)+1` hexagons. Evidence: an 11-year-old record, untouched by 2026 global optimization; the area bound `sqrt7` is far below.
* **PS6 (FINITE-EXACT; a closed question).** The `Out(S_6)` twist fixes exactly the Rédei–Berge functions with `c3 = 2 d33`, and maps every other `U_T` outside the family.

## 7. Reproduction

```bash
cd <worktree>
python3 04-computation/experiments/procgen_sources_20260923_run.py   # writes 05-knowledge/results/procgen_sources_20260923.out
# about 25 s in total; the scripts run one at a time; peak RSS about 540 MB (the base-2, 20-round log table)
```

**Sections of the `.out`:**
* M1–M6: moonshine.
* H1–H6: hexagons.
* L1–L5: the `{log n}` construction.
* T1–T7: interplay tests.

**Methods:**
* Group orders come from explicit closures, cross-checked with sympy's Schreier–Sims for `M_24`.
* Codes are enumerated completely.
* 2-adic identities are integer congruences mod `2^N`.
* The hexagon checks use exact `a + b sqrt3` arithmetic with rational `a` and `b`.
* The ILP uses `scipy.optimize.milp` with pairwise fork constraints.
* The Radchenko–Wheeler check is at 40 digits (mpmath).
* Downloaded sources are in `scratch/procgen_sources/` (not for commit).

## Sources

`[P]` means the primary text was read in this lane, `[R]` that a named secondary source was used, and UNVERIFIED that neither was possible.

* **Primary items.**
  * D. Radchenko, C. Wheeler, *Real quadratic fields and finite quantum dilogarithms I*, arXiv 2609.21892v1. [P: abstract; §1 Theorems 1–2; §2.1; §4.1–4.4 (Definition 4, Theorem 5, Proposition 10, Theorem 7, Remark 6); §5 (Theorems 8–9); §6; App. A.1.]
  * S. R. Karingula, S. Lovett, *An elementary proof of the Komlós conjecture*, arXiv 2609.20979v1. [P: abstract; §1 (Theorem 1.2, Lemmas 1.4–1.5, Conjecture 1.6); §2.1.]
  * E. Ezeunala, H. Jiang, *Rank-one matrix discrepancy and algorithmic Kadison–Singer*, arXiv 2609.17266v1. [P: abstract; §1 (Theorem 1.1, overview); §2.]
  * L. Kruer, J. Kohlmeyer, *Erdős Problem 859: a refutation of the logarithmic density law* and *Erdős Problem 1062(ii): irrationality of a divisibility density* (conjectures.io expositions, 21 and 22 Sep 2026). [P: full texts. The Lean sources (36,197 lines, and at least 74,198 lines) were not read; the statement of `rational_three_place_subspace_cover` is UNVERIFIED.] Also the conjectures.io results page (accessed 2026-09-23). [P]
  * `long-mathematics/deterministic-von-neumann-fair-extractor`, README and repository metadata. [P, cross-reference only]
* **Transcendence.**
  * M. Waldschmidt, *Sur la nature arithmétique des valeurs de fonctions modulaires*, Sém. Bourbaki exp. 824, Astérisque 245 (1997) 105–140. [P: introduction, Théorèmes 2–4, bibliography [1].]
  * K. Barré-Sirieix, G. Diaz, F. Gramain, G. Philibert, Invent. Math. 124 (1996) 1–9, DOI 10.1007/s002220050044. [R: statement via Waldschmidt; metadata via Springer]
  * Nesterenko 1996. [R, via Waldschmidt]
  * The Subspace Theorem in Schlickewei's p-adic form (1977) and Schmidt's LNM 1467. [UNVERIFIED, not read; standard]
  * K. Mahler, *On the fractional parts of the powers of a rational number II*, Mathematika 4 (1957), using Ridout's theorem. [UNVERIFIED, from memory]
* **Moonshine.**
  * G. Mason, *M24 and certain automorphic forms*, Contemp. Math. 45 (1985) 223–244; D. Dummit, H. Kisilevsky, J. McKay, *Multiplicative products of eta-functions*, Contemp. Math. 45 (1985). [R: search-result snippets; the DKM page numbers are UNVERIFIED]
  * Wikipedia, *Monstrous moonshine*: McKay 1978, FLM 1988, Borcherds 1992, the decompositions. [R]
  * Conway–Sloane SPLAG (ch. 4 §11, ch. 10–11), Conway 1985, Griess 1982, FLM 1988, Conway–Norton 1979, Thompson 1979. [UNVERIFIED in this lane; every number used from them is recomputed in M1–M6]
* **Packing.**
  * E. Friedman, *Hexagons in Hexagons*, Erich's Packing Center, erich-friedman.github.io/packing/hexinhex (accessed 2026-09-23), with the picture `7.gif`. [P]
  * T. Berthold, D. Kamp, G. Mexi, S. Pokutta, I. Polik, arXiv 2605.04850. [P: §1, §6.2, Table 2]
* **Repository.**
  * [hard class](collatz_procgen_20260922_hard_class.md) (Theorems D and Y, §1.5, Proposition M, HC1–HC6); [synthesis](collatz_procgen_20260922_synthesis.md) (§2c, §4).
  * HYP-9123, HYP-9127; THM-481; THM-3027; THM-467 and THM-002.
  * [level11_short](level11_short_20260922.md); HYP-3706, HYP-3726; HYP-2456; THM-4035 and its sixty-clocks reflection; the [p-adic zeta audit](../reference/p-adic-zeta-irrationality-source-audit-20260825.md).
