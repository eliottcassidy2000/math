# Independent adversarial audit of opus-2026-10-14-S35

**Scope and method.** An independent subagent audited the session's files (2026-10-14).
- **Files audited:**
  - `paley_least_resistance_runs_one_modular_form_20261014.md`
  - `one_twenty_eight_atlas_balance_gold_20261014.md`
  - HYP-9275, HYP-9277, HYP-9278
  - the cited scripts and outputs
- **Not audited:** the session filled the `[N8]` sweep result after the audit began.
- **Method.** All code in this folder was written independently; none of the session's scripts were imported. The subagent could not save report files, so the session wrote this report from its final message, keeping the content.

Every recomputable number checked out. The real problems are credit and typing:
- The wheel formula and the `T_6` law are classical: Greene's hypergeometric definition, and Koike's conjecture as proved by Frechette–Ono–Papanikolas.
- Planar duality is classical: Kramers–Wannier / Wu–Wang.
- The Catalan mechanism repeats one that THM-438 already records as false, and THM-438 is not cited.
- Several statements about `b(q)` beyond `q = 700` rested on a perfect-square test until this audit closed the gap.

## What was recomputed independently

**(a) `A_10`, `A_12`** (`audit_a_saw_dp.py`). A signed subset DP over self-avoiding walks, checked against brute force at `q = 7, 11, 13`.
- `A_12(19) = 3,706,601,472` and `A_10(19) = 226,059,264`, matching the closed forms; `q = 23` matches as well.
- Both closed forms also hold at `q = 27 = 3³`, a prime power outside the 14 interpolation primes, reading `a_27` as the Frobenius trace 260. This is the first genuinely out-of-sample test.
- `A_10(11) = 478,720` is confirmed by brute force.

**(b) `K_6 − P_6` and the newform** (`audit_b_k6p6_sums.py`, `audit_b_newform.py`, `audit_b_twoadic.py`).
- **Sums.** Brute force and an independent exact `O(q²)` evaluator, built from Legendre-curve traces, were used.
  - The square law and `W − Y = 7 + 14χ − 2χq` hold at every prime `3..47`, and HYP-9276's `W` values reproduce.
  - The evaluator matches the session's table at all 124 primes, and the square law holds at 25 random larger primes.
  - Exact values: `M_8191 = −2^26`, `M_131071 = −63·2^28`, `M_213247 = −213248²`.
- **Newform.** `b` was rebuilt from theta products, which span the full 12-dimensional `M_3(Γ_0(32), χ₋₄)`.
  - The `4i`-eigenvector matches the session's `b` and audit C's closed form well past the Sturm bound.
  - It was tabulated to `n = 524,288` (`audit_b_newform_bq.json`), and `|b(q)|` equals the session's `c_q` at all 9,005 primes.
- **2-adic statements.** All of the following hold to 524,288:
  - the depth table;
  - the minima 8/10/7/11;
  - "`v_2 = 7` exactly at every `q ≡ 11 (mod 16)`";
  - the law `gcd(q+1, 2^7) | b(q)`, with exponent 7 sharp.
- **Zeros.** The zeros below 200,000 are 127 and 8191. A further zero is **`b(213247) = 0`, which is not a Mersenne prime**.

**(c) Wheels** (`audit_c_wheels.py`, `audit_c_greene_fop.py`).
- The proof of `T_n = Σ_ψ (χ(−1)J(ψ, χ))^n` is correct. Numerical checks pass in both classes, and `U(W_n) = q·T_n` holds.
- The independent level-8 eigenforms of weights 4, 6 and 8 span the full spaces. The weight-6 newform `b6` and both weight-8 newforms agree with Frechette–Ono–Papanikolas' published `q`-expansions.
- `T_4` and `T_6` hold at all 39 primes. The general form `T_6 = (q−1)(−b6_q − q·a_q − χ(−1)q²)` holds in both classes, at 44 primes.
- `T_8` is not in the span of the level-2/4/8 forms: tested at 39 primes against all 15 basis functions.

**(d) Planar duality** (`audit_d_duality.py`).
- The proof, orientation and sign are correct.
- Plane embeddings were computed independently, by brute force over rotation systems. Duality holds for every graph tested, including the sign when `|E|/2` is odd and graphs with several embeddings.
- `K_6 − P_6` is self-dual.
- `U(Q_3) = q·U(K_(2,2,2))` at 35 primes from 5 to 157.

**(e) Tournaments** (`audit_e_tournaments.py`, `audit_e_closed_family.py`, `audit_e_strong8_probe.py`). Brute force over permutations for `n <= 6`, and a full DFS path enumeration at `n = 7`.
- These reproduce, for `n <= 7`:
  - the minimum `1/3` and its equality characterisation;
  - the counts 2/16/120/1040/10080;
  - the least strong balance 2/5 then 4/11;
  - the Claim A minima;
  - the strong `δ` spectra.
- The `n = 8` Claim A counterexample checks over all `8!` orders: `H = 45`, ratio `13/45`, `δ = 4/9`, and the tournament is strong. The `n = 9, 10` witnesses also check.
- The closed-family results reproduce: minimum `1/4`, 176 gold failures at `n = 6`, and the semicomplete `n = 4` case. The 4-path witness is exactly `HP(T) ∩ Ext(Q)`.
- `H(P_19) = 1,172,695,746,915` is confirmed.
- **New:** every strong tournament with `n <= 7` has a balanced arc, and so do the `n = 8, 9, 10` counterexamples.

## Errors

- **E1 (main note §3, HYP-9275).** The Catalan explanation ("closed walks folded onto trees … the Catalan structures") is refuted in THM-438's addenda (MISTAKE-060/061): trees give 1, 3, 13, 69, … and the Catalan numbers appear only after signed cancellation. THM-438 already states the Catalan law and is not cited.
  - **Fix:** cite THM-438, use its mechanism, and record that `A_10` and `A_12` add the cases 42 and 132.
- **E2 (atlas §0.3, §3b).** "The three perfect nontrivial codes" is false: all Hamming codes are perfect, plus the two Golay codes.
  - **Fix:** say "the binary `[7,4,3]` Hamming code and the two Golay codes are QR codes".
- **E3 (atlas §3d).** "The limma `2⁸/3⁵` gives the 1/13 cycle" is false. The 1/13 cycle has one odd step and four halvings, so `1/13 = 1/(2⁴ − 3)`. That `2⁸ − 3⁵` also equals 13 is a coincidence.
- **E4 (atlas title, HYP-9277 title).** Both call an open conjecture a "theorem", the MISTAKE-611 pattern.
  - **Fix:** say "conjecture (verified `n <= 7`)".
- **E5 (main note contents list).**
  - "Rigidity circuits are exactly the graphs whose double covers are Calabi–Yau" is false. The Calabi–Yau condition is only the count `|E| = 2|V| − 2`, and the note's own table lists "`K_4` plus a degree-2 vertex", which meets the count but is not a circuit.
  - "The `K_6 − P_6` terminal is the symmetric square …" states the open HYP-9276 as fact.
- **E6 (main §4).** "The congruence forces `|b(127)| ≡ ±128 (mod 2^7)`" does not follow: for `q ≡ 7 (mod 8)` the depth table gives only mod 64.
  - **Fix:** derive `{0, 128}` from the valuation law plus `|b| <= 254`.
- **E7 (atlas §0.6, §6).** "All known Borsuk counterexamples are two-distance sets from strongly regular graphs" is false for the Kahn–Kalai line (dimensions 1325 down to 298).
  - **Fix:** restrict the statement to the lowest-dimensional ones (65, 64, and the unrefereed 63).
- **E8 (main §4b, HYP-9278: `T_8`).** The claim is true, but the cited fit output is underdetermined (13 primes, 14 functions). It prints a meaningless "exact fit" and lacks the `a_3 = −84` newform, and no session output shows the full weight-8 space.
  - **Fix:** cite an overdetermined test, like this audit's.
- **E9 (atlas §3d, minor).** The doublings list omits `k = 3`. The doubling positions are a complementary Beatty sequence, not convergent denominators of `log₂ 3`.
- **E10 (minor).** There are seven Mersenne primes up to `2^19 − 1`, not six: `q = 3` has `|b(3)| = q + 1`.
- **E11 (trivial).** The Sato–Tate mean 0.231 is over 123 primes; over the stated 124 it is `0.2326 ± 0.021`.

## Imprecisions

- **I1 (major).** The wheel formula, `T_n = (q−1)q^{n−1}·ₙF_{n−1}(1)`, is Greene's definition, not new. `T_6` is a theorem (Koike / Frechette–Ono–Papanikolas, IMRN 2004, Cor. 1.2), so it should be typed CITED rather than FINITE-EXACT. This also answers HYP-9278's open item on why the wheel ladder stays at level 8.
- **I2 (major).** Planar Fourier duality is classical (Kramers–Wannier / Wu–Wang 1976).
  - **Fix:** "PROVED (specialisation of a classical duality)".
- **I3 (major).** Between 700 and 200,000 the session checked only that `M_q + (q+1)²` is a perfect square. Its claims about `b(q)` there (the Mersenne values, the zeros, the 2-adic laws) were conditional until this audit's newform comparison.
  - **Fix:** state what was computed and cite this audit.
- **I4.** HYP-9276 is described as covering primes to 30,000. As audited, it covers every prime to 100,000, 1,000 random primes to 299,983, and 22 prime powers. "Beyond HYP-9276's range" is wrong for `[30,000, 100,000]`.
- **I5.** The FFT values are rounded floating-point results and should be typed VERIFIED, not FINITE-EXACT. Exactness now comes from this audit's exact evaluator and the newform match.
- **I6.**
  - The proposed 2-adic mechanism cites McCarthy–Papanikolas Theorem 1.8, which covers `p ≡ 1 (mod 4)`, while the congruence concerns `q ≡ 3 (mod 4)`.
  - The `q ≡ 1 (mod 4)` remark is vacuous; the right factorisation there uses `q − 1`.
  - HYP-9276 (e) now says "partially explained", not "unexplained".
- **I7 (atlas §0.2, §3a′).** "2-adic continuity of Frobenius traces", "Mersenne primes are where the factor is smallest", and "127 enters exactly through `127 ≡ −1 (mod 2^7)`" are ANALOGY or NUMEROLOGY presented as structure. The non-Mersenne zero at 213247 counts against them.
- **I8 (atlas §0.4, §4).** "Gold partition is a convexity phenomenon of order polytopes" and "pure combinatorics cannot prove it" overclaim what the closed-family counterexamples show.
- **I9 (HYP-9278 status).** Only the canonical-class formula is proved. "Calabi–Yau" for these singular double covers needs a resolution.
- **I10.** The cube is not "the first" pure core in the Fano regime: all five connected cubic graphs on 8 vertices qualify.
- **I11.** The Claim B proof says "transitive" where the argument actually rules out directed 3-cycles. The claim itself is true.
- **I12 (atlas §3c).** THM-868's shell index is `x = Σ(2s_i − 7)²`, not `Σ(s_i − 7/2)²`.
- **I13.** The "octahedron minus two edges" row was tested by the session only for `q ≡ 3 (mod 4)`. It holds in both classes (20 primes).
- **I14 (additions the session can use).**
  - the `q = 27` out-of-sample check of `A_10` and `A_12`;
  - balanced arcs exist exhaustively for `n <= 7`;
  - the zero `b(213247) = 0`;
  - the DP computes `A_14` at `q = 19, 23, 27` cheaply.

## Literature checks

- **Confirmed:**
  - Akiyama–Frougny–Sakarovitch Corollary 4, verbatim from the authors' preprint, including the Dubickas comparison;
  - Flatto–Lagarias–Pollington Theorem 1.4;
  - Brightwell–Felsner–Trotter `(5 − √5)/10`;
  - Peczarski's weak `t_0 >= t_1 + t_2`, matching the session's code;
  - OEIS A000930;
  - OEIS A076227, the same sequence and indexing as the Collatz sieve;
  - A038375, known to `n = 11`;
  - Bondarenko 65 and Jenrich–Brouwer 64;
  - `e8 = so(16) ⊕ 128`.
- **Not found:**
  - no prior tournament analogue of 1/3–2/3;
  - no modular evaluation of `₈F₇(1)`, so `T_8` stays open.
- **Credit:**
  - Credit to mac-mini (THM-4653/HYP-9276) and tiling-modular-atoms is present, apart from the out-of-date ranges. The session's `K_6 − P_6` orientation is identical to tiling's.
  - Missing credit: THM-438, Greene, Frechette–Ono–Papanikolas, Wu–Wang.

## Verdicts

| Claim | Verdict |
|---|---|
| `A_2..A_12` closed forms | CONFIRMED, including out of sample at `q = 27` |
| Octahedron identity | CONFIRMED (FINITE-EXACT) |
| Catalan coefficients | Numbers CONFIRMED; mechanism REFUTED (E1) |
| `K_6 − P_6` square law | CONFIRMED WITH FIXES (I3–I5) |
| `W − Y` normalisation link | CONFIRMED |
| Newform `b` | CONFIRMED |
| 2-adic depth table and valuation law | CONFIRMED to 524,288 (fix E6, I7) |
| Mersenne values, `M_8191 = −2^26` | CONFIRMED |
| Wheel Jacobi-sum formula | CONFIRMED, but standard |
| `T_4`, `T_6` laws | CONFIRMED; `T_6` is a published theorem |
| `T_8` not in the level-8 span | CONFIRMED with new evidence; the session's evidence was inadequate (E8) |
| Planar duality and `U(Q_3) = q·U(K_(2,2,2))` | CONFIRMED; the duality is classical |
| Rigidity/Calabi–Yau dictionary | OPEN; typing fixes E5, I9, I10 |
| Tournament results `n <= 7`, Claim A, closed families, Paley `P_19` | CONFIRMED |
| Proved lemmas | CONFIRMED (I11 wording) |
| Collatz sieve | CONFIRMED |
| `n = 8` sweep; unrefereed Gupta and hexagonmath inputs; Mahler and Borsuk cross-references | NOT CHECKED |

## Sources consulted by the audit

- Akiyama–Frougny–Sakarovitch preprint: https://www.math.tsukuba.ac.jp/~akiyama/papers/3half_H65fullTH-060406.pdf
- Frechette–Ono–Papanikolas, IMRN 2004: https://uva.theopenscholar.com/files/ken-ono/files/084_8.pdf
- Frechette–Ono–Papanikolas, combinatorics of traces: https://uva.theopenscholar.com/files/ken-ono/files/088_8.pdf
- OEIS A076227: https://oeis.org/A076227
- OEIS A038375: https://oeis.org/A038375
- Flatto–Lagarias–Pollington (IMPAN listing): https://impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/70/2/108521/on-the-range-of-fractional-parts-p-q
- Brightwell–Felsner–Trotter: https://trotter.math.gatech.edu/papers/97.pdf
- Abelian spin models survey (duality): https://arxiv.org/pdf/1112.4404
- Planar duality for spin models (J. Algebraic Combin.): https://emis.univie.ac.at/journals/JACO/Volume4_2/h3w2062p556h3521.html
- Dawsey–McCarthy survey: https://arxiv.org/pdf/2101.06303
