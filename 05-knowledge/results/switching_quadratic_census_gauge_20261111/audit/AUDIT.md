# Independent audit: switching_quadratic_census_gauge_20261111

Auditor: independent agent, 2026-11-11.

Scope:
- README.md (all of it), PREREGISTRATION.md, THM-5139, THM-5140, HYP-9669, HYP-9670 and out/;
- the cited prior art.

Method: own programs only, all in this directory; nothing outside `audit/` was modified.

## Verdicts

| item | statement | verdict |
|---|---|---|
| A | Theorem G (Lemmas 1, 2, induction, row sums, symmetry, (*), (**), base cases) | **PASS** |
| B | Corollary 1, Corollary 2 (incl. triple sums of H = 0 mod 4) | **PASS** (one citation to add, m5) |
| C | Corollaries 3, 4, 5, 6 | **PASS** |
| D | Section 2.3 (reduction, digit lemma, equivalent form); HYP-9669 | **PASS WITH CORRECTIONS** (m1-m4); no counterexample found |
| E | Proposition B1 | **PASS** (novelty wording, m6) |
| F | Section 3.3 enclosures | **PASS** (digits reproduced by an independent code path) |
| G | Propositions B2, B3 | **PASS WITH CORRECTIONS** (credits, m7, m8) |
| H | Section 3.5 table, s(n)/a(n), section 4 arithmetic | **PASS WITH CORRECTIONS** (m9) |
| I | Section 5 (Yang-Mills reading) | **PASS WITH CORRECTIONS** (m10, m11) |

**No MAJOR findings.** No statement typed PROVED contains a mathematical error. The minor findings below concern:
- the precision of side claims;
- missing credits;
- the wording of novelty claims.

## Proof checks by hand (A-C)

**Lemma 1.**
- `s_T(u) = sum_{v in N+(u)} s_{T-u}(v)`.
- `t_T(u) = H(T-u) - sum_{v in N+(u)} t_{T-u}(v)`. This uses `sum_v t_R(v) = H(R)`, which is valid at every order, including order 1.
- `H(T-u)` is odd (Redei).
- Correct for `n >= 2`.

**Lemma 2.** Switching at `u` swaps `N^+` and `N^-` at `u`. So `s_T(u) + s_{T^u}(u) = H(T-u)`, and likewise for `t`. Correct.

**Induction.**
- If `u` is not in `X`, then `N^+_{T^X}(u) = N^+_T(u) + X` and `T^X - u = R^X`.
- The expansion gives (*). It uses the symmetry and the zero diagonal of `Lambda_R`, which are part of the induction hypothesis at order `n-1`.
- Row sums: `T^{V-u} = T^{u}` together with Lemma 2. Correct.
- Extension to sets `X` that contain `u`: pass to `V - X`. Correct.
- Class invariance holds because the affine law determines the matrix uniquely.

**Symmetry.**
- Lemma 1 in `T-u` gives `e_{T-u}(w) = 1 + <b, e_Q>`.
- Applying (*) to `T-u` at row `w` needs order `n-1 >= 3`.
- Summing over `a` gives (**).
- (**) is symmetric under `(u,a) <-> (w,b)` by the symmetry of `Lambda_Q` at order `n-2`.
- The induction hypotheses are used in a consistent order.

**Base cases.**
- `n <= 2`: `e` is `(0)` or `(1,1)`.
- `n = 3`: `Lambda = J - I`.

**Corollary 1.** Over the Hamiltonian paths `P` of `T-a`:
`hc(T) + hc(T^a) = sum_P [alpha(1-beta) + (1-alpha)beta] = sum_P (alpha+beta) (mod 2) = 1 + e_T(a)`.
The induction on `|X|` is correct.

**Corollary 2.**
- The cycle-length parities follow from `(T^X)[W] = T[W]^(X cap W)`.
- `H = 1 + 2 c_odd (mod 4)`.
- A quadratic function has vanishing third differences, so `sum H(T^S) = 8 H(T) = 0 (mod 4)`.
- Correct.

**Corollary 3.**
- There are `2^(n-3)` odd sets `Q` inside `R`. This number is even, so the constant 1 cancels.
- The linear part sums to `E_R(v)`.
- The quadratic part sums to `M_R(x,y)`; the case `x = y` contributes nothing because `M_R` has zero diagonal.
- Base values: `M = 0` at `n = 4`, and `M = J - I + Lambda` at `n = 5`.

**Corollary 4.**
- `E_R(v)` counts the transitive triples in which `v` is the source or the sink. That number is `C(d,2) + C(3-d,2)`, which is odd.
- Hence `M(a,b) = 1 + d_a + d_b`. Correct.

**Corollary 5.**
- The logic is sound: a function of `F_2`-degree at most 2 that takes one value at `0`, at every `e_i` and at every `e_i + e_j` is constant.
- Coverage was checked in `../character_of_two_records_20261111/cycle_parities.py`, mode `pairs`:
  - the loop `for u in range(N): for v in range(u, N)` visits all singles and all pairs (300 at `N = 24`, 210 at `N = 20`);
  - the normalised member is computed separately;
  - `switch()` is correct;
  - every odd length is recorded.
- Each output file shows one parity vector per class.
- `100000011 = C(9,k) mod 2`, so `H = 3 (mod 4)`.

**Corollary 6.**
- This is Dickson's character sum on `F_2^V/<1>`.
- The all-ones vector lies in the radical because the row sums are zero.
- The histograms were reproduced (next section).

**Remark (1).** A class and its converse class are distinct for `n >= 3` and have the same `Lambda`, so every fibre is even. Correct.

## Independent computation (own encodings and DPs)

### `audit_G.c`

Outputs: `audit_G_n3_6.out`, `audit_G_n7.out`.

**Orders 3-6, every labelled tournament.** The program uses a start-indexed subset DP. Every check below gave 0 violations:
- Lemmas 1 and 2;
- the affine law for **all** `X`;
- symmetry, zero diagonal and zero row sums of `Lambda`;
- class invariance for all `X`;
- (*) and (**), with the `Lambda` of each subtournament computed by switching inside `T[W]`;
- Corollary 1 for all `X`;
- the `c_l` parity law with the stated polar form, for every `l` and every `X`;
- the `H mod 4` formula with `lambda_a` and `M`, for every `X`;
- triple sums of `H`;
- Corollaries 3 and 4.

At order 6 this is, for example, 2,097,152 checks of the affine law and 8,388,608 checks of the `c_l` law.

**Order 7.**
- Lemma 2, the affine law for all `X`, the symmetry/diagonal/row-sum checks, class invariance, Corollary 1 and the triple sums: 0 violations.
- Class histograms reproduced exactly:

| order | members with `H = 1 (mod 4)` per class (number of classes) | imbalance (`#H=1 - #H=3`) |
|---|---|---|
| 4 | 0 (2), 8 (6) | constant on each class |
| 5 | 10 (24), 12 (40) | +4, +8 |
| 6 | 20 (all 1,024) | +8 |
| 7 | 24 (1680), 28 (720), 32 (5040), 36 (19280), 40 (6048) | -16, -8, 0, +8, +16 |

- Distinct `Lambda` per order, orders 3-7: 1, 4, 22, 91, 2972. At order 6, `Lambda = 0` for 184 classes.

### HYP-9669

**Exhaustive at orders 3-7.** All tournaments and all non-empty `U` contained in `V`: 0 violations.

**`audit_2adic.c`** (output `audit_2adic.out`). All `2^n` switchings and all `U` contained in `V`, on:
- `P_7` and `(source) + P_7`;
- transitive and near-transitive tournaments at orders 8-12 (one arc reversed; transitive only at orders 8-11);
- circulants of orders 9 and 11, including `P_11`;
- `(source) + P_11` at order 12;
- 5, 20, 8 and 1 random tournaments at orders 8, 9, 10 and 11.

Result: no violation, and digits 1, 2, 3 of `H` have degree at most 2, 4, 8.

### `audit_census.py`

Output: `audit_census.out`. All results match the packet.
- A000568(1..12).
- Level bit-loss bound: 0 violations over 8,656 cycle types, `3 <= n <= 40`.
- `e_d(1..n-1) <= C(n,2)^d/d!`: 0 violations.
- `T_K <= B_K` at (12,0), (20,1), (30,2), and in all 437 cases with `3 <= n <= 40`.
- Log-convexity of `t_k`: holds.
- `B_0(n) < 1` for `n >= 8`; `B_0(7) = 1.0175`.
- `N_odd = A049313` for odd `n <= 13`, and the clock-shift identity holds exactly.
- Leading digits of `a(n)` were recomputed with mpmath `loggamma` and agree in all 31 digits and in the exponent for every `n` in 3.3, up to `10^1000`.
- B2 moments: confirmed.
- `s(n)` ratios: confirmed.

### `census_forecast.py cmd_giant` (F)

The ball provably contains `log10 a(n)`:
- monotonicity together with the arb `union`;
- `unique_fmpz` as the ambiguity guard;
- a conservative tail;
- a safe fallback for `arb.max`;
- the log-convexity count;
- arb rounding.

No defect.

## MINOR findings

**m1. "Attained at every `|U| <= 6`" is stated as a general fact but holds only for the R1 runs.**
- Location: README 2.3, last sentence; HYP-9669 status; PREREGISTRATION outcomes.
- Problem: the sentence reads "for odd `n` the bound is attained at every `|U| <= 6`", which is false in general.
  - `P_11`: the least `v_2` at `|U| = 5, 6` is 4.
  - Circulants of order 9 with connection sets {1,2,3,4} and {1,2,4,6}: the least `v_2` at `|U| = 5` is 4.
- Fix: restrict the claim to "the random tournaments of R1".

**m2. The orders 6-8 checks are not evidence for the `H` part of HYP-9669.**
- Location: README 2.3 evidence; PREREGISTRATION.
- Problem: for `n <= 8`, `2 floor(n/3) <= 4`, so `c_U(H) = 0` for `|U| > 4` by THM-4524. Beyond the proved range the `H` part is vacuous at these orders. The first order with content is 9.
- Fix: say so.

**m3. R1 tests a restricted family of sets `U` without saying why that suffices.**
- Location: R1 method.
- Problem: R1 tests only the sets `U` contained in `V` minus one vertex, and only at one base member.
- This is without loss of generality, but the note gives no argument:
  - change of base: the coefficients at another base are signed sums of the `c_W` with `W ⊇ U`;
  - removing the fixed vertex: use `e_{v0} = sum_{w != v0} e_w` modulo `<1>`, together with `Δ_w Δ_w = -2 Δ_w`.
- Fix: add one sentence giving this argument.

**m4. The digit lemma and the equivalent form need precise hypotheses.**
- Location: README 2.3.
- (a) The digit lemma's sketch assumes non-negative coefficients. Fix: first reduce each coefficient mod `2^(k+1)`.
- (b) The "equivalent form" presupposes `v_2(c_U) >= k` for `|U| = 2k`. Fix: state that hypothesis explicitly.
- The phrase "for odd `n` the two are equivalent" should read "at all odd orders `<= n`".

**m5. Missing credits for Corollary 2 and a possible prior source for Corollary 1.**
- Location: THM-5139 credits; README Corollary 2.
- Credit Grinberg-Stanley, arXiv:2307.05569, Thm 7.1 (`H = 1 + 2 c_odd mod 4`); THM-478 already cites it.
- Forcade, "Parity of paths and circuits in tournaments", Discrete Math. 6 (1973), could not be read. It must be checked before Corollary 1 is typed NEW.
- Repository search found no prior statement of `e`, of `Lambda`, or of a quadratic law on the switching cube.
  - Patterns searched: "endpoint", "start or end", "starting or ending", `s(v)+t(v)`.
  - Files read: THM-030, THM-466, THM-478, THM-4524, THM-4876, THM-5127, HYP-9667, THM-1405, THM-1460, THM-1465, THM-1470.
- THM-478 (the tiling/arc cube) and THM-5127 (C) (the vertex cube) concern different cubes.

**m6. The note overstates what Session271 left open.**
- Location: README section 1, "sections 3.2-3.3 answer those two points"; THM-5140 NEW list.
- Problem: Session271 `forecast/PROOF.md` section 2 already proves a remainder bound by support, `U_M(n)`, valid for all `n` and `M`, and the sandwich `b(1+C_s) <= a <= b(1+C_s+U_{s+1})`. What it leaves open is a remainder that is uniform by level and outward-rounded digits.
- Fix: present B1 as the level-wise form of that bound.
- B1 itself is correct.

**m7. Parts of Proposition B3 are already in THM-479.**
- Location: README 3.6 (ii); THM-5140 B3.
- "`N_lev = 0` for odd `n`" and the values of `N_lev` up to `n = 16` are in THM-479 (3). Credit them.
- `N_lev <= 2^(n^2/4)` follows at once from THM-479's exact formula `N_lev = sum 2^(o2-k)/z`.
- The clock shift is a two-line consequence of THM-479.
- The first correction `(64/3)(n)_3 4^-n` was not found in the repository.
- The leading term `2^C(n-1,2)/n!` is probably folklore; check Babai-Cameron 2000.
- Fix: type only the correction term as NEW.

**m8. The NEW claim in Proposition B2 is too broad.**
- Location: README 3.4; THM-5140 B2 NEW list.
- Problem: the indeterminacy of `q^(-C(M,2))`-type moment sequences and their discrete theta solutions on `q^Z` are classical (Stieltjes; the log-normal examples).
- Fix: narrow NEW to the application to the A000568 cycle index, that is, the mixed moments of `H_p`.
- The proof was checked and is correct: completing the square, `o = C(M_1,2) + (1/2) sum_{d>=3} φ(d) M_d^2`, mutual singularity.

**m9. The order-18 density in the section 4 heuristic does not follow from the stated rule.**
- Location: README section 4, heuristic paragraph.
- Problem: "about `9e19` at order 18" does not follow. With `a(18) = 1.78e30` and `n!/2^(n-1) = 4.9e10`, the stated rule gives at most about `7.3e19`. Using half of max `H` gives less.
- Fix: write "about `10^19`-`10^20`".
- The rest of this section was re-derived and is correct:
  - the densities at orders 5-8;
  - the Poisson figures;
  - the arithmetic of HYP-9670;
  - the first moments and `I(mu,nu)` in the 3.5 table, by hand;
  - `1 - s/a ~ n 2^(2-n)`, numerically.

**m10. `Winter.convergence` is reported without its hypothesis.**
- Location: README 5.3, items 1(a) and 2(d).
- Problem: `Winter.convergence` (`YM/FullVolume/Glide/Vac.lean`, line 46) assumes `Tendsto (Z i).bonusB atTop (nhds b)`. The note states it without this hypothesis.
- Fix: add the hypothesis.

**m11. The growth exponent quoted for `AllCopy` is not in the cited file.**
- Location: README 5.3, item 2(a).
- Problem: in `AllCopy` (`YM/OS/Slabs.lean`, lines 109-113) the bound is `C^|w| (|w|/r)^(A|w|)`, with `A` existentially quantified. The exponent `4(D+1)n` was not found there.
- Fix: quote the bound as it stands in the file.

**Wording.** README 2.2 says "`n >= 4`" before (*), which holds for `n >= 3`. Attach the condition to (**) instead.

## Section 5 spot-checks

Read-only; nothing in the clone was built, run or imported. Pinned commit `3456e5e7` confirmed. Checked:

| item | where | what was found |
|---|---|---|
| `Winter`, `Winter.convergence` | Vac.lean, lines 31 and 46 | present; see m10 |
| `Winter.mk` | `verification/declaration-usage.tsv` | recorded as used |
| `reaches` | `YM/RG/History.lean`, line 239 | `β N p 0 = clock.bar` |
| `Transfer` | `YM/RG/History.lean`, line 223 | `rep m .eval = push (rep (m+1)).eval` |
| `path_shared` | | present |
| `FixedBounds` | `YM/FullVolume/Discrete.lean`, line 50 | present |
| `first_grade` | `YM/Trace2/Transfer.lean`, line 32 | for `d < 4`, every field of weight `d` is a scalar multiple of the unit |
| `AllCopy` | `YM/OS/Slabs.lean` | present; see m11 |
| `local_action` | `YM/EveryPhase/Charge.lean`, line 143 | present |
| "Not `Tendsto` along a single diagonal" | `YM/Extraction/Orders.lean`, line 8 | present |
| `forest` | `YM/Ac/Tree.lean`, line 89 | hypothesis `(M+M+1)^4 γ <= 1` |
| `pieces` | Sprig.lean | present |
| `first_surface`, `Completed`, `HighEnergy`, `smallField` | | present |
| `b0 = 11 C_A/(48 π^2)` | `YM/Gauge/Colors.lean`, line 184 | `slopeConstant = 11 C_A/(24 π^2) = 2 b0` |
| README disclaimers | README, line 26 | present |

**Analogies Y1-Y14.** They are typed honestly: Y9 KNOWN, Y13 WEAK, Y14 OPEN.

**SUGGESTION paragraphs 2(d) and 4(d).** Both carry adequate conditional caveats. Add m10.

## Not checked

- Forcade 1973 and Babai-Cameron 2000, full texts.
- Grinberg-Stanley 2023, beyond Thm 7.1.
- The order-20 and order-24 cycle counts themselves; only the coverage and the logic of Corollary 5 were checked.
- The order-12 forecasts (not scorable).
- The (S) survey claims beyond the names spot-checked above.
- `rank_stats` was not rerun.
- The R1 rows at orders 12-14 were not rerun; analogues at orders 9-12 were.
