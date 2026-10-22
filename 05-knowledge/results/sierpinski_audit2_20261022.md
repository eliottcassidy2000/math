# Independent audit 2: THM-4757 (loop expansion) and HYP-9365 (two-loop law), 2026-10-22

**Scope.** The two pre-audit results of session mac-mini-2026-10-22-sierp:
* `01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`;
* `05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md`;
* with their context (HYP-9364, THM-4755, §11 of `sierpinski_steiner_collatz_tournaments_20261022.md`).

**Method.** All code is new and independent of the session's scripts: `04-computation/experiments/sierpinski_ocf_20261022/audit2/` (Python with exact rationals, python-flint, numpy/mpmath; C for Hamiltonian counts; nauty `gentourng`/`labelg`). The session's scripts were read only to locate claims.

**Scope note.** The audited texts are THM-4757 and HYP-9365 as of about 03:00 (unchanged since), and §11 of the results file. Not audited: material the session added later, namely §12 of the results file (commits `6d93cbaf5c`, `5c505b1be5`), and HYP-9367 with its loop engine (`2bf99b9fff`, `9dd30913b2`, `f40f0ffabd`).

**Status words** as in the repo: PROVED, CITED, FINITE-EXACT, DERIVED, NUMEROLOGY, CONJECTURE. **Severity** (as requested): high = a false statement labelled PROVED, or a wrong number; medium = overstated evidence or missing credit; low = wording.

## Verdict

**No stated result is wrong.** THM-4757 (1)–(4) are correct; this audit confirms them with exact brute force. Every exact number in HYP-9365 that could be recomputed was reproduced.

**High, by the letter of the scale ("a wrong number"), but immaterial.** Three digit errors; none affects a conclusion:
* THM-4757 §4 quotes the brute-force values as −0.19910, −0.16270 and −0.0068140. The true values are −0.1991111, −0.1627059 and −0.0068139 (edit T6).
* HYP-9365 (c): `n·Θ → −0.23703` should be −0.237054 (edit H7).
* HYP-9365 (e): the ratio range 1.62–1.72 is 1.62–1.69 with the session's own fit (edit H10).

**Medium.**
* THM-4757's own "independent check" compared 4 printed digits, on vertex-transitive tournaments, at orders where the dumbbell cannot appear (edit T2).
* A side remark in THM-4757 is false: Paley character sums do not stop collapsing at three loops (edits T10, R4).
* HYP-9365 (b) "settles" an approach that is conditional on the conjecture, and its unspecified `c` is `1 + y_∞` (edit H6).
* HYP-9365 (d) claims more scaling evidence than two orders give (edit H9).
* Two credits are missing: the random-regular value `e` (Wormald's Lemma 3 with McKay's method; Alon's observation, recorded by Yuster), and Bhattacharjee's expansion (edits T3, H3, H6).

**Low.** The `e ≥ 0` display, "vertex-transitive" for "circulant", the `Σ x²` sentence, the choice of extrapolations, and `y` against `y_pred` notation.

| item | claim | verdict | severity |
|---|---|---|---|
| A | forest identity `κ = Σ_L N_L/(n)_L` | PROVED; verified by direct forest enumeration, `n ≤ 8` | — |
| B | cluster expansion, excess grading, no tree level for regular `T` | correct; verified exactly against `log Φ_T` | low (display `e ≥ 0` before regularity is assumed) |
| B | one loop `X_0 = ½ Σ log((1+z²x²)/(1−z²x²))` | PROVED; proof checked; exact through `z^16` | — |
| C | two-loop kernels, symmetry factors, periodicity, resummed `X_1` | correct; exact through `z^8`–`z^20` on 30+ matrices, dumbbell included | high-immaterial (quoted digits); medium (weak own check; false remark on Paley character sums) |
| D | HYP-9365 numbers, correlations, limits | all reproduced; limits rederived | high-immaterial (`n·Θ` digits); medium ((b) "settles") |
| E | doubly regular tournaments | data confirmed; mechanism correct, and a theorem for `X_0, X_1, X_2` | medium (`O(n^(−3))` "as required" overstated) |
| F | wobble vs `F(R_n,−1)/n!` | values and signs confirmed; ratio fit-dependent | high-immaterial (1.72); low (labels) |
| G | prior art | nothing found for the loop expansion, two-loop formula or doubly regular universality | medium (credit for the random-regular `e`) |

## A. The forest identity (THM-4757 (1))

**Proof.** Correct. `Π_i(1 + S_(σ_i σ_(i+1)))` is `2^(n−1)·[σ is a Hamiltonian path]`. A directed linear forest with `L` arcs has `n − L` components, so it lies in `(n − L)!` orderings.

**Check** (`a_forest_identity.py/.out`). Direct enumeration of all Lah structures, i.e. every directed linear forest of `K_n`, with exact integers. Tested on:
* all 74 labelled tournaments with `n ≤ 4`;
* 40 random tournaments for each `n = 5, 6, 7`, and 12 for `n = 8` (394,353 directed linear forests each);
* `R_5`, `R_7`, `P_7`, and the non-circulant regular 7-tournament (`H = 171`).

Results:
* `2^(n−1)H = Σ_F (n−|F|)! Π S_e`: exact in every case.
* The top term `N_(n−1) = (−1)^(n−1) F(T,−1)`: exact in every case.

**New observation (small lemma).** Reversing one odd path is a sign-reversing involution: choose, for instance, the odd path through the least vertex that lies on an odd path. So only forests whose paths all have even length contribute.
* In particular `N_L = 0` for every odd `L`.
* `N_(n−1) = 0` for even `n`.
* The same involution applied to a single run kills every cluster that contains an odd run, which is why the kernels vanish for odd runs.
* Verified for every tournament above.

**Link with THM-4756** (`a_forward_step_link.py/.out`, exact on random `n ≤ 7`). Since `Σ_σ Π_i(1 + tS_i) = Σ_σ (1+t)^fwd (1−t)^bwd`,

    Σ_L (n − L)! N_L t^L = (1 − t)^(n−1) F(T, (1 + t)/(1 − t)).

* At `t = 1` this is (1), so (1) is a reparametrisation of the forward-step polynomial.
* `N_odd = 0` is the palindromy `F(T, x) = x^(n−1) F(T, 1/x)`.
* Worth one line in THM-4757 (optional, not among the required edits).

**Classical one-sided form.** Inclusion–exclusion over the steps that use an arc of the converse gives `H(T) = Σ_F (−1)^|F| (n − |F|)!`, over the linear forests `F` of `T^op`. This is Berge's inclusion–exclusion, the route to Rédei's theorem recalled by Grinberg–Stanley (arXiv:2307.05569).
* (1) is its symmetrised form.
* THM-4757's "classical" attribution is fair.

## B. Clusters, excess, one loop (THM-4757 (2), (3))

**Machinery** (`b_partitions.py`). Brute force over (run-length tuple, set partition of run positions) with Möbius weights and generic `einsum` contraction, in exact rational arithmetic.

**Validation** (`b_validate.out`):
* The sum over all excess levels equals `[z^L] log Φ_T`, computed directly from the `N_L`. Tested on 8 random non-regular tournaments (`n = 5, 6, 7`) and on `R_7`, `P_7`, for `L ≤ 6`.
* Tuples with an odd run cancel, excess level by excess level.
* For regular `T`, pruning the degree-1 blocks changes nothing.
* **Tree terms.** Non-regular `T` have nonzero excess `−1` terms (for example `[z²]X_(−1) = −|Ŝ1|²`; −0.32 for a random 5-tournament). Regular `T` have none, as claimed.
* So the display `log Φ_T(z) = Σ_(e ≥ 0) X_e(z)`, which in §2 comes before regularity is assumed, should read `e ≥ −1` (edit T4/T5, low).

**One loop.** The proof is correct in every step:
* excess 0 with minimum degree 2 means a single cycle;
* interior positions are singletons and run ends are paired, so `μ = (−1)^p`;
* the single-cycle matchings number `(p−1)! 2^(p−1)`, with every orientation pattern of runs `2..p` equally often (enumerated for `p ≤ 7`, `b_cycle_matchings.out`);
* backward runs give `(−Ŝ)^ℓ`;
* `A = y/(1−y)`, `B = 2y²/(1−y²)`, and `X_0 = −Σ_k (A/B) log(1+B)`;
* the `±ix` pairing replaces `A` by its even part `B/2`, and the zero eigenvalue contributes 0.

Brute force agrees with `−Σ_(m odd) tr(Ŝ^(2m))/m` exactly:
* through `z^8`–`z^20` on every matrix of §C (orders as in its table);
* through `z^8` on `R_15` and `P_19`.

Remark: the convergence radius is better than stated. Since `2σ_max² ≤ Σ_k σ_k² = n(n−1)`, every tournament has `|x_k| < 1/√2`, so all series here converge for `|z| < √2`.

## C. Two loops (THM-4757 (4))

**The resummed formula is exact.** Two independent brute-force methods were compared with the closed formula, coefficient by coefficient, as exact rationals. The formula's `z`-series was built directly from `G, H, O`, not from the kernel table.

* **P (set partitions)**, `b_partitions.py`.
* **M (blow-up scaling)**, `b_mscaling.py`. For `M' = M ⊗ J_m/m` every excess-`e` cluster scales as `m^(−e)`, so `log Ψ_(M')(z) = Σ_e m^(−e) X_e(M)` exactly. `Ψ_(M')` is computed with no cluster theory, as `Λ[exp W]`: `W` is the walk generating polynomial in visit variables and `Λ: t^c ↦ Π (m)_(c_v)`.

| matrices | methods, orders | result |
|---|---|---|
| all 3 regular 7-tournaments | P to `z^10`, M to `z^12` | exact |
| all 15 regular 9-tournaments | P to `z^8`, M to `z^10` | exact |
| 5 non-vertex-transitive 11-tournaments, `R_11`, `P_11` | P and M to `z^8` | exact |
| regular 9-tournament, class 3 (non-vertex-transitive) | P to `z^12` (2.19M partitions) | exact |
| `R_15`, `P_19` (the session's cases) | P to `z^8` | exact |
| 2 random regular 15-tournaments (3-cycle-reversal chain; 15 and 11 distinct values of `diag S^6`) | P to `z^8` (one to `z^10`) | exact |
| non-circulant regular 7-tournament | M to `z^16` | exact, including the dumbbell |
| 4 generic Eulerian skew matrices (entries `0, ±1`, zero row sums, unequal degrees) | P to `z^8`, M to `z^10`/`z^12` | exact |
| one of them (`n = 7`) | M to `z^20` | exact; dumbbell nonzero at every order from `z^6` |

Sources: `c_two_loop_check.out` (`OVERALL: ALL EXACT`), `c_partition_L12.out`, `c_session_cases.out`, `c_dumbbell_n7_L16.out`, `c_eulerian_L20.out`.

For odd `L` every cluster sum vanishes, at every excess level, by the reversal involution of §A. This was checked with all run lengths at `L = 3, 5` in `b_validate.out`, and the odd walk polynomials of method M vanish identically. So `L = 5, 7` carry no information.

**The dumbbell.** For every regular `T`, the first term of `d = diag(H − G)` is the constant `2 diag(Ŝ²)`, and `G·1 = 0`. So the dumbbell has no `z^L` term with `L < 14`; its first term has loops `6, 6` and bridge 2, i.e. `(a, b; c) = (6, 6; 2)`.
* The session's check (`z⁴`, `z⁶`) therefore could not test it.
* Here it is tested in two ways. At `z^14` and `z^16` on the non-circulant 7-tournament, the exact dumbbell is `−221184/13841287201` at `z^14`, and the totals agree exactly. On Eulerian skew matrices it enters from `z^6`. Everything agrees exactly.
* The derivation uses only `Mᵀ = −M` and `M1 = 0`, and the formula indeed also holds on the four generic Eulerian skew matrices tested.

**The session's own check was weak** (medium). `two_loop_series_check.py` compares against brute-force values hard-coded to 4 significant digits. The quoted brute-force numbers are rounding artefacts:

| | z^4 | z^6 |
|---|---|---|
| `R_15`, quoted | −0.19910 | −0.060250 |
| `R_15`, exact | −0.1991111 | −0.0602495 |
| `P_19`, quoted | −0.16270 | −0.0068140 |
| `P_19`, exact | −0.1627059 | −0.0068139 |

The identity holds exactly (edit T2, T6).

**Kernel table** (`d_kernels.py`, `d_kernels_L12.out`, `d_kernels_L14.out`). The combinatorial coefficient of every (core type, sorted lengths) class was extracted from the brute-force enumeration of excess-1 partitions. Each partition was reduced to its core with its orientation sign; no kernel code was used.
* Compared with `K·#orderings/|Aut|` from the table, with `|Aut| = 8, 12, 8`.
* Every class with total length `≤ 14` matches: 79 classes with `Z_core ≢ 0`, from 42,964,024 excess-1 partitions at `L = 14` alone (`KERNEL TABLE: ALL MATCH`). Every nonzero table class occurs.
* This includes the dumbbell classes `(2,2;2) = 4`, `(2,2;4) = −4`, `(2,6;2) = 8`, `(2,6;4) = −8`, `(2,4;2) = (4,4;2) = (2,4;4) = (2,8;2) = 0`, and the theta classes with a length-1 path.
* Periodicity is visible directly: figure-eight `(2,10) = (2,6) = −7` and `(6,6) = (2,2) = −7/2`.
* Classes with an odd loop have nonzero coefficients but `Z_core ≡ 0`.
* Periodicity in a loop length also shows at `L = 14`: dumbbell `(6,6;2) = (2,2;2) = 4`, `(2,10;2) = (2,6;2) = 8`, `(2,6;6) = 8`.

**Four-periodicity.** The generating-function argument is right:
* one segment: `sz/(1−sz)`;
* two or more segments: `−s_1s_2 z²(1−z²)/((1−s_1z)(1−s_2z)(1+z²))`. This reduces to simple poles when `s_1 = ±s_2`, and has degree ≤ 0 at infinity.
* The kernel is a finite sum, over core configurations, of products of such univariate factors. So it is 4-periodic in each length for every length `≥ 1`.

PROVED (computer-assisted kernel values) is justified, and the values are now confirmed independently.

**Minor.**
* The Fourier statement for the theta, `(1/n) Σ_(k_1+k_2+k_3≡0)`, needs a circulant (cyclic Cayley) `T`, not merely a vertex-transitive one (edit T8, low).
* The series–parallel remark is right that `K_4` (excess 2) is the first non-series–parallel core.
* But its next sentence, "That is where, for example, Paley character sums stop collapsing to polynomials in `q`", is false (medium).
  * By the mechanism of HYP-9365 (d), every three-loop term of a doubly regular tournament, `K_4` included, is fixed by the 4-vertex census, i.e. by `n`.
  * So for Paley tournaments the three-loop character sums still collapse. The individual arithmetic of a doubly regular tournament (its 5- and 6-vertex censuses) enters at four loops.
  * The results file repeats the claim (edits T10, R4).

## D. HYP-9365: numbers, correlations, limits

**Exact `H`, recomputed** with `hcount.c` (inclusion–exclusion with rotation orbits, `__int128`) and `hcount_group.c` (orbits of `(Z_3)^3`):
* `R_n` for every odd `n ≤ 31` (`d_rot_paley_H*.out`);
* `P_19`, `P_23`, `P_31`;
* `P_27` on `F_27 = F_3[t]/(t³−t+1)`, `H = 418425945688463432031` (`d_paley27.out`; not recomputed in audit 1);
* every circulant class in the repo files for `n = 9..25` (`e_populations_circ.out`, 0 mismatches);
* all 15 regular 9-tournaments and all 1223 regular 11-tournaments (Held–Karp);
* `R_33`: `H = 5394851076280218488551446699`, by four-way parallel inclusion–exclusion (`d_rot33_H_parallel.out`, 1860 s), which matches the session's `R33.out`. (Line filled by the session from the auditor's recorded output; the auditor's run stalled after writing it.)

All agree with the repo values.

**Population statistics.** These are an independent implementation (`e_predictor.py`; dense matrices, cross-checked against an FFT version to `1e-13`). Spearman uses average ranks.

| population | corr | Spearman | mean `y − y_pred` | sd |
|---|---|---|---|---|
| all regular `n = 9` (15) | 0.9700 | 0.9643 | −0.0528 | 0.0718 |
| all regular `n = 11` (1223) | 0.9777 | 0.9731 | −0.0608 | 0.0227 |
| circulants `n = 21` (88) | 0.9969 | 0.9920 | −0.0409 | 0.0142 |
| circulants `n = 23` (94) | 0.9971 | 0.9956 | −0.0398 | 0.0110 |
| circulants `n = 25` (208; all `H` recomputed) | 0.9984 | 0.9974 | −0.0374 | 0.0090 |

* At `n = 11`: Spearman(`κ`, `ρ`) = −0.7428 and Spearman(`κ`, `ρ(1 − y_pred/n)`) = +0.8956. The argmax is `P_11` (`H = 95095`) for the corrected predictor; for `ρ` alone it is `R_11`.
* All as claimed.
* **Uniformity** (`e_uniformity.out`): `n(y − y_pred)` lies in [−2.69, +0.07] over all 1223 regular 11-tournaments and in [−1.39, −0.39] over the 208 circulant classes at `n = 25`.
* The mean offset times `n` is −0.72, −0.79, −0.86, −0.86, −0.92, −0.94 for circulants with `n = 15..25`. This is the same slow drift as for `R_n`; HYP-9365's refutation target is the right one to watch.
* The `n = 9` ranking numbers in `two_loop_regular.out` (−0.714, +0.832) come from tie-breaking by `argsort`. With average ranks they are −0.700 and +0.842; HYP-9365 does not quote them.

**Families** (`f_families.out`, `f_tower_limits.out`).
* `y(R_n)` and `y_pred(R_n)` for `n = 9..33`, `y(P_q)` for `q = 19, 23, 27, 31`, and the towers to 729: all as in `two_loop_predictor.out`, `two_loop_towers.out` and `paley_y.out`.
* `n(y − y_pred)` for `R_25..R_33` is −0.860, −0.907, −0.948, −0.984, −1.017. This is the drift HYP-9365 records. Over this range the data cannot tell convergence (`a + b/n` with `a ≈ −1.5`) from slow `log n` growth.

**Closed-form limits, rederived.**

*Rotational `R_n`.*
* The eigenvalues are `x_k = (1/n)cot(πk/2n)` for odd `k` and `−(1/n)tan(πk/2n)` for even `k`. So the large values are `±2/(πj)`, `j` odd.
* From `tan z = Σ_(j odd) 8z/(π²j²−4z²)` and `tanh z = Σ_(j odd) 8z/(π²j²+4z²)`: `tr G → tan 1` and `tr H → −tanh 1`.
* Hence `nF8 → −2 tan1 tanh1 − (3/2)tan²1 = −6.0105035`.
* `D_1(z) → z tan z + z tanh z` and `D_2 = z∂_z D_1|_1`, so `E_0 → 1 + tan²1 + tan1 tanh1 = 4.6116314`.
* The theta vanishes: the large frequencies are small odd signed integers, so any three of them have an odd sum smaller than `n` in absolute value, never `≡ 0 (mod n)`. Numerically `n·nΘ → 1.5574`.
* So `y_∞(R) = tan1 tanh1 + ½tan²1 − 1 = 1.3988720` is correct. `y_pred(R_100001) = 1.3988197`.

*Doubly regular `T`.* With `G = (I − J/n)/(n−1)`, `H = −(I − J/n)/(n+1)` and `O = S/(n+1)`, the two-loop value is an exact rational function, checked against dense computation for `q = 19, 23, 31`. Its expansion `1/2 − 1/(3n) + O(n^(−2))` also appears in the session's §12 (added 03:06); the closed form is new here:

    y_pred(n) = (3n² − 8n + 1)/(6(n−1)²) = 1/2 − 1/(3n) − 1/n² + O(n^(−3)).

* The limit is 1/2, as claimed.
* The residual `q(y − y_pred)` is −1.553, −1.528, −1.496, −1.494 for `q = 19, 23, 27, 31`: a stable `O(1/q)` three-loop term.
* `κ(P_q)/(e(1 − 3/(2q))) − 1` is 0.0094, 0.0063, 0.0045, 0.0034, about `3.3/q²`.
* These are better evidence for (a) than the quoted "two-point extrapolation 0.502". The 2-, 3- and 4-point extrapolations give 0.502, 0.479 and 0.425 (edit H4, low).

*Ternary towers.* `y_∞ = 1.3457257` is confirmed (Richardson to `k = 11`, `n = 177147`).
* But `n·Θ → −0.237054`, not −0.23703 (edit H7; a wrong digit, immaterial).

**(b) is overstated** (medium).
* "This settles the shape of the approach to Wormald's constant" is conditional on the conjecture.
* The unspecified `c` is `1 + y_∞ = 2.398872`, not `y_∞`. The reason: `ρ(R_n) = (cosh1/cos1)(1 − 1/n + O(n^(−2)))`. Numerically `n(1 − ρ/(cosh1/cos1)) = 0.99278, 0.99767, 0.99977` at `n = 33, 101, 1001`. This follows from Bhattacharjee's expansion, in which `e^(1/n)ρ(R_n)` has no `1/n` term. HYP-9365 should credit it (edit H6).
* The quoted extrapolations 1.3975 and 1.406 are the mixed 3- and 4-point fits. The 5- and 6-point fits give 1.475 and 1.968.
* Fitting within each residue class mod 4 (4 to 7 points) gives 1.3971–1.4002, which is better support for 1.398872 (edits H5, H11, R1).
* **Out-of-sample point** (`f_R35_check.out`). The session's `H(R_35) = 1610822294345031078388661430315` (`R35.out`, untracked, written after HYP-9365; not recomputed here) gives:
  * `y(R_35) = 1.22049` against `y_pred = 1.25038`, so `n(y − y_pred) = −1.046` (the drift continues: −0.984, −1.017, −1.046);
  * residue-class fits for `n ≡ 3 (mod 4)` through 35 of 1.3969, 1.3990, 1.3983, 1.3981, 1.3986 (3 to 7 points), all within `2·10^(−3)` of 1.398872.

**Small wording points.**
* (a): "`Σ x_k² = (n−1)/n` for every regular tournament, `ρ = e(1 − 1/n + O(n^(−2)))`". The identity for `Σ x_k²` holds for every tournament. The conclusion needs `Σ x^(2m) = (n−1)/n^m`, which holds for doubly regular tournaments only (edit H2).
* The Statement section uses `y(T)` for the two-loop predictor, while the evidence uses `y` for the exact residual (edit H12).

## E. Doubly regular tournaments

All data confirmed (`e_drt.py`, `e_drt_noIE.out`, `e_drt23_ie.out`, `e_drt_census6.out`):
* All six codes of `drt_results.out` satisfy `SSᵀ = nI − J`.
* By nauty, the `P_19` and `P_23` codes are cyclotomic Paley tournaments.
* `D_19.1` is not Paley, and it is self-converse (claim confirmed).
* `D_23.1–3` are pairwise non-isomorphic, not Paley, not self-converse, and no one of them is the converse of another.
* Every `H` was recomputed and matches: Held–Karp at order 19, unreduced inclusion–exclusion at order 23 (about 2 minutes for all four).
* Census: the 4-vertex censuses are equal, `(1026, 570, 570, 1710)` at order 19 and `(2530, 1265, 1265, 3795)` at order 23. The 5-vertex censuses differ; for transitive 5-sets, `P_19` has 342 and `D_19` has 360.
* Every listed tournament has spectrum `x² = 1/n`, `ρ = ((n+1)/(n−1))^((n−1)/2)` and the same `y_pred`; the dumbbell vanishes.

**Mechanism: correct, and stronger than stated.**
* `Ŝ² = (J − nI)/n²`, `ŜJ = 0` and `Ŝ∘Ŝ = (J − I)/n²`. So `span{I, J, Ŝ}` is closed under matrix and entrywise products.
* Every excess-`e` core has at most `2e` vertices.
* In a doubly regular tournament the out-neighbourhood of each vertex is regular of degree `(n−3)/4`. Each transitive 4-set has a unique source, and the other three vertices form a transitive triple in its out-neighbourhood. So the number of transitive 4-sets is `n·m·C((n−3)/4, 2)`, `m = (n−1)/2`: 1026 at order 19 and 2530 at order 23. The other three 4-counts then follow from the score identities.

Hence `X_0`, `X_1`, `X_2` are the same formal power series for all doubly regular tournaments of one order. That is a theorem given THM-4757 (2), not a heuristic. Only the step to `H(D)/H(D') = 1 + O(n^(−3))` is conjectural (edit H8).

**Scaling evidence overstated** (medium). "They have the same size at both orders, as an `O(n^(−3))` (four-loop) difference requires."
* There are only two orders, the sign changes between them, and the tournaments are unrelated. So `n^(−3)` cannot be separated from `n^(−2)` or `n^(−4)`.
* A new observation is consistent with the mechanism. `D_23.2` and `D_23.3` have **identical** 5-vertex censuses, but their 6-vertex censuses differ. Their `H` differ by only `2.8·10^(−8)` relative. The other order-23 pairs differ by `3.3·10^(−6)` (`D_23.1` against `D_23.2/3`) and `1.6–1.9·10^(−5)` (against `P_23`) (edit H9, R2).
* "Regular 19-tournaments spread over several percent" (`e_random_regular19.out`).
  * The 30 circulant classes span 1.85%.
  * Twelve random regular 19-tournaments (long chains of 3-cycle reversals from `R_19`) all have `H/H(P_19) − 1` in [−0.55%, −0.45%].
  * So "several percent" describes the extremes. The fairer yardstick for the doubly regular gap of `3·10^(−5)` is the spread of about `10^(−3)` among typical regular tournaments, which is still 30 times larger (low).
* Not checked: the order-15 sentence (no data in the repo).

## F. The wobble (HYP-9365 (e))

**Confirmed (FINITE-EXACT; `f_fneg*.out`, `g_wobble.out`).**
* `F(R_n,−1)` recomputed for `n ≤ 21`: 6, 0, 112, −2304, 84480, −4499456, 330332160, −31989497856, 3951688548352, −606498931605504.
* The sign is `(−1)^((n+1)/2)` for every odd `n ≤ 21` except `n = 5`, which gives 0.
* `|F/n!|` shrinks by 0.534, 0.577, 0.584, 0.591, 0.597, 0.601, 0.6045 per vertex: "about 0.6", slowly rising.

**The ratio depends on the fit.**
* With the session's definition of `w_n` (cubic in `1/n` through the 4 nearest points of the other residue class), `w_n κ_n/(F/n!)` is 1.638, 1.663, 1.621, 1.693, 1.615 for `n = 13..21`. That is 1.62–1.69, not 1.62–1.72 (a wrong digit, immaterial).
* A quartic fit gives 1.60–1.70; a quadratic fit gives 0.40–2.76.
* The `w_n` themselves keep the sign pattern for the cubic and quartic fits through `n = 31`, and lose it for the quadratic.
* So "tracks" and the ratio are NUMEROLOGY (edit H10, R3). "Non-perturbative" is an interpretation; the loop expansion's terms are powers of `1/n`.

## G. Prior art (searched 2026-10-22)

**(i) Counts near `e·n!/2^(n−1)`.**
* For a *uniformly random regular* tournament, `E[H] ∼ e·n!/2^(n−1)` is known. Sources: Wormald, arXiv:2610.16245, the `d = 1` case of Theorems 1 and 5, via his Lemma 3 and McKay's method; and Alon's observation that the Adler–Alon–Ross bound holds for regular tournaments, recorded in Yuster, *On the maximum number of spanning copies of an orientation in a tournament* (https://math.haifa.ac.il/raphy/papers/container.pdf, arXiv:1511.07784).
* Wormald's Lemma 3 is the averaged counterpart of the forest identity: a fixed linear forest with `a` arcs and `s` nontrivial paths lies in a random regular tournament with probability `∼ 2^(−a) e^((a−s)/n)`.
* **Missing credit** (medium): neither THM-4757 nor HYP-9365 (a) cites these (edits T3, H3). Nothing was found on Paley, doubly regular or quasirandom tournaments individually.

**(ii) Loop or cluster expansions, and the determinant ratio for a single tournament.** None found.
* Wormald, https://arxiv.org/abs/2610.16245: blow-ups only, `K(y) = det(I − (I/2 + A)D_y)`, no skew-spectrum formula.
* Bhattacharjee, https://github.com/creelie/Breakthrough-problems-1 (release v1.0.0): blow-up gain only. Its relevant extra is `ρ_R(R_d) = (cosh1/cos1)(1 − (tan1 − tanh1)/(3d²) + O(d^(−4)))`, used above. It has no forests, loops, single-tournament `H(T)` or `y`.
* Adler–Alon–Ross, https://web.math.princeton.edu/~nalon/PDFS/aar4.pdf: triangle-correlated random tournaments, `e` from a Poisson(1) count of consecutive triangles. THM-4757's description is accurate.
* Cuckler, *Hamiltonian cycles in regular tournaments*, CPC 2007 (https://www.cambridge.org/core/journals/combinatorics-probability-and-computing/article/abs/hamiltonian-cycles-in-regular-tournaments/5CB605A35014A30F60F740D74199D3E7): lower bound `n!/(2+o(1))^n`.
* Not relevant: Calkin–Novick–Ushijima-Mwesigwa (arXiv:1506.00699), and Zhang's dense-digraph counting algorithm (arXiv:0812.1123).

**(iii) Near-equal `H` across non-isomorphic doubly regular tournaments, or skew Hadamard matrices.** None found.
* Classification counts of skew Hadamard matrices (orders 16, 20, 24: 2, 2, 16) are tabulated, unverified, in https://math.ipm.ac.ir/~tayfeh-r/papersandpreprints/SkewH36x36.pdf. No Hamiltonian data there.

So the single-tournament loop expansion, the two-loop formula and the doubly regular universality mechanism appear new, as far as these searches go.

## Scripts and outputs (`04-computation/experiments/sierpinski_ocf_20261022/audit2/`)

**Forests and clusters**
* `a_forest_identity.py/.out`, `a_forward_step_link.py/.out`: (A).
* `common.py`: tournaments and exact `X_0`/`X_1` series.
* `b_partitions.py`: set-partition brute force.
* `b_mscaling.py`: blow-up-scaling brute force.
* `b_validate.py/.out`, `b_cycle_matchings.py/.out`: (B).

**Two-loop checks**
* `c_two_loop_check.py/.out`, `c_partition_L12.py/.out`, `c_session_cases.py/.out`, `c_dumbbell_n7_L16.py/.out`, `c_eulerian.py` (+ `c_eulerian_L20.out`), `c_random_regular15.py/.out`, `c_random_regular15_L10.py/.out`: (C).
* `d_kernels.py`, `d_kernels_L12.out`, `d_kernels_L14.out`: kernel table.

**Exact counts**
* `hcount.c`: Held–Karp, inclusion–exclusion with orbits, `F(T,−1)`.
* `hcount_group.c`.
* `d_rot_paley_H.sh/.out`, `d_rot_paley_H_big.sh/.out`, `d_rot33_H.out`, `d_paley27.py/.out`.

**Predictor and families**
* `e_predictor.py`, `e_populations.py`, `e_populations_reg.out`, `e_populations_circ.out`, `e_uniformity.py/.out`.
* `f_families.py/.out`, `f_tower_limits.py/.out`, `f_R35_check.py/.out`.

**Doubly regular tournaments**
* `e_drt.py`, `e_drt_noIE.out`, `e_drt23_ie.out`, `e_drt_census6.py/.out`, `e_random_regular19.py/.out`, `drt_codes.txt`, `drt23_codes.txt`.

**Wobble**
* `f_fneg.sh`, `f_fneg.out`, `f_fneg_small.out`, `g_wobble.py/.out`.

**Edits**
* `h_required_edits.py`: the edits below as data. It checks that each old string occurs exactly once (all do).

## Required edits (exact old -> new)

Each old string occurs exactly once in its file, checked by `h_required_edits.py` on the current files. Apply in any order.

**T1** `01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`

old:
```
  PROVED, pre-audit (mac-mini-2026-10-22-sierp): (1)-(3) are elementary identities, (3) as a formal power series in z,
```
new:
```
  PROVED (mac-mini-2026-10-22-sierp; independently audited 2026-10-22, 05-knowledge/results/sierpinski_audit2_20261022.md:
  (1)-(4) confirmed by exact rational brute force; one side remark corrected): (1)-(3) are elementary identities, (3) as a formal power series in z,
```

**T2** `01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`

old:
```
  the closed formula's z^4 and z^6 coefficients equal the brute-force excess-1 cluster sums for R_15 and P_19, and the
  brute-force excess-0 sums equal the log rho coefficients for L <= 6 (cluster_check.py, two_loop_series_check.py).
```
new:
```
  the closed formula's z^4 and z^6 coefficients agree with the brute-force excess-1 cluster sums for R_15 and P_19 to the
  4 digits printed (vertex-transitive cases; the dumbbell has no term below z^14 for any regular T), and the brute-force
  excess-0 sums equal the log rho coefficients for L <= 6 (cluster_check.py, two_loop_series_check.py). AUDIT
  (sierpinski_ocf_20261022/audit2/): exact rational equality of the excess-0 and excess-1 cluster sums with (3) and (4),
  by set-partition enumeration and by the blow-up scaling log Psi(M (x) J_m/m) = sum_e m^(-e) X_e(M), for all regular
  tournaments with n = 7 (to z^12) and n = 9 (to z^10), R_11, P_11 and five non-vertex-transitive n = 11 classes (to
  z^8), R_15 and P_19 (to z^8), the non-circulant regular 7-tournament to z^16 (dumbbell nonzero at z^14, z^16) and
  Eulerian skew matrices (one to z^20; dumbbell from z^6); the kernel table equals the brute-force coefficient of every core class
  of total length <= 14.
```

**T3** `01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`

old:
```
  walks; Wormald's generating functions use (I - AD)^-1 for blow-ups.
```
new:
```
  walks; Wormald's generating functions use (I - AD)^-1 for blow-ups. Wormald's Lemma 3 (a fixed linear forest with a arcs
  and s nontrivial paths lies in a uniformly random regular tournament with probability ~ 2^(-a) e^((a-s)/n), by McKay's
  method) is the averaged counterpart of the forest identity (1).
```

**T4** `01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`

old:
```
      log Φ_T(z) = Σ_(e ≥ 0) X_e(z),
```
new:
```
      log Φ_T(z) = Σ_(e ≥ −1) X_e(z),
```

**T5** `01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`

old:
```
of **excess** `e = #arcs − #blocks`, i.e. cyclomatic number `e + 1`.
```
new:
```
of **excess** `e = #arcs − #blocks`, i.e. cyclomatic number `e + 1`. Trees (`e = −1`) occur exactly when `T` is not regular (`audit2/b_validate.out`).
```

**T6** `01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`

old:
```
* The `z⁴` and `z⁶` coefficients match brute force: −0.19911 against −0.19910 and −0.060250 against −0.060250 for `R_15`; −0.16271 against −0.16270 and −0.0068139 against −0.0068140 for `P_19` (`two_loop_series_check.py`).
```
new:
```
* The `z⁴`, `z⁶`, `z⁸` coefficients equal the brute-force cluster sums exactly: −0.1991111, −0.0602495, −0.0526505 for `R_15`; −0.1627059, −0.0068139, −0.0012067 for `P_19` (`audit2/c_session_cases.out`; `two_loop_series_check.py` compared 4 printed digits). Non-vertex-transitive tournaments and the dumbbell: `audit2/c_two_loop_check.out`, `audit2/c_dumbbell_n7_L16.out`.
```

**T7** `01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`

old:
```
* For vertex-transitive `T`, `d` is constant and `G·1 = 0`, so the dumbbell vanishes. The figure-eight becomes
```
new:
```
* For vertex-transitive `T`, `d` is constant and `G·1 = 0`, so the dumbbell vanishes. For every regular `T` the first term of `d` is the constant `2 diag(Ŝ²)`, so the dumbbell has no `z^L` term with `L < 14` (first term: loops `a = b = 6`, bridge `c = 2`, i.e. `(a, b; c) = (6, 6; 2)`). The figure-eight becomes
```

**T8** `01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`

old:
```
a function of the spectrum alone. The theta becomes `(1/n) Σ_(k_1+k_2+k_3 ≡ 0)`
```
new:
```
a function of the spectrum alone. For circulant `T` the theta becomes `(1/n) Σ_(k_1+k_2+k_3 ≡ 0)`
```

**T9** `01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`

old:
```
The top term is `N_(n−1) = Σ_σ (−1)^(#backward steps) = (−1)^(n−1) F(T, −1)`, the alternating value of the forward-step polynomial of THM-4756.
```
new:
```
The top term is `N_(n−1) = Σ_σ (−1)^(#backward steps) = (−1)^(n−1) F(T, −1)`, the alternating value of the forward-step polynomial of THM-4756.

Reversing one odd path (say the one through the least vertex lying on an odd path) is a sign-reversing involution, so only forests whose paths all have even length contribute to `N_L`; in particular `N_L = 0` for odd `L` (audit, `audit2/a_forest_identity.out`). The same involution on a single run kills every cluster that contains an odd run.
```

**T10** `01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`

old:
```
That is where, for example, Paley character sums stop collapsing to polynomials in `q`.
```
new:
```
For doubly regular tournaments, Paley included, the three-loop terms (`K_4` among them) still reduce to functions of `n`, because they are fixed by the 4-vertex census (HYP-9365 (d)); Paley character sums do not stop collapsing at three loops.
```

**H1** `05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md`

old:
```
  OPEN, pre-audit (mac-mini-2026-10-22-sierp). The expansion is derived (THM-4757: exact forest identity, proved one- and
```
new:
```
  OPEN (mac-mini-2026-10-22-sierp; independently audited 2026-10-22, 05-knowledge/results/sierpinski_audit2_20261022.md:
  every exact H, y and y_pred recomputed and confirmed; wording and credits corrected). The expansion is derived (THM-4757: exact forest identity, proved one- and
```

**H2** `05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md`

old:
```
* Since `Σ_k x_k² = (n−1)/n` for every regular tournament, `ρ = e(1 − 1/n + O(n^(−2)))`. Hence
```
new:
```
* Since `Σ_k x_k² = (n−1)/n` for every tournament, and here `Σ_k x_k^(2m) = (n−1)/n^m`, `ρ = ((n+1)/(n−1))^((n−1)/2) = e(1 − 1/n + O(n^(−2)))`. Hence
```

**H3** `05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md`

old:
```
* This sharpens HYP-9274 (A) with an explicit rate. The `1/n` coefficient
```
new:
```
* This sharpens HYP-9274 (A) with an explicit rate. The constant `e` itself is the classical random-regular value: a uniformly random regular tournament has `E[H] ∼ e·n!/2^(n−1)` (Wormald, arXiv:2610.16245: the case `d = 1` of Theorems 1 and 5, via Lemma 3 and McKay's method; Alon's observation that the Adler–Alon–Ross bound holds for regular tournaments, recorded by Yuster, arXiv:1511.07784). New here are the single-tournament statement and the rate `3/(2n)`. The `1/n` coefficient
```

**H4** `05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md`

old:
```
* Data: `y(P_q) = 0.398, 0.417, 0.431, 0.440` at `q = 19, 23, 27, 31`; the two-point extrapolation gives 0.502.
```
new:
```
* Data: `y(P_q) = 0.398, 0.417, 0.431, 0.440` at `q = 19, 23, 27, 31`. Polynomial extrapolation in `1/q` gives 0.502, 0.479, 0.425 from the last 2, 3, 4 points, so it fixes the limit only loosely. For every doubly regular `T` the two-loop value is exactly `y_pred = (3n² − 8n + 1)/(6(n − 1)²) = 1/2 − 1/(3n) − 1/n² + O(n^(−3))`, and the residual `q(y − y_pred) = −1.55, −1.53, −1.50, −1.49` is stable, as an `O(1/q)` three-loop term should be; `κ(P_q)/(e(1 − 3/(2q))) − 1 = 0.0094, 0.0063, 0.0045, 0.0034`, about `3.3/q²` (`audit2/f_families.out`; `H(P_27)` recomputed in `audit2/d_paley27.out`).
```

**H5** `05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md`

old:
```
* Exact data: `y(R_33) = 1.2106`, with extrapolations 1.3975 (3 points) and 1.406 (4 points).
```
new:
```
* Exact data: `y(R_33) = 1.2106`. Mixed polynomial fits in `1/n` are unstable (1.381, 1.3975, 1.406, 1.475, 1.968 from the last 2–6 points) because of the `n mod 4` term of (e); fitted within each residue class mod 4, the 4- to 7-point fits give 1.3971–1.4002 (`audit2/f_families.out`).
```

**H6** `05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md`

old:
```
* Together with `ρ(R_n) → cosh 1/cos 1`, this settles the shape of the approach to Wormald's constant: `κ(R_n) = (cosh 1/cos 1)(1 − c/n + …)`.
```
new:
```
* Since `ρ(R_n) = (cosh 1/cos 1)(1 − 1/n + O(n^(−2)))` (Bhattacharjee's expansion of the regular-block gain `e^(1/n) ρ(R_n)` has no `1/n` term), the conjecture predicts `κ(R_n) = (cosh 1/cos 1)(1 − (1 + y_∞)/n + O(n^(−2)))`, with `1 + y_∞ = 2.398872`. This is conditional on HYP-9365; it does not yet settle the approach to Wormald's constant.
```

**H7** `05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md`

old:
```
Here the theta term survives, `n·Θ → −0.23703`, because
```
new:
```
Here the theta term survives, `n·Θ → −0.237054`, because
```

**H8** `05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md`

old:
```
* Hence non-isomorphic doubly regular tournaments of one order should agree to relative `O(n^(−3))`.
```
new:
```
* Hence `X_0`, `X_1`, `X_2` are the same formal power series for all doubly regular tournaments of one order (a theorem, given THM-4757 (2)). That non-isomorphic ones then agree to relative `O(n^(−3))` in `H` is the conjectural part.
```

**H9** `05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md`

old:
```
They have the same size at both orders, as an `O(n^(−3))` (four-loop) difference requires.
```
new:
```
They have the same size at both orders, which is consistent with an `O(n^(−3))` (four-loop) difference; two orders, with a sign change between them, cannot separate `n^(−3)` from `n^(−2)` or `n^(−4)`. At order 23, `D_23.2` and `D_23.3` have identical 5-vertex censuses (their 6-vertex censuses differ), and their `H` differ by only `2.8·10^(−8)` relative (`audit2/e_drt_noIE.out`, `audit2/e_drt_census6.out`).
```

**H10** `05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md`

old:
```
and `w_n κ_n / (F(R_n,−1)/n!) ≈ 1.62–1.72` for `n = 13–21` (`fneg1.out`).
```
new:
```
and `w_n κ_n / (F(R_n,−1)/n!) ≈ 1.62–1.69` for `n = 13–21` with the cubic smooth fit of `rot_limit_fit.py` (`fneg1.out`). The ratio depends on the fit: 1.60–1.70 with a quartic, 0.40–2.76 with a quadratic (`audit2/g_wobble.out`). The values and signs of `F(R_n, −1)` are FINITE-EXACT (`n ≤ 21`); the ratio is NUMEROLOGY.
```

**H11** `05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md`

old:
```
| `R_n`, `n ≤ 33` | 1.2106 (33) | 1.2415 | 1.398872 | 1.3975 / 1.406 |
```
new:
```
| `R_n`, `n ≤ 33` | 1.2106 (33) | 1.2415 | 1.398872 | 1.3971–1.4002 (within residue classes mod 4) |
```

**H12** `05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md`

old:
```
**Conjecture.** `κ(T) = H(T) 2^(n−1)/n! = ρ(T) (1 − y(T)/n + O(n^(−2)))`, uniformly over regular tournaments.
```
new:
```
**Conjecture.** `κ(T) = H(T) 2^(n−1)/n! = ρ(T) (1 − y(T)/n + O(n^(−2)))`, uniformly over regular tournaments.

*Notation below:* `y` is the exact residual `n(1 − κ/ρ)` and `y_pred` the two-loop value `y(T)` defined above; the conjecture says `y = y_pred + O(1/n)`.
```

**R1** `05-knowledge/results/sierpinski_steiner_collatz_tournaments_20261022.md`

old:
```
| rotational | `y → tan 1 tanh 1 + ½ tan² 1 − 1 = 1.398872`; exact-data extrapolation 1.3975–1.406 |
```
new:
```
| rotational | `y → tan 1 tanh 1 + ½ tan² 1 − 1 = 1.398872`; exact-data extrapolation within residue classes mod 4: 1.3971–1.4002 |
```

**R2** `05-knowledge/results/sierpinski_steiner_collatz_tournaments_20261022.md`

old:
```
The scaled gap `n³(H/H(P) − 1)` is ≈ 0.2 at both orders 19 and 23, the `O(n^(−3))` four-loop scaling (`drt_results.out`).
```
new:
```
The scaled gap `n³(H/H(P) − 1)` is ≈ 0.2 at both orders 19 and 23, consistent with the `O(n^(−3))` four-loop scaling (two orders only; `drt_results.out`).
```

**R3** `05-knowledge/results/sierpinski_steiner_collatz_tournaments_20261022.md`

old:
```
with ratio ≈ 1.65 for `n = 13–21`.
```
new:
```
with ratio ≈ 1.65 for `n = 13–21` (fit-dependent, NUMEROLOGY; `audit2/g_wobble.out`).
```

**R4** `05-knowledge/results/sierpinski_steiner_collatz_tournaments_20261022.md`

old:
```
There the arithmetic of the tournament (Paley character sums, subtournament censuses) appears for the first time.
```
new:
```
There the 4-vertex subtournament census enters as a genuine 4-index sum; for doubly regular tournaments, Paley included, that census is fixed by `n`, so their individual arithmetic (5- and 6-vertex censuses) enters only at four loops.
```
