# Independent audit of HYP-9391 (the seven-vertex Paley census)

2026-10-26.

The auditing agent could not write this file; the session saved its final report text.

**Method.**
- All checks used the auditor's own code (C, Python, PARI 2.17.4).
- No existing file was edited and nothing was committed.
- The auditor's one binary intermediate, `mE_fit307_9999.npy` (fitted residuals), was moved out of the repository by the session. It can be regenerated with `part_c_analysis.py`.

**Bottom line.**
- **Exact parts.** The TT7 formula (a) is exact at all 618 primes `q ≡ 3 (mod 4)` up to 9967, against the auditor's own counts. Both fourth-moment identities (b) are true.
- **Overclaims in (c).** There are three:
  - the stated projection-error bound is wrong by a factor of about 3;
  - the statistics cannot pick out USp(4), or rank 4;
  - (b) is a known theorem, but it was presented as FINITE-EXACT and uncredited.

## Verdicts

| Claim | Verdict |
|---|---|
| (a) O(q²) closed formula for TT7 | **Confirmed.** Exact at all 618 primes `q ≡ 3 (mod 4)` in 7..9967, with independent TT7 counts and independent Y sums. The counts also match all 69 lane data-file values (19..811). |
| (a) 37 identities `S_H = known ± Y` | **Confirmed for all 37**, not just a sample: brute-force sums at 18 primes 7..223. q = 211 and 223 are new for H23, H34, H36 and H37. `Σ m_H S_H = g(q)` and `S(K7 − C7) = (q−6)S(K6 − P6)` also hold at all 18. |
| (a) Y closed forms and cone reductions | **Confirmed** at 618 primes. Literal O(q³) cone sums with a different fixed pair agree at the 14 primes below 110. |
| `Y04 = −Y02`, `Y06 = Y05` | **True, and PROVED** here for every `q ≡ 3 (mod 4)`; the HYP gave no proof. In `cheapY_5000_9999.txt` these two columns are imposed by construction, not data (`run_cheap_Y4.py` writes them). |
| (a) "Exactly four new directions" | **Confirmed.** The base has rank 9; adding the four Y gives exact rank 13; the 16 "new" sums span the same space. |
| (b) Twisted and untwisted fourth moments | **True, and KNOWN.** Both are the k = 6 case of Frechette–Ono–Papanikolas, IMRN 2004, Prop. 2.1. The status should be PROVED with credit, not FINITE-EXACT. Re-checked at all 428 odd primes 5..2999. |
| (b) Monodromy / Euler–Poincaré count | **Correct**, and consistent with the constants −1 and −3. |
| (c) Irreducible, weight 5, multiplicity one | **Plausible heuristic.** It assumes integral multiplicities and orthogonality on the class `q ≡ 3 (mod 4)`. The Gram matrix is within 0.047 of the identity on 587 primes. |
| (c) Not GL(2) | **Supported**, but the stated "< 0.05" error bound is false. |
| (c) Rank ≥ 4 | Rank ≥ 3 is supported. Rank ≥ 4 depends on the conjecture that odd-weight self-dual motives are symplectic. |
| (c) "USp(4) Sato–Tate statistics" / "genus-two-type" | **Not established.** USp(6), USp(8) and a Gaussian fit equally well or better. |
| (c) Not quadratic twists, not induced | **Supported** for the characters tested. The induced-motive test as worded is flawed. |
| (c) Siegel weights (4) or (3, 2) | Incomplete as stated. The mod-q test supports excluding the omitted Tate-twisted alternatives. |
| (c) Mechanism (Katz convolution) | **Misapplied**; at most heuristic. |
| Results §4b, k = 7 table row, HYP-9390 | Overclaims ("rank-4", "resolved structurally", "match USp(4)") and minor inaccuracies. |
| Numbers vs outputs | All reproduced, except the `P(|x| < 0.2)` row and the "3–5%" figure, which came from a different sample. |

## What was run

### A. Graph sums

`brute_sums.c`:
- brute force over all injective maps, fixing different vertex pairs from the lane's (0, 1);
- all 38 sums equal the lane's at every shared entry (38 sums at q ≤ 131; 33 at q = 211 and 223).

Side finding: the 37 genuine sums take only 19 distinct values up to sign. For example `S_H36 = S_H05 = S_H06 = −S_H21`, which is why H36, with no degree-3 vertex, satisfies the Y05 identity.

### B. The Y sums

`ysums.c` computes the closed forms directly in O(q²), together with independently hand-derived closed forms for the H04 and H06 cone reductions. All 618 lane `cheapY` rows match.

**Proof of `Y04 = −Y02`.**
- Let `A(u) = Σ_{s≠0,1} φ(s)φ(u−s)L(s)²`.
- Then `−Y02 = Σ_u φ(u)L(u)A(u)` and `Y04 = Σ_u φ(u)φ(u−1)L(u)A(u)`.
- On `{u : φ(u−1) = −1}`, the map `u ↦ u/(u−1)` is a fixed-point-free involution. It fixes `L(u)` and `A(u)`, and flips the sign of `φ(u)`, so that part of the sum is zero.
- Replacing `φ(s)` by `φ(s)φ(s−1)` gives `Y06 = Y05` in the same way.
- Every step was checked termwise (`y_symmetry_proof_check.out`).

### C. TT7 counts

`tt7_count.c`:
- `x ↦ a²x + b` acts simply transitively on arcs, so `TT7 = (q(q−1)/2)·T5(D)` with `D = {x : φ(x) = φ(x−1) = 1}`;
- this is O(q⁴), not the HYP's "O(q⁵) or worse"; q = 4999 takes under a second;
- `TT7 = 0` exactly at q = 7, 11, 19, 23.

### D. Moments

The twisted and untwisted second and fourth moments, and the Sym⁴ bookkeeping `Σ φ·Sym⁴ = −(1 + c + d)`, `Σ Sym⁴ = −(3 + d)`, hold at all 428 odd primes 5..2999 (211 of them `≡ 1 (mod 4)`).

### E. Part (c)

Every lane coefficient and statistic was reproduced to the printed digits.

**Projection error.**
- At primes ≥ 2000, the fitted known part has a per-prime standard deviation of 0.14 (mean) to 0.19 (max).
- Coefficient standard errors: about 21 for `q²`, 2.3 for `|b|²`, 100 for `a_q`, and `2·10⁴` for `q`.
- The fitted low-order terms themselves contribute 0.085–0.305 there.
- Refitting with the 280 primes in 5000..9999 moves the `q²` coefficients by up to 17; for Y09 it goes from −13.1 to +4.2.
- Rounding the near-integer coefficients changes individual values by up to 0.27.

**Not GL(2).** On 464 primes ≥ 2000 (refit on 587 primes):
- the counts of `|x| > 2.5` are 5, 7, 8 and 3, against an expected 0.02–0.5 for a GL(2) motive plus this error;
- the maxima are 3.17, 3.18, 3.07 and 2.97.

**Statistics on 587 primes.**

| | m1 | m2 | m3 | m4 | USp(4) |
|---|---|---|---|---|---|
| `E[x⁴]` | 3.13 | 2.96 | 2.88 | 2.70 | 3 |
| `E[x⁶]` | 16.0 | 13.0 | 13.6 | 10.5 | 14 |
| share `|x| > 2` (q ≥ 2000) | 5.2% | 5.0% | 4.5% | 4.5% | 4.95% |

Kolmogorov–Smirnov distances (5% critical value 0.056):
- rejected: SU(2)×SU(2) 0.062–0.088; Sym³ SU(2) 0.085–0.103;
- not rejected: USp(4) 0.021–0.032; USp(6) 0.018–0.031; USp(8) 0.017–0.030; N(0,1) 0.018–0.031.

Since `E[x⁴] = 3` for every USp(2n) with n ≥ 2, the rank is not determined.

The reference values are exact, not Monte Carlo:
- USp(4): moments 1, 3, 14, 84; `P(|x| < 0.2) = 0.1698`; `P(|x| > 2) = 0.0495`;
- SU(2): moments 1, 2, 5; `P(|x| < 0.2) = 0.1262`.

**Mod-q test (CRT + LLL lattice).**
- No small-integer congruence `Y ≡ c0 + c1 a_q + c2 c_q + c3 d_q + c4 |b|² (mod q)` holds for Y02, Y05, Y09, Y18 or Y05 + Y18.
- This stays true after adding `(2/q)`, `(3/q)`, `h(−q)`, `(2/q)h(−q)` and `h(−8q)`. Positive controls are recovered.
- So each new trace is nonzero mod q: each motive has a (5,0) Hodge part. This heuristically excludes the "q × weight-3 motive" and "q² × abelian surface" alternatives.

**Twists and identification re-runs.**
- Over 10 quadratic characters, `|E[χ m_i m_j]|` is at most 0.174 on 307 primes and 0.136 on 587.
- Re-running the lane's identification scripts (43 and 71 primes, 260 odd-level orbits, 329 products): none reduces a new direction. The HYP had no `.out` receipts for them; they are now `rerun_identify_*.out`.

### F. Literature

**(b) is Frechette–Ono–Papanikolas, IMRN 2004, no. 60, 3233–3262, Prop. 2.1.**
- Part (1) is Ahlgren, Finite Fields Appl. 8 (2002) (the `η(2z)¹²` fivefold).
- Part (2) is Ahlgren–Ono, Monatsh. Math. 129 (2000).
- The proposition is stated for all even k ≥ 4. Both identities follow, using:
  - `Σ_{λ=2}^{p−2} f(λ²) = Σ_μ (1 + φ(μ)) f(μ)`;
  - `Tr₆(Γ₀(8)) = c + 2d` (the space has dimension 3);
  - `Tr₆(Γ₀(4)) = d`.

**Katz's concatenation theorem** (Exponential Sums and Differential Equations, ch. 8) is stated correctly but misapplied.
- The true s-integrands are `φ(s)L(s)²`, `φ(s−1)L(s)²`, `φ(s)φ(s−1)L(s)` and `φ(s−1)L(s)`.
- `L(s)²` has rank 4, and the `φ(s−1)` twist gives local monodromy `−J₂` at `s = 1`, which is not a pseudo-reflection. So none of these is a hypergeometric convolution; Katz's middle convolution is the relevant operation.
- By the auditor's Euler–Poincaré count, the "rank-8" local system is right for Y02 only; Y05, Y09 and Y18 give rank 12.

**Prior work worth citing** (motivic weight 3, not 5): rank-4 weight-3 Calabi–Yau motives and weight-3 paramodular forms.
- Golyshev–van Straten, "Congruences via fibered motives", PAMQ 2023.
- Dummigan–Tornaría, arXiv:2412.14289.
- Poor–Shurman–Yuen, arXiv:2302.00764.

No prior art was found for these weight-5 triple products.

## Errors by severity

**Major.** None. Parts (a) and (b) are mathematically correct.

**Moderate.**
1. The bound "low-order terms contribute less than about 0.05" is false. The real error is about 0.14–0.19 per prime, up to about 0.3 under refits. The not-GL(2) conclusion survives; small-value frequencies become unreliable.
2. "USp(4) Sato–Tate statistics", "genus-two-type", "rank-4" and "four rank-4 Siegel-type motives" are not supported: USp(4) cannot be separated from USp(6), USp(8) or a Gaussian.
3. (b) is a known theorem presented as uncredited and FINITE-EXACT ("may be a known case; not checked").
4. Numerical inferences from (c) are stated as facts in the HYP title, results §4b, the k = 7 table row, and the HYP-9390 title and bullet ("Resolved structurally").

**Minor.**

5. The Katz mechanism is misapplied, and "rank-8" is wrong for Y05, Y09 and Y18.
6. The `P(|x| < 0.2)` row (0.178, 0.153, 0.143, 0.153) and the "3–5%" figure:
   - they come from `induced_test.py` run on 314 primes: the 307 plus 7 rows written meanwhile by the then-running job;
   - so they are not reproducible from the committed data;
   - on the stated 307 primes the values are 0.179, 0.150, 0.169 and 0.156, and the share is 2–6%.
7. The induced-motive test, "3–5% against about 50%", ignores that the projection noise is three times the 0.05 window. An induced motive would show about 15%.
8. `Y04 = −Y02` and `Y06 = Y05` were asserted "exactly" without proof, and some data columns are imposed by construction.
9. "O(q⁵) or worse for counting" is wrong: counting is O(q⁴).
10. The Siegel weight list omits Tate twists; the mod-q test supports excluding them.
11. Receipts and script list:
    - the `identify_*.py` outputs are missing;
    - `run_cheap_Y4.py` and `untwisted_moment.py` are not listed;
    - `run_fast.py` is described as "all survivor sums" with g "checked at every prime", but it computes 33 of them and cannot check g.
12. Results §4b:
    - it says "43 primes", but the four sums without a degree-3 vertex have only 23;
    - it says the identities come from cone vertices, which is not true for H36.
13. Wording: "odd weight forces symplectic" is a conjecture; and "compared with Monte Carlo values" is wrong, since the reference values are exact.

## Files

Everything is in this directory (`audit_hyp9391/`).

**Code:**
- `tt7_count.c`, `ysums.c`, `brute_sums.c`;
- `forms_dump.gp`, `modq_hodge_test.gp`;
- `literal_cones.py`, `tt7_formula_check.py`, `identity_check.py`, `rank_check.py`, `moments_check.py`;
- `y_symmetry_proof_check.py`, `part_c_analysis.py`, `lowblock_contrib.py`, `usp_rank_references.py`, `gl2_noise_test.py`;
- `build_required_edits.py`.

**Outputs:** each script's `.out`, plus `rerun_identify_{new,odd_levels,Y,products}.out`.

**Data:**
- `forms_12000.txt`, `graphs38.txt`, `brute_sums_all.txt`;
- `tt7_T5_7_4999.txt`, `tt7_T5_5000_9999.txt`;
- `ysums_q3mod4_7_4999.txt`, `ysums_q3mod4_5000_9999.txt`, `ysums_oddprimes_3_3000.txt`;
- `modq_data.gp`.

**Edits:** `required_edits.json` (29 edits) — 19 to HYP-9391, 7 to the results note, 3 to HYP-9390. Every `old` string occurs exactly once, and they apply cleanly in sequence.

Compiled binaries were left in the session scratchpad, not in the repository.
