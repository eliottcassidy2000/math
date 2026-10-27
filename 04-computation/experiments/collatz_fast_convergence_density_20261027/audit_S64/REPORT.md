# Independent audit of opus-2026-10-27-S64 (three-pillars note, HYP-9425, scripts, ledger rows)

Auditor: fable (adversarial, independent), 2026-10-27. Own scripts in this directory, written before the
session's scripts were read: `audit_tau_census.py` (+ `.out`, n < 2^24), `audit_tau_census_28.py` (+ `.out`,
n < 2^28, adds skewness/median of z, quarter-dyadic, 4-adic, 3-adic, finite-X model sums),
`audit_model_checks.py` (+ `.out`: sigma, E(C), C_F, entropy bound, 3x-1 seeds, Lagarias-Weiss constant).
No repository file was edited.

## Verdict: SOUND WITH CORRECTIONS

Every number recomputed reproduces (all of 2^24 exactly; the 2^28 tables, sawtooth figures, maxima and the
P1-P5 outcomes as well). Proposition F is correct. The sigma derivation, E(C), the seed lemma and the trunk
identity are correct. The corrections concern typing consistency, two stale references, one overstated word in
the title, one inequality weak by a factor 4 with a sign slip, and a few unverifiable or mis-attributed sentences.

## Corrections (major first)

1. **The pillar typing contradicts the note's own definitions and the atlas's recorded blindness** (note
   section 2 definitions; section 3 rows M, K-L, S6; section 0 items 1-2). Sentence: "*Positivity* (P): the
   statement uses `n > 0` and nonnegative weights, so that it changes under `n -> -n`"; row M: "K (weights; max
   over mean)"; row K-L: "K (counts)"; 0.1: "its SHEET wall (the sign) is the owner's positivity". Problem: by
   that (conclusion-based) definition a statement whose truth does not change under `n -> -n` drops P. The
   barrier atlas types Mazur SHEET B and INTEGRAL B (section 8) and Krasikov-Lagarias SHEET-blind and
   INTEGRAL-blind (section 1), and the note itself proves the seed lemma for `3x-1`, says "the mechanism cannot
   tell the sheets apart" (5.7) and "The mechanism is sign-blind" (0.6). So "Mazur's argument is the one row
   carrying all three pillars" is inconsistent with the note's own P (and I: the conclusion holds unchanged for
   the `Z_(2)` control). Likewise S6 has "O = K" while its wall column says "`3n-1`, `5n+1` rigid too", a drop
   under the O definition ("changes when the multiplier changes"). Fix: define K/D explicitly as what the proof
   USES, delete "so that it changes under ..." from the three definitions, and replace 0.1's identification of
   pillars with walls by "pillars are proof inputs, walls are conclusion controls; a proof can use positivity
   (Mazur's weights) and still be sheet-blind in conclusion". Or keep the conclusion rule and retype M and K-L
   as P = D, I = D, which removes the headline claim.
2. **Stale lattice-proximity test** (withdrawn in commit `9c744f56e2`, HYP-9425 (4) body). Still present in:
   the HYP-9425 **title** ("the floor set's `log_2 n` histogram at resolution 1/16 peaks within 1/8 of lattice
   points"); note 5.3 ("The clean test is to histogram `log_2 n` ... at resolution `1/16` ...; HYP-9425 (4)
   pre-registers it"); note section 1 table ("the lattice test (HYP-9425 (4))"); note 9 obligation 2 ("and the
   lattice test (HYP-9425)"); INDEX line (4). Fix: replace by the word-model test of the HYP body and record
   its reported PASS (not audited here) in the HYP status, the INDEX line and the note.
3. **Title overstates: "a no-parameter fluctuation law verified to `2^28`"** (title; 0.4 "holds to 5%"). What
   holds: block sd of `z` `0.949 -> 0.963`, mean `-0.038 -> -0.034` (reproduced). What does not: the Gaussian
   SHAPE. `audit_tau_census_28.out`: block skewness of `z` is `+0.60 -> +0.66` (k = 20..28, rising), the block
   median of `z` is `-0.168 -> -0.146` against the mean `-0.04`; the model's `G_k(8)` is off by `0.066`
   absolute at k = 20 (0.3597 vs 0.2934) and `G_k(10)` by `0.044`. Fix: "agrees in scale to within 5% to
   `2^28`; the shape is right-skewed (skewness `0.6`, OBSERVED), so the Gaussian profile is not verified".
   Consequently 5.2's "the median ... rising towards `C*` at about `0.025` per doubling (the model's
   finite-size bias ... consistent with a constant offset `b/(sigma sqrt(ln n))`)" is wrong as an explanation:
   a constant offset of `-2.4` steps puts the block median at `C* - 0.13` at k = 28, not at the observed
   `9.88 = C* - 0.55`; the median deficit is the skewness. HYP-9425 (1)'s "`G_k(C) = Phi(...) + O(1/sqrt(k))`"
   is unfalsifiable at these k (`1/sqrt(20) = 0.22`); tighten or drop.
4. **Trunk inequality: weak form and a sign slip** (5.7, 0.7, PROBLEM-LEDGER). The chain
   "`max_(j >= 2) H_n(R_j) >= (H_(n+1)(1) - (3/4) H_n(1))/(1/4) >= H_(n+1)(1) - (3/4) H_n(1)`" is false in its
   second step whenever `H_(n+1)(1) < (3/4) H_n(1)` (then `4x < x`); the final inequality still holds,
   trivially, because the left side is nonnegative. The sharp statement is
   `max_(j >= 2) H_n(R_j) >= 4 H_(n+1)(1) - 3 H_n(1)`; state that. (Weights `3*4^(-j)`, `j >= 2`, sum to `1/4`:
   verified; `R_1 = 1`: verified; the digest's list is `rho_18` and the `3/2` conversion is right: 0.4245,
   0.3405, 1.740, 0.0705, 1.4985; no mislabel.)
5. **Proposition F: one implicit step** (5.4). "the orbit reaches `m_0` after exactly `k_0 + d_0` ordinary
   steps" then "`tau(n) >= k_0 + d_0 + log_2 m_0`" needs that the orbit has not hit `1` before T-step `k_0`;
   true since `T^j(n) >= n/2^j >= 2^(k_0 - j) > 1` for `j < k_0`, but unstated. Add it. Every other step checked
   (two per class; Terras; `sum_(d <= pk) C(k,d) <= 2^(k h_2(p))` for `p <= 1/2`, verified to k = 400;
   `m_0 >= 3^(d_0)`; `2^A >= m`; `C_F = 3.307351`; `2.585 = 1 + log_2 3`; `k_0 >= log_2 n - 2`; the dyadic sum);
   `O(X^(1-c))` and "density zero for `C < 3.307`" are justified. Prior art: I do not know this statement in
   print. Known to me: the all-`n` bound `sigma_oo(n) >= log_2 n` (Lagarias 1985; Lagarias-Weiss 1992);
   Applegate-Lagarias 2003 (infinitely many `n` with `sigma_oo(n) > 6.14316 log n`, T-steps: a limsup lower
   bound, not almost-all); Lagarias-Weiss 1992 (renewal and large-deviation statements for the stochastic
   models only); Korec 1994 (the same LLN at horizon `log_2 y`, used for the upper bound `T^n(y) < y^c`). The
   two-sided LLN behind Korec gives the constant implicitly, so folklore is plausible; the flag is right and
   the statement should not be filed as a theorem without reading Lagarias 1985 section 2.4 and the
   Applegate-Lagarias introduction.
6. **Sign-blind transfer: one unverifiable sentence and one direction change** (5.7, HYP (3)). "the
   time-constant layer uses `S(y) <= 2y`, which `3x-1` satisfies": `S(y) <= 2y` holds for both maps
   (verified), but the paper's Section 7 text is not in the repository, so I cannot confirm it is the only
   sign-dependent input. Checkable now: for `3x-1` the size identity is `x = M 2^A 3^(-d) / prod(1 - 1/(3 x_j))`
   (factors `< 1`) against `prod(1 + 1/(3 x_j))` for `3x+1`; the inequality between `log(x/M)` and
   `A log 2 - d log 3` reverses (it favours the `3x-1` time constant), so "hold verbatim" should read "hold with
   the size inequality reversed, which only helps". (Side remark, outside scope: the digest's P7 row prints
   `A log 2 - d log 3 <= log(x/M)`, the wrong direction for `3x+1`; the identity gives
   `log(x/M) = A log 2 - d log 3 - log Pi <= A log 2 - d log 3`.)
7. **`c = beta/(32M)`** (section 4, "what this changes" (i)) conflicts with digest audit item 80
   (`c^(-1) = 256 M^2 m/3`): `1/M` against `1/M^2`. Cite one formula with its equation number or neither.
8. **Quotation marks on paraphrases** (0.1 and 2(a)): "without adding invariant vectors" and "the maximal
   characteristic-two depth does not produce `n` legal toggles" are not sentences of `cancellation/RESULT.md`
   (it has "Higher root order supplies no additional XOR directions" and "Root depth25 is not a
   25-dimensional legal toggle group"). Theorem 1(3) itself is cited faithfully and the file is "PROVED,
   independently audited"; kind-pasteur's Lemma 1 is cited faithfully and credited. Remove the quotation marks.
9. **Small wording slips.** (a) 5.1: THM-4504's `41.677` is in T-steps (recomputed: `41.6776` at odd share
   `0.609`; ordinary-step analogue `67.06`); say so. (b) Section 1: "every parity word is the word of a rational
   cycle point" -> "every periodic parity word". (c) Section 3, S2 wall: "the negative cycles' complements" is
   unclear; THM-4476's status says "every divergent positive orbit of `3n+1` or `3n-1` has summable
   reciprocals", so "holds on both sheets" is right and the wall should say that. (d) 0.2: the identity
   `1/x = 3^d 2^(-A) prod(1 + 1/(3x_j))` is Steiner-Eliahou (digest Theorem C (c), THM-4633), not THM-4813 (A).
   (e) 5.3: "`n = (2^k - c_w)/3^d` lies within a factor `2` of `2^k/3^d`" is false for general words (`c_w` can
   exceed `2^(k-1)`); say "for the admissible low-weight words of the floor set".
10. **Script-level remark** (`fast_convergence_density_20261027.py`): densities come from a histogram of
    `tau/ln n` with bin width `0.01`, counted by bin centre `<= C`. For integer `C`, `10.46` and `10.48` this
    equals the exact count (verified to 4 decimals at every k), but the `3.307` column is `r < 3.31` (4718 at
    k = 20 against the exact floor count 4699, computed exactly in its own column) and medians have
    resolution `0.01`. Harmless for every number quoted; state it in the docstring.
11. **Ledger tier definition** (TRUSTED-ADVANCES, "Trust tiers"): the Mazur rows are T2 with the right caveats
    (no replay, no named-expert reading; the natural-density row says "axiom profile not shown"; the `x^0.90`
    row names `native_decide`), but the T2 definition ("Lean-checked ... with an open dispute ..., or a human
    preprint not yet refereed") does not cover "platform-kernel-checked, no dispute, no expert reading". Add
    that clause. Nowhere does the note call a Mazur theorem PROVED or VERIFIED (grep: PROVED/VERIFIED attach
    only to the note's own items). P3 and P4 are reported FAILED in note 0.5, 5.2, 7, HYP status, INDEX and
    PROBLEM-LEDGER (all checked).

## Checks run and results

1. tau census, own code, n < 2^24 (0.7 s) and n < 2^28 (2.9 s): tau(27) = 111, tau(97) = 118, tau(871) = 178,
   tau(6171) = 261, tau(77031) = 350, tau(837799) = 524 OK; max tau 956 at 226,588,897, max orbit value
   6,404,797,161,121,264 = 2^52.51 REPRODUCED. F_X(523/50) = 0.555047 / 0.542628 / 0.540373 at 2^16 / 2^20 /
   2^24 and 0.5550 ... 0.5420 for k = 16..28 REPRODUCED (`<=`/`<`, n >= 1 / n >= 2 variants agree to 1e-4).
   Block sd of z 0.94874, 0.95062, 0.95273, 0.95425, 0.95637 (k = 20..24) -> 0.963 at 28; means -0.0379,
   -0.0381, -0.0377, -0.0375, -0.0369 -> -0.034 REPRODUCED. G_k(8) = 0.3597 / 0.3249 / 0.2988, G_k(13) =
   0.7425 / 0.7651 / 0.7796 at k = 20 / 24 / 28 REPRODUCED. Floor counts 3489, 4699, 8153, 12233, 24130, 32643
   (k = 19..24), 67829, 91956, 167377, 239065 (k = 25..28), ratios 1.347 ... 1.428 REPRODUCED exactly.
   Quarter-dyadic densities (24.1, 19.8, 8.5, 18.8; 30.4, 28.7, 13.2, 27.5; 32.8, 40.7, 31.1, 57.1), 4-adic
   ratios 2.773, 2.490, 2.785, 2.814, 2.544, 3-adic ratios 1.947 ... 2.395, two-block exponents (mean 0.703,
   sd 0.052) REPRODUCED. Digest's "54.22% of n < 10^6": 542194/999999 = 54.2195% VERIFIED.
2. sigma: E[a] = 2, Var a = 2, C* = 10.428178, Var(c - C* delta) = 2(1 - C* ln 2)^2 = 77.5825, / E[delta] =
   0.287682 gives 269.68, sigma = 16.4220 VERIFIED (direct numeric variance agrees). The formula is the
   renewal-reward CLT: i.i.d. cycles (delta_i, c_i), N(t) = #cycles with cumulative drift <= t, R(t) = sum c_i,
   (R(t) - r t)/sqrt(t) -> N(0, Var(c - r delta)/E delta), r = E c/E delta (Smith 1955; Asmussen, Applied
   Probability and Queues, ch. VI), with t = ln n: the right form. 523/50 is 0.305% above C* VERIFIED.
3. Proposition F: every step checked (correction 5); entropy bound verified numerically; C_F = 3.307351.
4. E(C): rederived; unconstrained max 1 at delta/kappa = 1/2, kappa = 6.9521 = 1/(ln 2 - ln 3/2), tau/ln X =
   10.4282 (a scan over p confirms the maximum); boundary closed form p = (C ln 2 - 1)/(1 + C ln 3);
   E(2.5) = 0.5841, E(3.307) = 0.7383, E(8) = 0.9861, E(C*) = 1 VERIFIED. Finite-X sums of the same model:
   ln S/ln X = 0.63 at 2^16, 0.68 at 2^40, 0.716 at 2^140, local exponents 0.72-0.76 at 2^60-2^100: the
   asymptotic 0.738 is approached slowly, so the measured 0.70 +- 0.05 should be compared with the finite-X
   model value as well (not wrong as written, but the fair comparison is missing).
5. Trunk inequality: identity, weights, R_1 = 1, arithmetic VERIFIED; correction 4. The rho_18 values at R_j
   were not recomputed (level-18 FFT law).
6. 3x-1 seeds: r_j odd integer, 3 r_j - 1 = 2^(2j+1), v_3(2^(2j+1)+1) = 1 + v_3(2j+1) for j < 2000; r_0..r_(3^q-1)
   permute Z/3^q for q = 1..6; ord_(3^(q+1))(4) = 3^q; S(y) <= 2y for both maps VERIFIED.
7. Typing table against the atlas and coalescence rows: K-L "keeps positivity and integrality" contradicts the
   atlas (SHEET-blind, INTEGRAL-blind) under the note's definitions; "truncates at depth 11" is a fair
   description of the NLP (gamma = 0.84175 at k = 11), not of the difference inequalities; THM-4476 "holds on
   both sheets" CONFIRMED from its status; Tao sheet-blind via GGM covering 3N-1 CONFIRMED (atlas Tao and GGM
   rows). Correction 1.
8. Concurrent citations: faithful and credited (correction 8 for the quotation marks).
9. Ledger tiers: consistent with the caveats (correction 11). Web evidence: not re-fetched.
10. Wording: corrections 3, 9; P3/P4 reported failed everywhere.
11. Spot-check of `word_model_fine_resolution_20261027.py` (coordinator's late addition): weight
    `comb(k, d) 3^(-d) (k-d)/k = C(k-1, d) 3^(-d)` CORRECT for words with an even last step; carry
    `c_w = sum_i 3^(d-1-i) 2^(j_i)` (i-th odd step at 0-based position j_i) CORRECT (`2^k T^k(n) = 3^d n + c_w`).
    The rest of that script (cell assignment, `n = num // 3^d`, sampling) and its reported PASS (correlation
    0.992, sign agreement 0.889) NOT checked.

## Not checked

ProofAtlas pages and Tao's blog comments (not re-fetched); the Mazur paper's text (Section 7, Lemma 6.1, the
`c` formula, the "`S(y) <= 2y`" claim); the level-18 values `rho_18(R_j)`; section 10a's referee record (its
scratch scripts are not committed); the word-model test beyond the two formulas above; the literature question
beyond my own knowledge.
