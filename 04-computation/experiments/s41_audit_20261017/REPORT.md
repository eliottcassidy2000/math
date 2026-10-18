# S41 independent audit (opus-2026-10-17-S41 claims)

Independent referee (general-purpose subagent), 2026-10-17. Audited state: commit e160c7c72b. The checkpoint commits
e32b58e6fb and e160c7c72b swept the referee's untracked files a1-a13 into Git; the committed copies match the final
versions. Line numbers refer to e160c7c72b. The referee could not write this file itself; it was saved verbatim from the
referee's final message by the session.

## Summary

Three errors need fixing:
1. **Calibration slope.** Logarithmic budgets "with b of order 10^3" are only certified if the logarithm is counted from
   the census size, i.e. `log2(X/2^37)`, which is how HYP-9283 writes them. In HYP-9288's own form `a + b log2 X`, the
   2^37 certificate accepts only `b <= 45.7` and the 2^40 certificate only `b <= 338`. The chain of implications itself
   is sound.
2. **Credit for γ_4-γ_6.** THM-4667 called them new, but clock_duality published γ_4-γ_10 at 19:31, before S41's first
   commit at 20:07.
3. **Two shell inequalities were conflated.** The note said the max-inside inequality is a necessary condition for
   HYP-9288 (B). It is not: (B) implies only the weaker sup-outside form, which the files also mislabelled "maximal
   inequality".

Everything else checked out apart from the overstatements and minor points below. That covers statements 1-4 and 8-10,
the incremental scan, the ideal ladder to d = 12, the operator norms, and all six-colour computations. clock_duality's
earlier repairs (the unrooted energy term in statement 2 and the Syracuse-offset identification in statement 10) are
correct.

---

## Errors (must fix)

**1. Calibration slope.**
- **Where:** THM-4667 lines 67-69; HYP-9288 line 36 and status lines 5-7; root_transport_energy_20261017.md lines 259-263.
- **Why:** HYP-9283(B) budgets are `a + b j` with `j = log2(X/N)` and `N = 2^37` (or 2^40). So HYP-9288's `a + b log2 X`
  equals `(a + 37b) + b j`.
- **Evidence (a7):**
  - The budget enters the certificate only through one sum `S(a, b)`; the certificate passes iff `S < S*_N`.
  - The three published margins of each certificate lie on one level set (spread about 10^-5), with `S*_37 = 145.555`
    and `S*_40 = 410.195`.
  - In the `log2 X` form the largest certified b is 45.7 (with `a = 51/7`) at `N = 2^37`, and 338 at `N = 2^40`.
  - At `b = 1000`, `S/S*` is 4.67 (2^37) and 1.72 (2^40), so neither certificate passes.
- **Fix:** restate (B) as `X^2 (a + b log2(X/2^37))` for `X >= 2^37`, or quote `b <= 45` (`N = 2^37`) and `b <= 338`
  (`N = 2^40`) for the `log2 X` form.

**2. Credit for γ_4-γ_6.**
- **Where:** THM-4667 lines 76-77, "the rest are new".
- **Why:** root_positive_transport_20261017g.md section 6 gives γ_4-γ_10 as rational enclosures, committed before S41's
  first commit.
- **Fix:** only γ_11, γ_12 and the energies E_d are new. The note's section 1 (lines 49-50) should cite the same table.

**3. Shell inequalities conflated.**
- **Where:** note lines 213-214 call `Σ_N max_d Vroot_d^2` "a necessary condition for HYP-9288 (B)"; HYP-9288 lines 3,
  38-39 and 62, and note line 254, call the sup-outside form a "maximal inequality".
- **Why:** (B) bounds `sup_d Σ_N E_d >= sup_d Σ_N Vroot_d^2`, so it implies only the sup-outside form. The max-inside form
  is strictly stronger and is not implied by (B). The data (at most 5.20 X^2) are fine either way.
- **Fix:** call `sup_d Σ_N Vroot_d^2` the uniform shell bound (implied by (B)) and `Σ_N sup_d Vroot_d^2` the shell
  maximal inequality (observed, not implied).

## Overstatements (should soften)

4. **Typing labels.** THM-4667 lines 7-10, 80, 88; note lines 22, 25, 124, 136. The d = 7-12 ladder is computed in
   float64, which AGENTS.md types VERIFIED, not FINITE-EXACT. "NUMERICAL" is not a repository label. The values are right
   to every printed digit.
5. **"Beyond Krasikov-Lagarias" drops its qualifier.** THM-4667 line 11; HYP-9288 lines 6-7; note lines 262-263. The
   conclusion is a limsup, so it holds only "along a sequence of scales", as HYP-9283(C) says. It is also conditional on
   the corrected constants from finding 1.
6. **Six colours: the larger-source counts include an unfinished run.** Golden note lines 35 ("five larger") and 144
   ("six larger"), table line 161, HYP-9289 line 10. contact2's additive pre-screen never produced a result; both .out
   files stop after building the edges. Four larger sources completed, all with 0 blocked patterns. contact2's edge
   counts are correct.
7. **Six colours: the v^2 result is stated too broadly.** Golden note lines 96 ("cannot force at any size") and 174
   ("closed, as v^2 now is"). Only the source with these coefficient sets (|A| = 54, |C| = 10) was tested, and for that
   source the conclusion is correct.

## Minor

8. THM-4667 line 94: the actual E_6/N octave means are 1.69138, 1.69902, 1.69322, 1.69333, 1.69271 against the ideal
   1.693330. That is agreement within 0.35%, not to three decimals.
9. THM-4667 line 103: at octave 2^10, d = 2, INTRA is 91.0% of the change; the range is 91-100%.
10. Note line 207: at d = 40 the participation ratio is 1.766 for N = 66463 (not 1.0-1.3).
11. Note line 194: the octave profiles already range over 1.39-1.54 at d = 15 (not "from about d = 25").
12. The band-tail code in root_transport_energy_ladder_banded_20261017.py (lines 51-59, 69):
    - It starts from the initial band distribution truncated at b < 200, and cuts the convolution at 1000 bands.
    - The written argument is rigorous, and correcting both truncations gives the same bound to 4 digits
      (2.296e-25 at d = 12).
    - Float rounding is not enclosed.
13. The "brute force on all seeds below 2·10^6" check for statement 8 has no committed script, and the sibling-anatomy
    docstring claims an on-the-fly pair check the code does not do. The claim is true: the referee checked every seed
    below 4·10^6 (a2).
14. THM-4667 lines 72-73 ("equals the actual E_d(N)/N to leading order") are asserted without proof.
    - This is the l2 analogue of clock_duality's R8.3.
    - a14 confirms the new matching bullet: for N = 2^16, 2^18 and d <= 6, actual and ideal coefficients differ only at
      ROOT, at 3^d, and at 15 or 91.
15. The three-entry incremental update (THM-4667 line 165) is clock_duality's equation (3) in
    marked_contacts_root_current_20261017e.md and should be credited.
16. The golden note (lines 6-7, 42, 140, 146) and HYP-9289 cite "G5" for the golden-contact classification.
    - G5 is the three-contact theorem for different contacts.
    - The right citation is G1 of marked_contacts_root_current_20261017e.md section 5, a proved lemma whose hypotheses
      are certified in joint_contact_field_20261017f.md section 5.
    - The mixed-edge rule used (Δa = dq, Δb = Δc = q, |q| = 1) is exactly G1's family with e = 1, so it rests on a proved
      result.
17. Golden note line 63: "cyclic-shift correspondence (DP) colouring" holds only when each forbidden set S_q is a
    singleton; in general it is a Z/6 colouring with forbidden difference sets.
18. Golden note lines 18-19 present the additive extension lemma as new. Additive (coordinate-sum) colouring CSPs already
    appear in clock_duality's graded_tail_carrier_20261017b.md (committed 14:10, lines 32 and 585), which should be cited.
19. Golden note lines 102 and 119: with |D(C) ∪ D(L)| = 130, the e and product censuses cover 1.436 x 10^9 differences,
    not 9.1 x 10^8.
20. root_transport_crossterm_anatomy_20261017.py prints a d = 1 diagonal of 4.856988 instead of E_0 = 4.857143, because
    it applies a periodic sum to a non-periodic quantity. The true diagonal equals E_(d-1) to 10^-12, so statement 6 is
    right.
21. THM-4667 line 26 lists files that do not exist: symbol_norm has no .out and fourier_check has no .json.

## Issues already repaired by clock_duality (458c4b89d0); repairs confirmed

- **R1:** the unrooted energy complement in statement 2. The repaired identity holds, and so do `R(N)^2 <= liminf E_d(N)`
  and `D_X <= liminf T_d(X)`. So the calibration no longer assumes universal convergence.
- **R2:** statement 10 had equated the unfolded phases with ordinary Syracuse offsets. That was false: for letters 1-8,
  only 2/8, 4/64 and 10/512 words agree, while the factorisation `C S^(d-1) A` matches 8/8, 64/64 and 512/512. The
  P_3 = 1 vs S_3 = 8 (mod 9) example and the 11 -> 17 -> 13 clock example check out.
- **R3:** HYP-9288 (C) originally used sibling covariance only, without INTER; the repair adds INTER.
- **R4:** the Rudin-Shapiro normalisation and the probabilistic scope of Doob's inequality are now stated correctly. The
  unproved "sub-Poissonian number variance" reading of Tao, previously under a PROVED heading, has been removed.

## Verified as correct

- **Ledger (a1, a2):**
  - Statements 1, 2, 3, 4, 8 and 9 hold exactly for 464 values of N (28,835 depth checks), including the INTRA/INTER
    split and the (s, 2s+1) pair shape.
  - The smallest scale ratio between different blocks of a fibre is 32.0013.
  - For every seed below 4·10^6:
    - the class partition and the residue claims for x ≡ 7 (mod 8) and x ≡ 11 (mod 16) hold;
    - every fibre matches formula (10);
    - the spacing is s' = 64s + 49;
    - the singleton (2u - 1)/3 appears exactly when u ≡ 2 (mod 3).
- **Incremental scan (a3, a10, a11, a12):**
  - The referee's own incremental code reproduces the committed octave traces exactly for X = 2^1-2^13, including
    final = D_X, the maximum at d = 0, and T_5(2) = 4.
  - The referee's own direct transport agrees at every depth for 770 values of N, including octaves 2^12, 2^13, 2^16
    and 2^17.
  - The session code reproduces every committed octave statistic for 2^1-2^21.
  - The global maximum over N < 2^22 is 10.062705 at N = 7910; N = 1 gives 10.0.
  - At N = 6904411 the direct transport gives sup 6.269811 at d = 172 and R = 6466, which matches H(4N) - H(N) from
    independent clock parities.
  - The ROOT shares, the sibling-anatomy table and S40's bound D_X <= 4.5367 X^2 all check out.
- **Ideal ladder (a4):**
  - An exact forward push gives γ_0 = 16/7, γ_1 = 155896/87381, the claimed γ_2, E_0 = 34/7 and E_1 = 32898/9709.
  - It confirms R8.1 for d <= 5 and R8.2 for d <= 6, and the diagonal part equals E_(d-1) to 10^-12.
  - E_0-E_12, γ_0-γ_12, the cross terms and all eleven ratios match to every printed digit.
- **Fourier form (a5):** the recursion matches the referee's code to 1.6e-13 relative. The branch map, μ_o, the sign and
  the push-forward are correct. The Tao caveat is scoped correctly.
- **Operator norm (a6):**
  - Balanced norms 1.283080, 1.323878, 1.335294, 1.346489; full space 1.447009.
  - Truncated operators converge from below, which confirms the Toeplitz reduction, the fold, and the balanced subspace.
  - Explicit real balanced inputs grow in energy by factors 1.646-1.813, so "not a contraction" follows.
- **Calibration chain:** THM-4658's descent applies with the same D_X, and constant budgets up to 1806.6 are correct.
- **Six colours (a8, a9, a13):**
  - From the receipt alone, the axis families are Cartesian, and the 570 patterns are exactly the canonical
    6-colourings of the port graph.
  - All 570 have an additive colouring with c_1 = c_2, found with a different encoding and checked on all 30,856 edges.
  - The pentagon ports span K_5, and 465 patterns fail the two-distance test (105 pass).
  - A geometry built from first principles reproduces the receipt's edge list exactly, matching the "axis + (dq, q, q)"
    rule.
  - All six recipe counts are reproduced.
  - The v^2 sieve is sound. An independent float census of 905,531,494 candidates finds exactly 22 units, all pure
    v^2-axis; the next nearest is 3.4e-9 off.
  - The e-direction census has 26 pure-axis units.
- **Credits.** R5, formula (10), the energy kernel G_X, P1-P4 and equation (5) are credited. The referee found no earlier
  explicit statement of the full fibre decomposition, the l2 energy ledger, or the ROOT/shell identity.

## Reproduction

All scripts run in under about 2.5 minutes, with peak memory about 3.5 GB (a4 at d = 12):

    cd 04-computation/experiments/s41_audit_20261017
    python a1_ledger_bruteforce.py
    python a2_fibre_classes.py 4000000
    python a3_direct_vs_incremental.py
    python a4_ideal_ladder.py 12 200
    python a5_fourier.py
    python a6_operator_norm.py
    python a7_calibration.py
    python a8_golden_receipt.py
    python a9_golden_numeric.py
    python a10_sibling_anatomy.py
    python a11_session_run_reproduce.py 21
    python a12_direct_large_N.py 2109637 6904411
    python a13_golden_recipes_and_e.py
    python a14_l2_matching_bonus.py

## Note added by the session after the report

THM-4668 and the product-direction censuses (ue, ve, uv, e^2) were committed after the audited state. They were not
audited here, except that a13 independently confirmed the 26 pure-axis units of the e-direction census.
