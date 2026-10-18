# Audit F3 — merge_18_19_11_22 (second brief), THM-4690 "Why the h ≡ 7 primes", HYP-9325 consequence 4

Independent referee, 2026-10-18. Audited HEAD `ebb9216f49`; the session's commits `7e07b63f27` (q = 31 Hoffman margin,
Weil justification) and `ebb9216f49` (depth 1–6 clock-graded table) were checked too. Nothing was committed or pushed,
no existing file was edited, and only the referee's own code and PARI were run. The subagent could not write this
file, so the session saved it verbatim from the referee's final message.

## Verdicts

- **(a) §2.1, fixed places of F: CORRECT WITH FIXES (minor).**
  - PARI on the absolute degree-48 field, with an explicit automorphism c applied to `idealprimedec` for every
    p <= 1300, and the residue action checked.
  - c-fixed primes occur only above 29, 149, 239, 359, 479, 569, 659, 809, 1019, 1229 and 1289. All have e = 1,
    f = 2 or 6, and c acts as Frob^(f/2).
  - No prime above 2, 3, 5, 11 or 37 is fixed.
  - At the 29-place, Lemma R's hypotheses hold: residue F_841, c = Frobenius, and the prime is inert over F⁺ because it
    is unramified and c-fixed.
  - However, each listed p has **8** fixed primes; `fixed_places_golden_field.py` prints "2".
- **(b) §2.2/2.2b, κ: CORRECT WITH FIXES (2 MAJOR).**
  - The four colouring certificates are PROPER. They were checked with μ_(q+1) built from a primitive element, not
    from the norm equation.
  - Hoffman: 4.0148 (q = 23), 4.3082 (29), 4.0000439 (31; interval arithmetic gives 3λ_min + 32 = 4.685e-4 > 0).
  - The two-point LP equals Hoffman (195.2069 at q = 29).
  - λ = −Kl(N(a)) holds to 1e-14 for q = 7…101. Weil then gives |λ| <= 2√q, the floor exceeds 6 exactly for
    q > 97.99, and it is 8.76 at q = 239.
- **(c) §2.3 and the THM-4690 remark: CORRECT WITH FIXES (minor).**
  - In K' = Q(ζ5, √-11), all 4 primes above 19 are c-fixed with c = Frobenius. b and c have valuation −1 at all of
    them, and b/c is a unit.
  - In F' there are 12 primes above 19. Exactly the 4 with w ≡ 4 (h ≡ 7) are fixed, with c = Frob. At these,
    v(u) = v(u2) = 0.
  - In F, 24 primes lie above 19 and none is fixed. √-3 is not in F' (nor in K').
  - G2 of `rhombus_overlap_elliptic_20261011.md` is described correctly.
- **(d) §3.1: CORRECT WITH FIXES (minor).**
  - The denominators check out. The pattern also holds at d = 3: S41 a4's exact E_3 has denominator (2^162−1)/3^5.
  - The lag constants were derived exactly by closed form: M_1 = −5/7, M_5 = M_7 = 1/112, M_6 = −17/448,
    M_11 = 1/1792, M_12 = −17/28672, M_13 = 1/28672.
  - The next ones are M_17 = −5/458752, M_18 = 17/917504, M_19 = −5/1835008, and M_(D+18) = 2^-18 M_D for D >= 1.
  - The total reproduces E_1 = 32898/9709 exactly. The l1 mass 155896/87381 and E_1 also come out exactly from the
    period-18 band coefficients.
- **(e) §3.1b: CORRECT (minor wording).** The referee's own forward code at N = 2^20 reproduces every entry of the twin
  table exactly.
  - Twin pairs 35002 / 96656 / 201843 / 386023; shares 0.1204 / 0.2379 / 0.3189 / 0.3295; all the sums.
  - The ideal depth-2 share is 0.120306, from residues mod 27.
  - The exact depth-2 meet sums equal (E_2 − E_1)/2 to 7.6e-26, with lag 0 = −0.0467312620. E_2 was recomputed
    exactly.
  - "Children up to scale 6" is right: a cutoff at 4 changes the total to −0.2123; a cutoff at 12 changes nothing.
  - An independent clock-graded recursion reproduces the new d <= 6 table and E_3…E_6 to every printed digit.
- **(f) §3.2 / HYP-9325 consequence 4: CORRECT WITH FIXES (minor).**
  - Δ(N) = R(N); checked R(N) = Σ c_0 t = H(4N) − H(N) for N <= 1500 against brute-force parities.
  - THM-4658's normalized descent is quoted exactly.
  - The chains 4^j·2^40 and 4^j·2^41 cover every dyadic X >= 2^40, and the fixed-point argument gives
    e_X <= 2445.3 + o(1).
  - HYP-9283(B) at N = 2^40 certifies a <= 14348.
  - Calibration thresholds: 15.26, 25.51, 129.26.
- **(g) Prior art and typing.** Missed prior art: decalion89 (pinned commit cb660f3), THM-4628,
  `ei21_kappa_q_20261012.md`, and the oai2 note §4.4. PROVED/FINITE-EXACT typing is otherwise sound. The NUMEROLOGY
  labels are fine, but two "explains why" phrases overreach (fixes 7 and 8).

## Required fixes

All refer to `05-knowledge/results/merge_18_19_11_22_collatz_hn_20261018.md` unless another file is named.

1. **MAJOR — prior art on κ(29).**
   - decalion89 `notes/local_colourings.md` §4/§14 certifies **κ(29) >= 6** with a three-point bound,
     α(G_29) <= 163.25 < 841/5. It certifies κ(37), κ(41) >= 6 the same way.
   - The repo already records "κ(29), κ(41) >= 6 (decalion89, three-point bound)" in
     `oai2_openai_math_second_reading_20261006.md` §4.4 (and audit A).
   - Change "Lower bound (PROVED). `κ(29) >= 5`" to "κ(29) >= 6 (CITED, decalion89); Hoffman alone gives 5".
   - Change `κ(29) ∈ {5, 6, 7}` to **{6, 7}**.
2. **MAJOR — prior art in the κ table.**
   - decalion89 already has κ(23), κ(31) ∈ [5, 8], with 8-colourings "by intervals" and κ(23) "likely >= 7".
   - So "proven here <= 8" at q = 23 and 31 reproduces KNOWN results. Only κ(29) <= 7 appears new.
   - The "known (THM-4558)" column needs decalion89 added.
3. **MAJOR — wrong row and weak evidence.**
   - Row "37, 41 | 9 colours not reached within the budget" is wrong. Explicit line-projection colourings, verified by
     the referee, give **κ(37) <= 8**, κ(41) <= 9 and κ(43) <= 8.
   - Both TabuCol implementations stall at 37–79 conflicts with 9 colours on G_37. Tabu stalls are therefore weak
     evidence.
   - q = 31, "so 8" → "7 or 8". The heuristic independence number is about 146 >= 138, so 7 colours are not excluded.
   - For κ(29) = 7, add stronger evidence: independent sets of size 127 are found, 128 never in 7 runs of up to 3·10^7
     iterations, while a 6-colouring needs a class of 141.
   - Controls for that search: it finds α(G_13) = 36, decalion89's exact value, and α(G_17) of about 56–57 < 57.8.
   - "needs α <= 140 or UNSAT" → "could come from α <= 140 (three-point gives 163.25) or UNSAT".
4. **Cite the Hoffman–Weil source.** THM-4628(b),(c) already contain both the eigenvalue = −Kloosterman identity and
   "floor exceeds 6 for q >= 101", as does decalion89 Prop. B.
5. **The third condition of §2.1 is vacuous**, since a cubic over F_p has a root or is irreducible.
   - Reword as two conditions; the factorisation of the cubic only decides whether the residue field is F_(p²) or
     F_(p^6).
   - State that there are 8 fixed primes per such p.
   - Fix the counts printed by `fixed_places_golden_field.py` and commit a `.out` for it.
6. **Scope the claim.** "no residue colouring caps G(F) below 7" → "no single-prime (Lemma R) residue colouring".
   Colourings through products of the eight 29-places are not analysed.
7. **"in G(F′) the copies lie in different components" is ill-typed.**
   - The copies are not subsets of F′: M and C contain η = 1 + ζ3.
   - Correct wording: they lie in different cosets of Γ(F′), the group generated by the unit vectors of F′.
   - THM-4690 remark, "So statement 5's separation is THM-4558's Lemma R.1 in F′": statement 5 covers every step
     integral at the h ≡ 7 primes, including ±ζ15^l, which are not in F′. Lemma R.1 explains why u and u_2 are
     integral there and extends the separation to steps in μ_30·U(F′).
   - The same wording appears in `results/INDEX.md`.
8. **Why 19 is in the denominators.**
   - 19 enters b and c through N(h ± j) = −76 = −4·19, i.e. Res(p, x²+11) = 2^4·19^2.
   - The c-fixedness of the 19-places (19 ≡ −1 mod 5, (−11/19) = −1) is a separate fact. Reword §1 row 19 ("§2.3
     explains why") and the end of §2.3.
9. **The lag table and the 19 sentence.**
   - −99820/67963 is unreduced; it equals −14260/9709.
   - "is new": the finite-N split with 1/112 and −17/448 is already in THM-4690's remark. What is new is the
     ideal-model evaluation.
   - Replace "about 1e-5" with the exact values above.
   - Replace "lags >= 17 … 4^9 = 2^18 apart" with: every M_D has denominator 7·2^k. Lags 1–13 sum to 2Σ = −47/32
     exactly. The tail sums to 3/310688 = 3/(2^5·9709), so 19 and 73 enter only through the geometric tail.
10. **Mechanism wording.** The period proof (R8) also needs the three-lift cancellation, not only the order of B. The
    pattern can be stated FINITE-EXACT for d <= 3.
11. **§3.1b wording.**
    - Say the shares count only sources with c_0 ≠ 0. Over all seeds they are 0.107, 0.216, 0.292, 0.302.
    - Label the per-pair model ("+E[c_0²] / −E[c_0²]/2, cancel exactly") HEURISTIC.
    - "HYP-9325(B1) … is a statement about twin residues" is a reading. (B1) concerns first arrivals into V_1000,
      including same-clock pairs of unequal depth.
12. **The clock-graded subsection and §4.**
    - "Depths 1-5" → 1-6; "mod 3^(e+1)" → 3^(d+1).
    - Lag 3 is positive only for d >= 2.
    - The explanation for d = 7–12 extrapolates beyond the d <= 6 split; label it HEURISTIC.
    - §4 "Not proved … lag constants beyond lag 1" and "Next (Collatz)" are now stale.
13. **Wrong maximum.** "observed maxima are about 1.4 (V_1000 from 2^14 on)" is wrong. HYP-9325's own table gives
    1.790 at 2^14, and at most 1.557 from 2^15 to 2^27.
14. **§3.2 and HYP-9325 consequence 4 — two hidden assumptions.**
    - State the off-basin convention: E_V, h_t and G are taken over sources that reach V, with 1 ∈ V (Lemma F's
      convention). Otherwise the hypothesis presupposes convergence and the B(X) conclusion is vacuous.
    - HYP-9283(B) needs the budget at every dyadic X >= 2^40, so the O(log N) and o(1) terms must be made explicit. The
      margin 14348/2445 ≈ 5.9 makes any explicit bound sufficient.
15. **Receipts.**
    - `tabucol_runs.out` lacks maxiter and seeds, so "2·10^7 iterations" appears in no receipt.
    - `nineteen_place_subfield.py` tests `c == frob(t) or c == t`; the "or" branch should go.
    - §5 Files omits about a dozen used files, e.g. `verify_colouring.py`, `twin_residues.py`,
      `depth2_twins_ideal.py`, `hoffman_table.out` and `clock_graded_ladder.py`.

## Optional improvements

- State χ(G(F)) <= 7 and χ(G(F′)) <= 5 explicitly as PROVED.
- Upgrade the depth-1 lag constants from "recognized" to "derived" (closed form in `f3_band_exact.py`).
- THM-4690's remark "lags 11-13 give about 10^-3 in total": in the ideal model they sum to exactly 0.
- Cite `ei21_kappa_q_20261012.md` (θ(G_q) table) and THM-4663.

## Reproduction

Run from `audit_F3/`:

```
gp -q f3_fixed_places.gp > f3_fixed_places.out      # ~10 min
gp -q f3_nineteen.gp > f3_nineteen.out
python3 f3_kappa_check.py; python3 f3_lp.py 19 23 29 31 37 41
python3 f3_projection.py; python3 f3_projection_certs.py 19 37 41 43
cc -O2 -o /tmp/f3_tabu f3_tabu.c; /tmp/f3_tabu 37 9 20000000 1; /tmp/f3_tabu 29 6 30000000 7
cc -O2 -o /tmp/f3_mis f3_mis.c; /tmp/f3_mis 29 127 10000000 2; /tmp/f3_mis 29 128 30000000 11
python3 f3_band_exact.py 90; python3 f3_clock_graded.py 6 64
python3 f3_twins.py 20 5; python3 f3_twins_ideal.py 60; python3 f3_root_invariant.py 1500
```

Files: `f3_*.py` / `.gp` / `.c` with matching `.out` files, and `kappa_{19,37,41,43}_colouring_k*_projection.txt`
(new certificates).
