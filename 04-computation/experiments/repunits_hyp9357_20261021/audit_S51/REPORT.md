# Independent audit of the S51 packet (THM-4748, THM-4749, HYP-9359, note)

The auditor's report, saved verbatim by the session (audit subagents cannot write report files). Line numbers refer to the pre-audit files of commit 4712463a7b.

I found no numerical error anywhere in the packet. Every formula, count, coefficient and pre-registered verdict reproduces exactly with code written independently of the session's, including the strong inventory and canonical form. What needs fixing is two false sentences in a THM-4748 remark, a typing error in the THM-4749 title, overstated pre-registration claims in HYP-9359, and prior art and credits that were left out. I did not commit or modify anything; everything I wrote is in `audit_S51/`. The parent session added `fan_uniqueness.py` (commit 416b090462) while I was working; I checked its output too.

## ERRORS (must fix)

1. **THM-4748, line 124 (remark "Why the asymmetry").** "Subtractions come only from arcs sharing an endpoint" is false in general. The crossing formula subtracts `[q=0]2^p` and `[p=0]2^q`, and the nested formula subtracts `2^(x+y+1)`, for arcs that share no endpoint. The sentence is true only inside one-apex tournaments.
2. **THM-4748, line 124.** "Removing everything but 2^(k−1) and the final −2" is off by one. The surviving term is `2^(a+b−1) = 2^k`.
3. **THM-4749, line 3 (title).** "(all odd numbers >= 3 except 7 and 21; KNOWN, Session24)" is wrong. That value set is a conjecture, certified only below 80,265, and Session24 explicitly does not claim it.
4. **THM-4749, line 79 (proof of 2).** "When 2^k − 1 is prime ... w = 1" fails at k = 3: 7 is not a value, so `N(n,7) = 0`. The inequality still holds; only the sentence is wrong.
5. **HYP-9359, line 7, and PREREGISTRATION line 3.** "pre-registered ... before the data existed or were examined" and "Nothing in 1100 < d < 2100 has been looked at" are overstated (`c7_chronology.out`):
   - The order-9 slice of the test range had been examined. It is 95,847 of 587,031 test classes (16.3%). The scratch script `rep1.py` (18:54) correlated `N(9,h) = S_9(h)` over `h ≤ 2349` with 3|h, 5|h and primality, and the S50 census was already public.
   - The orders ≥ 10 data did not exist until 19:29:15.6. That is 9.6 s after the R1–R7 commit object, by local timestamps only. GitHub's push events only show 917fef912a on the remote by 19:35:48, which is after the data existed.
   - G1–G4 were written after the test-range S(d) and the R verdicts were known. Only the decomposition was blind: `repunit_rarity.py` was created 72 s after the G commit.
6. **HYP-9359, line 80.** The residuals "63 (−1.68), 127 (−0.55)" come from the full model refitted on [25, 1098], not from the [301, 1098] model in the table. No committed script prints them.
7. **Note, line 86.** "Labelled share ... 0.25, 0.80, 0.81, 0.94, 0.993" are densities relative to 1/4. The actual shares are 0.0625 up to 0.248.
8. **Note, lines 104–108 (the "partial mechanism").** Only orders 7–8 are reported. At order 6 the fa = 2 strong classes have a 7 (mod 8) share of 0.462 (6 of 13). At order 9, fa = 2 gives 0.082 and fa = 3 gives 0.332. There is no consistent depletion; this belongs in §8 as refuted.
9. **Note, line 141.** "Residues mix in one step" is not supported by any committed output: `refuted_leads.py` computes the transition matrix but never prints it. My check confirms it at orders 6→7 and 7→8 (distance from uniform 0.091 and 0.023), but only roughly at 5→6 (0.234, with 6 parents).

## FIXES (retyping and wording)

**THM-4748**
- Statement 6 (lines 62–64) says the 7 (mod 8) values are "no longer on finitely many lines". That is an all-m claim backed only by m ≤ 13 data. It can be proved:
  - Take a crossing pair with p = q = 0 disjoint from a C_3, so `H = 3(2^(r+2) + 2·3^r − 1)` with r odd.
  - This is ≡ 7 (mod 8), and the odd part of H + 1 is at least `3^(r+1)/(4(r+1))` for r ≥ 5.
  - Add this proof or retype the clause.
- Line 70: the "c + 1 admissible runs" rule needs "excluding endpoints of the used arcs".
- Line 40: touching intervals give one strong tournament, a one-vertex chain gluing, not "two strong blocks".
- Line 47: 23 = 3·2³ − 1 is already on the first line, so listing it separately is redundant.
- Cite THM-1370 by its slug (`THM-1370-h-spectrum-omits-7-21-all-n`); the ID is used by two files.

**THM-4749**
- Line 7: only statements 2–4 are PROVED; statement 5 is FINITE-EXACT, as its own heading says.
- Statement 3 (lines 43–44): w(63) = 1 and `N(n,63) = 2(n−7)` are KNOWN from Session24 PROOF.md §3 and the §5 table.
- Proof of 4 (line 87): the restriction to i ≤ 4 is unnecessary. `F_i = H(NT_(2^i+2))` is a strong value for every i, so `w(2^(2^j) − 1) ≥ j` holds unconditionally.
- Line 62: "conditional for k ≥ 17" is over-conservative. The w values at k = 18, 20, 21, 22, 24, 25, 26, 30, 36 are certified by factors below 80,265.

**HYP-9359**
- State that 0.237 → 0.139 and 0.217 → 0.135 use a different model (no prime or blocked terms) from the table's 0.239 and 0.219.
- After both removals the test β₃ is still +0.010 ± 0.005, about 2 SE. "Only" and "entirely" should become "almost". "Primality neutral throughout" is untested in the decomposition fits, which have no prime term.
- The prereg baseline "+0.049 ± 0.016" does not reproduce (I get +0.046 ± 0.012).
- The prereg says "trivial Aut" where G2 and the script use 3 ∤ |Aut|. G2 and G3 pass under both readings (`c5b`).

**Note**
- Line 25: k = 12 (4095) is not a Fermat product.
- Line 71: cite the proof of unbounded odd parts.
- HYP-1754 and HYP-1762 have no standalone files; they live only in `INDEX-HISTORICAL-THROUGH-2026-07-21.md`, and HYP-1762 is a collided ID.

**Credits**
- mac-mini's MSG-2400, sent 19:50, before THM-4749's first commit at 20:10, contains item 8 (7 | 2^k − 1 iff 3 | k, never 2^k + 1) and item 2 (the Fermat monoid). THM-4749 and the note credit neither.
- THM-4748 statements 2–3 (19:48) precede MSG-2400; the fan (19:56) follows it by 6 minutes. The existing credit wording is fair.

## PRIOR ART

**Literature**
- **Moon 1972** (Publ. Math. Debrecen 19, 101–104, read in full by my literature sub-agent): the one-arc count `1 + 2^(n−2)` and the touching-chain product. These are statement 5 and the disjoint and touching cases of statement 1.
- **bof, MathOverflow 234029 (2016-03-19)**, which I read: reversing v1v3 and v2vn gives `3·2^(n−4) + 3`. This is the crossing formula at p = r = 0, q = n − 4 (I checked it for n = 5..16). The sub-agent reports Math.SE 1688155, which gives the converse case.
- **El Sahili and Abi Aad**, "About the number of directed paths in tournaments", DAM 277 (2020) 115–126, doi 10.1016/j.dam.2019.09.005. The 2020 paper is joint with El Sahili, which THM-4748 omits. Its abstract says it computes exact H for tournaments made from the nearly transitive tournament by reversing extra arcs. That plausibly covers the nested and shared-endpoint cases. It also studies tournaments whose H is a power of 3.
- **Abi Aad**, "About the number of directed Hamiltonian paths in special tournaments", DAM 394 (2026) 220–228, doi 10.1016/j.dam.2026.07.017, SSRN 5804403. The repo's Session16 packet (`feedback_word_hamilton_spectrum_20261020/README.md`) reports that its Proposition 1 overlaps the one-apex formula. THM-4748 should carry that pointer.
- Neither DAM paper's full text was reachable (ScienceDirect CAPTCHA, SSRN 403). Someone should open both in a normal browser.
- **bof and Royle, MathOverflow 232751 (2016, edited 2022)**: origin of the "7 and 21 only" conjecture; all odd values below 80,557 except 7 and 21 occur on at most 11 vertices. THM-4749 calls this "Chan–Pak"; it should credit the 2016 origin.
- **Bar-Noy and Naor**, SIAM J. Discrete Math. 3 (1990) 7–20: H(T) equals the number of minimal feedback arc sets. This is useful context, not prior art for the formulas.
- Nothing in these sources gives the general crossing formula, the fa ≤ 2 mod-8 lines, the three-arc Mersenne theorem, the Fibonacci fan, or the N(n,h) degree law:
  - El Sahili and Ghazo Hanna (arXiv 2101.00713 and 2207.11510), read in full, and the papers citing them;
  - Moon's *Topics on Tournaments*;
  - Busch 2006.

**Repo (not cited in the packet)**
- THM-250 (kind-pasteur, 2026-03-17): the one-arc formula with the same prefix/suffix proof.
- `05-knowledge/results/deep_extension_s20al.out` §2 (2026-03-22), "TWO-TILE FLIPS: H FOR TWO REVERSED ARCS": all 45 two-arc values at n = 6, with interactions by geometry but no closed form.
- `arithmetic_seams_20260921_tournaments.md` §5, eq. (12): `H(E_n) = F_r` at n = 2^r + 2.
- Session24 §3 and §5 for 63 (above).
- `single-core-cycle-count.md` already lists the THM-344 H = 63 signatures 1001100 and 1100110. These are the k = 6 fan and its converse, so the k = 6 fan was already in the repo.
- The THM-4749 degree law, its citation of Session24, the 7·3^a loss and the Φ₃(2) reading of 7 and 21 are correctly sourced.

## CONFIRMED

All checks below come from my own helper library (`alib51.py`, self-tested in `c0`).
- **THM-4748 (1):** all 18,564 two-arc configurations with spans ≥ 2 and m ≤ 16 match the DP, boundary cases included.
- **THM-4748 (2):** every tournament with fa ≤ 2 on m ≤ 16 vertices (every reversal of at most two arcs, any spans) has 7 (mod 8) values only on the listed set. Every line value occurs, first at m = j + 3 and j + 4. Each step of the written mod-8 case analysis checks.
- **Minimum orderings:** no minimum ordering has a span-1 back arc (2,251 orderings of 1,044 classes, n ≤ 8).
- **THM-4748 (4):** for k ≤ 18 the fan is strong with 3 back arcs and fa = 3 for k ≥ 5, checked by DP. H = 2^k − 1 by DP for k ≤ 18 and by a path-decomposition count independent of OCF up to k = 3000. The one-apex formula matches the DP for m ≤ 13 and is off by exactly +1 if 1 ∈ A. The fan is the unique 3-in-arc one-apex solution for m ≤ 18.
- **THM-4748 (6):** the three-arc census is reproduced exactly, and k = 12, 13, 14 are not reached.
- **Class census** n ≤ 9 (191,536 classes; 178,133 strong): least fa is 3 for 31, 63, 127, 255 and 1 for 33, 65.
- **Strong inventory to 2100:** my own ear census and canonical form give 723,544 classes with identical per-order counts. The H multiset per order equals the session's `inv_2100.pkl` (read with a restricted unpickler). The counts below 625 are 13,439 / 8,769 / 3,536 / 658, and every odd 3 ≤ d < 2100 except 7, 21 occurs.
- **Session24 formula:** reproduced for every odd h < 2100 and n ≤ 30, including degree, onset and leading coefficient. S(11) = 1, S(23) = 4, S(31) = 4, S(41) = 6, S(89) = 23, S(683) = 258, giving 184, 516, 24 and 18. The w table and the n = 8, 10, 14 count tables match.
- **HYP-9359:** the regression coefficients and standard errors match to 0.0005. The decomposition and cover flags match, with the cover flags also cross-checked by a slow re-implementation. R1–R7 and G1–G4 all PASS, and R5 and R6 are robust to the choice of model.
- **HYP-9357:** H(R_27) and H(R_31) equal mac-mini's committed `04-computation/experiments/paley_extremal_20261021/R27.out` and `R31.out` (commit 269c9ec1f1).

## VERDICT

THM-4748 statements 1–5 may stay PROVED. The proofs are correct and complete as written, and they are confirmed far beyond the session's range. Statement 1 should state that priority is unsettled (El Sahili–Abi Aad 2020, Abi Aad 2026, bof 2016). The remark errors must be corrected and the last clause of statement 6 proved or retyped. In THM-4749, statement 1 is a faithful citation of the reviewed Session24 packet. Statements 2–4 are PROVED given the certified value set below 80,265, except that the 63 facts are KNOWN from Session24. Statements 5–6 are verified FINITE-EXACT. The title's "KNOWN" typing of the value set must be removed and the credits added. HYP-9359's numerical evidence is sound: everything reproduces exactly from an independent inventory. Its pre-registration claims, however, must be downgraded. The R1–R7 timing is supported only by local timestamps with a 9.6 s margin, the order-9 slice was not blind, and G1–G4 were blind only to the decomposition.
