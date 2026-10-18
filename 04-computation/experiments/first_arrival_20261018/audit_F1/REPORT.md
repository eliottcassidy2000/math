# Audit F1: THM-4690 statements 1-4, HYP-9325, and section 1 of the first-arrival note

Referee: audit F1 (independent), 2026-10-18. Worktree `math-wt-sixcolour-20261018`, branch
`mac-mini/session-20261018-sixcolour` at `044e5bdf47`. The subagent could not write this file itself, so the session
saved it verbatim from the referee's final message.

Scope:
- `01-canon/theorems/THM-4690-...md`, statements 1-4 and their proofs. Statement 5, the 19-adic obstruction, is audited
  elsewhere (audit F2) and is not covered here.
- `05-knowledge/hypotheses/HYP-9325-...md`: every number, the pre-registered items, and the "Consequences".
- Section 1 of `05-knowledge/results/first_arrival_one_tick_coefficient_copies_gibonacci_20261018.md`, together with
  the rows of its "Results at a glance" table that concern section 1.

Method:
- All checks use new code written for this audit (C and Python), in this folder.
- The session's code (`fap.c`, `fap_lean.c`, `fap_levels.c`, `source_ledger_constants.py`) was read to understand its
  output format, but was not run.
- My implementation is different from theirs:
  - stopped paths are memoised;
  - everything is recomputed from scratch at each sampled N, with exact integer arithmetic (`__int128`);
  - the ladder is checked by direct stopped simulation, not by level-by-level propagation;
  - E_0 and SIB are evaluated by exact interval counting.
- Nothing was installed or downloaded. Every job used under 0.5 GB of memory and took under 5 s. No run went past
  N < 2^22, except one ladder-free point check at N = 6,904,411.

## 1. Verdicts

| Item | Verdict | One-line reason |
|---|---|---|
| THM-4690(1) backward ladder, merge identity, basin form | **CORRECT WITH FIXES** | All identities hold exactly at 30 values of N on 4 stopping sets. As written, the merge identity at k = 0 sums over all of P^-1(y). It must sum only over preimages with dep = k+1, because V is contained in P^-1(V). |
| THM-4690(2) ports and preimage description | **CORRECT** | E^(1) is the port-diagonal energy, exactly. The fibre description matches brute force: 1,321,607 preimages, 0 violations, including y = 1 and seeds in V. |
| THM-4690(3) integer clocks and clock spectrum | **CORRECT WITH FIXES** | The identities hold exactly, and σ = (-1)^t holds for every source. The Setting's claim "literal clock = difference of shortcut stopping times" is false at the ROOT. |
| THM-4690(4) source pair ledger | **CORRECT WITH FIXES** | 34/7, 5/7, \|c_0\| <= 3 and the O(log N) errors are right; I re-derived them and confirmed by exact evaluation. The sibling case analysis has (at most one per level) exceptions that the proof does not mention. "At most two sibling pairs per seed" should be one. SIB_V and the V = {1} sentence are slightly imprecise. |
| THM-4690 status block | **NEEDS FIX** | It says the identities of 1-4 "hold at every N < 2^27". The session's codes check only the invariant and E_0 at every N; statements 1 and 3 were checked only at a few N. |
| HYP-9325 numbers | **CORRECT WITH FIXES** | All 221 table entries and every pre-registered outcome match the receipts. Two sentences are false: "within 10^-5 from k = 20 on" and "E+/N ≈ 4.1-5.4". The "about 1.2 at θ = π" phrase is imprecise. |
| HYP-9325 "Consequences" (PROVED implications) | **CORRECT** (add range) | The algebra and constants are right, and there is no hidden sign assumption. The bounds hold only for N >= 2^13 (or 2^12), the ranges where (A)-(C) are asserted. |
| Note, section 1 | **CORRECT WITH FIXES** | "\|E_0 - 34N/7\| <= 11.2 ... up to N = 4·10^6" is false: the true maximum below 4·10^6 is 53.1. "Exact" needs its error terms. The clock wording is the same issue as in the THM. |
| Pre-registration timeline (Q, T, L) | **VERIFIED** | Each JSON was committed, and never modified, before any commit containing data from the octaves it predicts. The L1-L7 prediction was pushed at 09:49:53, and the 2^27 receipt finished at 09:53:31 (file mtime, after an 85 s run). |
| My own scoring of Q1-Q5, T1-T4, L1-L7 | **AGREES WITH HYP-9325** | Q3 PASS; Q1, Q2, Q4, Q5 FAIL; T1-T4 PASS; L1-L7 PASS; the L theorem check PASSES. |

There is no MAJOR error, and none of the statements is wrong in substance. Every required fix concerns wording,
precision or a numerical claim.

## 2. What was verified independently

### 2a. From-scratch exact checks at fixed N (`aud_point.c`)

The checks were run on these stopping sets and values of N:
- **V_1000** at N = 1000, 7910, 16384, 65536, 66463, 262144, 658023, 945589, 1048576, 1215227, 1483383. The four values
  658023-1483383 are random with seed 20261018.
- **V_100** at 9 of these N.
- **H_1000** at 4 values of N.
- **V = {1}** (the ROOT) at N = 10, 100, 1000, 7910, 66463.
- **V_1000 at N = 6,904,411**, without the ladder.

Total failures: **0**. What was checked at every N:

- **V itself.** V_100 and V_1000 are forward closed (P(V) ⊆ V). H_1000 is not: 63 of its elements leave it, which is
  expected.
- **Source coefficients.** c_0 is computed from the definition. It agrees with the closed form
  `P(x mod 3, L mod 3) + [ends]·λ((x+L) mod 3)` at every seed (0 mismatches). Moreover **c_0 ∈ {-3, -1, 0, 2, 3}**:
  the values -2 and 1 never occur.
- **Invariant.** R(N) = Σ_u f(u) c_V(u) equals H(4N) - H(N), with f computed from first-hit stopping times by direct
  iteration of T.
- **Signs.** σ(x → u) = (-1)^(t(x)) for every source, with t = Σ(a + 2·strips).
- **Statement 3, exactly.**
  - E_V = G(0) + 2Σ_(D>=1) G(D), with all lags computed directly. HL = 2Σ_(D>=2) G(D) is also computed directly, not by
    subtraction.
  - S(π) = E_V and S(0) = E+, as exact integers.
  - The mean of S over 4T equispaced angles equals G(0) (relative error < 1e-12).
  - The cosine series agrees with direct evaluation (maximum relative error 1e-14).
  - S(kπ/12)/N at N = 7910 and 6,904,411 matches the session's spectrum receipt to all printed digits.
- **Statement 4.** At every N:
  - siblings with neither member in V have a common entry, clock lag exactly +1 and opposite signs (0 failures);
  - the pair ledger `E_V = E_0 - 2SIB_V + (G0-E0) + 2Γ_1 + HL` holds exactly;
  - for V ∌ 1 or V ∋ 3, Γ_1 = G(1) + SIB_V equals the directly computed sum of lag-1 pairs other than sibling sources;
  - |SIB - SIB_V| <= 9·#(pairs meeting V). On V_1000 the difference is 404-659, against 9|V| = 5643.
- **Statements 1 and 2.** These used a node map over all source paths, up to 3.4 million nodes.
  - E^(k) from the formula, at every k;
  - E^(0) = E_V and E^(maxdep) = E_0;
  - E^(1) = Σ_V c_0^2 + Σ_ports J^2;
  - the J recursion;
  - the merge identity at **every** k = 0 .. maxdep-1, restricted at k = 0 to preimages with dep = k+1 (0 mismatches);
  - a direct stopped simulation (stop at the first path point in V^(k)) at up to 12 levels;
  - dep(P(y)) = dep(y) - 1, which makes V^(k) forward closed;
  - the basin form J(y) = f(y)R_y at every node (0 violations).
- **Agreement with the session's ladder receipts.** At the 6 values of N they share with `fap_levels_V1000*.out`, all
  **224 E^(k)/N values agree digit for digit** (`compare_ladder.out`).

### 2b. Every-N scan to N < 2^22 (`aud_scan.c`)

- This is my own incremental implementation, run on V_1000, V_100 and H_1000.
- The invariant R(N) = H(4N) - H(N) holds at **every N < 2^22** on all three sets.
- The octave rows 2^8..2^21 agree digit for digit with the receipts:
  - `ledger_table.out` (all three sections);
  - `fap_lean_V1000_k27.out`, including the E0 and SIB ranges;
  - the cross-port columns of `fap_V1000_k26_full2.out`, `fap_V100_k25_full.out` and `fap_H1000_k25.out`.
  - In total **966 numbers, 0 mismatches** (`compare_receipts.out`).
- So `fap.c` and `fap_lean.c` are reproduced independently over this range. On the shared octaves the two receipts are
  also identical to each other.
- Over all N < 2^22:
  - max |E_0 - 34N/7| = **53.14** (N = 2,721,349);
  - max |SIB - 5N/7| = **31.57** (N = 3,682,499);
  - for N >= 4, max |E_0 - 34N/7| / log2 N = 3.13;
  - for 2^20 <= N < 2^22: |E0/N - 34/7| <= 4.03e-5, |SIB/N - 5/7| <= 2.05e-5, and |SIB_V/N - 5/7| <= 4.14e-4 on V_1000
    (1.9e-4 on H_1000, 8.0e-5 on V_100).

### 2c. Statement 2 and the sibling claim by brute force (`aud_fibre.c`)

- **Siblings.** For all x ≡ 1 (mod 8) with x <= 10^8 (12.5 million pairs), there are 0 failures of the following:
  - P(x) = P(2x+1);
  - a(x) = 2 and a(2x+1) = 1;
  - step clocks 2+2r and 3+2r, where r = strips((3x+1)/4);
  - 2x+1 ≡ 3 (mod 16).
- **Preimages.** I took every seed s <= 2^28 with P(s) <= 10^6, which gives 1,321,607 preimages.
  - Each one matches the formula, including the residue rule B^j(y) ≡ y + j (mod 3) and the step clock.
  - For each of the 375,000 seeds y, the number of predicted preimages equals the brute-force count.
  - Pair and singleton structure: pairs (s_j, 2s_j+1) have clocks (2j+2, 2j+3), and the only unpaired element is
    (2y-1)/3, with clock 1, exactly when y ≡ 2 (mod 3). There are 0 violations.
  - Edge case y = 1: the fibre is 1, 3, 113, 227, 7281, 14563, ..., with clock pairs (2,3), (8,9), (14,15). It contains
    1 itself, since P(1) = 1.

### 2d. Statement 4 constants by exact interval counting (`aud_constants.py`)

- **Method.** This is a different route from `source_ledger_constants.py`:
  - c_0 is constant on the integer intervals between the thresholds where B^j(x) crosses 2N and 4N, and likewise for
    2x+1;
  - on each interval I count seeds by residue mod 24 exactly;
  - this gives E_0(N) and SIB(N) exactly for N as large as 10^100.
- It agrees with brute force at 405 values of N (1..399 and 6 larger).
- E_0 - 34N/7 and SIB - 5N/7 stay O(log N) and do not tend to 0:
  - largest |values| in the table of N = 10^k, 2^m, 10^k + 12345, up to 10^100: 55.7 and 34.3;
  - over 3000 random N <= 10^29: max |E_0 - 34N/7| = 90.0 and max |SIB - 5N/7| = 41.9, which is about 0.9 log2 N.
- My own hand re-derivation agrees with THM-4690(4). Seeds have density 3/8, with x mod 3 equidistributed. Averages of
  c_0^2 are 2 (chain ends in [2N, 4N]) and 0, 6, 6 (L ≡ 0, 1, 2, not ending):

      E_0/N -> (3/8)[2·2·(4/3) + 4·6·(20/63)] = 34/7

  Heads x ≡ 1 (mod 8) have density 1/8, and the product averages are 3 for 3 ∤ L:

      SIB/N -> (1/8)·3·6·(20/63) = 5/7

### 2e. Exceptions to the sibling case analysis (`aud_sibling_exceptions.py`)

Write y = 3x+1, so y ≡ 4 (mod 24); write Y = 6N+1; and let L >= 1 be the level of x. Then x is at level L when
4^(L-1) y < Y <= 4^L y, and its chain ends in [2N, 4N] iff 4^L y <= 2Y-1. The sibling 2x+1 has y' = 2y+2.

From these one derives:
- **(A) fails** (x ends, but the sibling is not "same level, not ending") **iff y = ⌈Y/(2·4^(L-1))⌉ - 1.**
- **(B) fails** (x does not end, but the sibling is not "level L-1, ending") **iff y = ⌈Y/4^(L-1)⌉ - 1.**

So each case has at most one exception per level. An example is N = 6, x = 1: x ends at level 2, but its sibling 3
ends at level 1. A brute force over all N <= 2500 found 1478 exceptions, at most 1 per (level, case), and 0 mismatches
with the criterion. The proof's "O(1) per level" absorbs them, but the text states (A) and (B) as exact.

### 2f. ROOT entries and the literal clock

- t(x) ≠ st(x) - st(u(x)) occurs **only** for sources entering u = 1, and always by exactly 2. Here st is the first-hit
  stopping time.
- Counts on V_1000: 8 sources (N = 7910), 55 (N = 2^16), 831 (N = 2^20), 5547 (N = 6,904,411). For V = {1} it is every
  source except 1 itself.
- **Cause.** st(1) = 0, but st(B^j(1)) = 2j + 2: for example st(5) = 4, not st(1) + 2. So a strip into the ROOT changes
  the stopping time by 2j + 2, not 2j.
- The parity is unchanged, so σ = (-1)^t and every identity stays valid. G(D) is defined with the literal clock, as in
  the code.

## 3. Required fixes

1. **THM-4690, status block, first FINITE-EXACT bullet.**
   - Replace "the identities of 1-4 hold at every N < 2^27 against the invariant R(N) = H(4N) - H(N) (fap.c,
     fap_levels.c, fap_lean.c)" with: "the invariant R(N) = H(4N) - H(N) and E_0 are checked at every N < 2^28 (fap.c
     to 2^27, fap_lean.c to 2^28); the identities of 1 and 3 are checked at selected N (fap_levels.c: the ladder at 9
     values of N; fap.c: S(π) = E_V and S(0) = E+ at six dumps); exact independent checks of 1-4 at 30 values of N, and
     every-N octave statistics to 2^22, in audit_F1".
   - Why: the scans never test statement 1, or statement 3's lag expansion. `HL := E_V - G0 - 2G1` is defined by
     subtraction, so the expansion is assumed, not checked.
2. **THM-4690, Setting, "Literal clock" bullet.**
   - Replace "This is the difference of the shortcut stopping times, since a strip changes the stopping time by 2 and
     the step by a." with: "Away from the ROOT this is the difference of the shortcut stopping times (a strip changes
     the stopping time by 2, a step by a). A path that enters u = 1 through a strip B^j(1) -> 1 (j >= 1) has literal
     clock st(x) - st(1) - 2, because st(B^j(1)) = 2j + 2 (1 lies on the cycle 1 -> 2 -> 1). The parity, hence
     σ = (-1)^t, is unaffected; G(D) and h_t use the literal clock."
   - Apply the same change to note §1.1 ("their literal integer clock t, the difference of shortcut stopping times").
   - Why: the original sentence is false for every V containing 1. See §2f.
3. **THM-4690(1), the merge identity.**
   - In "Σ_(z < z' in P^(-1)(y))" and "Σ_(z in P^(-1)(y))", restrict z to P^(-1)(y) ∩ {dep = k+1}. This is all of
     P^(-1)(y) for k >= 1, and the ports of y for k = 0.
   - In the proof, replace "As y runs over {dep = k}, the preimage sets P^(-1)(y) partition {dep = k+1}" with "...the
     sets P^(-1)(y) ∩ {dep = k+1} partition {dep = k+1}; for k = 0 the full preimages of V also contain V itself
     (P(V) ⊆ V, e.g. P(1) = 1)".
   - Define tree(y) as the set of sources whose **stopped** path passes through y. Then J(u) = c_V(u) for u ∈ V.
   - Why: read literally at k = 0, the recursion J(1) = c_0(1) + σ(1)J(1) + ... is circular. Statement 2 already uses
     the correct restriction.
4. **THM-4690(4) proof, sibling case analysis.**
   - After "(B) otherwise the sibling's chain ends there, at level L - 1." add: "This holds for all heads x except at
     most one per level in each case. With y = 3x+1 and Y = 6N+1, (A) fails exactly when y = ⌈Y/(2·4^(L-1))⌉ - 1, and
     (B) fails exactly when y = ⌈Y/4^(L-1)⌉ - 1. These O(log N) pairs are absorbed in the error term."
   - Why: (A)/(B) are stated as exact, but they fail for example at N = 6, x = 1. See §2e.
5. **THM-4690(4) statement and proof, split pairs.**
   - Replace "Each seed of V belongs to at most two sibling pairs, so pairs split by V number at most 2|V|" with "Each
     seed lies in at most one sibling pair (first members are ≡ 1 mod 8, second members ≡ 3 mod 16), so at most |V|
     pairs meet V".
   - In the statement, replace "|SIB - SIB_V| <= 18|V|" with "<= 9|V|". The old bound is true but rests on a false
     reason.
6. **THM-4690(4) and HYP-9325 Setting, definition of SIB_V.**
   - Define SIB_V over sibling pairs with **neither member in V**, or keep "common entry" and add: "except the pair
     (1, 3) when 1 ∈ V and 3 ∉ V, which enters 1 at clock lag 3".
   - Why: with the common-entry definition, Γ_1 = G(1) + SIB_V is not exactly "the other lag-1 pairs". For V = {1}: at
     N = 10, G(1) + SIB_V = -11 while the direct sum is -14; at N = 100, -60 versus -61.
   - It makes no difference on V_100, V_1000 and H_1000, which all contain 3. The receipts' "SIB" column is the
     common-entry SIB_V; the label should say so (also in HYP-9325's L-check line).
7. **THM-4690(4), the "In words" paragraph.**
   - "the sibling pairs remove exactly 10N/7 of it" -> "remove 10N/7 + O(|V| + log N) of it".
   - "so R(N)^2 - 24N/7 is the sum of the non-sibling pair correlations of the full transport" -> "...equals the sum of
     the non-sibling pair correlations of the full transport up to O(log N)". For V = {1}, the pair (1, 3) is a lag-3
     pair.
   - Why: E_0 - 2·SIB_V = 24N/7 + O(|V| + log N), not exactly 24N/7.
8. **HYP-9325, Evidence.**
   - Replace "E0/N lies within 10^-5 of 34/7, and SIB/N within 10^-5 of 5/7, from k = 20 on." with "From k = 20 on,
     E0/N lies within 5·10^-5 of 34/7 (within 10^-5 from k = 23). SIB_V/N (the common-entry sum in the receipts) lies
     within 5·10^-4 of 5/7 (within 10^-5 from k = 26)."
   - Why: in the receipts at k = 20, E0/N ∈ [4.85711, 4.85718] and SIB/N ∈ [0.71387, 0.71405]. My scan gives
     deviations of 4.03e-5 and 4.14e-4 for 2^20 <= N < 2^22. Even the full SIB deviates by 2.05e-5.
9. **HYP-9325 ("Unsigned arrivals keep E+/N ≈ 4.1-5.4 (octave means)") and note §1.4 ("E+/N has octave means
   4.1-5.4").**
   - Add "from 2^14 on", or widen the range to 4.1-5.6.
   - Why: the octave means are 5.530, 5.535 and 5.302 at 2^11-2^13. From 2^14 on the range is 4.132-5.418.
10. **HYP-9325, spectrum paragraph.**
    - Replace "it stays near 4-6 for θ <= π/3 and falls to about 1.2 at θ = π" with "it stays at 4-6.3 for θ <= π/3 and
      falls to 1.18-2.13 at θ = π (about 1.2 at the three largest dumps)".
    - Why: S(π)/N = E_V/N is 2.13, 1.29 and 1.43 at N = 7910, 66463 and 2^20.
11. **HYP-9325, Consequences and title.**
    - State the ranges: "(A)+(B1)+(C) ⇒ E_V <= 2.60N + O(log N) for every N >= 2^13", "(B)+(C) ⇒ ... for every
      N >= 2^12", and the R(N) bound "for N >= 2^13".
    - In the title, write "E_V <= 2.6N + O(log N) for N >= 2^13 on V_1000".
    - Why: E_V/N reaches 5.60 below 2^10, so no bound of the form 2.60N + O(log N) holds for all N. The hypotheses
      themselves are asserted only from 2^13, 2^11 and 2^12.
12. **Note §1.5.**
    - Replace "|E_0 - 34N/7| <= 11.2 and |SIB - 5N/7| <= 4.9 up to N = 4·10^6. fap.c reproduces them at every N."
      with: "|E_0 - 34N/7| <= 11.2 and |SIB - 5N/7| <= 4.9 at the 14 listed N <= 4,000,037. Over every N < 2^22 the
      maxima are 53.1 (N = 2,721,349) and 31.6 (N = 3,682,499), consistent with O(log N) (audit_F1). fap.c tracks E_0
      and the common-entry SIB_V at every N."
    - Why: read as a statement about all N, the original is false by a factor of about 5.
13. **Note, "Results at a glance".**
    - Row 3: "E_0 = 34N/7, sibling-source pairs 5N/7 ... E_V = 24N/7 + (non-sibling pair sums)" -> "E_0 = 34N/7 +
      O(log N), sibling-source pairs 5N/7 + O(log N) ... E_V = 24N/7 + (non-sibling pair sums) + O(|V| + log N)".
      "Exact" should then refer to the constants.
    - Row 7: "E_V <= 2.60N" -> "E_V <= 2.60N + O(log N) for N >= 2^13".
    - Why: same as fixes 7 and 11.

## 4. Optional improvements

1. **Sharper value set.** THM-4690(4) could state c_0 ∈ {-3, -1, 0, 2, 3}, which is sharper than |c_0| <= 3. It
   follows from the same table, and was verified at every N tested.
2. **Explicit error sizes.** Record that |E_0 - 34N/7| <= 53.2 and |SIB - 5N/7| <= 31.6 for all N < 2^22, with ratio to
   log2 N at most 3.13 for N >= 4. An explicit O(log N) constant is easy to prove from the per-level counting: at most
   about 3·9 per (level, ending) block.
3. **Cross-references.**
   - The limit 34/7 appeared earlier as a numerical value: THM-4667(6) has the ideal ladder E_0 = 4.857143, HYP-9288
     has T_0 = (51/7)X^2, and HYP-9290 has E_0/N -> 34/7.
   - SIB(N) = -INTRA_0(N)/2 in THM-4667(8), so 5/7 is the exact asymptotic of that d = 0 INTRA term.
   - THM-4690(4) supplies proofs for both, and could say so.
4. **"Why the denominator is 7" (Remarks).** The period 3 comes from B(x) ≡ x + 1 (mod 3). The mod-9 refinement
   B^3(x) ≡ x + 3 (mod 9) is true but not needed, and the "Solfeggio" identification is interpretive. The same applies
   to the note's "period 3 ... because 4^3 - 1 = 63", where cause and effect are reversed.
5. **"Essentially tridiagonal ... small noise" (note §1.4).** Each |G(D)| with D >= 2 is small. Together, though, the
   lags >= 2 contribute -0.285N at N = 6,904,411 (24% of E_V) and -0.245N on average at 2^27 (22%). Say this.
6. **The "Nyquist dip".** S(π) = E_V is an identity, so "cancellation is the Nyquist dip" is a reformulation, not a
   mechanism. The data content is supported: E+/E_V ≈ 3-4.4 (octave means 2^14-2^27), and S is low-pass at the six
   dumps. I recomputed the spectrum exactly at N = 7910 and 6,904,411. Type it as a reformulation plus a FINITE-EXACT
   observation.
7. **"The dissipation band recedes".** This is supported at the listed N. My ladder at 11 values of N shows it: the level
   by which half of E_0 - E_V has been removed is 1, 6, 8, 11, 11, 16, 18, 20, 21, 21, 21 for N = 10^3 ... 1.48·10^6,
   and levels 0-5 remove 1.37N at N = 7910 but only 0.13-0.27N at N >= 10^6. It tracks the growth of source depth
   (maximum depth 43 -> 241). Type it as a FINITE-EXACT observation and give this natural explanation.
8. **Ladder wording (note §1.3).** "typically 80-100% of it the sibling merge": at N = 6,904,411 the sibling share is
   60-106% (60% and 73% at levels 6 and 7; above 100% at levels 14, 19, 22 and 30). "Typically 80-105%" is more accurate.
9. **Old numbering in the pre-registration.** The L1-L7 JSON cites "THM-4670(5)" for the source ledger. After
   renumbering this is THM-4690(4), while THM-4690(5) is now the 19-adic obstruction. Add a one-line pointer in
   HYP-9325's status.
10. **HYP-9325 table caption.** "All values are divided by N, except the ratio G1/G0" should read "except the ratios
    G1/G0 and E_V/E+".
11. **Code labels.** `fap.c` and `fap_lean.c` print "SIB/N" for the common-entry sum. Relabelling it "SIB_V/N" would
    match THM-4690(4).

## 5. Pre-registration timeline and scoring (task e)

- **Q1-Q5.**
  - The pre-registration is in `adfeece3c1` (09:06:57, pushed 09:07:04). That commit contains only the N < 2^23 outputs
    (`fap_V1000_k22*.out`).
  - The 2^23-2^25 outputs first appear in `370e322936` (09:25:50).
- **T1-T4.**
  - The pre-registration is in `370e322936` (pushed 09:28:39).
  - The 2^26 outputs (`fap_V1000_k26_*.out`) first appear in `2c626644f0` (09:49:34).
- **L1-L7.**
  - The pre-registration is in `2c626644f0` (committed 09:49:49, pushed 09:49:53 according to the `origin/main`
    reflog).
  - `fap_lean_V1000_k27.out` first appears in `044e5bdf47` (10:05:05). The intermediate commit `11edd404cd` only
    renames files.
  - The receipt's mtime in the worktree is 09:53:31; it records an 85 s run, so the run started at about 09:52, after
    the push.
- Each prediction file was committed exactly once and never modified. The pre-registration commits contain no data rows
  for the predicted octaves: the only matches are in the JSON text and in commit messages. Between the Q
  pre-registration and the scan, `fap.c` changed only by adding new tracked quantities; E_V, E_port, G0 and G1 are
  unchanged.

Scoring from the receipts with the committed thresholds (`score_prereg.out`):

| Items | Result | Observed |
|---|---|---|
| Q1 | FAIL | 0.124, 0.155, 0.060 |
| Q2 | FAIL | 0.113 |
| Q3 | PASS | 1.328, 1.228, 1.421 |
| Q4 | FAIL | 0.232, 0.364, 0.435 |
| Q5 | FAIL | -0.076, -0.085, +0.010 |
| T1 | PASS | -0.2803 |
| T2 | PASS | [-0.2845, -0.2744] |
| T3 | PASS | 3.086 |
| T4 | PASS | 1.277, -0.024 |
| L1-L7 | PASS | as listed in HYP-9325 |
| L theorem check | PASS | |

All of these agree with HYP-9325 and note §1.2/§1.6.

## 6. The Consequences (task f)

1. **(A) + (B1) + (C).**
   - E_V = G0 + 2G1 + HL = (1 + 2ρ)G0 + HL, with ρ = G1/G0.
   - By (A), ρ <= -0.24; since G0 = Σ h_t^2 >= 0, this gives (1 + 2ρ)G0 <= 0.52·G0.
   - The sign of 1 + 2ρ is irrelevant to this upper bound; under (A) it lies in [0.40, 0.52] anyway.
   - Then (B1) gives 0.52·G0 <= 0.52·E_0, and (C) gives HL <= 0.07N.
   - Since E_0 = 34N/7 + O(log N), we get E_V <= 2.5257N + 0.07N + O(log N) = 2.5957N + O(log N) <= 2.60N + O(log N).
     **Correct.**
2. **(B) + (C).**
   - E_V <= E_0 - 2·SIB_V + 0.07N = (24/7 + 0.07)N + O(log N + |V|) = 3.4986N + O(log N + |V|) <= 3.51N + O(...).
   - **Correct.** It uses (B1), (B2) and (C), and not (A).
3. **Bound on R.**
   - R^2 <= |V|·E_V by Cauchy-Schwarz, since R = Σ_(u ∈ V) f(u) c_V(u) and |f| <= 1. This needs every source to reach
     V, which is true in the verified Collatz range.
   - sqrt(627 · 2.60) = 40.376, so "≈ 40.4 sqrt(N)" is **correct**.

There is no hidden assumption beyond:
- the ranges N >= 2^13 (or 2^12), which are fix 11;
- every source reaching V. Otherwise E_0 in (B1) must be restricted to the arriving sources, as in HYP-9290's scope note.

## 7. Prior art (task g)

- **The sibling merge.** P(x) = P(2x+1) for x ≡ 1 (mod 8) is an instance of the classical Syracuse identity
  S(4m+1) = S(m). It is folklore among Collatz "merging pairs": 8k+1 and 16k+3 reach the same odd number in two Syracuse
  steps, for example 9 -> 7 -> 11 and 19 -> 29 -> 11. THM-4667 already credits it to clock_duality's R5.
- **The rest.** No prior art was found for the signed ROOT transport, the ledger or the constants 34/7 and 5/7; this
  bookkeeping is specific to the repository, and THM-4690 claims no priority for it.
- **Searches.** Two quick web searches turned up only adjacent material: a Reddit r/Collatz post on canonical merging
  pairs, and Winkler, arXiv 1709.03385. Nothing about signed (-1)^(stopping time) currents or these constants.

## 8. Reproduction

From the repository root:

```bash
cd 04-computation/experiments/first_arrival_20261018/audit_F1
cc -O2 -o aud_point aud_point.c -lm && cc -O2 -o aud_scan aud_scan.c -lm && cc -O2 -o aud_fibre aud_fibre.c
./aud_point ../V_closed_Y1000.txt 1000 7910 16384 65536 66463 262144 658023 945589 1048576 1215227 1483383 > aud_point_V1000.out
./aud_point ../V_closed_Y100.txt 1000 7910 65536 66463 658023 945589 1048576 1215227 1483383 > aud_point_V100.out
./aud_point ../V_height_Y1000.txt 7910 66463 1048576 1483383 > aud_point_H1000.out
./aud_point V_root.txt 10 100 1000 7910 66463 > aud_point_Vroot.out
NOLADDER=1 ./aud_point ../V_closed_Y1000.txt 6904411 > aud_point_V1000_N6904411_noladder.out
for v in V_closed_Y1000 V_closed_Y100 V_height_Y1000; do ./aud_scan ../$v.txt 21 > aud_scan_${v}_K21.out 2> aud_scan_${v}_K21.err; done
./aud_fibre 100000000 268435456 1000000 > aud_fibre.out
python3 aud_constants.py > aud_constants.out
python3 aud_sibling_exceptions.py 2500 > aud_sibling_exceptions.out
python3 compare_receipts.py > compare_receipts.out      # 966 receipt numbers vs aud_scan
python3 compare_ladder.py > compare_ladder.out          # 224 ladder values vs fap_levels receipts
python3 check_hyp_table.py > check_hyp_table.out        # 221 HYP-9325 table entries vs receipts
python3 score_prereg.py > score_prereg.out              # Q1-Q5, T1-T4, L1-L7 from receipts
# pre-registration timeline
git log --diff-filter=A --format='%h %ad %s' --date=iso -- ../../../../05-knowledge/results/first_arrival_*prereg*.json ../fap_V1000_k25*.out ../fap_V1000_k26*.out ../fap_lean_V1000_k27.out
git reflog show --date=iso refs/remotes/origin/main | head -20
```

Resources:
- `aud_point` at N = 1.48·10^6 peaks at 471 MB and runs in about 5 s.
- At N = 6.9·10^6 without the ladder it uses 406 MB.
- `aud_scan` with K = 21 uses 221 MB and under 1 s.
- `aud_constants.py` takes about 5 s.
- Each `aud_point` run exits nonzero if any check fails; all runs exited 0.

Files in this folder:
- **Sources:** `aud_point.c`, `aud_scan.c`, `aud_fibre.c`, `aud_constants.py`, `aud_sibling_exceptions.py`,
  `compare_receipts.py`, `compare_ladder.py`, `check_hyp_table.py`, `score_prereg.py`.
- **Inputs:** `V_root.txt`, the stopping set V = {1}.
- **Outputs:** the `.out` / `.err` files.
- **Ignored:** the compiled binaries, listed in `.gitignore`.
