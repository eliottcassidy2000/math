        # Message: mac-mini-2026-10-06-oaimath2: second openai/math reading (12 manuscripts) - Heegner rungs 3-chromatic (roadmap refuted), algebraic spherical sets Ramsey (cond. #172), gap-determined Erdos-592 witnesses die (THM-521 D for gap rules), Littlewood tower, slice lemma; two audits applied

        **From:** mac-mini-2026-10-07-S?
        **To:** all
        **Sent:** 2026-10-07 00:30

        ---

        mac-mini-2026-10-06-oaimath2, close-out.

**Owner prompt:** "spend another similar session picking out another handful or two of the papers from that repo [openai/math] and looking for connections and extensions and solutions to problems".

**What was read.** Twelve manuscripts, in six read-only lanes (no OpenAI code run):
* #165 crossing numbers;
* #090 triangular optimality and #022 Duffin–Schaeffer;
* #076 Littlewood (+2) and #155 monotile;
* #158 plane not 5-colourable and #172 Euclidean Ramsey;
* #164 Hindman FS∪FP, #189 cycle–clique, and the sharp log-exponent pair;
* #180 Barnette and #047 affine cancellation.

Two independent adversarial audits ran after the "audit in progress" checkpoint 9258bf9d52. All their corrections are applied (MISTAKE-577 sharpened, MISTAKE-579).

**Records.**
* Note: 05-knowledge/results/oai2_openai_math_second_reading_20261006.md.
* New canon: THM-4558 to THM-4562.
* New hypotheses: HYP-9215, HYP-9216.
* Mistakes: MISTAKE-574 to MISTAKE-577 and MISTAKE-579 (578 is opus S19's).
* Updates: THM-913, THM-922, THM-431 (both files), HYP-2267, THM-418, THM-470, THM-521; PROBLEM-LEDGER; the HYP index (also indexing HYP-9212 to HYP-9214).
* Scripts: oai2_20261006_{crossing_books, lattice_u21, littlewood_tower, fields_ramsey_checks}.py, each with .out.
* Reader scripts in oai2_20261006_readers/; audit reports and scripts in oai2_20261006_audits/.

**What changed.**

1. **Crossing (#165, Lean claimed).**
   * THM-922 (III) holds for every `n`.
   * The parity sum-class drawing of `K_(m,m)` has exactly `Z(m,m)` crossings (telescoping proof; the upper bound itself is classical).
   * The class-colouring minimum is unconditional for `m ≤ 8`.
   * The 2-page Zarankiewicz conjecture holds CONDITIONAL on #165.
   * THM-913 is the classical DDS construction (MISTAKE-574).
2. **Lattice (#090).** The triangular lattice attains `u(21) = 57` (MISTAKE-575).
3. **Littlewood (#076): THM-4561.**
   * The skew tower's rows are Thue–Morse-like (`F < 2`; `max F → 0` up to a sketch).
   * Rudin–Shapiro switching of the same tournament gives every row `sup ≤ 8.24√N`.
4. **Plane colouring (#158): THM-4558.**
   * Residue colourings at conjugation-stable primes are a KNOWN technique: Woodall, Fischer 1990, Madore, plus the MildlyMeticulous hn-2adic and decalion89 repositories from 2026, which the audit found.
   * New here: the Heegner rungs `Q(√−3, √−(4N−1))`, `N ≡ 2 mod 3`, are 3-chromatic; the `N = 3n` family is 4-chromatic; the Polymath field is exactly 5-chromatic.
   * So our Heegner roadmap is refuted. Class number does not govern the 5-chromatic step (Exoo–Ismailescu's `√−247`). "Measure cannot reach 5" was also false (MISTAKE-576).
5. **Ramsey (#172, Lean claimed): THM-4559.** Every algebraic spherical set is Ramsey, CONDITIONAL on #172 (separability idempotent; audit PASS).
6. **Erdős 592: THM-4560.**
   * The Galvin–Glazer idempotent-ultrafilter method shows that every gap-determined triangle-free graph on `N^n` leaves a binary subgrid independent.
   * So `t_dead(Finv) < ∞` for every `n`, and THM-521 D is unconditional for gap-determined witnesses. Row-invariant witnesses are not covered.
   * THM-470's Finv label is corrected: the `n = 2` cutoffs are 4 (gap-determined) and 5 (row-invariant) (MISTAKE-577).
7. **Barnette (#180): HYP-9215.** A knight-torus prescribed-path blocking law for paths of at most 4 moves, FINITE-EXACT for `n = 5..8`. Knight tori are P7-Hamiltonian.
8. **Cancellation (#047): THM-4562.**
   * A slice lemma, with the iff corrected by the audit.
   * THM-1300 is fully non-slice: fibre counts are PROVED for all `q` prime to 6.
   * #047's base is our Russell cylinder.
   * HYP-9216 asks whether `D(X) ≅ A_4`.
9. **Cross-link.** The rational plane reduces mod 7 onto the 7 × 7 knight torus (per coset, `z ↦ (1+2i)(z − z₀)`).

**Proof status.** No major open problem is closed. THM-4560 settles the gap-determined case of the repo's Erdős 592 strong-witness question. The rest are exact values, corrections, conditional extensions and new hypotheses. `κ(13)` and `κ(17)` stay in `[5, 6]`; our 25-minute 5-colouring runs did not finish, and decalion89 reports `κ(13) = 6`.

**Next obligations.**
* Decide `Finv (3,7)`.
* A finite (Milliken–Taylor) bound for THM-4560.
* Row-invariant strong witnesses (HYP-2558).
* `κ(17) = 5`? It would make `Q(√−3, √−7, √−11)` exactly 5-chromatic.
* Re-audit THM-431 (C) and HYP-2301 with product sets.
* `u_A(N)` for `13 ≤ N ≤ 20`.
* The tower's all-rows-flat switching question.
* HYP-9215's merging step.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
