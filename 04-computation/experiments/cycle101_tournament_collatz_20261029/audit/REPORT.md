# Independent audit: cycle101 (THM-4849, THM-4848, HYP-9446, results note)

The independent auditor (a fresh agent, 2026-10-29) returned this report as text, because its harness does not let
subagents write report files. The session saved it here verbatim. The auditor's scripts and outputs are the `a*`
files in this directory. The corrections were then applied by the session (MISTAKE-696).

**Method.** Every claim was rechecked with my own code in `audit/`. The session's scripts were read but never
imported. The only session file read as data is `out/circulant_collatz.out`, which `a16_tables.py` compares against
values I computed independently.

**Version audited.** Commit `75bd8b27b6`, which the session pushed after its checkpoint (`a551aed5d0`). It adds the
P7, k = 4 data and refutes HYP-9446 (3), so the audit covers that too.

## Overall verdict: SOUND WITH CORRECTIONS

- **Numbers.** Every number I checked reproduces exactly: all tables, every β, the norms, det B₀, the Smith data,
  the locking and interval tables, the HYP-9446 counts, the cycle spectra and the inversion pairs. I recomputed P19
  and the 840 × 840 P7 four-run observer from scratch.
- **Proofs.** The proofs of THM-4849 (A), (B′), (C)–(G) and THM-4848 (i)–(v) are correct apart from the points
  below.
- **One literal ERROR.** "Two base points iff achiral" is stated for every cycle of the shape, but holds only for
  tournament cycles.
- **One GAP.** The "at most two" step treats only half of the multipliers; the fix is one line.
- **HYP-9446 (4) is not open.** It is proved below for every odd q, together with a closed count of tournament
  cycles.
- **Credit.** C3's pair equals the 5x+1 pair by construction, not by discovery. Reframes C2 and five-mirrors
  Proposition 5 need crediting.
- **Everything else is minor wording.**

## Findings

**1. ERROR — "two base points ⟺ achiral" is false for non-tournament cycles.**
- **THM-4849 (B), lines 96–98:** "A shape-`(k, q)` cycle has at most two tournament base points. It has two iff its
  connection sets are `m·{1, …, k}` for some `m`, i.e. it is a cycle of the rotational tournament `R_q`.
  Equivalently, iff it is reversal-symmetric (achiral)."
- **THM-4849 title:** "(B) Two-base-point lemma: a cycle carries at most two tournament base points, and exactly two
  iff it is a rotational-tournament (interval) cycle, iff it is reversal-symmetric."
- **Results note §2 table:** "| rotational `R_q` | the only cycles with two base points = the achiral ones |"
- **Problem.** The proof only covers cycles that already have a base point. At q = 11 the word (1,1,1,4,4) is achiral
  (its reversal is a rotation), yet no translate of {0,1,2,3,7} is a connection set.
- **Evidence (`a14_achiral.out`).** Achiral necklaces with no base point number 5 at q = 11, 14 at q = 13, 62 at
  q = 17 and 117 at q = 19. Among tournament cycles the equivalence holds for every prime q ≤ 31
  (`a1_dictionary.out`).
- **Correction (lemma):** "A shape-`(k, q)` cycle has at most two tournament base points. A tournament cycle (one with
  a base point) has two iff its connection sets are `m·{1, …, k}` for some `m`, i.e. it is a cycle of `R_q`;
  equivalently, iff it is reversal-symmetric (achiral). Cycles without a base point can be achiral from `q = 11` on:
  `(1,1,1,4,4)`; 5 such necklaces at `q = 11`, 14 at `q = 13`."
- **Correction (title):** "…, and a tournament cycle carries exactly two iff it is a rotational-tournament (interval)
  cycle, iff it is reversal-symmetric."
- **Correction (results-note row):** "| rotational `R_q` | the only tournament cycles with two base points = the
  achiral tournament cycles (cycles with no base point can also be achiral from `q = 11` on) |"

**2. GAP — the "at most two" step skips the multipliers k < m ≤ 2k.**
- **THM-4849, line 234:** "None does: for `2 ≤ m ≤ k`, the element `x = ⌊k/m⌋ + 1 ∈ I` has `k < mx ≤ 2k`."
- **Problem.** No argument is given for k < m ≤ 2k. The case is trivial, since m·1 = m ∉ I, but it is missing.
- **Correction:** append "; for `k < m ≤ 2k`, already `m·1 = m ∉ I`."

**3. MINOR (status) — HYP-9446 (4) is PROVED, for every odd q.**
- **Quoted text:**
  - HYP title: "(4) OPEN: at most two tournament base points per cycle for every odd q (holds for q <= 27)."
  - HYP §4 heading: "## 4. At most two base points for every odd `q` (OPEN)".
  - HYP §4: "Whether that is impossible is the open step."
  - THM-4849 line 123: "- *The bound.* No cycle has more than two base points; this is conjectured for all odd `q`
    in HYP-9446 (4)."
  - Results note line 188: "No cycle with three base points was found (HYP-9446 (4))."
  - Results note line 287: "4. At most two base points for every odd `q`."
- **Theorem.** For odd q ≥ 3, a connection set S ⊂ Z_q has at most one t ≠ 0 with S + t a connection set.
- **Proof.** Suppose S + t and S + u are both connection sets, with 0, t, u distinct. Put f = 1_S, A = ⟨t⟩ and
  B = ⟨u⟩, of odd orders n_A, n_B ≥ 3. By (B′), f is constant on every coset of A other than A, and in the coordinate
  j ↦ j·(−2t), f|_A is the indicator of I_A = {1, …, (n_A − 1)/2}. Likewise for B.
  1. **A ≠ B.** Otherwise u = mt with m ≢ 1 (mod n_A), and m·I_A = I_A. As in (B), that is impossible: for
     2 ≤ m ≤ (n_A − 1)/2 take x = ⌊(n_A − 1)/(2m)⌋ + 1, and for larger m take x = 1.
  2. **A ∩ B = 0.** Let E = A ∩ B, of order d, and suppose 0 ≠ E ≠ A.
     - Each x ∈ A∖E lies outside B, so f is constant on x + B ⊇ x + E. So f|_A is constant on every coset of E in A
       other than E.
     - In the coordinate of A, E = (n_A/d)Z with 3 ≤ d < n_A. The coset 1 + E contains 1 ∈ I_A and also
       n_A + 1 − n_A/d, which lies strictly between (n_A − 1)/2 and n_A and so is not in I_A. That is a
       contradiction.
     - Hence E = 0 or E = A. Symmetrically E = 0 or E = B, and by step 1, E = 0.
  3. **Contradiction.** Take a ∈ A∖0 and b ∈ B∖0. Then f(a + b) = f(a), by constancy on a + B, and f(a + b) = f(b),
     by constancy on b + A. So f is constant on A∖0, which contradicts f(a) + f(−a) = 1. ∎
- **Corollary (count).** The number of ordered pairs (S, t) is Σ_{n|q, n>1} φ(n)·2^((q/n−1)/2). So the tournament
  cycles number 2^((q−1)/2) − ½·Σ_{n|q, n>1} φ(n)·2^((q/n−1)/2). Of the doubly-based cycles, φ(q)/2 are rotational.
  For prime q the count is 2^k − k.
- **Evidence (`a1_dictionary.py`).** Exhaustive over every S and every t at q = 9, 15, 21, 25, 27, 33, 35 and 39:
  - the (B′) iff holds, and no cycle has more than 2 base points;
  - the pair counts are 10, 24, 40, 36, 62, 104, 80 and 176, matching the formula;
  - the HYP-9446 (4) table is reproduced.
- **Corrections.**
  - **HYP title:** "(4) PROVED (independent audit): at most two tournament base points per cycle for every odd q; the
    tournament cycles number 2^((q-1)/2) - (1/2) sum_{n | q, n > 1} phi(n) 2^((q/n-1)/2)."
  - **HYP status:** "apart from (1) and (3), which are REFUTED, (4), which is PROVED (independent audit,
    2026-10-29), and the proved parts marked PROVED".
  - **HYP §4:** change the heading to "— PROVED (audit)". Replace the Conjecture paragraph and the last paragraph with
    the theorem, proof and corollary above.
  - **THM-4849 (B′):** "- *The bound.* No cycle has more than two base points, for every odd `q` (proved in HYP-9446
    (4)): base points `S + t`, `S + u` with `⟨t⟩ ≠ ⟨u⟩` would force `⟨t⟩ ∩ ⟨u⟩ = 0` and then `f` constant on
    `⟨t⟩ ∖ 0`."
  - **Results note §2:** "No cycle has three base points, for any odd `q` (HYP-9446 (4), proved by the audit)."
  - **Results note §7 item 4:** "4. PROVED (audit): at most two base points for every odd `q`."

**4. MINOR — the locking theorem's title and summary drop the hypothesis p ∤ 2^v − 3.**
- **THM-4849 title:** "(E) Locking theorem: every marking cancels p^e iff every letter is congruent to v modulo
  ord_(p^e)(2)."
- **Results note lines 205–206:** "Every marking is ≡ `ξ_v (mod p^e)` iff every letter is ≡ `v` (mod
  `ord_(p^e)(2)`)."
- **Problem.** Without the hypothesis the title statement is false.
  - w = (3), v = 3, p = 5: D(wv) = 55, the letter satisfies 3 ≡ v (mod 4), but C(w·3) = 11, so nothing cancels.
    The same happens for (4) with v = 4, p = 13 and for (5) with v = 5, p = 29 (`a12_edge.out`).
  - `a2_locking.out` has 965 tests with p | 2^v − 3; in 428 of them (2) holds and (1) fails.
  - The body of (E) states the hypotheses correctly.
- **Correction (title):** "(E) Locking theorem: for p^e | D(wv) with p not dividing 6(2^v - 3), every marking
  cancels p^e iff …"
- **Correction (results note):** prefix the statement with "For `p^e | D(wv)` with `p ∤ 6(2^v − 3)`:" and add "(the
  hypothesis is needed: `w = (3)`, `v = 3`, `p = 5`)".

**5. OVERCLAIM/CREDIT — the C3 pair equals the positive 5x+1 three-cycles by construction.**
- **Quoted text:**
  - THM-4849 (G), lines 205–206: "**The owner's C3 pair is exactly the pair of positive 5x+1 three-cycles, and both
    are `R7`.**"
  - Results note §0, item 4: "4. **The same pair is the pair of positive 5x+1 three-cycles.**"
  - Results note §5 row: "| C3 pair = the two positive 5x+1 three-cycles | STRUCTURAL (cited transport plus a new
    criterion) | …"
- **Problem.** collatz_fermat_reframes built C3 from its own C1 table. The table rows are `(1,1,5)` ↔ `13,33,83` and
  `(1,3,3)` ↔ `17,43,27`. C3 says "The pair `u=(1,1,5)`, `u'=(1,3,3)`", and the packet lists "the positive `5x+1`
  cycle `13,33,83`" among its hostiles.
- **What is genuinely new.**
  - The all-letters-odd criterion.
  - A structural reason that both are R7, which the files do not give. Three distinct odd parts sum to at least
    9 > 7, so an all-odd word of shape (3,7) repeats a letter. It is then a rotation of its reverse, hence achiral,
    hence R7 by (B).
- **Correction (G):** "The owner's C3 pair was taken by collatz_fermat_reframes from its own C1 table, so it is the
  pair of positive 5x+1 three-cycles by construction (inherited). What (G) adds is that both are `R7`, since at
  `q = 7` an all-odd word repeats a letter and is therefore achiral, and the criterion below, which says why exactly
  these two of the five (3,7) words are integral for 5x+1."
- **Correction (item 4):** "4. **The same pair is, by construction, the pair of positive 5x+1 three-cycles**
  (reframes took C3's words from its C1 table)."
- **Correction (§5 row):** "INHERITED (C3 used C1's 5x+1 words); the criterion all letters odd ⟺ `3 | C_5` is new and
  STRUCTURAL".

**6. CREDIT — reframes C2 and THM-4484 already have the five 3x+101 cycles.**
- **Where.** The THM-4849 PRIORITY list credits reframes C1 and C3 only. Its (C) and its title ("The five 3x+101
  cycles are …") rest on the count of five.
- **Prior statement.** Reframes C2: "At `(k,H)=(3,7)`, all five primitive word necklaces for `a=3` have denominator
  101. Consequently `Z_(3,b)(3,7)=5` if `101|b`, and zero otherwise." THM-4484's shift criterion gives the same.
- **Correction (add to PRIORITY):** "- that the five (3,7) necklaces all have denominator 101, so that 3x+b has exactly
  five (3,7) cycles iff 101 | b, is C2 of collatz_fermat_reframes_20261029 (and the shift criterion of THM-4484);
  (C) adds their identification with the two circulant tournaments on Z_7."

**7. CREDIT — five-mirrors already defines reversal on rational cycles.**
- **THM-4849 PRIORITY list:** "- reversal invariance of the cycle trace is Proposition 6 of
  collatz_five_mirrors_20260929;"
- **Problem.** Five-mirrors Proposition 5 and its reversal survey (§5) already have:
  - reversal as an involution on rational cycles;
  - its fixed points, "exactly the cycles whose word is a rotation of its reverse";
  - genuine pairs, e.g. the 3x+37 cycles {23,53,49} ↔ {29,31,65}.
- **Correction:** "- reversal as an involution on rational cycles, with fixed points the cycles whose word is a
  rotation of its reverse and genuine reversal pairs, is Proposition 5 and the reversal survey of
  collatz_five_mirrors_20260929, and reversal invariance of the cycle trace is its Proposition 6; (B) and (F) add that
  on tournament cycles reversal is the converse and the achiral tournament cycles are exactly the R_q cycles;"

**8. CREDIT (minor) — related packets not cited.**
- **clock101_joint_ports_20261029.** It was committed at 11:50 (`c7ca37bdc0`), before the 12:04 checkpoint, and works
  on the same numbers 101, 19, 79 and 169. Its RECOVERY.md maps older sources on each; its P1 has the 7x+169 cycle
  79/169, 361/169, 337/169. It does not overlap THM-4848/4849.
  - **Correction:** add it to results note §1.
- **ramanujan_heegner7_trivial_cycle_20261006.** Its Proposition 1, that the code's Gauss period ζ + ζ² + ζ⁴ equals
  τ₇ − 1, is the q = 7 precedent of (F).
  - **Correction:** in (F), add "(the q = 7 Gauss period of the trivial cycle's code is Proposition 1 of
    ramanujan_heegner7_trivial_cycle_20261006)".

**9. MINOR — a vacuous case in (D).**
- **THM-4849, line 148:** "- If `p | 2^v − 3`, nothing cancels."
- **Problem.** If 2^v ≡ 3 (mod p), then D(wv) ≡ 3·D(w) (mod p). So p | D(wv) forces p | D(w), which (D) has
  excluded, and the case never occurs. 2224 instances confirm the implication, and the case never arose in 6582 law
  tests (`a2_locking.out`). For the same reason, the p-integrality clause in (E)(1) is automatic.
- **Correction:** "- The case `p | 2^v − 3` does not occur here: then `D(wv) ≡ 3·D(w) (mod p)`, so `p | D(wv)` would
  force `p | D(w)`."

**10. MINOR — Turner's theorem is stated too strongly in (A).**
- **THM-4849, line 86:** "By Turner's theorem these are all the isomorphisms between circulants of prime order."
- **Problem.** The isomorphisms are the affine maps x ↦ m·h·x + b with h ∈ H. Since the automorphism group has odd
  order, Burnside's theorem puts it inside AGL(1,q). Turner 1967 covers graphs; the directed case is Babai's 1977
  theorem. The frontmatter's CITED item is correct.
- **Correction:** "By Turner's theorem (graphs) and Babai's (digraphs), circulants of prime order are isomorphic only
  when their connection sets are multiplier-related, and every isomorphism is affine, `x ↦ m′x + b`."

**11. MINOR — THM-4848 table: two signs missing on det B₀.**
- **Lines 152 and 155:** `3^3·7·13^3` and `3·11^2·19^3·23·43·241^2·22291^2`.
- **Problem.** det B₀ does not depend on the choice of representatives, so its sign is meaningful, and both values
  are negative:
  - P7, k = 2: det B₀ = −415233;
  - P11, k = 2: det B₀ = −71065187852341089324093.
  The table signs the P7, k = 3 and P19 entries.
- **Correction:** `−3^3·7·13^3` and `−3·11^2·…`.

**12. MINOR — THM-4848, line 163: wrong column.**
- **Quoted text:** "is the second column's entry".
- **Problem.** The quadratic norm is the fourth column.
- **Correction:** "the entry in the column 'quadratic norm of β'".

**13. MINOR — THM-4848, Z₁₃ row: one claim is unverified.**
- **Quoted text:** "none (all have `3 | f`)".
- **Problem.** The 170-digit block is unfactored. Each of its primes appears as (Z/p)³, which is consistent with f = 3
  and also with f = 1 (p ≡ 1 mod 13). In the f = 1 case it would be a non-trivial instance of (iii). Data:
  N(det B₁) = 3^15·13^6·1907³·2011³·c170³ (`a8_z13.out`).
- **Correction:** "none at 3, 1907, 2011 (`f = 3`); the 170-digit block is unfactored, so its residue degree (1 or 3)
  is unknown".

**14. MINOR — THM-4848 controls note.**
- **Quoted text:** "That every odd prime in the controls' `X_p` is `≡ 1 (mod q)` is the generic behaviour … It is
  not a theorem here."
- **Problem.** By (iii), every multiplicity is a multiple of f, and in R7 and R11 each odd prime occurs exactly once.
  So f = 1 is forced, including for R11's unfactored 196-digit block.
- **Correction:** "Each odd prime of the controls' `X_p` occurs exactly once; since every multiplicity is a multiple
  of `f` by (iii), this forces `f = 1`, i.e. `p ≡ 1 (mod q)` (also for the prime factors of R11's unfactored
  196-digit block). What is generic, not a theorem, is that the odd part of `N(det B_1)` is squarefree."

**15. MINOR — THM-4849, line 292.**
- **Quoted text:** "101 is `D_7`, the last prime at which every cycle of the shape is a tournament."
- **Problem.** It reads as if 101 were the last prime; the last prime q is 7.
- **Correction:** "101 is `D_7`, and 7 is the last prime `q` at which every cycle of the shape is a tournament."

**16. MINOR — the "19² only" justifications are wrong.**
- **Quoted text:** THM-4849 line 298, "These follow from `C mod 19²` only."; results note §5 row, "divisibility of
  `C` by `19²` only".
- **Problem.** The 361 comes from 5 | C. The 3x+19 cycle needs 5·19 | C with 19² ∤ C. The GENERIC typing itself is
  fine: the reduced denominators at q = 11 are 1805 (22 cycles), 361 (3), 95 (1) and 19 (1).
- **Correction:** "These follow from the residues of `C` modulo 5 and modulo `19²` (i.e. modulo `D_11`) only."

**17. MINOR — results note §2 table, decimation row.**
- **Quoted text:** "decimation of the word (moves between cycles with the same denominator)".
- **Problem.** Decimation preserves D_q, but not the reduced denominator: P11's QR cycle has 361 and its NQR cycle
  1805.
- **Correction:** "(moves between cycles with the same `D_q`; reduced denominators can differ, e.g. P11's QR and NQR
  cycles, 361 and 1805)".

**18. MINOR — results note line 70.**
- **Quoted text:** "but 101 is never named there."
- **Problem.** THM-4484 does contain 101 ("d = 14303, least element 101"), just not as the clock 2^7 − 3^3.
- **Correction:** "but the clock `101 = 2^7 − 3^3` is never named there (101 appears only as an unrelated least
  element, of the `d = 14303` cycle)."

**19. MINOR — HYP-9446 (2) counts odd primes only.**
- **Quoted text:** "299 of the 782 primes below 6000" and "… primes `q < 6000`".
- **Problem.** π(6000) = 783, so the 782 are the odd primes. The 299 is correct.
- **Correction:** "782 odd primes".

**20. MINOR — HYP-9446 (1), description of the test.**
- **Quoted text:** "every prime `p` with `ord_p(2) ≥ 3` dividing a difference of letters".
- **Problem.** The condition is that ord_p(2) divides every letter difference. Only the wording is affected: my
  definition-based recount reproduces the table.
- **Correction:** "every prime `p ≥ 5` whose order `ord_p(2)` divides every difference of letters (by THM-4849 (E)
  the only candidates)".

**21. MINOR — HYP-9446 (3), "Open residue" is ill-posed.**
- **Quoted text:** "…always divisible only by 2, `q`, inert primes and a bounded set of split primes?"
- **Problem.** "Bounded set" is undefined.
- **Suggested text:** "For fixed `q`, is the set of split primes dividing the content `gcd(a, b)` of `β_k(P_q)`, over
  all `k`, finite? (Data: none for `k ≤ 3`; `{71}` for P7, `k = 4`.)"

**22. MINOR — THM-4849's PROVED list is incomplete.**
- **Quoted text:** "PROVED (…): (A) given the cited realisation theorem; (B); (D); (E) in both directions, with the
  interval corollary; the Gauss-period remark in (F)."
- **Problem.** Section 2 also proves (B′), the trace relation in (F), and (G).
- **Correction:** "(A) given the cited realisation theorem; (B); (B′) in both directions; (D); (E) in both
  directions, with the interval corollary; (F), the trace relation and the Gauss-period remark; (G)."

## Verified and found correct (scripts in audit/)

**THM-4849 (`a1_dictionary.py`, `a2_locking.py`, `a6_misc.py`, `a16_tables.py`)**
- **Conventions.**
  - The cycle formula, exact valuations and T-word odd positions, on 400 random words.
  - For every connection set at q = 3–13: T(z_v) = z_(v−1); z_v odd iff v → 0; exact period q; U-word = gap
    sequence; decimation; reversal = converse.
  - The values D₃ through D₁₃ as stated.
- **(B).**
  - 2^k − k cycles and at most two base points for every prime q ≤ 31.
  - Cat_k necklace counts for q ≤ 23.
  - Class tables at q = 11 (5, 10, 10, 2) and q = 13 (6, 12, 12, 12, 12, 4).
  - At q = 7 a brute force over all 7! relabellings gives exactly 2 isomorphism classes. At q = 11 and 13, invariants
    separate the multiplier orbits.
- **(B′).** Proof re-derived, including the converse; the iff checked exhaustively for q = 9 to 39.
- **(C).**
  - The whole (3,7) table, with the per-marking denominators.
  - The fibre {7:1, 25:1, 35:5, 175:8}.
  - Reduced denominators and traces of all 90 cycles at q = 7, 11, 13 match the session file.
- **(D).** The identity, checked symbolically; the attachment law on 6582 tests.
- **(E).**
  - Both directions re-derived.
  - 12003 tests, 532 of them with e ≥ 2: 711 locked cases, all agreeing with the theorem.
  - The interval table for q < 100 is identical to `locking_iff.out`.
- **(F).** The trace relation on 500 random cycles. The Gauss-period identity for all 31 primes q ≡ 3 (mod 4) with
  q < 300.
- **(G).** The 5x+1 cycles and their words; the all-odd criterion; C(1^(k−1), k+2) = 3^k − 2^k.
- **Hamiltonian paths.** H(R7) = 175 and H(P7) = 189.
- **Citations** are plausible and correct: Turner 1967 JCT 3; Babai 1977 Acta Math. Hungar. 29; Böhm–Sontacchi 1978;
  Lagarias 1990 Acta Arith. 56; Dvoretzky–Motzkin 1947.

**THM-4848 (`a3_runobs.py`, `a4_smith*.gp`, `a5_p7k4.py`, `a8_z13.py`, `a13_p19*.{c,py}`, `a15_monomial.py`,
`a6_misc.py`)**
- **Validation of my observer code.** My Smith form of M₃(P7) reproduces THREE_RUN_GLUE (R9) exactly, and its
  cokernel order matches (R12).
- **Proofs.** All parts re-derived, including the Z_(p) splitting and the orbit count |H⟨p⟩|/f.
  - The explicit monomial conjugation holds numerically.
  - Exact σ_m-invariance holds for all m in H.
  - det M = det B₀·N(det B₁) modulo three primes, in every case.
- **P7, k = 2.**
  - β = 11904 − 23232η; quadratic norm 2^14·3²·7·1451.
  - Smith 2-part {1:6, 2:9, 6:3}; X₃ = (Z/3)^6; (Z/1451)^3.
  - Converse: β(NQR) = 35136 + 23232η, the conjugate.
- **P7, k = 3.**
  - det B₀ = −2^14·3^18·17^4.
  - β and the quadratic norm 2^22·5^4·11²·79·5861²·392809, as stated.
  - Nontrivial 2-part (Z/2)^21 ⊕ (Z/8)^3 ⊕ (Z/16)^9.
  - At 79 the translation clock's cubic is X³+13X²+12X−1 for QR and X³+67X²+66X−1 for NQR, so (v) holds.
- **P11, k = 2.**
  - (Z/89)^5 and (Z/62987)^5, absent from B₀.
  - 2-part (Z/2)^30 ⊕ (Z/8)^10.
  - Every odd-exponent prime is split.
- **Controls.**
  - R7: the primes 29, 71, 337, 303073 and 2277647 each occur once.
  - R11: nine primes ≡ 1 mod 11 each occur once, and X₂ = (Z/2)^60.
- **P19.**
  - det B₀ as stated; σ₄-invariance and the ninth-power law hold.
  - The 310-digit quadratic norm = 2^8·5²·11·47·61·263·c299.
  - Content 2^4; (0, 2) at 5.
- **Z₁₃.** Invariance holds, and two labellings give identical data.
- **P7, k = 4.**
  - β matches digit for digit, with content 2^41·5³·7³·71².
  - At 71, (v_𝔭, v_𝔭̄) = (2, 3): the refutation of HYP-9446 (3) stands.
  - Local Smith forms of the 840 × 840 M₄: (Z/71)^15, (Z/2843)^9, (Z/79²)^3, (Z/151)^9 and (Z/7529)^6, with no
    contribution from B₀.
- **Credits.** CLOCK_MINIMAL_PAIRING C5–C8 and THREE_RUN_GLUE R9/R12 contain exactly what THM-4848 says it inherits.

**HYP-9446 (`a2_locking.py`, `a3`/`a5`/`a13`, `a6_misc.py`)**
- **(1).** The table q = 5–23 is reproduced from the definition: 18 locked cycles at q = 23, 14 of them
  non-rotational. The example word checks out.
- **(2).**
  - 299 of 782; the first 15 values; the pattern first fails at q = 109 with gcd 23.
  - The p = 5, 7, 11, 13 congruences hold for every q < 6000.
  - The Ailon–Rudnick (Acta Arith. 113, 2004) and Bugeaud–Corvaja–Zannier (Math. Z. 243, 2003) attributions are
    correct. The restricted question implies the (2,3) Ailon–Rudnick conjecture, so it is hard.
- **(3).** Every chirality entry, and min(v_𝔭, v_𝔭̄) = v_p(gcd).
- **(4).** The table, plus the proof in finding 3.
- **(5).** 5 | D_q iff q ≡ 3 (mod 8), checked for q < 2000. The Paley-pair entries at q = 11, 19 and 43.

**Priority.** Repo-wide greps found no earlier statement of:
- the odd-positions circulant dictionary;
- the two-base-point lemma;
- the locking theorem.

The session101 credits are accurate. Their centre C ≡ −2^H/3 is equivalent to x ≡ 1/(2^v − 3).

**Results note**
- **Cycle spectra.** My own C search to 2·10^7 for b = 19, 23, 35, 79, 101, 169 and 175 gives exactly the cycles in
  `cycles_3xb.out`, so every §6 claim holds (`a7_cycles.c`).
- **Inversion pairs.** 2805 periodic rationals and 192 mutual pairs, recomputed (`a11_inversion.py`).
- **§1 spot checks.**
  - The 23→35→53→5 path.
  - The two shapes of the carry 19, and the ladder 19/(2^q − 27).
  - The 7x+169 cycle {79, 361, 337}.
  - ord₇₉(2) = 39, ord₁₀₁(2) = 100, and (2/101) = −1.
  - 101 | s(19) − 2⁹ (`a9_camion.py`).
  - 169 passes the strong Camion test while 1313, 1501, 1919 and 3211 fail, by direct recurrence; 13² ∥ D(13)
    (`a9_camion.py`).
  - det(S + J) = 169 for R₁₃.
  - H = 169 is absent at n = 7 and present at n = 8 (gentourng, `a10_hspectrum.py`).
  - The R7 skew polynomial x(x⁶+21x⁴+35x²+7).
- **Not re-verified** (cited prior work only): THM-853/854, HYP-9080, the session54 ranks, THM-4566, THM-4820,
  THM-4751, THM-4656, and the "pushed 11:08–11:11" timing.

**Reproduction.** Run the `a1_…` to `a16_…` scripts in `audit/`. Each one's output is saved beside it as `.out`, and
the gp Smith-form scripts read the matrices in `mats/`. The C programs `a7_cycles.c`, `a11_cycles_all.c` and
`a13_p19.c` compile with `cc -O2`.

## Application by the session (2026-10-29)

All 22 findings were applied to THM-4848, THM-4849, HYP-9446, the results note and both INDEX entries.
- The proof in finding 3 is now THM-4849 (B''), credited to this audit.
- HYP-9446 (4) is marked PROVED.
- The ledger entry is MISTAKE-696.
