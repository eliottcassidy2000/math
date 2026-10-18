S43b REFEREE REPORT -- THM-4685, HYP-9310, THM-4686, HYP-9311 (opus-2026-10-18-S43)
=====================================================================================
Referee: independent adversarial audit subagent, 2026-10-18.
Worktree: C:\Users\Eliott\AppData\Local\Temp\math-wt-explicit-modal-38 (branch codex/session-explicit-modal-38).

Scope and state
  - Audit started at HEAD 1d7fd58ab0. The session kept committing during the audit: 7f331c46e5 (HYP-9311),
    eb53cfd571 / 91c22f341d (HYP-9311 part (O)), e8bc51fca4 (THM-4686 incl. the "KLMSY Conjecture 9 at s = 1"
    section), 5b770a4571, dfce3de3c6 (probes; not audited).
  - THM-4685, HYP-9310 and the results note did not change after 1d7fd58ab0. THM-4686 is audited as of e8bc51fca4,
    HYP-9311 as of 91c22f341d (both were added to the scope by the coordinator mid-audit).
  - Nothing was committed, pushed or modified. All referee code and outputs are in
    04-computation/experiments/s43b_audit_20261018/ (written from scratch; no session module imported).
  - Literature PDFs (KLMSY, Beukers-Vlasenko, Z.-W. Sun) were downloaded to a scratch directory outside the repository
    and text-extracted in isolated mode.

Files audited
  01-canon/theorems/THM-4685-half-period-twins-from-dixon-and-reflection-and-the-cm-harmonic-law-at-the-self-dual-points.md
  05-knowledge/hypotheses/HYP-9310-cm-harmonic-log-law-for-sum-binom-2k-k-cubed-over-m-to-the-k.md
  05-knowledge/results/halfperiod_twin_classes_cm_log_law_20261018.md
  01-canon/theorems/THM-4686-the-cm-harmonic-law-follows-from-the-cm-supercongruences-at-both-ends-of-a-2-isogeny.md
  05-knowledge/hypotheses/HYP-9311-cm-harmonic-law-for-the-four-ramanujan-families-and-the-truncated-canonical-coordinate.md
  and the scripts/outputs they cite: halfperiod_twin_classes, halfperiod_twin_class_forms, halfperiod_twin_mechanisms,
  cm_harmonic_log_law(_extended), cm_points_forms_legendre, cm_law_lucas_block_form, cm_detector_lucas_block_scan,
  cm_law_from_known_supercongruences, klmsy_conjecture9_s1_check, cm_harmonic_law_four_families,
  pshift_lemma_four_families, cm_lucas_block_mod_p3_probe, cm_lucas_block_ordinary_mod_p3 (all *_20261018).


SUMMARY VERDICT
---------------
* THM-4685: MATHEMATICALLY CORRECT. No gap, sign error or missing hypothesis was found.
  - Every step was re-derived by hand: Lemma 1 (S1-S6, M1-M4 and the consequences), Lemma 2 (a)-(d), part (i)
    (i-a)-(i-d) with the exact factors (-1)^n, 4^n, 16^n, Lemma 3, and part (ii) (the m <-> 4096/m equivalence,
    m = -64, m = 64 in both parities, the twin form).
  - All of it holds in exact arithmetic for every prime 5 <= p <= 3000 (mod p^2) and 5 <= p <= 6000 (mod p).
  - p = 3 must be excluded and is: (i-a) fails mod 9.
  - The calibration of part (i) (corollary of the classical zeta(2) supercongruences) is adequate, and could be
    stronger: (i-c) follows termwise from b(n) (n1).
  - The census statement "exactly seven" is right under the stated convention (eps in {1, (-1)^n}). It was reproduced
    by a full pairwise method on 76 primes p <= 400 (m7).
  - The modular identifications are all correct for p <= 1000.
  - BUT the theorem's remarks treat the CM harmonic law and (H) as open (HYP-9310). They are published: see MAJOR M1.
* HYP-9310: every numerical claim reproduced, and extended.
  - The law (L_p) has 0 failures at all 8 CM points for all 9590 primes p <= 10^5 and for the 40 primes in
    [199500, 200000]; (K) likewise.
  - The law fails at generic m (195 non-CM values, pass rate <= 4.9%).
  - The Lucas-block lemma, the equivalence (L_p) <=> (B1), and the F_(p^2) identity all hold.
  - The detector scan was reproduced exactly; its survivor lists are complete in range and every survivor is CM.
  - BUT the main statement is Kibelbek-Long-Moss-Sheller-Yuan (KLMSY) Corollary 3, (H) is KLMSY Corollary 4, and the
    "Lucas-block form" is KLMSY Lemma 10 with their proof of Corollary 3 (M1). It is not an open conjecture, and the
    "pre-registered prediction" predicted a published theorem.
* THM-4686: the proof is CORRECT.
  - Algebra verified; z^n == (z/p)(1 + p q_p(z)/2); eps_1 = (-1)^n (z/p) eps_2.
  - alpha^2 agrees at lambda and 1/lambda, because E_lambda is 2-isogenous to the -2 twist of E_(1/lambda).
  - The p = 7 caveat is right.
  - It is a NEW ROUTE to the published KLMSY Corollary 3 (Theorem 2 at both ends + the reflection), not an upgrade of
    an open hypothesis (M1).
  - It copies a sign misprint of KLMSY Theorem 2: the character ((lambda-1)/p) must be ((1-lambda)/p). Harmless for the
    proof (m1).
  - The Conjecture-9 (s = 1) corollary is correct; digits a >= p and carries are handled. But it is a two-line
    consequence of KLMSY's own Lemma 10 (stated for any n) plus their Corollary 3. No later proof of KLMSY Conjecture 9
    was found (M1).
* HYP-9311: all numbers reproduced, including (O) exactly: mod p^3 the congruence holds at 908/908 ordinary
  prime-instances and only sporadically (16/1038) at supersingular ones.
  - New observation: the p-shift relation also holds at carry indices, so (L_p^A) <=> (B1') for all four families.
    The "pre-registered" law test at the 33 detector survivors was therefore largely implied by the detector's own data
    (m3).
  - The heuristic attributes to Beukers-Vlasenko more than they state (m5).
* Pre-registration: git cannot order docstring against run for the HYP-9310 law, the detector (P1, P2) or the HYP-9311
  law (same commits). HYP-9311(O) has a verifiable prediction-before-output commit order (E).


FINDINGS
--------
MAJOR

M1. The CM harmonic law, (H), the Lucas-block equivalence and the Conjecture-9 (s = 1) corollary are all in KLMSY 2016,
    a paper the session cites but did not read.
    Where:
    - HYP-9310: status "OPEN, strongly supported", "Pre-registered prediction PASSED", "The S43 reading of the paper
      found no statement of (H) itself", and the "Lucas-block form (PROVED equivalent)" section without attribution.
    - THM-4685: remarks table (classes 3 and 7 "<=> L(1), HYP-9310"); "Novelty not established".
    - Results note: sections 3-5 ("Literature", "Open: prove (L_p) at m = 1, -8, 16 ...").
    - THM-4686: status ("This upgrades HYP-9310's main statement from OPEN to PROVED ... The contribution is the
      identity ..."), and its Conjecture-9 "Novelty" paragraph.
    - HYP-9311: the A1 row and "OPEN".
    What is wrong: J. Kibelbek, L. Long, K. Moss, B. Sheller, H. Yuan, "Supercongruences and complex multiplication",
    J. Number Theory 164 (2016) 166-178, arXiv:1210.4489. Numbering below is from arXiv v3; the published PDF was not
    accessible (HTTP 403).
    - Corollary 3. If E_lambda is CM, then for almost all primes p with lambda in Z_p,
          sum_{i<=(p-1)/2} C(2i,i)^3 (lambda/64)^i ( 6(H_2i - H_i) + ((lambda/64)^(p-1) - 1)/p ) == 0  (mod p).
      With lambda/64 = 1/m and q_p(1/m) == -q_p(m), this is exactly (L_p): 6 h_m == q_p(m) F_m.
      - It covers every CM lambda != 1, i.e. m in {1, 4096, -8, -512, 16, 256, -64}.
      - Only m = 64 (lambda = 1, where E_lambda degenerates) is outside it.
    - Corollary 4. "For all primes p > 3", sum_{i=1}^{(p-1)/2} C(2i,i)^3 sum_{j=1}^{i} 1/(i+j) == 0 (mod p). Since
      sum_{j<=i} 1/(i+j) = H_2i - H_i, this is the census twin (H), i.e. classes 3 and 7.
    - Lemma 10. For (p-1)/2 <= k < p and ANY n,
          a_(k+pn) == a_k a_n + p b_k sum_{i<=n} i C(2i,i)^3 (lambda/64)^i  (mod p^2),
      with a_n = sum_{i<=n} C(2i,i)^3 (lambda/64)^i and b_n the Corollary-3 sum.
      - This is HYP-9310's block form.
      - k = (p-1)/2, n = 1 gives "(B1) <=> (L_p)". KLMSY Section 4 states that the supercongruence is equivalent to the
        auxiliary congruence.
      - k = n = p-1 gives HYP-9310's identity F_(p^2) - F_p^2 == p(6h - qF) sum a t_a.
      - k = (p-1)/2 and n = (M-1)/2, with Corollary 3 (b_((p-1)/2) == 0 mod p), give THM-4686's (D):
        F_((Mp-1)/2) == F_((p-1)/2) F_((M-1)/2) (mod p^2) for every odd M.
      - So KLMSY Conjecture 9 at s = 1, all odd m, follows in two lines from KLMSY Lemma 10 + Corollary 3 + Theorem 2.
        KLMSY do not state this.
    - KLMSY derive Corollary 3 from Lemma 10 together with the Coster-van Hamme supercongruences; their proof uses the
      Legendre value at (3p-1)/2 (see n4). KLMSY is in the session's own references (THM-4685 status; HYP-9310
      "related"), and the previous audit (S43 REPORT n9) flagged KLMSY next to the m = 1 case.
    Repair:
    - HYP-9310: retype the main statement and (H) as PROVED/CITED (KLMSY 2016 Cor. 3, Cor. 4). Credit the Lucas-block
      lemma, the block equivalence and the F_(p^2) identity to KLMSY Lemma 10 / Section 4 (the mechanism goes back to
      Gessel 1982 for Apery numbers, as KLMSY say). Keep the runs as FINITE-EXACT records, not as evidence of novelty.
    - THM-4686: retitle and retype as "an alternative proof of KLMSY Corollary 3, from KLMSY Theorem 2 at both ends of
      the 2-isogeny plus the mod-p^2 reflection". That route, and the explicit prime coverage (all p >= 5; p = 7
      direct), is the contribution, and may be new. The Conjecture-9 corollary should read "immediate from KLMSY
      Lemma 10 + Cor. 3 (+ Thm 2); not stated by KLMSY; no later proof found".
    - THM-4685: update the remarks table (classes 3 and 7: PROVED via KLMSY Cor. 4). Note that m = -64 is also
      KLMSY Cor. 3 (lambda = -1), proved here without CM input; m = 64 (lambda = 1) is the case outside KLMSY.
    - Results note: fix sections 3-5 accordingly.
    - Log a MISTAKE entry: cited paper not read; prior audit n9 not followed up.
    - Genuinely open after this: (K); the "only if" of Conjecture (c); the A2-A4 laws (HYP-9311); HYP-9311(O); the
      heuristic.

MINOR

m1. KLMSY's sign misprint is propagated into THM-4686.
    Where: THM-4686 Step 2 (eps_1 = ((z-1)/p), eps_2 = ((1/z - 1)/p)); the Conjecture-9 statement and "Conclusion"
    (((lambda-1)/p) alpha^2).
    What is wrong:
    - KLMSY v3 prints ((lambda-1)/p) in Theorems 1 and 2, Proposition 8 and Conjecture 9. Their own proof gives
      eps^2 = ((1-lambda)/p), and their text compares with Sun's conjectures using ((1-lambda)/p), consistent with
      Ahlgren-Ono-Penniston.
    - Exact check on their model E_lambda: y^2 = (x-1)(x^2 - 1/(1-lambda)), with alpha^2 == a_p^2 - 2p (mod p^2):
      "F(lambda) == ((lambda-1)/p) alpha^2" fails at 217 of 1148 (lambda, p) pairs, 7 < p < 1000. The version with
      ((1-lambda)/p) has 0 failures.
    - Example: lambda = 64, p = 11: a_p = 4, alpha^2 == -6, F == 115 == -6 (mod 121), but ((63)/11) = -1.
    - THM-4686's proof uses only eps_1/eps_2 = (-1)^n (z/p), which is invariant under the correction, so the theorem
      stands.
    Repair: use ((1-z)/p), ((1-1/z)/p), and ((1-lambda)/p) in the Conjecture-9 quotation, with a note "(KLMSY v3
    misprints ((lambda-1)/p))".

m2. THM-4686 Step 2: "The curves at z and 1/z are 2-isogenous, so they have the same Frobenius" is imprecise.
    What is wrong:
    - With x -> x+1, E_z: y^2 = x(x^2 + 2x + 1 - c), c = 1/(1-z). The 2-isogeny with kernel <(0,0)> gives
      Y^2 = X(X^2 - 4X + 4c). X = -2x turns this into -8 * [x(x^2 + 2x + 1/(1-z))], i.e. the twist by -2 of E_(1/z).
    - So a_p(E_z) = (-2/p) a_p(E_(1/z)) (verified for all 7 < p < 1000 and all three pairs). The traces differ by
      (-2/p); alpha^2, which is all the proof needs, agrees.
    Repair: "E_z is 2-isogenous over Q to the quadratic twist by -2 of E_(1/z); hence the same alpha^2".

m3. HYP-9311: the p-shift relation extends to carry indices, so the law is equivalent to the detector test.
    What is wrong:
    - Exact check: A(p+j) == A(1)(A(j) + p[A(j) w_j]) (mod p^2) for ALL 0 <= j < p, carry indices included, with
      [A(j) w_j] the p-integral product. 491 checks per family, p <= 61, 0 failures.
    - Factorial bookkeeping proves it. The extra factor prod C(c_i + e_i, e_i)^(+-1) equals
      1 + sum +-c_i H_(e_i) (mod p) whenever v_p(A(j)) = 1 (checked case by case for the four families), and both
      sides vanish mod p^2 when v_p(A(j)) >= 2.
    - Hence (L_p^A) <=> (B1') whenever p !| A(1) m: truth table over CM and generic m, p <= 200, 2790 rows,
      0 disagreements.
    - Consequently the "pre-registered" law test at the 33 detector survivors was implied, for 7 <= p <= 200, by the
      detector's own confirmation data. Only the primes 211..300 were new.
    Repair:
    - State the equivalence. It strengthens HYP-9311, since the law and the detector are one statement.
    - Reduce the evidential claim of "pre-registered PASSED".
    - Note that the p-shift lemma holds without the carry-free hypothesis at a = 1.

m4. HYP-9311(O): a wrong expected count, and missing context.
    What is wrong:
    - "16 sporadic ... about 21 would be expected": the sum of 1/p over the 1038 supersingular instances is 17.6. 16 is
      consistent with chance, but the stated expectation is wrong.
    - KLMSY already remark that at rational CM lambda their Theorem 2 congruences "appear to hold modulo p^3 for
      ordinary primes".
    - Long-Ramakrishna (Adv. Math. 2016, arXiv:1403.5232, Thm 3) prove the lambda = 1 case mod p^3, which they call a
      conjecture of Kibelbek.
    Repair: correct the number and cite both as prior context. The block statement (O) is different but parallel.

m5. The Beukers-Vlasenko attribution in the HYP-9310 route and the HYP-9311 heuristic goes beyond the paper.
    Where: HYP-9310 "Corollary 8.2 gives ... at fixed points t_0^sigma = t_0, which are the CM points"; results note "the
    CM points are the fixed points"; HYP-9311 "The CM points are these fixed points ... Their Corollary 8.2".
    What is wrong:
    - BV identify CM fixed points of the excellent lift only for the n = 2 hypercubic family, where t(q) is modular.
      For n >= 3 they only conjecture t^sigma(1/2^n) = 1/2^n, the cusp-type point; for n = 3 that is t = 1/8, m = 64.
    - "CM points = fixed points" for the n = 3 family is the session's heuristic, not BV's.
    - In the hypercubic coordinate m = 1/t^2, so t is in Q_p only when m is a square in Q_p. For m = -8, -64, -512 and
      half of the primes the route needs a quadratic extension.
    - HYP-9311's reformulation x_0^(p-1)(1 + p Phi_p(x_0)) == 1 needs F_p to be a unit (ordinary p). At supersingular
      primes Phi_p is not p-integral, although the law still holds.
    - Everything else that was checked in BV is accurately summarised: Def. 7.14 (q^sigma = gamma^(p-1) q^p),
      Thm 7.3, Cor. 6.4 (the log(t^sigma/t^p) theta-term), Cor. 8.2, Conj. 7.5 (mod p^(2s), "unable to prove this for
      any example"), and the n = 3 family 3F2(1/2,1/2,1/2;1,1;64t^2). Prop. 7.17 (a proved mod-p^(2s) supercongruence
      for the polynomials P_Q at the excellent lift) is not cited by the session and does not give the law.
    Repair: attribute precisely, and keep the "fixed point" step labelled HEURISTIC.

m6. Ranges and provenance.
    - THM-4685 status: "every statement for 5 <= p <= 300 (census and forms scripts)". The census output says 5..200;
      only the class-forms run reached 300.
    - cm_law_from_known_supercongruences docstring promises "(2) the derived value of G(1/z) ... from the two
      supercongruence values alone", but the code checks (R) and the law directly.
    - cm_points_forms_legendre docstring says p < 300; the committed run is p <= 200.
    Repair: align the text with the runs.

m7. Census conventions.
    - The script confirms each bucket only against its first member, a possible source of false negatives. The full
      pairwise union-find recomputation found none: identical 7 classes on all 76 primes p <= 400.
    - The docstring promises flags for identical and vanishing classes that are never printed. No product vanishes at
      every half-period, so nothing is lost.
    - "Seven" depends on eps in {1, (-1)^n}:
      - with eps = (2/p) allowed, sum C(n,k)C(n+k,k)C(2k,n)C(n,2k) = [n even] C(n,n/2)C(3n/2,n/2) == (2/p) a_4C(p)
        (a CDE corollary) joins class 1;
      - classes 1, 2 and 4 are one class up to the scalings 4^n, 16^n (THM-4685 (i) itself);
      - allowing 4^(+-n), 16^(+-n) gives 56 classes, mostly through the termwise identity of n1.
    Repair: state the convention in THM-4685 and the results note.

m8. HYP-9310 Lemma: the hypothesis "2a < p" is unnecessary. C(2N,N) = C(2a,a) * (p-adic unit) exactly for every a >= 0
    and j <= n (verified for 0 <= a <= 2p). THM-4686's corollary needs, and correctly states, the general version.

m9. The partial converse is weaker than it reads. In "the converse holds whenever sum_(a<=n) a t_a !== 0 (mod p)", the
    condition fails at every supersingular prime tested at the CM points except m = 64. Example, m = 1: exactly the
    17 primes <= 151 inert in Q(sqrt-7). So F_(p^2) == F_p^2 implies the law only at ordinary primes. At generic m,
    5 of 268 (m,p) rows have F_(p^2) == F_p^2 with the law false. State this.

NOTES
n1. THM-4685 (i-c) needs no Dixon input. Termwise, C(n+k,k) C(2n,n+k) == (-1)^n 16^n C(n,k) (mod p^2) for 0 <= k <= n,
    so (i-c) == (-1)^n 16^n b(n). (i-b) and (i-a) do need Lemma 2(a)/(b). This sharpens the calibration ("what is new
    in (i)" is smaller still).
n2. THM-4685 (ii): the m = -64 law is KLMSY Cor. 3 at lambda = -1, but THM-4685 proves it by pure reflection with no CM
    input, which is worth saying. m = 64 is outside KLMSY.
n3. THM-4686 status lists van Hamme/Ishikawa/Ahlgren as the "input at lambda = 1", but its proof takes m = 64 from
    THM-4685(ii), which does not use them. In the corollary, lambda = 1 is not a CM value of E_lambda (degenerate),
    so Conjecture 9 is not formulated there; (D) still holds (checked).
n4. KLMSY's written proof of Cor. 3 "rewrites a_((3p-1)/2) mod p^2 as P_((3p-1)/2)(sqrt(1-lambda))^2". Taken
    literally this is off: for p <= k <= (3p-1)/2 the pair 2j-1 = p contributes (9p^2 - p^2)/4, a factor -8. The
    combination a_((3p-1)/2) - a_((p-1)/2) a_1 equals -1/8 times the Legendre combination mod p^2, so their conclusion
    stands. THM-4686 supplies a clean derivation of Cor. 3 from Thm 2. Relevant if HYP-9310's "route" is pursued.
n5. The Legendre-square identity is Z.-H. Sun's (1.7) (arXiv:1012.3898), cited by KLMSY. The congruence
    F == P_((p-1)/2)(sqrt(1-lambda))^2 (mod p^2) is in KLMSY Section 3. Cite both. Re-verified: the identity for
    N <= 40, and the congruence at all 8 CM points for p <= 300.
n6. Twin T_m = sum_j C(2j,j) C(p-1-j,j)^2 m^-j: no prior statement found (limited search). Its m = 1 case is
    equivalent to KLMSY Cor. 4 through the session's (S5).
n7. Sun arXiv:2210.07238 is accurately summarised: 17(iii) is mod p^5 with right side -112 p H_(p-1), so zero mod p;
    29(i)-(iv) give the weights 6(H_2k - H_k), 3H_3k + 2H_2k - 5H_k, 4(H_4k - H_k), 3(2H_6k - H_3k - H_k); 36 is
    m = 4096, 42k+5, mod p^3. Sun does not cite KLMSY, which explains the miss.
n8. THM-4685 cites "Sun's conjecture, proved by Z.-H. Sun and by KLMSY". KLMSY Theorem 2 proves all cases of Sun's
    Conjecture 5.2 mod p^2 (modulo the misprint of m1). Z.-H. Sun's abstract covers only some cases mod p^2.
n9. Gukov-Vafa analogy (HYP-9311): RCFT for T^2 needs both tau and the Kaehler modulus rho in the same imaginary
    quadratic field, so "rational exactly at CM moduli" is imprecise. It is labelled ANALOGY.
n10. HYP-9310 "archimedean analogue REFUTED" refutes the session's own naive guess. That is fine as typed; the numbers
     match (|h/F| = 0.001-0.017).
n11. THM-4686's p = 7 caveat is correct.
     - On KLMSY's model Delta = 64 lambda^2/(1-lambda)^3.
     - lambda in {64, 1/64}: bad primes {2,3,7}. lambda in {4, 1/4, -8, -1/8}: {2,3}. lambda = -1: {2}.
     - Any curve over Q with CM by an order in Q(sqrt-7) is bad at 7.
     - The law at p = 7 for m = 1, 4096 holds (direct check).
n12. The law has a crisp equivalent, implicit in THM-4686. For every m with p !| m,
         L(m)  <=>  F_m == (-m/p) F_(4096/m)   (mod p^2).
     Proof: by (R) and z^n == (z/p)(1 + p q_p(z)/2), with (z/p) = (m/p).
     Verified as a truth table: 1074 (m,p) rows, 0 disagreements, p <= 300 (audit_reflection_equiv.out).
     Stating it would make THM-4686's logic transparent: the law is exactly the agreement of the two CM
     supercongruence values across the 2-isogeny.


COORDINATOR'S SPECIFIC QUESTIONS (THM-4686, added mid-audit)
-----------------------------------------------------------
(a) Every algebra step is correct.
    - (R): reindex k -> n-k. a_(n-k)^3 == (-1)^n a_k^3 (1 + 3p Delta_k) mod p^2 gives
      F(z) == (-1)^n z^n [F(1/z) + 3p G(1/z)].
    - z^n == (z/p)(1 + p q_p(z)/2) (mod p^2): put z^n = (z/p)(1 + pu); squaring, 1 + 2pu == z^(p-1) = 1 + p q_p(z),
      so u == q_p(z)/2.
    - Step 3: G(1/z) == -(q_p(z)/6) eps_2 A (mod p), using eps_1 = (-1)^n (z/p) eps_2. This holds because
      ((1/z - 1)/p) = ((1-z)/p)(z/p) and ((1-z)/p) = (-1)^n ((z-1)/p). It holds for any common rescaling of eps_1 and
      eps_2, in particular for the corrected characters of m1.
    - Step 4: H_2k - H_k == -Delta_k + q_p(2) (mod p) by M1-M3, so h(1/z) == (q_p(64z)/6) F(1/z), i.e. L(4096/m).
    - Conversely, given (R), the law at 4096/m is EQUIVALENT to F(z) == (-1)^n (z/p) F(1/z) (mod p^2). So exactly
      "the same alpha^2 at both ends" is needed.
(b) KLMSY Theorem 2 covers lambda in {64, -8, 4, -1, 1/4, -1/8, 1/64}.
    - KLMSY list "lambda = -8, 1, -1/8, 4, 1/4, 64, 1/64, -1 ... all of the CM values for E_lambda over Q", with
      lambda = 1 degenerate.
    - The hypotheses are lambda != 1 CM and good reduction of E_lambda over Z_p. At supersingular p, alpha = 0.
    - The printed right-hand side ((lambda-1)/p) alpha^2 is a misprint for ((1-lambda)/p) alpha^2 (m1).
    - alpha^2 coincides at lambda and 1/lambda: E_lambda is 2-isogenous over Q to the -2 twist of E_(1/lambda) (m2),
      so alpha_lambda = +-alpha_(1/lambda).
(c) The p = 7 caveat is correct and sufficient (n11). Every other prime p >= 5 is good for all seven lambda.
(d) Yes, HYP-9310's main statement is PROVED given the citations.
    - It was in fact already a published statement: KLMSY Cor. 3 for lambda != 1, Cor. 4 for (H), and THM-4685(ii)
      for m = 64.
    - THM-4686 gives an independent, cleaner derivation from KLMSY Thm 2.
    - Remaining open:
      - (K), the class-6 congruence: verified p <= 10^5, no literature found;
      - the "only if" of the CM detector;
      - the A2-A4 laws, which are equivalent to the detector congruence (m3) and have no reflection;
      - HYP-9311(O) mod p^3; prior context in m4;
      - the canonical-coordinate heuristic.
Conjecture 9 at s = 1:
  (i) Wording (KLMSY v3): "For CM values lambda of the family E_lambda: y^2 = (x-1)(x^2 - 1/(1-lambda)), such that
      lambda in Z_p and p is ordinary, for all positive integers m and s with m odd, 3F2(1/2,1/2,1/2;1,1;lambda)
      truncated at (mp^s - 1)/2 == ((lambda-1)/p) alpha_(p,lambda)^2 times the same truncated at (mp^(s-1) - 1)/2
      (mod p^(2s))".
      - THM-4686 quotes it faithfully, including the sign misprint (m1).
      - "The eight rational CM values" should be seven: at lambda = 1 the conjecture is not formulated, though (D)
        holds there.
      - Ordinary p is a hypothesis of the conjecture. (D) also holds at supersingular p, where both sides are
        0 mod p^2.
  (ii) Digits are handled correctly.
      - The p-shift lemma is an exact factorization C(2N,N) = C(2a,a) * unit for every a >= 0 (2a >= p and a >= p
        included). When p | C(2a,a), block a vanishes mod p^3.
      - Terms with last digit n < j < p carry, so their cubes are 0 mod p^3.
      - m^(-a(p-1)) == 1 - apq for every a.
      - Verified for all odd M <= 2p+3 (digits up to p+1), p <= 61, and for M <= 11, p <= 150.
  (iii) Nothing found in later literature.
      - No proof of KLMSY Conjecture 9 (any s) was found.
      - s = 1, all odd m, is a two-line consequence of KLMSY's own Lemma 10 (any n) and Cor. 3 (M1).
      - KLMSY's Prop. 8 (via Coster-van Hamme) proves the Legendre-polynomial analogue for all s; BV Prop. 7.17
        proves a mod-p^(2s) analogue for the hypercubic P_Q polynomials. Neither is the truncated-sum statement.


INDEPENDENT CHECKS (referee code in 04-computation/experiments/s43b_audit_20261018/)
-----------------------------------------------------------------------------------
audit_thm4685.py -> audit_thm4685.out
  [A1] The following, at all 428 primes 5 <= p <= 3000: all hold.
       - S1-S6 mod p^2 for every k (S4: every j <= 2n);
       - M1-M4 and the three consequences mod p;
       - T_n = (-1)^n S_n and THM-4670's odd/even cases;
       - Lemma 2 (a), (b), (c); (d) both as the identity 2F == -3p sum a^3 Delta and as F == 0;
       - (i-a)-(i-d) with their intermediate congruences;
       - the class-2 partner C(n,n/2)^2 [n even];
       - the twin identity T_m == F_m - 2p h(64/m) mod p^2 at m = +-64 and two random m.
       F == 0 mod p^2 at all 217 primes p == 3 (4), and F != 0 mod p^2 at all 211 primes p == 1 (4).
  [A2] The following, at all 781 primes 5 <= p <= 6000: all pass.
       - L(64) and L(-64);
       - Lemma 3 (F and h) at 3 random z per prime;
       - the defect identity delta(m) == -(-1)^n z^n delta(4096/m) mod p at random integer m;
       - q_p(4096/m) == 12 q_p(2) - q_p(m).
  [A3] p = 3: (i-a) == 6, (-1)^n 4^n F == 0 (mod 9). The hypothesis p >= 5 is needed.
audit_census.py -> audit_census.out
  - 935 products (220 + 715). 648 distinct nonzero sequences by identity on n <= 40, the same as by n <= 12. No product
    is identically zero, and no sequences split between n <= 12 and n <= 40.
  - Twin classes by union-find over the full relation, eps in {1, (-1)^n}, all 76 primes 5..400: exactly 7, with the
    same members as the session's.
  - The trivial class (zero at every half-period) is empty.
  - eps in {+-1, +-(-1)^n, +-(2/p), +-(-2/p)}: still 7 classes, class 1 gains [n even]C(n,n/2)C(3n/2,n/2).
  - With the scalings 4^(+-n), 16^(+-n): 56 classes.
audit_forms.py -> audit_forms.out
  - The weight-3 eta-products 4C, 8A, 7AB and eta(2t)^3 eta(6t)^3: a(p) = +(4x^2 - 2p) at every split p
    (p = x^2 + d y^2; x odd for 4C) and 0 at every inert p, 5 <= p <= 1000.
  - F_m == chi(p) a_f(p) mod p^2 at all 8 CM points with the stated forms and characters, every 5 <= p <= 1000:
    1: 7AB; 4096: (-1/p)7AB; -8, 64: 4C; -512: (2/p)4C; -64: (-1/p)8A; 16: E12; 256: (-1/p)E12.
  - Census class values, all p <= 1000:
    - mod p^2: b and sum C(n+k,k)^3 == a_4C; (i-a) == 4^n a_4C; (i-b), (i-c) == 16^n a_4C; sum C(2k,k)^3 and
      T_1 == a_7AB; Apery A == a_4A;
    - mod p: classes 6 and 7;
    - mod p^2 the class-6/7 values are only sporadic (class 7a: p = 7, 53, 89; class 6: p = 97).
  - The j = (m-16)^3/m values and their discriminants (-7, -28, -4, -4, -16, -8, -3, -12) are as stated.
  - 2^3 6^3 is not among the 21 M24 cycle shapes.
audit_cm_law.py -> audit_cm_law.out (method: discrete logarithms, vectorized; different from the session's)
  - (L_p) at all 8 CM points:
    - all 9590 primes 5..10^5: 0 failures, F == 0 mod p at 50.0-50.2% of primes;
    - the 40 primes in [199500, 200000]: 0 failures.
  - (K): 0 failures on the same primes.
  - Generic m: all 195 non-CM 0 < |m| <= 100, p <= 3000: maximum pass fraction 0.049, median 0.028.
  - Direct mod-p^2 twins, all 194 primes 5..1200: classes 3, 6 and 7 hold. The referee's own reduction
    u6a - (-1)^n u6b == p sum_j a_j^3 a_2j K_j (mod p^2) holds, so class 6 <=> (K) with constant 1.
audit_lucas_klmsy.py -> audit_lucas_klmsy.out
  [1] p-shift lemma for 0 <= a <= 2p (2a >= p and a >= p included), 5 <= p <= 31: 3582 checks, 0 failures.
  [2] block_a == t_a (F + ap(6h - qF)) for 0 <= a <= 2p+1, p <= 61, 16 values of m: 0 failures.
      (B1) <=> (L_p): 252 rows, 0 disagreements.
  [3] KLMSY Lemma 10, (p-1)/2 <= k < p, n <= 2p, p <= 23: 8352 checks, 0 failures.
  [4] F_(p^2) - F_p^2 == p(6h - qF) sum_(a<=n) a t_a, p <= 151: 0 failures.
      F_(p^2) == F_p^2 at all 8 CM points, 34/34 primes each.
      sum a t_a == 0 mod p exactly at the supersingular primes (m != 64).
  [5] (D) for all odd M <= 2p+3, 5 <= p <= 61: 4224 CM checks, 0 failures. Generic m: 1581 of 2028 fail.
  [6] THM-4686's reflection (R) mod p^2 at the 7 lambda values, p < 1000: 0 failures.
      KLMSY Theorem 2 as printed: 217 failures out of 1148 pairs. With ((1-lambda)/p): 0 failures.
      a_p(E_lambda) = (-2/p) a_p(E_(1/lambda)), 7 < p < 1000: no failures.
audit_conj9_range.py -> audit_conj9_range.out
  - (D) at the 8 CM points, M = 3..11, all 33 primes 5..150: 0 failures.
  - Generic m = 3, 5, 7, 10 fail at 25-30 of 32-33 primes. The counts per M are identical to the session's.
audit_detector.py -> audit_detector.out (truth tables mod p^2 plus residue lookup, a different method)
  Scan: all 0 < |m| <= 10^5, primes 7..200, then re-confirmation on all primes 7..600. The survivors are identical to
  the session's:
    A1: -512, -64, -8, 1, 16, 64, 256, 4096
    A2: -8640, -1728, -192, -27, 8, 64, 108, 216, 1458, 3375
    A3: -82944, -12288, -3969, -1024, -144, 81, 256, 648, 2304, 20736
    A4: -32768, -3375, 1728, 8000, 54000
  Independent CM signature for EVERY m (supersingular density over the primes 101..1200):
    - density >= 0.30 occurs exactly at the survivors (0.49-0.54);
    - the maximum among non-survivors is 0.082-0.112.
    So the lists are complete in range relative to the CM signature.
  j-invariant CM test: tau = j^-1(j) is imaginary quadratic for every survivor. Discriminants:
    A1: -16, -8, -4, -7, -3, -4, -12, -28
    A2: -75, -51, -27/-3, -15, -8, -11, -12, -24, -48, -60
    A3: -52, -36, -28/-7, -20, -12/-3, -7, -8, -16/-4, -24, -40
    A4: -11, -7, -4, -8, -12, which are exactly the rational CM j with 0 < |j| <= 10^5
  A1 completeness is also algebraic: an integer root of (m-16)^3 = j m divides 4096, and checking +-2^e (e <= 12)
  leaves exactly Sun's 8 values.
check_p5.py -> check_p5.out: (B1') holds at p = 5 for every applicable survivor in all four families.
audit_hyp9311.py -> audit_hyp9311.out
  - Law at the 33 parameters, primes 7..400: A1 600, A2 750, A3 749, A4 375 prime-instances, 0 failures.
  - No non-p-integral term.
  - The p-shift relation holds with carries (m3).
  - (B1') <=> law: 0 disagreements.
  - (O) mod p^3, p <= 300: ordinary 221+275+278+134 = 908, all hold. Supersingular 4+6+4+2 = 16 of 1038, at
    p in {7, 11, 13, 17, 23, 41, 109}.
expected_sporadic.py -> expected_sporadic.out: the sum of 1/p over the 1038 supersingular instances is 17.56.
audit_legendre_square.py -> audit_legendre_square.out: the identity holds for N <= 40, and F_m == P_n(x)^2 mod p^2 at
  the 8 CM points for all p <= 300.
audit_reflection_equiv.py -> audit_reflection_equiv.out: L(m) <=> F_m == (-m/p) F_(4096/m) mod p^2. 1074 rows over
  18 values of m, p <= 300: 0 disagreements. It holds in all 480 CM rows.

E. Pre-registration (git, read-only)
  - HYP-9310 law (cm_harmonic_log_law_20261018.py/.out): docstring, output, HYP-9310, THM-4685 and the (H) results
    were all added in the single commit 8056a055ae (10:07:12). The order of docstring and run is not verifiable.
  - Detector P1/P2: script and output were added together in cf98ae5bce (10:20:34); not verifiable.
    - P1 restates behaviour already observed at the 8 points in 524379c0ed. A1 completeness is also forced
      algebraically (above).
  - HYP-9311 law: script, output and HYP were added together in 7f331c46e5 (10:40:56); not verifiable. It was largely
    implied by the detector data (m3).
  - HYP-9311(O): the prediction-bearing script was committed without output in eb53cfd571 (10:45:22); the output
    followed in 91c22f341d (10:47:18).
    - This is the only verifiable ordering, though it cannot exclude an uncommitted local run.
    - The prediction was informed by the 8-parameter probe committed in the same commit.
  - THM-4686 Conjecture-9 check (e8bc51fca4): no pre-registration was claimed.


Literature consulted
--------------------
- J. Kibelbek, L. Long, K. Moss, B. Sheller, H. Yuan, Supercongruences and complex multiplication, arXiv:1210.4489v3
  (full text read): Thm 1, Thm 2, Cor. 3, Cor. 4, Prop. 5, Thm 7 (Coster-van Hamme), Prop. 8, Conj. 9, Lemma 10,
  Section 4. Published as J. Number Theory 164 (2016) 166-178 (LSU repository record; published PDF not accessible).
- F. Beukers, M. Vlasenko, Dwork crystals III, arXiv:2105.14841v3 (IMRN 2023), full text read: Def. 7.14, Thm 7.3,
  Rem. 7.4, (19), Conj. 7.5, Prop. 7.17, Cor. 6.4, Section 8 (Cor. 8.2, n = 1, 2, 3 remarks), the hypercubic family.
- Z.-W. Sun, Series with summands involving harmonic numbers, arXiv:2210.07238v9 (full text read): Conj. 17, 29, 36.
- Z.-H. Sun, Congruences concerning Legendre polynomials II, arXiv:1012.3898 (abstract; identity (1.7) as quoted by
  KLMSY).
- L. Long, R. Ramakrishna, Some supercongruences occurring in truncated hypergeometric series, arXiv:1403.5232 (ar5iv:
  Thm 3 and the "conjecture of Kibelbek").
- Coster-van Hamme, J. Number Theory 38 (1991) 265-286 (as stated in KLMSY Thm 7).
- Chowla-Dwork-Evans 1986 (used via the formula C((p-1)/2,(p-1)/4)^2 == 2^(p-1)(4a^2 - 2p); not re-read).
- Web searches for later proofs of KLMSY Conjecture 9, for direct statements of (H)/(L_p), and for the twin
  sum C(2j,j)C(p-1-j,j)^2 m^-j found nothing beyond the above.
