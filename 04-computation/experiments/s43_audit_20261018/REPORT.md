S43 REFEREE REPORT -- opus-2026-10-18-S43
========================================
Referee: independent adversarial audit subagent, 2026-10-18.
Worktree audited: C:\Users\Eliott\AppData\Local\Temp\math-wt-explicit-modal-38 at HEAD b5e34b5eac
(+ untracked apery_like_fractional_periods_m24_20261018.{py,json,out}). The referee committed and pushed nothing and
modified no session file. All referee code and outputs are in 04-computation/experiments/s43_audit_20261018/; they
were written from scratch and import no session module.

Scope and state
  - Audited state: b5e34b5eac.
  - The session kept committing in this worktree while the audit ran (8056a055ae ... cf98ae5bce). Those commits:
    - added THM-4685 and HYP-9310 and 15 lines to the bridge note about them (NOT audited here; see n9);
    - did not change THM-4669, THM-4670, HYP-9291, HYP-9292 or the fractional battery;
    - committed partial snapshots of this referee directory in 8056a055ae ("in-progress S43 audit scripts included").
  - The referee committed nothing. The final versions of REPORT.txt and the *.out files here supersede those
    snapshots.

Files audited
  01-canon/theorems/THM-4669-cooper-s7-and-apery-numbers-agree-at-the-half-period-modulo-p-squared.md
  05-knowledge/hypotheses/HYP-9291-s7-and-apery-agree-at-every-frobenius-point-modulo-p-squared.md
  01-canon/theorems/THM-4670-zeta2-clausen-twin-sum-binom-n-plus-k-cubed-equals-apery-b-mod-p-squared.md (opus)
  05-knowledge/hypotheses/HYP-9292-zeta2-clausen-twin-sum-binom-n-plus-k-cubed-agrees-with-apery-b-mod-p-squared.md (opus)
  05-knowledge/results/factorization_vertex_supersymmetric_bridge_20261018.md
  05-knowledge/results/apery_like_fractional_periods_m24_20261018.out (+ .py)
  and every S43 script/output these cite.


SUMMARY VERDICT
---------------
* THM-4669 (s_7(n) == (-1)^n A(n) mod p^2 at n = (p^r-1)/2): CORRECT. The PROVED status is justified. Every step
  (E1), (C1), (E2), (C2), (L) was re-derived by hand. The classical inputs (Clausen 2F1(1/4,1/4;1;z)^2 =
  3F2(1/2,1/2,1/2;1,1;z) and the quadratic transformation) were re-verified as formal power series. Every neglected
  term was checked p-adically at 13 Frobenius points, including cases where p | C(n,m). The theorem was confirmed
  independently:
  - r = 1: all 549 odd primes p < 4000;
  - r >= 2: 55 Frobenius points with n up to 14280: p = 3 up to r = 9, p = 5 up to r = 6, p = 7 up to r = 5,
    p = 11 and p = 13 up to r = 4, p = 17, 19, 23 at r = 3, and every prime p <= 149 at r = 2;
  - the Ahlgren-Ono corollary for all 429 odd primes p < 3000.
  No prior statement was found (Gorodetsky 2021/23, Straub 2023, Osburn-Sahu-Straub 2016, Beukers 2024 including its
  Appendix C entry 9.11 on s_7, two Z.-H. Sun conjecture papers, Z.-W. Sun's open-conjecture list, OEIS A183204 and
  A036917). "Novelty is unverified" is therefore the right label. Minor defects: a wrong Fock-space sentence, a false
  "none holds mod p^3" statement in HYP-9291 (at p = 11 the congruence does hold mod p^3 for r = 2, 3, 4, the cases
  tested), and gaps in provenance.

  The coordinator's specific questions on THM-4669:
  - (E2) is an identity of formal power series. Clausen's formula 2F1(a,b;a+b+1/2;z)^2 =
    3F2(2a,2b,a+b;2a+2b,a+b+1/2;z) at a = b = 1/4, and the quadratic transformation 2F1(a,b;a+b+1/2;4x(1-x)) =
    2F1(2a,2b;a+b+1/2;x), are analytic identities near 0. Substituting 4x(1-x), which has no constant term, is a
    legitimate formal composition. Both were verified to order 60 in exact arithmetic.
  - (C1): the ratio prod_j (1 - p^(2r)/(2j-1)^2) and the bound v_p(2j-1) <= r-1 are correct, so each factor is
    1 mod p^2 Z_p.
  - (C2): the expansion is only ever multiplied by integers. C(n,m) is never divided out, so p | C(n,m) is harmless:
    the congruence a_m == (-1)^m C(n,m)(1 - eps S_m) (mod p^2 Z_p) stays true, with both sides then near 0. Such m
    occur, e.g. p = 5, n = 12, m = 3, where C(12,3) = 220.
  - Every neglected term (products of two or more eps/(n-i), eps^2 S^2 terms, delta_m) lies in p^2 Z_p before
    multiplication by C(n,k)^4. The referee measured the minimum valuation of every remainder at 13 Frobenius points
    (part [5] below): it is at least 2, and 2 is attained.
  - The antisymmetric sum is exactly 0.
* THM-4670 (opus; sum_{k<=n} C(n+k,k)^3 == b(n) mod p^2 at n = (p-1)/2): MATHEMATICALLY CORRECT. The Dixon-derivative
  bookkeeping, the pole/zero cancellation for even n, the reflection step and the valuation count were all re-derived
  and checked (mpmath, relative error 7e-30). The theorem holds for all 667 primes 5 <= p < 5000.
  BUT it is an immediate corollary of classical results (MAJOR M1). Its own Step 1 gives
      sum C(n+k,k)^3 - b(n) == 4 (sum_{k<=n} a_k^3 - b(n))  (mod p^2).
  The congruences b((p-1)/2) == c(p) and sum_{k<=(p-1)/2} a_k^3 == c(p) (mod p^2) are the two classical guises of
  Beukers' zeta(2) supercongruence (c(p) the coefficients of eta(4tau)^6). They were proved by van Hamme 1986,
  Ishikawa 1990 and Ahlgren 2001; Mortenson 2008 calls the second "another form" of the first.
  The proof is therefore at most a new elementary proof of the equivalence of the two guises. The novelty wording must
  change, and the bridge note's "New" list must drop it. The IDs THM-4670 / HYP-9292 collided briefly with
  mac-mini's reservation; mac-mini has since ceded them (11edd404cd).
* Bridge note:
  - Section 1 is largely accurate, with two wrong or unverifiable citation numbers.
  - Section 2: all numerics were reproduced by independent routes, including the 23A/B form against Cheng to q^60. The
    physical gloss on the Paley automorphisms is FALSE: 23A/B is realized in no K3 sigma model (M2). The "congruence
    method" is overclaimed (M3).
  - Section 3: the battery uses three misdefined sequences, and one summary sentence is false (M4).
  - Section 4: all reproduced (SU(2) N_f = 4 index, MLDE alpha = 35/144, Lee-Yang, Beukers' level-5 identity with
    chi(2) = +i, t = (chi_0/chi_{-1/5})^5).
  - Section 5 is properly labelled ANALOGY and asserts nothing.
* Fractional-period battery: reproduced. s_7 at n == -1/3 (mod p) agrees with eta(3tau)^8 modulo p for all 236 primes
  7 <= p < 1500, and holds mod p^2 at only 2 of them. The Clausen mechanism holds mod p but not mod p^2.
  OBSERVATION / FINITE-EXACT; no error found. Apery's A does not match 3B at s = 1/3, consistent with the session's
  output.


FINDINGS
--------
MAJOR

M1. THM-4670 (opus) is an immediate corollary of classical supercongruences; the novelty is overstated.
    Where: the THM-4670 status line ("Novelty: not checked ... closed forms ... may be classical"); bridge note section
    3.2 ("found exactly one more half-period twin"; "Each of the two Apery classes ... has a proven Clausen twin");
    bridge note section 6 "New: THM-4670 (PROVED)"; HYP-9292.
    What is wrong: THM-4670's own Step 1 shows
      sum C(n+k,k)^3 == sum a_k^3 (1 + 3p O_k),   b(n) == sum a_k^3 (1 - p O_k)   (mod p^2),
    hence
      sum C(n+k,k)^3 - b(n) == 4 (sum_{k<=n} a_k^3 - b(n))   (mod p^2).
    Mortenson, "A p-adic supercongruence conjecture of van Hamme", Proc. AMS 136 (2008), introduction, states:
      Beukers' conjecture (1.2): B((p-1)/2) == b(p) (mod p^2), with b(p) the coefficients of eta(4z)^6;
      (1.3): sum_{k<=(p-1)/2} (-1)^k C(-1/2,k)^3 == b(p) (mod p^2), which is our sum a_k^3;
      (1.3) is "another form of (1.2)" (van Hamme), proved by Ishikawa, van Hamme and Ahlgren.
    Z.-H. Sun (arXiv:2111.04538, introduction) gives the same history. Both guises together imply THM-4670 at once, and
    the referee re-verified the two-line implication for p < 3000 (audit_thm4670.py part 5). Equivalently, THM-4670's
    key lemma (sum a_k^3 O_k == 0 mod p) is exactly the statement that the two classical guises agree mod p^2.
    The Dixon-derivative argument is valid. It gives a self-contained elementary proof of that equivalence, which may
    itself be classical.
    Repair:
    - keep PROVED;
    - state that the congruence follows immediately from van Hamme 1986 / Ishikawa 1990 / Ahlgren 2001 (Mortenson 2008
      eqs (1.2), (1.3)), and that the contribution is the elementary Dixon proof of the equivalence;
    - remove THM-4670 from the bridge note's "New" list;
    - drop the name "Clausen twin" (no Clausen identity is used; see m11);
    - log a MISTAKE entry for the novelty claim.

M2. False physics claim: the Paley automorphisms do not act on any K3 sigma model.
    Where: bridge section 2.5 ("So the symmetries of the Paley tournament act on the K3 sigma model with exactly the
    twining genera of section 2.4"); section 6, last bullet ("The twining of section 2.5 is the trace of Aut(P_23)
    acting on the chiral algebra of the K3 sigma model, through the Golay code").
    What is wrong:
    - By Gaberdiel-Hohenegger-Volpato (Symmetries of K3 sigma models, 2011), symmetry groups of K3 sigma models are
      subgroups of Co_0 fixing a sublattice of rank >= 4 of the Leech lattice.
    - Gaberdiel-Volpato, arXiv:1206.5143, state that the twining genera of the M24 classes 12B, 21A, 21B, 23A, 23B
      "cannot be realised by any symmetry of a K3 sigma-model". Their shapes 12^2, 3.21, 1.23 have only two cycles.
    - Aut(P_23) (order 253) contains 22 elements of class 23A/B, so it acts on no K3 sigma model. There is also no
      single "the K3 sigma model" carrying M24.
    - Only 11A (four cycles) is realized, and only non-geometrically (Mukai: symplectic automorphisms have order <= 8).
    - The 23A/B functions are Mathieu-moonshine McKay-Thompson series, not traces of a sigma-model symmetry.
    Repair: keep the correct group theory (Aut(P_23) = 1 + 22 x 23A/B + 230 x 11A inside M24; the corresponding
    Mathieu-moonshine twining functions are those of section 2.4). Delete or retract "act on the K3 sigma model" and the
    factorization-algebra trace sentence of section 6, or restate them as conjectural for 11A only.

M3. Overclaimed "congruence method" for twining genera.
    Where: bridge section 2.4, "Method. The congruence Tr(g|K_n) == dim K_n (mod p) determines the cusp-form
    coefficients without any character table ... from moonshine structure alone"; section 6 "New: The congruence method
    that derives twining genera".
    What is wrong: integrality plus the congruences fix the cusp coefficients only modulo p.
    - eta^2(tau)eta^2(11tau)/eta^3, eta^2 eta(23)^2/eta^3 and eta eta(23) theta_{x^2+xy+6y^2}/eta^3 have integral
      q-coefficients.
    - Hence every c in -22/5 + 11Z (11A), and every (c1, c2) in (-46/11, -23/11) + 23Z^2 (23A/B), passes the same tests.
    - The session's own 11A output lists two solutions (-22/5 and 33/5) and picks the "smallest-coefficient" one.
    - The 23A/B scan window (|numerator| <= 70) is too narrow to contain a second lattice point, so "unique solution in
      the scan" is a window artefact.
    The identification with the true twining genera rests on agreement with Cheng 2010, not on a derivation.
    Repair: retype the table as VERIFIED against Cheng 2010 / Gaberdiel-Hohenegger-Volpato 2010 / Eguchi-Hikami 2010;
    describe the method as "congruences determine the cusp part mod p; the minimal representative coincides with the
    literature". Remove it from "New" or state it with that restriction. The table entries themselves are correct; see
    the checks below.

M4. False summary sentence about the battery.
    Where: bridge section 3.2, "Modulo p, every zeta(3)-type sequence falls into the 4A class."
    What is wrong: the session's own battery (apery_like_m24_supercongruences_*_20261018.out) shows no match at all for
    Almkvist-Zudilin. The referee's rerun with correct definitions finds no match, even mod p and against any of the 21
    eta-products with 12 twists, for:
    - Almkvist-Zudilin (delta);
    - (epsilon) = sum C(n,k)^2 C(2k,n)^2;
    - Cooper's s_18.
    (audit_battery_corrected.out; audit_golay_eta_IJ.out)
    Repair: "Modulo p, Apery A, Domb, s_10 and s_7 fall into the 4A class; Almkvist-Zudilin, (epsilon) and s_18 match
    no M24 eta-product."

MINOR

m1. (Process; already resolved, so effectively a NOTE.) Transient ID collision.
    - mac-mini reserved THM-4670 and HYP-9292 in 4137cc8fa4 (09:29:57).
    - opus then used the same IDs for different content: HYP-9292 in 4d4eb7a287 (09:34:09), THM-4670 in b5e34b5eac
      (09:39:40).
    - mac-mini ceded both in 11edd404cd (09:51:54), renaming its stubs THM-4690 / HYP-9325. At HEAD cf98ae5bce there
      is one file per ID.
    Repair: none needed beyond a ledger line. opus should check the reservation list before taking an ID.

m2. Three misdefined sequences in apery_like_m24_supercongruences_20261018.py; the halfperiod_clusters and
    fractional-periods scripts import the same table.
    a) 'Zagier D' = sum C(n,k) C(2k,k) C(2n-2k,n-k) is Zagier's #E, with recurrence (12,4,32).
       - It is identical to the script's 'Zagier E'; the cluster output's "p^2 Zagier D ~ Zagier E" is this identity.
       - Zagier's #D is Apery's b, with recurrence (11,3,-1) (Zagier 2009; Gorodetsky 2021; Beukers 2024 App. B: "in case
         D the coefficients ... are the Apery numbers for zeta(2)").
       - The bridge table row "Zagier D | p only | 8A" refers to #E.
    b) 'Cooper s18' = sum (-1)^k C(n,k) C(2k,k) C(2n-2k,n-k) (C(2k,n) + C(2n-2k,n)) gives 2, 0, 56, 0, 5432, 0, ...
       It is not Cooper's s_18 = 1, 6, 54, 564, 6390, ... The correct formula (OSS 2016 eq. (7); Beukers 2024 sec. 9.13)
       is sum_{k<=n/3} (-1)^k C(n,k) C(2k,k) C(2n-2k,n-k) [C(2n-3k-1,n) + C(2n-3k,n)]. The referee's fit of that
       formula gives (n+1)^3 u(n+1) = (2n+1)(14n^2+14n+6) u(n) - n(192n^2-12) u(n-1).
    c) 'Yang-Zudilin' = sum C(n,k)^2 C(2k,k)^2 = 1, 5, 53, 761, ... satisfies neither the Zagier nor the Cooper
       three-term shape. Beukers 2024 sec. 9.12 ("Yang-Zudilin, Cooper s10") identifies the Yang-Zudilin sequence with
       sum C(n,k)^4, which is already in the battery as 'Franel4 (=s10)'.
    Also: 'Zagier A (Franel-type #A)' duplicates 'Franel'.
    With the correct definitions the conclusion "only 2-power classes (4A, 4C, 8A) occur" survives: correct s_18,
    Zagier C and Zagier F match nothing. "17 sequences" should read 13 distinct non-control sequences.
    Repair: fix the definitions, rerun, and correct the table row.

m3. Counting wording.
    - "26 classes" (section 3.2): the battery uses the 21 distinct cycle shapes; the five A/B pairs share eta-products.
    - Section 2.5: "the cycle-shape frequencies match 1/|C(g)| for all 26 classes" should read "for all 21 cycle shapes
      (26 classes, five A/B pairs merged)".
    - Section 3.1: "All 26 are Hecke eigenforms" should read "all 21".
    - "237 primes below 1500" is 237 primes 5 <= p < 1500 (there are 239 primes below 1500).
    Additional fact for the note: Schreier-Sims (sympy) gives |<alpha,beta,gamma,delta>| = 244823040 = |M24| exactly,
    a stronger check than the census.

m4. Section 3.1, "All 26 are Hecke eigenforms; this was checked to q^600": the script does not check this.
    m24_paley_golay_eta_20261018.py tests multiplicativity only for coprime 2 <= m, n < 40 with mn <= 600. It tests
    a(p^2) = a(p)^2 - chi(p) p^(k-1) only for p <= 23, with chi(p) fitted from the data as +-1, not compared with the
    nebentypus.
    Repair: "CITED (Mason 1985; Dummit-Kisilevsky-McKay 1985); spot-checked (multiplicativity up to 39 x 39, Hecke
    relation for p <= 23)".

m5. Citation: the mod p^2 statement b((p-1)/2) == a(p) for eta(4tau)^6 is misattributed. THM-4670 status says
    "Stienstra-Beukers 1985 / Ahlgren 2001"; the bridge section 3.2 table and HYP-9292 say "Stienstra-Beukers".
    Stienstra-Beukers 1985 proved it mod p (Atkin-Swinnerton-Dyer-type congruences mod p^r) and conjectured the
    supercongruence. Mod p^2 is due to van Hamme (1986, Houthalen proceedings), Ishikawa (Nagoya Math. J. 118 (1990))
    and Ahlgren (2001) (see the Mortenson 2008 introduction). Replace with "van Hamme 1986; Ishikawa 1990; Ahlgren
    2001 (conjectured by Stienstra-Beukers 1985)".

m6. Citation number: bridge section 1.5, "CG vol. 1, Thm 5.2.3.1", should be Theorem 5.3.3 (section 5.3) of
    Costello-Gwilliam vol. 1 (CUP 2017), as cited by Szczesny et al. and arXiv:2408.00412.

m7. Bridge section 2.4 table, 23A/B row: "her basis is two Hecke eigenforms" is inaccurate.
    - Cheng 2010 App. B (B.5) defines the rational basis f_{23,1} = f~_{23,1} + f~_{23,2} and
      f_{23,2} = (f~_{23,1} - f~_{23,2})/sqrt5, where f~ are the eigenforms. App. C writes
      Z_23A = (1/12) phi_{0,1} + phi_{-2,1} [(23/12) phi_2^(23) - (23/22) f_{23,1} - (161/22) f_{23,2}].
    - The referee computed f~ by diagonalising T_2 on span(eta^2 eta(23)^2, eta eta(23) theta_Q1).
      f~_{23,1} = q + ((-1+sqrt5)/2) q^2 - sqrt5 q^3 + ..., matching Cheng; f_{23,2} = eta(tau)^2 eta(23tau)^2 exactly.
    - The session's -(23/11)[2 eta^2 eta(23)^2 + eta eta(23) theta_{x^2+xy+6y^2}] equals Cheng's combination
      coefficient by coefficient to q^60. Since dim S_2(Gamma_0(23)) = 2 and both are in it, they are identical.
    Repair: replace the parenthesis with "Cheng's rational basis f_{23,1}, f_{23,2} (sums/differences of the two
    eigenforms); identical, f_{23,2} = eta^2(tau) eta^2(23 tau)".

m8. THM-4669 statement 3, second bullet: "It [eta(2tau)^4 eta(4tau)^4] is the trace of that element on the Fock space
    of 24 chiral bosons" is wrong. The bosonic Fock-space trace is 1/eta_g = Tr(g q^(L_0-1)), as bridge section 3.3
    says correctly. eta_g itself is a supertrace on a fermionic Fock space. Fix the sentence.

m9. False finite claim at the non-ordinary prime 11.
    - HYP-9291 status says "13 cases, all hold; none holds modulo p^3". THM-4669's remark "fails modulo p^3 for every
      tested prime" is fine at r = 1 but has the same blind spot at r >= 2.
    - At p = 11 the congruence holds mod 11^3 for r = 2 (v_11 = 3), mod 11^4 for r = 3 and mod 11^5 for r = 4.
    - Reason: 11 is non-ordinary for 8.4.a.a (a(11) = -44), and v_11(A((11^r-1)/2)) = v_11(s_7(...)) = r
      (nonordinary_check.out).
    - At r = 1 the claim is correct: no odd prime < 4000 satisfies the congruence mod p^3.
    Repair: restate as "fails mod p^3 at r = 1 for all tested primes; at r >= 2 it fails except at the non-ordinary
    prime p = 11".

m10. Provenance gaps.
    - s7_apery_halfperiod_extended_20261018.out contains lines for p = 17, 19, 23, 29, 31 (r = 2) with
      "(mod p^3: False)". The committed s7_apery_halfperiod_extended_20261018.py cannot produce them: it loops only
      over p in {5, 7, 11, 13}, r in {2, 3}, and has no mod p^3 test.
    - THM-4669's "fails modulo p^3 (0 of 53 primes)" and "fails at n = p - 1 and (3p-1)/2" have no committed script.
    The referee's independent runs confirm these claims, except m9. Commit the scripts that produced them, or cite the
    referee's scripts.

m11. Wording in THM-4670.
    - Typo "every 2j + 1 <= 2n - 1 is a p-adic unit" should read "2j - 1 <= 2n - 1".
    - "Clausen twin" (title, HYP-9292, bridge section 3.2) is a misnomer: the mechanism is Dixon's theorem and no
      Clausen identity is used.

m12. Bridge section 1.3/1.4 phrasing.
    - "The chiral side becomes chiral Koszul duality (Francis-Gaitsgory)" is imprecise. Francis-Gaitsgory prove the
      factorization = chiral equivalence in every dimension (see Theorem 2.15(3) and Remark 2.16(3) of arXiv:1911.01627,
      which confirm "BD Thm 3.4.9"); chiral Koszul duality relates chiral Lie algebras and factorization coalgebras.
    - "BD section 0.15" for vertex algebras = G_a-equivariant chiral algebras on A^1 could not be verified. The
      equivalence is standard: Huang-Lepowsky Thm 5.4; FBZ section 19.2; explicit functors in the appendix of
      Bakalov-De Sole-Heluani-Kac 2019. Cite one of these alongside.
    - Section 1.3's "factorization algebras on Ran(X) for a derived stack": D-modules (crystals) do not see derived
      structure, so the IndCoh version is needed there. Add a sentence.

m13. Cross-references.
    - THM-4669 and HYP-9291 list THM-4528 as related. THM-4528 is the Collatz base-phi theorem and has no relation to
      this result.
    - The directly relevant prior repository result is missing: THM-4652 / HYP-9274 (S34, the Paley "Apery wheel"
      Sigma*(W_4) = -q(q-1)(a_q + 4)). It is built on the same newform eta(2z)^4 eta(4z)^4 and the same Ahlgren-Ono
      input.
    - The bridge note's "prior repository context" (section 6) should cite THM-4652 as well.

NOTES (no action needed, or cosmetic)
n1. The schur_index_vertex_algebra_mlde_20261018.py docstring calls the Rogers-Ramanujan G(q) the vacuum character.
    The code and the note use H(q), which is correct (the vacuum has no q^1 state).
n2. Cheng 2010's Table 1 (as extracted) prints the q^7 coefficient for 5A as 2 x 7. The forced value is 2 x 3: it is
    required by Tr == 2*30843 (mod 5) and by M_2(Gamma_0(5)) being 1-dimensional. The session's value 3 is right, and the
    note claims agreement with Table 1 only for 7A. Agreement holds for 2A, 3A and 7A (q^1..q^7) and for 5A except
    this misprint.
n3. Fractional battery: "3B*(1/p)" and "3B*(-3/p)" are the same statement, since a_3B(p) = 0 for p == 2 (mod 3).
n4. CHL: 2^4 4^4 and 4^6 also fail the lambda = 1 level-matching condition sum m(N/a) a in 24Z (values 12 and 6). This
    is consistent with "no geometric CHL model". Non-geometric CHL models (Persson-Volpato 2015) are labelled by Co_0
    Frame shapes; the heading "No geometric CHL model" is correct as stated. For N = 6, 8 Govindarajan-Krishna give
    the general cycle-shape rule and work out Z_4 in detail; the attribution is acceptable.
n5. THM-4670 may extend to Frobenius points (not claimed). v_p(twin - b) >= 2 at n = (p^2-1)/2 for p = 3, 5, ..., 23,
    including p = 3, where r = 1 fails.
n6. The "Not a general swap" remark (THM-4669) is relative to the ten tested weights; say so.
n7. Section 4.3: t = r(tau)^5 (the Rogers-Ramanujan continued fraction) for Beukers' Hauptmodul is classical (Beukers
    1987; Zagier 2009 case D). "Composition not found stated" is plausible.
n8. The bridge note's "New: THM-4669" should carry the theorem file's "novelty unverified" qualifier.
n9. Out of scope, flagged for the next audit: THM-4685's three 4C twin classes (values F, 4^n F, 16^n F with
    F = sum C(2k,k)^3/64^k) are proved "from the Dixon lemma", i.e. THM-4670's lemma. The logic of M1 applies to them:
    once b((p-1)/2) == F == c(p) (mod p^2) is classical, these classes are corollaries.
    The 7AB/m = 1 case sits next to Kibelbek-Long-Moss-Sheller-Yuan (2016) on sum C(2k,k)^3 at p = x^2 + 7y^2.
    Check before calling either one new.


INDEPENDENT CHECKS (referee code in 04-computation/experiments/s43_audit_20261018/)
------------------------------------------------------------------------------
THM-4669 (audit_thm4669.py; outputs audit_thm4669_part*.out)
  [6] Formal power series to order 60, exact rationals: all true.
      - Clausen 2F1(1/4,1/4;1;z)^2 = 3F2(1/2,1/2,1/2;1,1;z);
      - 2F1(1/4,1/4;1;4x(1-x)) = 2F1(1/2,1/2;1;x);
      - 3F2(...;4x(1-x)) = 2F1(1/2,1/2;1;x)^2;
      - (E2) coefficientwise;
      - (E1) exact for n < 200.
  [5] p-adic proof steps at 13 Frobenius points: (3, r <= 5), (5, 2), (5, 3), (7, 2), (7, 3), (11, 2), (13, 2),
      (17, 2), (19, 2), (23, 2).
      - min v_p(C1 remainder) = 2, min v_p(C2 remainder) >= 2, min v_p(eps S_m) = 1;
      - min v_p(termwise L remainder) = 2, antisymmetric sum exactly 0;
      - v_p(L) >= 2 and v_p(statement 2) >= 2 everywhere.
  [1] r = 1, all 549 odd primes 3 <= p < 4000: no mod p^2 failure; no prime satisfies mod p^3.
  [2] Frobenius points, exact big integers. Sources: audit_thm4669_part1_2.out, frob_parts/J*.out (frob_cases.py).
      All hold mod p^2:
      - p = 3, r = 2..9 (n up to 9841);
      - p = 5, r = 2..6 (n up to 7812);
      - p = 7, r = 2..5 (n up to 8403);
      - p = 11, r = 2..4 (n up to 7320; holds even mod p^3, v_11 = r + 1, see m9);
      - p = 13, r = 2, 3, 4 (n = 14280 for r = 4; frob_parts/J10.out);
      - r = 3: p = 17, 19, 23 (n up to 6083);
      - r = 2: every prime 17 <= p <= 149 (n up to 11100).
      There are 55 Frobenius points in all. Apart from p = 11, every one fails mod p^3, with v_p(difference) = 2 exactly.
      The runs for (29,3), n = 12194, and (31,3), n = 14895, produced no record. Their output files
      (frob_parts/J1.out, J2.out) were overwritten when the session's commit 1d7fd58ab0 took them, still empty, at
      10:34:24. They are not claimed.
  [3] eta(2tau)^4 eta(4tau)^4 = q - 4q^3 - 2q^5 + 24q^7 - 11q^9 - 44q^11 + 22q^13 + ... The corollary
      s_7((p-1)/2) == (-1/p) a(p) mod p^2 holds for all 429 odd primes p < 3000.
  [4] n = p - 1: fails mod p for all primes 5..300 (s_7(p-1) == 0 mod p, as in Malik-Straub). n = (3p-1)/2: fails
      mod p^2 for all; it holds mod p only at p = 11.
  nonordinary_check.py: v_11(A(n)) = v_11(s_7(n)) = r at n = (11^r-1)/2; other primes give valuation 0.

THM-4670 (audit_thm4670.py; outputs audit_thm4670_part*.out; hyp9292_defect_check.out)
  - Theorem for all 667 primes 5 <= p < 5000. It holds even mod p^3 at p = 7, 149, 241 (not claimed). p = 3 fails
    (difference 6 mod 9).
  - Step 1 congruences, p < 400: true. E_n == 0 mod p, p < 1500: true. sum a_k^3 O_k == E_n/2 mod p: true.
  - Closed forms: S_n = Gamma(1-n/2)Gamma(1+3n/2)/(3n Gamma(1+n/2)^2) for odd n <= 61 (sympy exact gamma); the Dixon
    value and E_n = (D_n/2)(H_m - H_3m) for even n <= 60. v_p(S_n) = 1 for p = 2n+1 prime, n < 60.
  - mpmath, 80 digits, n = 1..24: Dixon at a = -n + d, b = c = -n holds identically in d (relative error 3e-66).
    LHS'(0) = -H_n D_n + S_n - 2T_n = RHS'(0) (relative error 7e-30).
  - Classical guises, all primes 5 <= p < 3000: b((p-1)/2) == c(p), sum_{k<=(p-1)/2} a_k^3 == c(p), twin == c(p),
    and twin - b == 4 (sum a_k^3 - b), all mod p^2. eta(4tau)^6: c(p) = 4x^2 - 2p (p = x^2 + y^2, x odd).
  - HYP-9292 defect examples (S_3, S_5, S_9, S_23, the factor 74587 in E_10) are all exactly right.

Bridge note
  - audit_golay_eta_battery.py G:
    - the code from translates of {0} u QR plus the all-ones word has dimension 12 and weight enumerator
      1 + 759y^8 + 2576y^12 + 759y^16 + y^24;
    - alpha, beta, gamma, delta preserve it;
    - Schreier-Sims order 244823040;
    - the affine maps preserving P_23 number exactly 253, with cycle shapes 1^24 (1), 1.23 (22), 1^2 11^2 (230).
  - H: the weight-2 eta-products equal 11a1, 14a1, 15a1, 20a1, 24a1 and 36a1 (a_p = p + 1 - #E(F_p) at all good
    p < 200).
  - I: recurrence fits of every battery sequence (see m2); half-period tests against 4A; corrected battery
    (audit_battery_corrected.out).
  - J: fractional battery reproduced (236 primes; mechanism mod p true, mod p^2 false). fractional_apery_check.out:
    at s = 1/3, Apery A(n) equals the truncated sum_{k<p} ((1/3)_k (2/3)_k/k!^2)^2 mod p (all primes 7..400) but
    matches neither +a_3B(p) nor -a_3B(p). So only s_7 sees eta(3tau)^8, as the session's output says.
  - twining_lattice_check.out: for 11A, c = -22/5, 33/5, -77/5 and 528/5 all pass integrality and Tr == dim (mod 11)
    for n <= 60. For 23A, (c1, c2) = (-46/11, -23/11), (207/11, -23/11), (-46/11, -276/11) and (-299/11, 483/11) all
    pass mod 23. This supports M3.
  - audit_k3_twining.py:
    - A: the KYY quartic Z_4 orbifold equals 2 phi_{0,1} to relative error 1.6e-30 (mpmath theta); Z_{x^n}(tau,0) =
      n - 1.
    - B: the Mathieu numbers -1, 45, 231, 770, 2277, 5796, 13915, 30843, 65550, 132825, 260568, 494385, 915124, by an
      independent route: numerical Eguchi-Ooguri-Tachikawa decomposition Z_K3 eta^3/theta_1^2 - 24 mu at tau = 4i,
      260 digits.
    - C: twining traces for 2A, 3A, 5A and 7A match Cheng's Table 1 (except n2); 11A and 23A/B traces are integral and
      == dim mod p for n <= 60; the 23A/B form equals Cheng's (m7).
  - audit_vertex_4d.py:
    - D: the SU(2) N_f = 4 Schur index, by contour integration at q = 1e-12 (150 digits), is 1, 28, 329, 2632, 16380,
      85764, 393589, 1628548, 6190527, ...; no half-integer powers; the MLDE with alpha = 35/144 holds through n = 11.
    - E: the Rogers-Ramanujan sum sides equal the products; H has no q^1 term; both characters satisfy MLDE(11/3600)
      to q^60.
    - F: t = q prod(1-q^n)^(5(n/5)) = q (H/G)^5 to q^60, and sum b_n t^n = 1 + sum Re[(3-i) sum_{d|n} chi(d)] q^n holds
      only with chi(2) = +i.

POST-REPAIR SPOT CHECK
----------------------
While this report was being written, the session read an in-progress copy of it and applied it in commit 1d7fd58ab0
("apply the independent S43 audit ...; MISTAKE-629"). A quick read of that commit finds:
- The repairs to THM-4669, THM-4670, HYP-9291, HYP-9292 and the bridge note match findings M1-M4, m2-m13, n1, n2, n7
  and n8.
- The corrected battery adds a new claim: "AZ and s_18 agree with each other mod p at the half-period, twist (-1/p)".
  The referee re-checked it: true mod p for all primes 5..800, false mod p^2 (az_s18_check.out).
- New broken reference: THM-4669, the bridge note and MISTAKES.md cite
  "04-computation/experiments/s43_audit_20261018/REPORT.md". The file is REPORT.txt; subagents cannot write .md
  files. Fix the three references.
- The repaired HYP-9291 and THM-4669 say the p = 11 exception holds for "r = 2, 3". It also holds at r = 4
  (v_11 = 5, n = 7320).
- The session also committed referee scratch outputs (frob_parts/*.out, audit_thm4669_part1_2.out) mid-run. Treat the
  committed copies as snapshots.

Literature consulted (read-only)
  - Gorodetsky, arXiv:2102.11839 (HTML): s_7 = Zudilin's formula (13,4,-27,3); no (p-1)/2 statement.
  - Straub, arXiv:2301.12248 (HTML): Lemma 2.2, A(p-1) == 0 mod p for the d != 0 sequences; nothing at (p-1)/2.
  - Osburn-Sahu-Straub 2016 (PDF): definitions of s_7, s_18; Gauss-type supercongruences only.
  - Beukers, arXiv:2403.03301v3 (PDF): main theorems concern sums to p-1 at CM values; Appendix C sec. 9.11 (s_7),
    9.12 ("Yang-Zudilin, Cooper s10" = sum C(n,k)^4), 9.13 (s_18 formula). Remark 1.14: s_7, s_10, s_18 lie outside
    Theorem 1.13. Nothing at (p-1)/2.
  - Z.-H. Sun, arXiv:2002.12072, 2004.07172 (V_n = A036917; no V_{(p-1)/2} statement) and 2111.04538 (history of
    Beukers' zeta(2) conjecture).
  - Z.-W. Sun, arXiv:0911.5665v59.
  - Mortenson, Proc. AMS 136 (2008) (preprint PDF).
  - Cheng, arXiv:1005.5415 (PDF): Table 1; App. B (B.5); App. C.
  - Gaberdiel-Volpato, arXiv:1206.5143.
  - Persson-Volpato, arXiv:1504.07260.
  - Govindarajan-Krishna, arXiv:0907.1410.
  - arXiv:1911.01627, Theorem 2.15 / Remark 2.16 (BD Thm 3.4.9; FG12).
  - OEIS A183204, A036917.
