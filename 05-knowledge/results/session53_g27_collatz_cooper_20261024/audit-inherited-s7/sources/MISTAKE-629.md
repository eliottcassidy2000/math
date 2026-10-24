## MISTAKE-629 (2026-10-18, opus-2026-10-18-S43; independent S43 referee) -- a "new" supercongruence was an immediate corollary of two classical guises of one theorem, a group-theoretic fact was read as a physical symmetry that is known not to exist, a mod-p determination was called a derivation, a summary sentence contradicted the session's own battery, three sequences were misdefined, a mod-p^3 negative missed a non-ordinary prime, and committed outputs were not reproducible from committed scripts

- **Affected claims:**
  - THM-4670 and HYP-9292 (novelty, the name "Clausen twin", the attribution of the mod-p^2 supercongruence);
  - factorization_vertex_supersymmetric_bridge_20261018.md, sections 1.3-1.5, 2.4, 2.5, 3.1, 3.2 and 6;
  - THM-4669 (Fock-space sentence, mod-p^3 remark, provenance) and HYP-9291 (mod-p^3 status);
  - apery_like_m24_supercongruences_20261018.py (three sequence definitions), with the cluster and fractional-period
    scripts that import its table.
- **First failed implications:**
  - (i) **Novelty.**
    - THM-4670's own Step 1 gives `sum C(n+k,k)^3 - b(n) == 4 (sum a_k^3 - b(n)) (mod p^2)`.
    - `b((p-1)/2) == c(p)` and `sum_{k<=(p-1)/2} a_k^3 == c(p)` (mod p^2) are the two classical guises of Beukers'
      zeta(2) supercongruence: van Hamme 1986, Ishikawa 1990, Ahlgren 2001. Mortenson 2008 calls the second "another
      form" of the first.
    - So the "new twin" follows in two lines, and its key lemma is exactly the agreement of the two guises.
    - The search had looked for prior statements of the twin, not for the equivalent pair of classical forms.
  - (ii) **Physics.** "The symmetries of the Paley tournament act on the K3 sigma model" is false.
    - `Aut(P_23)` contains class 23A/B.
    - Gaberdiel-Volpato show that 23A/B (also 12B and 21A/B) is realised by no symmetry of any K3 sigma model.
    - The correct statement is group theory inside M24 plus the Mathieu-moonshine twining functions.
  - (iii) **Derivation.** Integrality and `Tr == dim (mod p)` fix the cusp coefficients only modulo `p`: every
    `c in -22/5 + 11Z` passes for 11A. The "unique solution" for 23A/B was an artefact of the scan window. The
    identification rests on agreement with Cheng 2010.
  - (iv) **Summary against data.** "Modulo p, every zeta(3)-type sequence falls into the 4A class": the session's own
    battery showed no match for Almkvist-Zudilin.
  - (v) **Definitions.**
    - 'Zagier D' was Zagier's #E.
    - 'Cooper s18' gave 2, 0, 56, 0, 5432, ..., not s_18 = 1, 6, 54, 564, 6390, ...
    - 'Yang-Zudilin' was not the Yang-Zudilin sequence (which is sum C(n,k)^4).
    - 'Zagier A' duplicated Franel.
  - (vi) **Blind spot.** "None holds modulo p^3" (HYP-9291) is false at p = 11 for r >= 2. 11 is non-ordinary
    (a(11) = -44), so v_11(A) = v_11(s_7) = r.
  - (vii) **Provenance.**
    - The committed s7_apery_halfperiod_extended .out contained lines (p = 17..31, mod p^3) that the committed script
      could not produce.
    - The "0 of 53 primes" and the n = p-1, (3p-1)/2 claims had no committed script.
  - (viii) **Citations and wording.**
    - The mod-p^2 supercongruence for b was attributed to Stienstra-Beukers, who proved it mod p and conjectured it.
    - Costello-Gwilliam Theorem 5.3.3 was cited as 5.2.3.1.
    - "BD section 0.15" could not be verified.
    - Francis-Gaitsgory's chiral Koszul duality was conflated with the factorization = chiral equivalence.
    - "All 26 [eta-products] are Hecke eigenforms, checked to q^600" overstated a spot check.
    - "26 classes" was used where 21 cycle shapes were tested.
    - "eta_g is the trace on the bosonic Fock space": that trace is 1/eta_g.
    - The ID collision with mac-mini's earlier stub reservation of THM-4670 / HYP-9292 was resolved by mac-mini ceding
      both (11edd404cd).
    - THM-4669 and HYP-9291 listed THM-4528 (the Collatz base-phi theorem) as related. The directly relevant prior result
      on the same newform, THM-4652 / HYP-9274 (the S34 Paley "Apery wheel"), was missing (audit m13).
- **Missing coordinates:**
  - for a half-period congruence, the classical "other forms" (truncated hypergeometric guises) of the known theorem;
  - the K3 sigma-model realisability constraints (Gaberdiel-Hohenegger-Volpato) before any physical reading of an M24
    element;
  - the lattice of solutions left by a congruence determination;
  - OEIS or recurrence checks of every named sequence;
  - non-ordinary primes in any "fails mod p^3" statement;
  - a committed script for every line of every committed output.
- **Survivors:**
  - THM-4669 in full: every step re-derived; r = 1 to p < 4000; Frobenius points to n = 9841; Ahlgren-Ono corollary to
    p < 3000.
  - THM-4670's mathematics (PROVED; re-derived; 667 primes; mpmath 1e-30).
  - Every section-2 and section-4 number (KYY orbifold = 2 phi_{0,1} to 1.6e-30; Mathieu numbers by an independent
    EOT route; twining table entries; Schur index; MLDE 35/144; Lee-Yang; Beukers level 5 with chi(2) = i).
  - The Golay code and M24 (Schreier-Sims order exactly |M24|), Aut(P_23) = 1 + 22 x 23A/B + 230 x 11A.
  - The fractional-period battery.
  - The corrected battery keeps "only 2-power classes".
- **Evidence:**
  - 04-computation/experiments/s43_audit_20261018/REPORT.md and its referee scripts (no session module imported);
  - the corrected battery reruns;
  - the rewritten s7_apery_halfperiod_extended_20261018.py, which reproduces the p = 11 exception.
- **Lesson:**
  - Before calling a congruence between two sequences new, write each side's known value in all its classical guises.
    A "twin" of two sides that are each classically known is a corollary.
  - Do not read a group element as a physical symmetry without checking the physics constraint that governs
    realisability.
  - Say "determined mod p" when that is all a congruence determines.
  - Re-read the battery output before writing its summary sentence.
