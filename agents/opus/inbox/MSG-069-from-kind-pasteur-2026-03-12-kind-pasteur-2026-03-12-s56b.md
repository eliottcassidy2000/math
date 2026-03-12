        # Message: kind-pasteur-2026-03-12-S56b: SPECTRAL FLATNESS REVERSAL — Paley vs Dirichlet kernel + D_{2p} representation theory

        **From:** kind-pasteur-2026-03-12-S?
        **To:** all
        **Sent:** 2026-03-12 13:34

        ---

        MAJOR SESSION — Deep dihedral group analysis of Paley H-maximization (user-directed).

BREAKTHROUGH FINDING — SPECTRAL FLATNESS REVERSES AT p=1 mod 4:

1. p=3 mod 4 (Paley exists): Flat spectrum = MAX H.
   - Paley achieves perfect flatness: all |lambda_k|^2 = (p+1)/4 (Gauss sum imaginary).
   - Correlation(H, spread) = -1.0 at p=7. Distance from spectral center = 0.
   - Paley orbit has size 2 in Z_p^* (stabilizer = QR_p).

2. p=1 mod 4 (NO Paley): CONCENTRATED spectrum = MAX H.
   - Cyclic interval has Dirichlet kernel spectrum: |lambda_k| = |sin(pi*k*m/p)/sin(pi*k/p)|.
   - One dominant eigenvalue (|lambda_1| = 4.15 vs others < 1.41 at p=13).
   - Correlation(H, dist_from_center) = +0.46 at p=13 (POSITIVE — farther = better!).
   - Cyclic interval orbit has size p-1 in Z_p^* (trivial stabilizer).

NEW THEOREMS/HYPOTHESES (HYP-457 through HYP-463):
- HYP-457: SPECTRAL DETERMINATION — H determined by eigenvalue multiset {|lambda_k|^2}
- HYP-458: SPECTRAL REVERSAL — flat wins at p=3(4), concentrated wins at p=1(4)
- HYP-459: ANTI-AUT PRESERVATION — reflection v->-v fixes all eigenvalues (Re(lambda_k)=-1/2)
- HYP-460: GAUSS SUM DICHOTOMY — g purely imaginary iff p=3 mod 4
- HYP-461: ORBIT SIZES — Paley (size 2), Satake (size 4), generic (size p-1)
- HYP-462: DIRICHLET KERNEL — cyclic interval = Dirichlet kernel (verified exactly at p=13)
- HYP-463: VARIATIONAL PRINCIPLE (CONJECTURE) — H max at center of spectral simplex

KEY QUANTITATIVE RESULTS:
- p=7: 2 Z_p^* orbits, 2 spectral classes, 2 H values (189, 175)
- p=11: 4 orbits, 4 spectral classes, 4 H values
- p=13: 6 orbits, 6 spectral classes, 6 H values (3711175 down to 3669497)
- Sum_y2 = p(p-1)/8 is CONSTANT across all circulant tournaments (Parseval)
- At p=7: 2 power sums determine H. At p=13: 4 power sums needed.
- H determined by (Sigma_2, Sigma_3, Sigma_4) at p=13

OPEN QUESTIONS FOR NEXT AGENT:
1. PROVE the Spectral Variational Principle (HYP-463) — or find counterexample at larger p
2. Test at p=17 (=1 mod 8) and p=19 (=3 mod 4) to confirm pattern
3. Find explicit polynomial H(Sigma_2, Sigma_3, ...) for small p
4. Why does the H-landscape flip? Algebraic explanation via Gauss/Jacobi sums?
5. Connection to representation ring of D_{2p} — is H a character?

NEW FILES: dihedral_spectral_analysis.py, spectral_orbit_analysis.py + outputs

        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
