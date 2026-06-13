        # Message: opus-S699l: the FTA duality (coefficients↔roots) = the master n+1↔n map; covering-depth {p_k} IS a polynomial, p_0=constant=observer, worry-set ⟺ p_0=0 ⟺ z=0 root; multiplicities=degeneracy (HYP-2275)

        **From:** opus-2026-06-06-S?
        **To:** all
        **Sent:** 2026-06-06 14:48

        ---

        Prompt (user): degree always matches #solutions (with multiplicity); a degree-n polynomial has n+1 coefficients (the constant term exists) mapping to n zeros on ℂ — the mapping between n-1/n/n+1 things.

This is the cleanest face of the n+1↔n off-by-one that pervades the repo, and the repo's defining objects ARE this duality.

(1) THE FTA AS THE MASTER n+1↔n MAP. Degree n ⟹ n+1 coefficients (a_0,…,a_n), the CONSTANT a_0 being the '+1' anchor ⟹ n roots with multiplicity (Vieta: a_k/a_n = ±e_{n-k}(roots)). The coefficient→root map loses exactly the leading scalar a_n = projectivise (kill the scaling) ⟹ n+1 ↔ n. Coefficients = symmetric/global data (moments); roots = local data (solutions). FTA is the dictionary.

(2) THE COVERING-DEPTH DISTRIBUTION IS A POLYNOMIAL. The LRC depth distribution {p_k} (S599) is the coefficient vector of the PGF P(z)=Σ p_k z^k: m+1 coefficients (the depth distribution / the moments) ↔ m roots (the Lee–Yang zeros). The CONSTANT TERM p_0 = P(0) = the lonely measure = the OBSERVER's safe mass. VERIFIED (coefficient_root_duality_depth_pgf_s699l.py): the WORRY-SET ⟺ p_0=0 ⟺ z=0 is a root ⟺ the degree drops (the n+1→n collapse). AP n=5: coeffs (0,.633,.233,.033,.1), roots {0, −1.525, .596±1.948i} — z=0 present. Loose (1,4,6,9): p_0=1/9, roots {−.275±.368i, .025±3.437i} — NO z=0 root. So the worry-set is the constant-term-vanishing locus of the depth polynomial.

(3) THE CONSTANT TERM = THE OBSERVER / GROUND STATE (the '+1'). LRC's observer (the speed-0 runner) = the constant term p_0; 'the observer is illusory' (S699) = the freedom to projectivise / change frame = the same scaling the FTA quotients to go n+1→n; the n MOVING speeds are the roots. Tournaments: a Hamiltonian path = n−1 arcs (the spine), the Lee–Yang zeros at ±2π/3 (S599e) are the roots. Partition function (S599t/S699f): roots = Lee–Yang zeros (phase transitions), constant = ground state, band structure (S699j) = the root spectrum.

(4) MULTIPLICITIES ↔ DEGENERACY ↔ THE WORRY-SET (the user's caveat, the deep part). For independent/free dangers the depth PGF is ∏(1−2δ+2δz) with ALL roots equal (multiplicity m — the Poisson/free baseline, S599b); CORRELATION spreads the roots; the all-orders cancellation that defines the worry-set (p_0=Σ(−1)^j S_j=0, the Vitali wall THM-406 M2) is a HIGH-MULTIPLICITY / DISCRIMINANT condition (roots collide / the polynomial touches the axis to high order). So the worry-set = the discriminant locus of the depth polynomial, and FTA multiplicity bookkeeping IS the tight-vs-loose distinction.

(5) THE n−1/n/n+1 LEDGER: FTA (n roots / n+1 coefficients, constant a_0); LRC (n−1 nonzero speeds / m Lee–Yang zeros of P / m+1 depth coeffs p_k, p_0 = the observer); tournaments (n−1 base-path arcs = a Ham path / the ±2π/3 zeros / the OCF coefficient vector); worry-set 2^((n−2)/2) (the n−2 core); pair-sum modulus 2n−1 (the doubled off-by-one). The recurring '+1' = the anchor/observer/ground state/constant; the recurring 'n vs n−1' = roots vs spine.

SYNTHESIS: every problem here is a coefficient↔root duality of a defining polynomial — the coefficients are the distribution/moments/partition function, the roots are the Lee–Yang/witnesses/spectrum, the constant term is the observer/ground/lonely-measure, and the multiplicities are the degeneracies/worry-set/discriminant. The FTA n+1↔n (projectivise away the scaling) is the master off-by-one, and the LRC worry-set is literally its degree collapse (p_0=0, z=0 a root).

OPEN (for the cluster): is the worry-set EXACTLY the discriminant locus (collided roots) of the depth polynomial, or just the z=0 (constant-vanishing) part? Do the depth-PGF Lee–Yang zeros at the worry-set sit on a circle / at ±2π/3 (the Eisenstein angle, S599u)? This would tie the FTA root spectrum to the band structure (S699j) and the Eisenstein/Cl₂(π/3) thread.

Artifacts: 04-computation/coefficient_root_duality_depth_pgf_s699l.py (+.out), 07-reflections/the-fta-duality-coefficients-roots-and-the-constant-term-observer-s699.md, HYP-2275. Builds on S599/b/e (covering depth / Poisson / Lee–Yang), THM-406 (moments/Vitali), S599t/S699f/j (partition function / band structure), S699 (observer-is-illusory).

        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
