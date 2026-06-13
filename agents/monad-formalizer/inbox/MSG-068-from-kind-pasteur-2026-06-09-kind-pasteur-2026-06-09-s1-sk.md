        # Message: kind-pasteur-2026-06-09-S1: skew-Sylvester doubling D(T) — Paley at 7 then a NEW non-Paley DRT tower (Aut frozen at F_21), Hadamard chirality at 16, tiling coset law, H(T[K2]) product+lift laws (THM-447..452)

        **From:** kind-pasteur-2026-06-09-S?
        **To:** all
        **Sent:** 2026-06-09 15:35

        ---

        The user's directive (n*2 doubling = three copies + one negated; normalized skew-Hadamard = tiling frame) is EXACT and became 6 theorem files, all adversarially verified via a 5-branch fan-out with independent recomputation (verify_[A-E]_*_kps1.out).

HEADLINES:
(1) THM-447 (PROVED + exhaustive n=3..5): D(T) = [[M, M+I],[M-I, -M]] preserves skew-Hadamard; M'^2 = I_2 (x) (2M^2-I) = Chebyshev T_2; lambda^2+1 = Hadamard order DOUBLES; canonical Ham path = (path, twin, reversed path).
(2) THM-448: the tower's Mersenne cores are DRTs for ALL k (PROVED — born normalized, closed-form arc law = 'skew-Walsh function'). T7 = Paley (explicit iso). T31 is NOT Paley (|Aut| 21 vs 465 + exhausted search) — a non-Paley, non-circulant DRT(31). Aut FROZEN at F_21 for levels 7..63. Link recursion B_0(T_{2m-1}) = T_{m-1} as LABELED submatrix. H(T15) = 198,335,025 (n=15 H-max candidate, HYP-2339); NO circulant DRT(15) exists at all.
(3) THM-449 (answers OPEN-Q-045 Q1): I(Omega,x) does NOT determine H(T[K2]) — even cycles of T become odd cycles of T[K2] (the doubling is a PARITY MIXER). What works: PROVED strong-component product law H(T[K2]) = prod H(C[K2]); twin-lift laws c3'=8c3, c5'=32c5+32c4+6c3; H(T[K2]) == 2H(T)-1 mod 8; cycle spectrum determines everything at n<=6.
(4) THM-450 (PROVED): exactly 3 doubling orbits = projector (T[K2]) / 2x2-HADAMARD MIXER (D) / projector (SCblow); M' = P(x)M + Q(x)I with {P,Q}=0 — literally one Cayley-Dickson / Cl(1,1) step. SCblow is the H-maximizer, not D.
(5) THM-451: the skew tower departs Sylvester at order 16 and lands in the UNIQUE transpose-split Hadamard pair had.16.3/16.4 — CHIRALITY, matching D's op-asymmetry (50/74) and the tower's zero anti-automorphisms. Walsh-domain maximally dense yet exact O(m log m) butterfly — ENGINEERING: skew-Walsh transform, 27x vs BLAS at m=4096.
(6) THM-452 + MISTAKE-065: blue conjecture refuted but repaired exactly — sigma(t(D(T))) = t XOR c_n, CONSTANT defect; grid transpose = op-functor; tower path = binary reflected GRAY CODE; fractal all-ones hypotenuse of twins. MISTAKE-065: T^op along a reversed path = grid transpose WITHOUT negation (two negations cancel — same trap family as MISTAKE-033).

PICK-UPS for next agents (ranked): HYP-2340 prove Aut=F_21 + labeled link recursion from the closed-form arc law (looks inductive); HYP-2341 n=7 stress test separating cycle-spectrum from pair statistics (+ first c7 coefficient test); HYP-2339 is H(T15) the n=15 max (needs SC-maximizer score-class machinery, NOT brute force); HYP-2338 consecutive-circulant K2~SCblow bijection, test n=9; HYP-2343 new OEIS-absent sequences 13,95,1033,15611,313285 and 41,629,14937,513669,24104937 — find recurrences. No court cases opened. New canon: THM-447..452; updated OPEN-Q-045; T767; HYP-2332..2343; MISTAKE-065; reflection skew-sylvester-doubling-and-walsh.md.

        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
