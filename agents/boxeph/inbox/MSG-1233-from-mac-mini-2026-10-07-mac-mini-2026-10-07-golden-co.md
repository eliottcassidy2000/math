        # Message: mac-mini-2026-10-07-golden: Collatz classes equidistributed (THM-4590); musical resonance spectrum (HYP-9230); cycles to period 301,993 = five Pythagorean (THM-4591); eleven squares = period-5 cycles (THM-4592); maximal sieve + 2-adic sign barrier (THM-4594); Collatz OPEN

        **From:** mac-mini-2026-10-07-S?
        **To:** all
        **Sent:** 2026-10-07 14:28

        ---

        mac-mini-2026-10-07-golden close-out.

Owner prompt: "finish a Collatz proof through creative means" — golden ratio, base φ, 223/233/332/425/105; {2,3,11}; the Platonic solids against Ellison's exceptions; mod 18/19 against the 18/19/101 counts. Inspiration: the eleven-squares PDF and the CrocSwap integer-multiplication result, treated as trusted.

Note: 05-knowledge/results/golden_collatz_resonance_20261007.md.

COLLATZ IS STILL OPEN. Results, after two independent audits (corrections in MISTAKE-584 and MISTAKE-585):
- THM-4590 (PROVED from THM-4581): every union of grand-orbit classes has o(x) cuts n|n+1. It is equidistributed mod every M, uniform across multiplicative windows, slowly varying, and almost affine-invariant.
  - So a counterexample's basin cannot hide in residues or leading digits.
  - Negative-integer test: the basins of -1, -5 and -17 have densities 0.327, 0.325 and 0.348.
- HYP-9230 (NUMERICAL): Collatz moves by fifths up and octaves down. Class densities have log-periodic modes at small ||f log2 3||:
  - 12-, 41- and 53-TET; the semiconvergents 94, 147, 200, 253, 359; harmonics such as 106 = 2*53, which is dominant for k >= 20.
  - Decay is the exact root of the fifth/octave walk. It matches every measurable mode, which resolves the factor-2 puzzle (audit A2).
  - Ellison's exceptions share the continued fraction: a common cause.
  - Prior art: Berg-Krueppel 1998 and Tavares 2026 (predecessor density).
- THM-4591 (FINITE-EXACT, two methods, independently re-run): the integer Terras cycles of period <= 301,993 are exactly the five Pythagorean ones.
  - Octave, fifth, fourth and whole tone (gap 1, THM-4484), plus the apotome 2187:2048 (the -17 cycle). Periods {1,2,3,11}.
- THM-4592: the eleven-squares Fibonacci torus is the period-5 Collatz cycle set: 0, the 1/13 cycle (QR mod 11) and the -5 cycle (non-residues).
  - G_5 = BS(1,3) mod 11.
  - For (mx+1)/2, a coalescence weight with s < 1 exists iff m < 4.
- THM-4593: witnessed finite lag sets (R <= 32; every nonempty set of odd lags D <= 61) have exponent 1/2, because partners coalesce first. All lags at once give <= (T+1) 2^-0.0346T.
  - HYP-9217's any-lag 0.69 was a coalescence transient.
- THM-4594: the maximal class-decided sieve leaves 0.493% mod 2^34 (descent: 0.884%), but beats Angeltveit's published 2-adic rules by only 4.5% (0.957x at 2^30).
  - Given Barina's 2^71, a minimal counterexample lies in 6,915,181 classes mod 2^30. The first exclusions beyond published rules are 11247, 12191, 12799, 23743 mod 2^15.
  - 2-ADIC SIGN BARRIER: for unrefined classes, -1 is never certified, and -5 and -17 not to depth 120. 3-adic refinement certifies all three at depth 0. No bounded-memory weight and no pathwise chain weight exists.
- HYP-9231: survivor lemma (FINITE-EXACT s <= 30).
- NUMEROLOGY verdicts:
  - Platonic solids (p = 0.78);
  - the golden features of 332 -> 233 -> 425 <- 223, which all lie on 27's trunk (the real content is the (334,335) coalescence);
  - 105's orbit;
  - base-phi binary readings (killed);
  - the 18/19 counts (the Borwein-Choi characterization is KNOWN, Borwein-Choi 2000);
  - 3^5 = 1 + 2*11^2 is real number theory (Ljunggren / Golay / Wieferich) but not Collatz.

REMAINING GAP: Collatz <=> every n > 2^71 has a class-decided certificate whose threshold is below n. ("U_inf ∩ Z>0 = ∅" is weaker, because cycles sit at thresholds.)

NEXT:
1. Pathwise, Baker-type control of k log2 3 mod 1 along one orbit (the musical bridge), combined with 3-adic data.
2. Prove HYP-9231.
3. Prove the HYP-9230 transfer law for class densities.
4. Unbounded-lag exponent.
5. Upgrade THM-4581 (4) from sketch level.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
