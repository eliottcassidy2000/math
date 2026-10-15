        # Message: opus-2026-10-14-S35: Paley run ladder to 12 steps; gcd(q+1,128) | b(q) for the K6-P6 newform; tournament 1/3-2/3 conjecture (verified n <= 8); audited MISTAKE-613

        **From:** opus-2026-10-14-S?
        **To:** all
        **Sent:** 2026-10-14 20:31

        ---

        opus-2026-10-14-S35: the Paley run ladder through 12 steps, the 128 inside the level-32 K6-P6 newform, and a tournament 1/3-2/3 conjecture (verified to n = 8). Audited, MISTAKE-613.

Owner's briefs:
1. Fractal least-resistance recursion behind 2q(q-1)(q-3)(q-7)(7q-55): triangular-Fibonacci, F=U+S, mod-8 Bott, Mahler 3/2, the two-point Chowla paper.
2. Look back at 128 across sessions and recurse a few levels.
3. Pasted with it: a genuine 1/3-2/3 and gold-partition attempt, Borsuk-4, hexagonmath as solid, and the golden-lane snippet.

Notes:
- 05-knowledge/results/paley_least_resistance_runs_one_modular_form_20261014.md
- 05-knowledge/results/one_twenty_eight_atlas_balance_gold_20261014.md
Hypotheses: HYP-9275 (updated), HYP-9277, HYP-9278.
Audit: 04-computation/experiments/paley_s35_audit_20261014/REPORT.md, with independent scripts.

RESULTS
- A_12 = 4q(q-1)(33q^5 - 1540q^4 + 28243q^3 - 249550q^2 + 1027794q - 1488872) + 28q(q-1)(17q - 251) a_q, with a = eta(2z)^4 eta(4z)^4. Confirmed by the audit's independent signed subset DP, and out of sample at q = 27.
  - Through 12 steps the only cores are the wheel and the octahedron.
  - The leading terms follow THM-438's Catalan law; this adds the exact cases k = 5, 6.
- K6-P6: mac-mini's HYP-9276/THM-4653 came first. I reached the same square law independently, in the unrestricted normalisation, from an exact table at 124 primes, with b rebuilt from eta quotients.
- 2-adic structure of b, FINITE-EXACT to 524,288 via the audit's exact theta-product tabulation of b:
  - |b(q)| = +-(q+1) mod 2^m, with m = 4, 5, 6, 6 for q = 11, 3, 7, 15 mod 16; this gives the 2^7/2^10 laws.
  - gcd(q+1, 128) | b(q), and the exponent 7 is sharp. It also holds at 580 primes q = -1 mod 512 up to 2e6, through the perfect-square test.
  - Zeros: 127 and 8191 (Mersenne) and 213247 (not Mersenne).
- Wheels: Tr((DS)^n) = sum_psi (chi(-1)J(psi,chi))^n is Greene's nF(n-1)(1) (standard). T_4 and T_6 are Ahlgren-Ono and Frechette-Ono-Papanikolas (CITED).
  - New: T_8 is NOT in the span of S_8(Gamma_0(32)) plus twisted lower weights (39 primes, rank 17 vs 18).
- Planar Fourier duality (a specialisation of Kramers-Wannier / Wu-Wang) gives U(cube) = q U(octahedron). This refutes my pre-registered "Laman-sparse cores are polynomial".
  - HYP-9278: new motives sit at 2D rigidity circuits, with weight |V| - 2 (an organising hypothesis).
- Tournament 1/3-2/3 (HYP-9277, a CONJECTURE).
  - Verified for all 2^28 labelled tournaments at n = 8 and below.
  - Equality holds exactly at transitive compositions of singletons and C_3. These are A000930(n) - 1 classes, 128 at n = 14, matching Gupta's 128 posets.
  - The least strong balance is 2/5 (n = 4), 4/11 (n = 5..7) and 62/157 (n = 8).
  - REFUTED: "every arc tail-first >= 1/3" (n = 8, 13/45). Also refuted: the gold partition on Hamiltonian paths with precedence answers (n = 6).
  - The audit found that every strong tournament with n <= 7 has a balanced arc.
- 128 atlas:
  - six corpus themes (2-adic order of 3, Collatz residues mod 128, Paley/127, the LRC14 F_2^7 cube, Fano/Hamming/E8, AMM R = 128) plus two new;
  - QR codes of the Paley primes 7, 11, 23;
  - the Collatz sieve, 128 at k = 11 (A076227);
  - Mahler's middle third (AFS Cor. 4);
  - Borsuk-4 status only;
  - cross-references to tiling-modular-atoms' same-day gold/Borsuk packets.

HONEST LIMITS
- No proof of 1/3-2/3, the gold partition, Borsuk-4 or Mahler.
- The 2-adic congruences of b are unexplained.
- The rigidity dictionary has already lost two pre-registered clauses: Fano => polynomial, and W_8 at level 16/32.

NEXT
- A_14 via a compiled signed subset DP at enough primes.
- Identify T_8 (nebentypus chi_8 / level 64?).
- Prove gcd(q+1,128) | b(q), perhaps via the descent form.
- HYP-9277 at n = 9, and the balanced-arc question.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
