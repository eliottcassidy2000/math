# Audit B report (mac-mini-2026-10-14-wilson, second half)

Independent adversarial audit of THM-4650, THM-4651, d19_parity_shift_synthesis_20261014.md and the HN additions
(findings 5-7 of hn_type_m_hinge_rotations_and_forced_radii_20261014.md, README items 8-10). Run by a separate agent;
this file condenses its returned report. The corrections are recorded as MISTAKE-610 and applied in the audited files.

## Verdicts

| File | Verdict | Notes |
|---|---|---|
| THM-4650 | CORRECT WITH FIXES | mathematics checks; prior art missing; one torsion group wrong |
| THM-4651 | CORRECT WITH FIXES | computations reproduce; one false statement (N(N+1) Giuga); main reading has prior art |
| synthesis note | CORRECT WITH FIXES | |
| HN additions | CORRECT WITH FIXES | the hinge-edge description was wrong |

## Major findings

- **M1.** THM-4650 statements 1-3 are Yang Ji, arXiv:2105.05250 (Fermat descent), plus Ieno, arXiv:2111.02888. Bremner-Ulas, arXiv:1502.07312, overlaps statement 5. MacLeod (SWJPAM 2 (2003) 9-17) is the base/height-ratio problem. New: the circumcircle (Ptolemy) and the explicit curve identifications.
- **M2.** Only N(N - 1) is a Giuga number; N(N + 1) is primary pseudoperfect (Giuga index of 42 = 20/21, of 1806 = 902/903).
- **M3.** "Giuga mod p^2 iff w_p == -1" and its data are known: primepuzzles Puzzle 998 (1277, 780887; Volpatti 2021 adds 56151923 and 11774118061 below 2*10^13); OEIS A306431, A268000, A079853.
- **M4.** Of the 120 new edges of W, 72 join radius-2 points; 48 join |z|^2 = 2/3 to |psi w|^2 = 8/3 (24 each way), same angular gap arccos(7/8). In M3 u psi P1 the 48 cross edges split 24 / 24 (`hn_cross_edges.out`).

## Minor findings

1. Torsion at 3/2, 5/6, 8/3 is Z/2 x Z/4.
2. `d19_isogeny.gp` aborted; the fact is true (`ellisogeny` with kernel <(-1,0)>).
3. The brute-force bound was |p| <= 3000, q <= 1500.
4. General type is Corvaja-Turchet-Zannier, arXiv:2403.02030 (K^2 = 8, p_g = 5). Bombieri-Lang is ineffective.
5. Odd Giuga numbers need >= 14 prime factors.
6. In general n' = a n + 1 (Grau-Oller-Marcen).
7. g13 + 1 is proved prime by gp's `isprime`.
8. Bounded-gap attributions (Polymath8b 246; Maynard 12 under EH; Polymath8b 6 under GEH).
9. Prime-merging corollary correct, periodicity mod 2^K confirmed exactly (n <= 2^18, K in {8, 12, 14}, h in {1, 2, 3}). The script uses THM-4627 (2)(a)'s event; over all n at 10^7 the fractions are 0.49647 (h=1), 0.40288 (h=2). The corollary is asymptotic only.
10. A divergent Collatz orbit is not a finite witness; the conclusion holds because Collatz is arithmetic.
11. THM-4628's eigenvalues are minus the Kloosterman sums at a conjugation-fixed prime. THM-4520 has no Riemann-type input and gives no mixing. The JC/Dixmier bullet is correct.
12. "Two oldest proofs" is hyperbole; Vojta for log2(3) is an analogy; the rational-cos link covers Parts' pieces only.
13. Only M3 u psi P1 was run; P1 u psi M3 is congruent by conjugation. "M3 is rigid" was not measured.
14. 1,441 is corroborated by Haugland, arXiv:2608.04542. The venue is unconfirmed, and 1441 = 1 + 24*60 is numerology.

## Independently verified

- **D19:**
  - Jacobians from the I, J invariants.
  - Conductors 32 and 40, j-invariants, torsion.
  - Rank 0 by ellrank and by analytic rank, L(E,1) = 0.6555..., 0.7422..., which makes it unconditional (Kolyvagin).
  - The 2-isogeny.
  - The j-formula.
  - All 46 ranks.
- **Wilson quotients** for all 17,983 primes below 2*10^5, by Legendre factorisation in C:
  - w_p == 0 at 5, 13, 563;
  - -1 only at 1277;
  - 1 at 3, 11, 107, 4931.
  - S(p) == p - 1 + p w_p directly for all 2,261 primes below 2*10^4.
- **Carmichael census to 10^8** via 2,057 base-2 pseudoprimes:
  - 255 Carmichael numbers, with an identical (G, omega) table;
  - 60 with c +- 2 prime.
- The 13 known Giuga numbers; the count 405,285; the sums 0.99896 and 1.03344; the merge data.
- **HN:**
  - Rebuilt in Q(zeta_24)(sqrt 5) with an all-pairs search: W = 64,513 / 542,472; M3 u psi P1 = 48,121 / 385,296.
  - The 336 near-miss pairs are correctly excluded.
  - L727 u rho S361 = 1,087 / 5,892 is not 4-colourable (glucose4, kissat404).
  - The four critical sets are all vertex-critical (union 784).

## Not checked

- The kissat404 run on M3 u psi P1 (48,121 vertices) had no verdict after about 2 h 10 min. Its output is `hn_sat_M3_u_psiP1_kissat404.out`.
- Forcing for M3 and P1.
- The Voronov-Neopryatnaya-Dergachev internal figures.
- Parts' 374 + 136 split.
- The full texts of Yang Ji and Ieno, and Heule's venue.
- Classical attributions, checked from memory only.
