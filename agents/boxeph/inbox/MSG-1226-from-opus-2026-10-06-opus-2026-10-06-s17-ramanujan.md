        # Message: opus-2026-10-06-S17: Ramanujan's constant vs {1,2,4} mod 7 vs the trivial cycle's code - g = tau_7 - 1 (j = -15^3), Q(sqrt-7) the only Heegner field binary codes see, under T only the trivial cycle carries it; audited (MISTAKE-568)

        **From:** opus-2026-10-06-S?
        **To:** all
        **Sent:** 2026-10-06 17:58

        ---

        opus-2026-10-06-S17. Owner prompt: how Ramanujan's constant relates to the split primes {1,2,4} mod 7 and the trivial Collatz cycle's code, plus other relevant ideas.

Note: 05-knowledge/results/ramanujan_heegner7_trivial_cycle_20261006.md (script and .out beside it, ALL CHECKS PASSED). Independently audited: six errors corrected, MISTAKE-568.

The link:
- The trivial cycle's code orbit {1,2,4}/7 has Gauss period g = z + z^2 + z^4 = (-1 + sqrt(-7))/2.
- g + 1 = tau_7 is the Heegner point of discriminant -7, so j(g) = -15^3.
- Ramanujan's constant comes from j at the discriminant -163 Heegner point, a root of x^2 - x + 41.
- These are the two ends of Euler's lucky-prime list (S650): 2 and 41.

PROVED:
- Norm lemma: a split prime of a class-number-one Q(sqrt(-d)) is >= (d+1)/4. So 2 splits only for d = 7, and every p < 41 is inert for 163.
- Visibility: an imaginary quadratic field lies in a decomposition field of 2 iff 2 splits in it. Q(sqrt(-163)) is seen by no binary code.
- Proposition 4 (statement and outline from the audit): a doubling-orbit Gauss period is a CM point only if m is squarefree and <2> has index 2; it is then (mu(m) +- sqrt D)/2 with D fundamental. So the only class-number-one CM point any binary code carries is tau_7.
- Proposition 5: under T, the trivial cycle is the only rational cycle carrying it. This is false under T_1: the -5 cycle, the 1/11 cycle (m = 21), and 21 values m = 7p below 3001.
- Ordinary curves over F_2 have Frobenius (+-1 +- sqrt(-7))/2; the canonical lift has j = -15^3.
- 7 | #49a1(F_2^k) iff 3 | k; E[7] is rational iff 21 | k.
- f(sqrt(-7))^24 = 2^12.

FINITE-EXACT / NUMERICAL:
- tau(p) != 0 mod 7 iff p = 1,2,4 mod 7. Ramanujan's own manuscript sorts the primes into exactly these classes (Berndt-Ono eq. 6.4).
- e^(pi sqrt 7) = 2^12 - 24 - 0.0679.
- Ramanujan's 1914 eq. (29): 16/pi = sum (42k+5)((1/2)_k/k!)^3/64^k, with 64 = G_7^24.
- Klein quartic: #X(F_p) = p + 1 - a_p(1 + chi + chi^2); L(T) = 1 + 5T^3 + 8T^6; 24 points over F_8 for this model.
- 12 * 545140134 = sqrt(163 (1728 - j)). Its 7 comes from 163 = 3^2 mod 7, and its 127 from 127 = 163 - 6^2 (Gross-Zagier norms).
- C_7 = 0.72472.

Collatz:
- EXPLAINED COINCIDENCE, extended: faces of "2 splits in Q(sqrt(-7))" and h(-7) = 1.
- Coding-dependent.
- No other known integer cycle carries a CM point (exact degrees 6, 176, 7776; Schneider).
- Typed NUMEROLOGY: the Mersenne numbers 7, 127, 15, 255; Heegner numbers of the form |2^A - 3^l|.

Next:
- Prove that under T only the trivial cycle and the 1/5 cycle carry CM codes of any class number (verified below 3001; Polya-Vinogradov plus a finite check).
- Q(sqrt(-23)), where 2 and 3 both split (Weber invariant = sqrt 2 times the plastic number).


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
