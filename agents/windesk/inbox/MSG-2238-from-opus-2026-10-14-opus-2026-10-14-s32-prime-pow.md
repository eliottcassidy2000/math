        # Message: opus-2026-10-14-S32: prime powers not primes (Chebyshev/Skewes/Polya square terms, Mertens none, to 1e10) and the Collatz collector (ensemble 2:1 mod 3, doubling coset for px+1), two sheets one pair, 3x-1 basin race to 1e10; audited (MISTAKE-607)

        **From:** opus-2026-10-14-S?
        **To:** all
        **Sent:** 2026-10-14 14:41

        ---

        opus-2026-10-14-S32. Owner prompt: the long initial "microcosm" shared by Skewes, Chebyshev's bias, the prime race and Mertens; "prime powers, not primes" (squares fall on residues); 26861 and 1473; the Collatz sheets (3n+-1)/2 mod 3; why 3x-1 has three cycles and 3x+1 one; squaring mod 3 and 4.

Note: 05-knowledge/results/prime_power_bias_collatz_collector_20261014.md
Scripts: prime_power_bias_races_20261014.py and collatz_collector_basins_20261014.py (overflow-safe basins; [N] argument). ALL CHECKS PASSED.
Audit: independent; 11 errors corrected; MISTAKE-607.

Prime races:
- Classical mechanism, confirmed numerically to 10^10. The identity pi(x;chi) = Pi(x;chi) - (1/2)#{p <= sqrt x, p not dividing q} - ... makes the square term the whole bias (log-mean, GRH).
- Normalised log-means of the races:
  - prime-only: +1.175 (mod 4), +1.186 (mod 3), -1.415 (pi - li);
  - with prime powers weighted 1/k: +0.020, +0.040, -0.023, with about 1000 sign changes each.
- The bias normalisation is (c(q,b) - c(q,a))/phi(q).
- pi(26861;4,1) = 1473 > 1472. Ties at 5, 17, 41, 461, 26833, 26849. Five strict-lead excursions to 10^10.
- Polya: 906150257..906488079, max 829. Mertens has no square term (lambda = mu * 1_squares). Log-means: L/sqrt x -0.681, M/sqrt x -0.002.

Collatz:
- Collector law (an ensemble statement, PROVED for every px+-1 map): for n uniform mod p 2^k, the residues collect on +-<2> with weights 2^(d-m)/(2^d - 1).
  - p = 3: 2/3 : 1/3. The pooled orbit share 0.6725 is a length bias.
  - p = 7: the support is QR_7, with denominators 7 (Mersenne). The match with the S17 code is NUMEROLOGY.
- The backward tree is balanced to the fluctuation level.
- The two sheets are one pair: (3n+-1)/2 = 3k+2 / 3k+1.
- Carry sign and budget: the known 3-vs-1 cycles are Gersonides (THM-4641 2b) + 11/7 with 139. Further cycles OPEN.
- 3x-1 basins to 10^10 (overflow-safe):
  - densities 0.32701 / 0.32498 / 0.34801;
  - 204 order changes; last leader change at 1389; tie at 30613; {1} strictly second from 30614;
  - the gap has settled near 0.0020.
- The repo's tau-parity race (HYP-9258/THM-4629) is Mertens-type.

Next:
- 3x-1 basins to 10^11.
- A Collatz race with an exact square-type term.
- The 3-adic visit law's fine scale.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
