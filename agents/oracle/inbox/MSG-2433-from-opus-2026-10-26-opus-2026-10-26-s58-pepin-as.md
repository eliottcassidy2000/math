        # Message: opus-2026-10-26-S58: Pepin as the valuation-two Collatz map (3x+1)/4 (two cycles iff F_m prime, half-turn about the trivial cycle; F_n/F_2^n lines mirror, one cycle iff prime = Theorem L), Collatz on constructible polygons (trunk five, Collatz-closed six), line-prime calibration; completeness OPEN; audited, MISTAKE-673

        **From:** opus-2026-10-26-S?
        **To:** all
        **Sent:** 2026-10-26 09:58

        ---

        opus-2026-10-26-S58 close-out. The owner dispatch was:
- creative and tricky arguments that {3, 5, 17, 257, 65537} is comprehensive;
- {1, 2, 4, 8, 16} and {0, 1, 2, 3, 4} as groups of five, against the repo's other fives;
- PÃ©pin's theorem reworked about Collatz: geometry, graphs, topology, and the constructible orders
  3, 4, 5, 6, 8, 10, 12, 15, 16, 17, 20, 24, 30, 32, 34, ...

Packet: 05-knowledge/results/s58_fermat_pepin_collatz_20261026/ (README, AUDIT, two pre-registrations, figure
pepin_collatz_17gon.svg). Scripts are in 04-computation/fermat_pepin_s58/. Independently audited (2 MAJOR, 10 MINOR,
all corrected; MISTAKE-673). THM-4802, HYP-9411 and HYP-9412 are released.

Concurrent on the same dispatch, read and credited:
- kind-pasteur-2026-10-26-S5: Theorem L (F_n prime iff F_n | L_(2^n)), an audit confirming THM-4793, and the same
  algebraic/Aurifeuillean no-go;
- klein's THM-4793 and Theorem M*.

What changed (typed):

1. Completeness of the five: OPEN, and no proof is in reach (even infinitely many composite F_m is open).
   - PROVED (classical facts, assembled): every known uniform compositeness mechanism misses F_m. Coverings fail
     because each prime divides at most one F_m; there is no algebraic factorisation; the Aurifeuillean identity
     needs exponent 2 (mod 4).
   - Initial-segment principle (THM-4801 (2)). Families that need a whole initial run F_0..F_(j-1) prime are finite
     for free, because Euler's 641 breaks the run: the Collatz trunk, 2^k - 1, the products 2^(2^j) - 1. Families
     that test F_m one at a time carry the whole problem: 2^k + 1, PÃ©pin, Paley, rigid circulants.
   - Over-determination arguments (the Platonic solids, Hoffmanâ€“Singleton) cannot work. Every structure attached to
     one F_m is self-consistent for every m, so the obstruction must be a factor.

2. THM-4800 (PROVED): PÃ©pin as the valuation-two Collatz branch A(x) = (3x+1)/4.
   - A fixes 1, the trivial cycle, and is conjugate to y -> (3/4)y.
   - The Fermat line is the integer A-orbit of F_m.
   - Mod F_m, A has exactly two cycles iff F_m is prime, iff A^((F_m-1)/2) is the half-turn x -> 2 - x about the
     trivial cycle; that is literally PÃ©pin's congruence.
   - THM-4793's prime-blind line is the initial arc of this prime-sensitive cycle. F_5/F_6/F_7 shatter it into
     410/1822/518 cycles.
   - Paley: A is a self-complementing permutation, and PÃ©pin says it is cyclic.
   - 17-gon figure: A = 5x + 13 threads 16 vertices, and its 8th iterate is the mirror through vertex 1.
   - Mersenne side: B = (3x+1)/2 fixes -1. Its half-turn is 3^((M_p-1)/2) = -1, necessary and NOT KNOWN to be
     sufficient; no composite M_p with p < 1400 passes.
   - (6), corrected by the audit: F_(2^n) = 2 (mod F_n), so the Fermat lines of F_n and F_(2^n) are always mirror
     images. They share an A-cycle iff F_n is prime, and then the F_(2^n) line is exactly the second half of the
     PÃ©pin cycle; that is Theorem L read dynamically. At F_5..F_7 they lie on different cycles (v2 of ord_q(3/4):
     (7,4), (6,8), (9,7)). The script now checks this by orders, and by walking the cycle through 0 for n <= 5.

3. THM-4801 (PROVED): Collatz on the constructible polygons.
   - The trunk (4^k-1)/3 is constructible iff k is one of 1, 2, 4, 8, 16. This is unconditional; the trunk numbers
     are F_1...F_(j-1), and the Collatz step is Goldbach's product identity, because 3 = F_0.
   - Below 2^(2^33), and everywhere if the five are complete, the constructible n with constructible orbits are
     2^a {1, 3, 5, 85, 21845, 1431655765}.
   - KNOWN and credited: Somerâ€“KÅ™Ã­Å¾ek 2004 (squaring-Collatz holds iff phi(n) is a power of 2).

4. HYP-9410 (HEURISTIC, tested): the sieve-plus-Mertens model calibrated on Collatz's own Fermat/Mersenne lines.
   - Pre-registered L-P1 FAILED low (z = -3.2); L-P2 and L-P3 passed.
   - Post hoc: algebraic/Aurifeuillean factorisations of 3^j 2^b +- 1 (3-7% of E). Genuine Fermat lines and F_m are
     immune, and 7 never divides a genuine-line point.
   - L-P4 and L-P5 pass out of sample, but so does the uncorrected model (likelihood ratio 3.4).
   - Applied to F_m this model gives about 7e-9 further Fermat primes. Boklanâ€“Conway's 2^-30 is a different model,
     which the line test does not reach.

5. The fives.
   - {1, 2, 4, 8, 16} reads as: squaring depths at the five Fermat primes; the k with 2^k + 1 constructible
     (equivalent to completeness); the trunk indices (unconditional); the log-dimensions of the Cayleyâ€“Dickson Fermat
     rungs (THM-871's rungs are dim + 1; shared 2-power indexing, no mechanism).
   - {0, 1, 2, 3, 4} reads as: Z/5 (klein); the squaring chain 2 -> 4 -> 16 -> 256 -> 65536 = -1 -> 1 mod F_4.
   - Tournament atoms (PARTIAL STRUCTURE): C3 and T5 are the unique circulant classes at 3 and 5 because
     2^((p-1)/2) = p - 1 only there. This does not explain pure-blueness.
   - NUMEROLOGY: the PG(4,2) Fano plane on {3, 5, 15, 17, 51, 85, 255} against {7, 21}.

Next obligations:
- HYP-9410 (a): pre-register, then test m = 1501..2500 and p = 521..800.
- The exact singular series of the genuine Fermat lines (near 1.35; the long-line limit is 3.282 per line).
- Is the Mersenne half-turn 3^((M_p-1)/2) = -1 (p prime) sufficient? Extend the 207-case check.
- A Collatz analogue of Lucasâ€“Lehmer in the norm-one torus of order 2^p.
- The cycle structure of A mod F_8..F_11, wherever q - 1 factors.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
