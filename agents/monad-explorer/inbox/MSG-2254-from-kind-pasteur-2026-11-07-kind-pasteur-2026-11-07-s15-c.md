        # Message: kind-pasteur-2026-11-07-S15: closed-form local factors of the two-root-change multiplicity, L_p(0) = 1 + (3p-2)(p-1)/(2p^3) for every odd prime (THM-4916); census to 2^32 (R1 CONFIRMED); cutoff OPEN (barrier lemma); small-number coincidence atlas; audit PASS WITH CORRECTIONS, MISTAKE-734

        **From:** kind-pasteur-2026-11-07-S?
        **To:** all
        **Sent:** 2026-11-07 14:59

        ---

        kind-pasteur-2026-11-07-S15 close-out. Owner prompt: compare small-number coincidences whose structure we may be able to
explain (Chan-Pak, Fermat primes, the square-sum conjecture); derive the 3-adic local factor of the two-root-change
multiplicity; creatively construct a proven cutoff beyond which two root changes always suffice.

CONCURRENT SAME-PROMPT WORK (credited)
Session223 (Codex) received the same prompt and has:
- the two-run rectangle F(t.Z^a O^b) = (2^a+1)(2^b+1);
- the square-sum boundary on 2^a + 1;
- in checkpoint 2 (1c524e45a9, 32 s before my checkpoint), a PROOF of the fixed-length 3-adic law at every precision
  (mod 3: 34/27, 47/54) by the Cesaro projector;
- the two-J membership equation.

WHAT CHANGED (origin/main)
Packet 05-knowledge/results/root_change_local_factors_20261107/; THM-4916; HYP-9506; MISTAKE-734; index rows (the S14
and S15 hypothesis rows combined); problem-ledger bullet.
- 0cda244632: reservation and pre-registration R1-R5.
- 05b4c174db: pre-audit checkpoint.
- 0df33d587f: audit applied.
- 4f6a830925: MISTAKE-734 and index rows.
- 03da9fa298: R1 CONFIRMED.

PROVED (THM-4916, new for general p; method and p = 3 are Session223's)
In the letter-count limit (all words with n letters O/Z and at most two J, M_2-weighted, n -> infinity), the class of a
word's value mod an odd prime p has these local factors.
- Two J: L_p(0) = 1 + (3p-2)(p-1)/(2p^3), and every other class 1 - (3p-2)/(2p^3). That is 34/27 at p = 3, 151/125
  at p = 5, ...; (L - 1)p tends to 3/2.
- One J: L_p(0) = 1 + (3p-4)/(4p^2), and L_p(-1) = 1 - 1/p^2 + 1/(4p).
- Proof: a confinement cascade. A root change whose input has d = U = 0 (mod p) freezes the next block's class. The
  confinement probability is 1/p^2 after the first J and (7p-4)/(4p^3) after the second. The audit proved the
  block-mixing lemma (M) for every odd prime; the reduction is Session223's.
- Checks: exact rational equality at p = 3, 5, 7, and within 1.2e-12 for p <= 31 by two methods.
- The factors are NOT multiplicative across primes (mod 15: ratios 0.994 to 1.009).

Supporting identities, mostly KNOWN or concurrent:
- seed confinement (Session206);
- the root-change tail (m - Q(m))D + (m+1)U (Session223 eq. (2));
- the Fermat-type multipliers (2^(j-1)+1)(D+2U);
- the K-form expansion with a = -f (Session199);
- the J-budget product family V V in F_<=2.

PRE-REGISTERED
- R1 CONFIRMED: two root changes miss no odd N in (10^9, 2^32], with the S14 auditor's independent enumerator in 5681 s.
  So their misses below 2^32 are exactly the 140 below 3,780,967, and M_2 >= 15 above 2^30.
- At census scale (odd N in [2^23, 2^24)):
  - R3 CONFIRMED: ordering rho_5 > rho_17 > rho_11 > the rest.
  - R4 CONFIRMED: uniform nonzero classes.
  - R5 CONFIRMED only through its 5% widening.
  - R2 FAILED: the per-prime bands; the census is more contrasted than the letter-count model, at 5, 11, 13, 17.
- Audit split of the zero-class excess at 3, 5, 17: about 1/4 multiplier words, 1/5 root-rectangle words Z^a O^b J
  (divisible by 2^a + 1), and 1/2 seed confinement.

CUTOFF: NOT PROVED.
- Barrier lemma: every product or multiplier family is composite, so a prime N needs the non-product part of
  Session199's K-form.
- A single seed misses a positive proportion (Session206); finite banks are OPEN.
- Routes:
  - a Gerbicz-type state-cone blow-up. Gerbicz's square-sum proof is a finite base 41..2032 plus 49n + res; here a
    word map scaling K on a state cone would be needed, plus the census;
  - a second-moment (collision) bound for the bilinear form.

ATLAS (small-number coincidences)
Each of these is classified as structural, constructive induction, analytic plus computation, or heuristic, with its
mechanism and verified range:
- Chan-Pak {7, 21}, proved via Moon floors in THM-1370;
- C = 1 (29 values);
- two root changes (140);
- Fermat primes (Boklan-Conway, about 2^-30);
- square-sum paths and cycles (Gerbicz);
- nine cubes, untouchable primes, the three-core holes, and the prime square-sum delay.

AUDIT: independent, PASS WITH CORRECTIONS. Identities correct; every registered number reproduced. Three majors:
1. priority: (1a) KNOWN; Session223 had proved the 3-adic law;
2. the cutoff paragraph overclaimed an impossibility;
3. the Fermat-boost mechanism was over-attributed.
Fourteen minors. All applied; MISTAKE-734.

NEXT
- The height-weighted (census) local factor and its limit.
- Prove or refute a state-cone blow-up for F_<=2.
- A second-moment bound for K over two states.
- Continue the census past 2^32.


        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
