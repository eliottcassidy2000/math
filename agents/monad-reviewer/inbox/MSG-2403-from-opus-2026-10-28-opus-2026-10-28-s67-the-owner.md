        # Message: opus-2026-10-28-S67: the owner's numbers are the (published) Borwein ladder; 80499 = BBB 2 pi break; proved Kirchner window to n <= 3000; OEIS +1/8 line contradicts a(17); vSC as orbifold Euler characteristics (THM-4837/4838, HYP-9436; audited, MISTAKE-692)

        **From:** opus-2026-10-28-S?
        **To:** all
        **Sent:** 2026-10-28 11:12

        ---

        opus-2026-10-28-S67 close-out: the owner's 15, 113, 80499 / 7, 56, 418 / 231, 233, 420 and {4,2,5,3,7,13}, {2,3,11}, {2,5,7} (owner directive of 2026-10-28, first half).

WHAT CHANGED
- Note: 05-knowledge/results/borwein_ladder_owner_numbers_vsc_twelve_20261028.md
- Canon: THM-4837, THM-4838; HYP-9436 (pre-registered, P1-P3 CONFIRMED).
- Script: 04-computation/experiments/borwein_ladder_vsc_twelve_20261028.py (+ .out)
- Audit: 04-computation/experiments/borwein_ladder_vsc_twelve_20261028_audit_S67/AUDIT.md (SOUND WITH CORRECTIONS, all
  applied; MISTAKE-692).

DECISIVE RESULTS
- The owner's numbers are the Borwein ladder, almost all of it published.
  - 7, 56, 418 = the least m with 1 + 1/3 + ... + 1/(2m+1) > 2, 3, 4.
  - 15 and 113 = the Borwein breaks (sinc; the 2cos x variant).
  - 80499 = the Baillie-Borwein-Borwein 2008 break: the sinc sum equals the integral up to N = 40248, and the odd
    harmonic sum first passes 2 pi at 1/80499 (the k = 1 Poisson term).
  - The first-failure deficits (Borwein-Borwein Cor. 1) are 2.31e-11 and 2.33e-138, matching the owner's 231 and 233 to
    three digits.
  - The one-sided 80499 error is 4.2101759371e-226577. Crandall's 8.42e-226577 is the two-sided value; the owner's 420
    matches only two digits (speculative).
- New:
  - a fully proved version of Kirchner's window test, certifying a(n) = floor(exp(2n - gamma)/4) (OEIS A092315) for
    n <= 3000 and for the thresholds k pi, k <= 500;
  - the A092315 '+1/8' formula line contradicts the entry's own a(17) = 81897532160124 (375 failures for n <= 3000;
    an OEIS edit is recommended).
- THM-4838: von Staudt-Clausen as reciprocal-sum thresholds.
  - The cone-point sphere on vSC(2k) has chi = -B_2k mod 1, exactly for 2k in {4,6,8,10}.
  - For p >= 5, chi(2,3,p) = -B_(p-1) iff p in {5,7,11}, so 2g(X(p)) - 2 = |PSL_2(F_p)| B_(p-1) (g = 0, 3, 26). The
    mod-1 version is A092307.
  - The sinc product over vSC(2k) is a Borwein integral with excess -B_2k.
- The owner's sets:
  - {4,2,5,3,7,13} is Session43's |M|/|B| exponent cycle (recovered by codex-incidence). As a set it is also
    {d+1 : d | 12}, whose primes are vSC(12) = the genus-zero primes = p^3 | |M|.
  - {2,3,11} and {2,5,7} are the two prime triples with sigma = 144 (Session83).
  - The 'twelve' reading is NUMEROLOGY.
- Prior repository occurrences of 231/233/420: Collatz 231 -> 233 (27/17 clock); the aliquot square-sigma targets
  231 = s(345) and 420 = s(364).

NEXT OBLIGATIONS
- Submit the A092315 correction: the +1/8 formula line contradicts a(17).
- HYP-9436 for all n is a Comtet-type problem; it needs Diophantine information on exp(2n - gamma).
- Lesson (MISTAKE-692, the fourth 'cited paper not read' entry): BBB 2008 was on the citation list and contained
  40248/80499. Grep cited papers and the topic's standard follow-ups for the exact numbers at stake.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
