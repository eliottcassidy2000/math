        # Message: opus-2026-10-26-S57: untouchables vs Pollack-Pomerance (model constant 0.17865689; class assumption fails 62-82% at computable sizes; residue signature), THM-4790 tc fibres 7/21 absent n<=11, THM-4792 merged sheet tables (Moree 1997), root-27 R-P1

        **From:** opus-2026-10-26-S?
        **To:** all
        **Sent:** 2026-10-26 00:16

        ---

        opus-2026-10-26-S57 close-out. The owner dispatch was:
- merge the tables mod odd n (3n-1 / 3n+1);
- a source-preserving rank on the root-27 orbit, with Frobenius-defect / clone-merger machinery;
- {2,3,5,7} = {3,5} | {2,7};
- a second odd untouchable or the exact untouchable share;
- the {7,21} tournaments as graph pairs;
- F1 with the aliquot fibres.

Packet: 05-knowledge/results/s57_untouchables_root27_20261026/ (README, AUDIT, two pre-registrations); scripts are
in 04-computation/untouchables_s57/. Independently audited by two referees. MISTAKE-672: MISTAKE-671 was taken by
mac-mini minutes earlier, and the entry was renumbered on rebase.

What changed (typed):

1. HYP-9400 (OPEN): the untouchable share against Pollackâ€“Pomerance.
   - VERIFIED: PP's model constant evaluates to Delta = 0.17865689, by window averages to 1e-8 and the audit's
     deterministic convolution (0.17865689139). The quoted "0.17" / 0.171822 are intermediate values:
     S(y)/log y = Delta - 0.1621075/log y reproduces PP's printed value and CGJM's Table 2 to 6 digits. The limit is
     not stated in PP, CGJM, OEIS, erdosproblems #418 or Booker.
   - FINITE-EXACT: exact census to 1e10 (= A057978).
   - At computable sizes PP's class assumption gcd(s(n),A_y) = gcd(n,A_y) fails jointly for 62% (y=13) to 82%
     (y=53) of inputs. Residue neutrality fails at each prime l for a share ~(log X)^(-1/(l-1)) (Ramanujanâ€“Watson).
   - Signature: the share varies across nonzero residues (spreads 1.45/1.29/1.28/1.20 mod 5/7/11/13 at 2^33,
     decaying like (log k)^(-0.75)). The class model cannot produce this; it comes from the classical kill/serve rule
     of the cores s(aP) = s(a)P + sigma(a) (PP sec. 2; Pomeranceâ€“Yang), shown by an exact core decomposition.
   - Unexplained: the deficient-end share (k/s(k) >= 1.4) is flat at 0.22 from 2^24 to 1e10, against PP's 0.14.
   - 13/13 pre-registered checks pass, 7 of them committed before the 1e10 array was analysed.

2. Odd untouchables. The Goldbach reduction is classical (kind-pasteur, session66). The Goldbach-dependent share of
   odd touchables is 11.25% -> 9.26% -> 8.24% (1e6 -> 1e8). The single-source odd numbers below 1e8 are
   3, 7, 9, 11, 39. No covering argument replaces Goldbach.

3. THM-4790 (PROVED + FINITE-EXACT, audited, no error): no tournament on n <= 11 vertices has tc = H/|Aut| in {7,21}.
   - Graph-pair chart: a class receives exactly tc labelled (backward, forward) graph pairs on [n-1].
   - Proof: reduction to strong; the Moon/Busch floor gives |Aut| >= 3,3,5,7,11; an odd-order cycle-type lemma
     (Cauchy, Z_3 wr Z_3, Sylow + Jordan); invariant enumeration.
   - The audit re-implemented it independently, and lemma-free runs at n = 9, 11 found no hit.
   - This extends HYP-9352 (2) from n <= 9.

4. THM-4792 (merge/reference; (c) KNOWN).
   - The centre labels -+1/2 of every odd multiplication table are the reflection centres of ceil/floor(3x/2).
   - The unit shift (the root-27 section change) conjugates the two one-step maps.
   - The halving-shift alternation of kind-pasteur's tables exists iff m | 2^t + 1, iff all ord_p(2) share one 2-adic
     valuation. That is Moree, Acta Arith. 80 (1997), Thm 1 (OEIS A014657). I first claimed it as new; the audit
     corrected this.

5. Root-27. The source-preserving rank was settled during the session by kind-pasteur S4 (Theorem B: rank =
   |w| - LCP with the carry word; Haar share sqrt5-2) and session66-root27 (switch cocycle, merger rule, unbounded
   refills). S57 adds:
   - R-P1, pre-registered: deficit excursion 20 at 1e5 returns, in the window [15, 35];
   - the phi^-d drop law on 2e6 exact symbols.
   The Frobenius defect does not apply (agreeing with session66). THM-4791, HYP-9401 and HYP-9402 are released.

6. F1 (with session66's monoid carrier). The fibres of s and tau - 1 coincide exactly at 0 and 1 (the owner's first
   two entries), and differ for every k >= 2. A hostile: q+1 is untouchable for 785 of 3661 odd prime powers q, so
   "2 = |P^1(F_1)|, 5 = |P^1(F_4)|" is numerology. sigma(14) = sigma(15) = 24.

Next obligations:
- HYP-9400 (4): a 2^40 census with residue and abundancy tallies, and a target-side heuristic built on the kill/serve
  rule.
- HYP-9352 (2) at n = 12 (order-27 subgroups of the Sylow 3-subgroup of S_12).
- The limit of the odd Goldbach-dependent share.
See OPEN-QUESTIONS (S57 update) and INVESTIGATION-BACKLOG (S57 entry).


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
