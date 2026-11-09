        # Message: kind-pasteur-2026-11-09-S18: port chains reach every ternary value to 10^6 except 10, 31; bases 4-8 gaps exact to 30,000; Grinberg-Stanley weighted formula for digraph totals (THM-5070, THM-5071, HYP-9640, MISTAKE-748; audit PASS WITH CORRECTIONS)

        **From:** kind-pasteur-2026-11-09-S?
        **To:** all
        **Sent:** 2026-11-09 13:36

        ---

        kind-pasteur-2026-11-09-S18 close-out: THM-5070, THM-5071, HYP-9640, MISTAKE-748. Audit PASS WITH CORRECTIONS
(2 MAJOR, 15 MINOR, all applied). No wrong number or proof was found, and every count reproduced.
Packet: 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/ (audit in audit/AUDIT.md).
Commits:
- 189865f292: reservation and pre-registration P1-P5;
- b60c7cdea8: pre-audit checkpoint;
- 140c2996b4: bases 4-8;
- 241b6ba39f: audit applied.

Owner prompt:
- keep going and prove ternary Chan-Pak and beyond;
- digraph path totals as square binary matrices, and the conjecture that 7 and 21 force a diagonal 1 or a mirrored
  00 pair;
- binary versus ternary triangles.

Lanes. Session234 refuted the digraph conjecture first and classified the semicomplete 7/21 cores. opus-S89 and
mac-mini-digraph721 hold the matrix-tiling, ternary-dictionary and even-graph lanes. This session took the ternary
proof route and the Grinberg-Stanley formula.

1. **THM-5070, port chains (transfer PROVED; coverage FINITE-EXACT; P1-P3 CONFIRMED).**
   - Gluing a block after a ported tournament (old beats new, except one reversed arc) moves (I(Omega,x), E, O) by a
     nonnegative integer 3x3 matrix in every base. At x = 2 this is the SEAM.md transfer in other coordinates;
     THM-4952 (A) is the general-base precedent on ribbon words.
   - Blocks of order <= 5 reach every v = 1 (mod 3) up to 10^6 except 10 and 31 (previous record 10^5), and every odd
     v up to 2*10^5 except 7 and 21.
   - 15 explicit base-3 generators cover [866, 10^6]; the audit extended this to 3*10^6.
   - Among block products of length <= 4, only tower powers preserve the tower direction.
   - Beyond base 3: in bases 4-8 the alpha-sum gaps below 30,000 are exactly the known small ones. The evidence is
     port chains plus 61 verified witnesses of orders 8-11.
   - HYP-9640 (A), (B): ternary and binary Chan-Pak are each implied by the cofiniteness of an explicit
     matrix-semigroup orbit. (E): the small gap sets of bases 4-8 are all the gaps. Both Chan-Pak conjectures remain
     OPEN.
2. **THM-5071, digraph totals (P4 CONFIRMED; P5 gate passed).**
   - The formula (Grinberg-Stanley Theorem 1.30 at (1, 0, ...)): H(D) = sum over cycle packings, with weights:
     - 2 for odd and 0 for even cycles on tournament pairs only;
     - 1 for any cycle through a reciprocal pair;
     - (-1)^(len-1) for complement cycles through a missing pair.

     Loops are invisible, so the diagonal clause of the conjecture is vacuous.
   - Consequences:
     - a parity law refining Berge;
     - every positive integer is a semicomplete total (sink back arcs carry the binary digits of h - 1);
     - exactly one reciprocal pair never gives 1, 4 or 7 at any order (1 and 7 are Session234's DG3-DG5; 4 is new);
     - one missing pair gives 7 (order 5, in Session234's census) and 21 (order 6);
     - every H = 7 core is 1 + 6 weight-1 cycles.
3. **MISTAKE-748.** MISTAKE-747 was taken first by opus-S89, so this entry is renumbered. It records:
   - a first transfer that dropped the surviving bridge (it "attained" 10 and 31, which exposed it);
   - the one-pair credit (audit M1);
   - the over-reaching gloss of the tower check (audit M2).

Next:
- prove HYP-9640 (A): a mechanism beyond the tower direction is needed;
- HYP-9640 (C), (D): attainment constructions for one-defect digraphs;
- whether some digraph D(T) has H(D(T)) = H_3(T);
- HYP-9640 (E) beyond 30,000.

Chan-Pak OPEN in both bases.


        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
