        # Message: collatz-procgen-20260922 wave 24: THM-4529..4534 + Lean (617 thms) - OPEN-Q-060 no natural bijection; Allikvere OP2 Theta(d^(1/3)) solved; all-odd tournaments to N=38; THM-4526 Thm A is classical (Grunbaum/Havet/El Zein)

        **From:** mac-mini-2026-10-01-S?
        **To:** all
        **Sent:** 2026-10-01 19:22

        ---

        collatz-procgen-20260922, wave 24 (2026-10-01): seven lanes on the owner's requests, all independently audited and in canon.

1. THM-4531 (OPEN-Q-060, bijective form): NO natural (relabelling-invariant) matching exists between switching classes of tournaments and even (untwisted) Euler graphs.
   - Classes to even Euler graphs: none bijective on types for 5 <= n <= 100.
   - Even Euler graphs to classes: no natural map at all for n >= 3.
   - Literature: THM-4524 E1 is the switching analogue of Royle-Praeger-Glasby-Freedman-Devillers 2023 ("tournaments and even graphs are equinumerous"). Both are instances of one twisted Mallows-Sloane principle. Their natural-bijection question is answered negatively, in the same equivariant sense, for 5 <= n <= 15 plus sporadic n.
   - Mallows-Sloane's even-n remark is explained.
   - Odd n: A049313 counts the tournaments with out = in (mod 4) at every vertex.
   - Byproduct: Higashitani-Ueyama's expected s_{4,n} = t_{4,n} (arXiv:2409.10904 Ex. 4.10) holds for every modulus.

2. THM-4534 (Allikvere arXiv:2608.09983, Open Problems 2-4):
   - OP2 solved in Theta form: 0.8146 <= liminf ln edim_m(Q_d)/d^(1/3) <= limsup <= 2.0526. C* = (3 sqrt2 ln2)^(2/3) is exactly the reach of the forest union bound; whether the limit exists is OPEN (HYP-9169).
   - OP3, partly: a closed form covers every d >= 17, and explicit sets cover 6..16 (new Q13..Q16). Finiteness for all d >= 6 therefore needs no union-bound computation.
   - OP4 answered.
   - 11 <= edim_m(Q7) <= 19.

3. THM-4529 (Petersen family, Paley coordinates):
   - Delta-Y on the 4 Fano lines avoiding 0 turns K6 = P7 - 0 into Petersen; on all 7 it turns K7 into Heawood. "7 translations <-> 7 graphs" is numerology.
   - Linear K6 has 3 linked pairs iff its Gale tournament is C3[TT2,TT2,TT2].
   - Anti-circulant antipodal parity.
   - QR127[mu14] is all-odd.
   - Collatz via the bridge: NO PROOF.

4. THM-4530 (Syracuse drops):
   - Every d in Z is a drop exactly 2 + N(6d+1) times over all odd integers (unit preimages 8d+1, -4d-1): the owner's "two copies", exact.
   - Exact densities.
   - Two consecutive drops determine the point.
   - Only q = 3 hits every integer.
   - Drop statistics cannot constrain cycles.
   - Builds on THM-4527 (S15).

5. THM-4533 (shaved tournaments). For opus S15:
   - THM-4526 Theorem A is classical (Grunbaum 1971 / Havet 2000 / El Zein arXiv:2204.11211).
   - Theorem C's constant is settled by Linial-Saks-Sos, Combinatorica 3 (1983): u(n) = n log n - O(n log log n).
   - A literature note was appended to THM-4526.
   - New results: a Redei-graph parity calculus; D69 settled (path + span-3 is Redei exactly for n = 4..6); the D_n family; u(9) = 14, kappa(9) = 22.

6. THM-4532 (the owner's Tournament Clock Prime Collatz):
   - Clock combination laws.
   - Two-sheet clocks D_r are all-odd for every odd r <= 19. That gives all-odd tournaments at every N = 2 mod 4 up to 38, including N = 34 and 38, where Paley cannot reach. C1 (all odd r) is in HYP-9170.
   - Syracuse mod-9 clock law (a reformulation).

7. Lean: 04-computation/lean/ProcgenSelfieEdim. 617 theorems, core Lean 4.30, no Mathlib, axioms propext/Quot.sound. Includes Redei's theorem from scratch, the selfie gauge theorem, the edim_m(Q6) bounds and lemmas, the drop theorems, anti-circulant parity for every odd m, and odd-abelian Cayley arc-evenness.

Process: a reboot wiped /tmp mid-wave. Lanes were recovered from subagent transcripts and relaunched in a persistent worktree.
Next obligations:
- HYP-9169: does the limit exist;
- HYP-9170 / C1: D_r all-odd for every odd r;
- HYP-9172: no natural bijection for every n >= 5;
- HYP-9173: Redei graph growth;
- N = 38 needs a second-engine confirmation;
- the Allikvere LRC(14) certificate audit is still owed and needs the owner's permission to download.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
