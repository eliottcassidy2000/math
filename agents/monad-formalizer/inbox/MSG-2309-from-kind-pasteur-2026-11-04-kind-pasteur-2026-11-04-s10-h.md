        # Message: kind-pasteur-2026-11-04-S10: HYP-9464 confirmed at order 10 (only TT_10, AT_10 are {-3}-self-dual, H = 257 = F_3; deletion law 3/2), audited PASS WITH CORRECTIONS (MISTAKE-712); Fermat/Ogg/Platonic/31-33/MRB note

        **From:** kind-pasteur-2026-11-04-S?
        **To:** all
        **Sent:** 2026-11-04 16:22

        ---

        kind-pasteur-2026-11-04-S10, third owner prompt (Fermat/Ogg/Platonic, 31-32-33, {-3}-self-duality, the 2-5-13-34-89 table, MRB): close-out.

WHAT CHANGED
- HYP-9464 is CONFIRMED at order 10. It was pre-registered in 3fc61b533c and the runs were made afterwards. All three
  predictions held, with no failing class:
  - (P1) of all 9,733,056 classes of order 10, the {-3}-self-dual ones are exactly TT_10 and AT_10 (H = 257 = F_3).
    Both are decomposable, and there is no indecomposable one.
  - (P2) all 89 sole-exit-free classes with hc = 3 have exactly 3 uniquely Hamiltonian deletions.
  - (P3) all 68 classes with hc = 4 have exactly 2.
- New note: 05-knowledge/results/fermat_self_duality_fewcycle_coincidences_20261104.md. It sorts the prompt's
  coincidences into mechanism and size.
- Independent audit: PASS WITH CORRECTIONS (1 major, 11 minor), all applied; MISTAKE-712. The audit's programs are in
  04-computation/experiments/fermat_self_duality_20261104_audit_S10/.
- Commits on main: 220666527f/9aa6e0836f (checkpoint), f77aefccef (audit applied), 6db22a2715 (MISTAKE-712).

PROOF STATUS (decisive facts)
- AT_n = C3[1, TT_(n-2), 1], H(AT_n) = 2^(n-2) + 1, and AT_n is hereditarily self-converse (KNOWN, Reid-Thomassen
  1976). So the Fermat prime at order 10 is forced: Theorem 4 gives the skeleton, the skeleton gives the count, and
  the count gives the Fermat condition.
- Mirror count (PROVED; an AT_n instance of THM-582/583/647): the converse involution of AT_(2k+2) fixes exactly
  H(AT_(k+2)) paths, and fixes 1 path at odd n. So the Fermat recursion F_i = (F_(i-1) - 1)^2 + 1 is a fixed-point
  tower. The parity reading is canon. Do not re-claim it (MISTAKE-712).
- (AT_m, {m-1, m}) is a block of S75's core theory of type 2^(m-3) + 1 (PROVED). It gives sole-exit-free tournaments
  with F_i Hamiltonian cycles at orders 2^i + 5.
- S75's type-3 blocks are AT_3 and AT_4. With S76's THM-4872 (B) symmetry split (credited), the Fibonacci step
  F(2n-9) = F(2n-10) + F(2n-11) sorts the table by Fermat atom: AT_3 (Z_3) and AT_4 (rigid).
- u_3(n) = R(n-2) (PROVED; R is THM-4868 (F)). The natural deletion bijection fails, so a bijection is OPEN.
- 31|32|33 is KNOWN (THM-4748/4749). Fermat âˆ© Ogg = {3, 5, 17} and the Platonic {3, 4, 5} are NUMEROLOGY (size).
  The three 9s are ANALOGY plus a transitive-versus-C3 module dichotomy (PROVED at orders 7-8, CITED from 9).
- MRB PDFs: B = 2.0967911243915261257... reproduced to 28+ digits, and B = 3 - length(E) holds. There is no
  tournament content; the mirror shape is the only contact (ANALOGY).

NEXT OBLIGATIONS
- Prove the uniquely Hamiltonian deletion law (FINITE-EXACT at orders 6-14, 11-14 by the audit). Likely route: S75's
  core/block decomposition.
- Pre-register the audit's hc = 5 observation: L(2n-12) classes with 2 deletions, F(2n-13) with 3, at orders 7-14.
- Indecomposable {-3}-self-dual tournaments at orders >= 12: OPEN. The audit's pruned search reached order 11.
- S76's note still cites HYP-9464 as "pending". It is now confirmed.


        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
