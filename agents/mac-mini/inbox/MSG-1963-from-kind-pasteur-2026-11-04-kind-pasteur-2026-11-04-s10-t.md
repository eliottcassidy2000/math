        # Message: kind-pasteur-2026-11-04-S10: THM-4874 PROVED - the uniquely Hamiltonian deletion law (core block decides hc(T - v); 3, 2, and 3 or 2 deletions); HYP-9467 confirmed at orders 10-15; audited PASS WITH CORRECTIONS (MISTAKE-714)

        **From:** kind-pasteur-2026-11-04-S?
        **To:** all
        **Sent:** 2026-11-04 17:36

        ---

        kind-pasteur-2026-11-04-S10, owner follow-up "prove it" (the uniquely Hamiltonian deletion law, HYP-9464 (B)): close-out.

WHAT CHANGED
- THM-4874 (canon, PROVED; audited PASS WITH CORRECTIONS, 0 major / 11 minor, MISTAKE-714). It is the uniquely
  Hamiltonian deletion law for sole-exit-free tournaments with k = 3 (n >= 6) or k = 4, 5 (n >= 7) Hamiltonian
  cycles.
- HYP-9467 was pre-registered in 671f7fb87e and is CONFIRMED. All 48 predicted joint counts at orders 10-15 are exact,
  with no other combination. Order 15 has 1,090,032 strong classes with hc <= 5, and the run took 16 min on the
  HYP-9464 auditor's independent census.
- HYP-9464 (B) is now PROVED in its range.
- Commits on main:
  - 671f7fb87e: reservation and pre-registration;
  - 460c918f9a: pre-audit checkpoint;
  - 6a32eea229: audit applied;
  - 7115e9704d: MISTAKE-714.

THE PROOF (decisive facts)
- By THM-4869 (C2), T = Q[z <- (D, u)]: a uniquely Hamiltonian prefix R, and a core D = chain => block (W, A) plus u,
  with u -> l the only arc out of D.
- The deletion lemma: for v not in {u, l}, hc(T - v) = e_u(D - v) * s_l(R - v), because every cycle crosses u -> l
  once.
- Hence hc(T - v) by position:

  | position of v | hc(T - v) |
  |---|---|
  | prefix | a multiple of k |
  | u | 0 |
  | l | 0, or k or h0(W, A) if R = {l} |
  | chain | k |
  | block | delta(v) = sum over x in A - v of end_x(W - v) |

- So the uniquely Hamiltonian deletions are exactly the block vertices with delta = 1. The eight blocks of type <= 5
  (S75's B5) give:
  - k = 3: 3 deletions;
  - k = 4: 2 deletions;
  - k = 5: 3 deletions iff the block is (AT_5, its two sole exits), else 2. That gives F(2n-13) classes with three
    and L(2n-12) with two.
- The block also fixes #{v : hc(T - v) = 2} for n >= 7. With S76's THM-4872 (B): a three-cycle sole-exit-free
  tournament has a nontrivial automorphism iff no deletion leaves exactly two cycles.
- Every one-arc join on (AT_m, {m-1, m}) with a uniquely Hamiltonian prefix has exactly three uniquely Hamiltonian
  deletions, at 1, m-1 and m.

CHECKS
- Block data re-derived over all labelled tournaments, by the session to order 6 and by the audit to order 7.
- Per-vertex table: 2,661 deletions on 270 built joins (session) and 9,536 deletions on all 917 in-range classes at
  orders 6-11 (audit), with 0 mismatches.
- The audit's own census reproduces the joint counts at orders 7-11.

NEXT OBLIGATIONS
- The law needs THM-4869's range. For hc >= 6, other cores appear (THM-4869 lists blocks only for type <= 5); a
  deletion law there is OPEN.
- Indecomposable {-3}-self-dual tournaments at orders >= 12 remain OPEN (HYP-9464 (A)).
- S76's note still cites HYP-9464 as "pending"; it is confirmed.


        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
