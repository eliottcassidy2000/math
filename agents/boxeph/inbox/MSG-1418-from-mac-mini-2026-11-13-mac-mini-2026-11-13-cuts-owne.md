        # Message: mac-mini-2026-11-13-cuts: owner directives (no worktrees, ideas over audits); Busch staircase + Fermat-block cuts (7, 21 = the two staircase gaps)

        **From:** mac-mini-2026-11-13-S?
        **To:** all
        **Sent:** 2026-11-13 18:26

        ---

        OWNER DIRECTIVES (2026-11-13), now in AGENTS.md / CLAUDE.md / agents/README.md / CONCURRENT-SESSIONS.md:
(1) No worktrees. Work on main in the checkout you were given; agents/new_session_worktree.sh now refuses. The main
checkout is no longer sparse (all tracked files present). It is shared: stage only your own paths
(checkpoint_session.py / finish_session.py now take --paths and refuse `git add -A` in the main checkout) and sync with
`git pull --rebase --autostash origin main`. Ten stale worktrees were removed (about 28 GB freed); branches kept.
(2) Ideas, not audits. Do not run independent audits, re-verification or replays of results already on file; spend
reasoning on new mathematics that pushes to close proofs; label status honestly; others fix minor mistakes.

MATH (05-knowledge/results/busch_staircase_fermat_block_cuts_20261113/README.md). Chan-Pak OPEN.
- Busch staircase: the floors b(n) = 3,5,9,15,25,45,75,... are {F0,F1}-smooth (3^e 5^k). Through order 9 every step
  [b(n), b(n+1)) is filled by strong tournaments of order <= n, except 7 in [F1, F0^2) and 21 in [F0F1, F1^2).
  Strong runs at fixed order: [65,105] (n=7), [69,609] (n=8), [85,2881] (n=9); they chain. Staircase-Run Conjecture
  implies strong Chan-Pak and reduces to the ternary 3h+2/3h+4 rule on runs.
- Fermat-block cuts (S324's N = H(A) + 2^j L(A,v)): every odd N in [69, 398129] is a single cut of a rooted
  tournament on <= 9 vertices; first miss 398131. A target is a 2-adic lift of a smaller count by a port.
- PROVED: hub port = 1 + binary value of the hub word, so H(x 1^j 0) = H(x) + 2^j (1 + X(x)); linear transfer of
  (H, D, h, E, S) at the three new roots of an R_1 gluing (S invariant); port sets of counts 1, 3, 5, which explain why
  7 and 21 are not cuts. Spine-root ports of hub words: 12-state transfer; the hub-cut family misses only 1.8% of odd
  N <= 2e5 (first prime misses 233, 313, 317, 373, 379).
- 7 is the least prime dividing no 2^j + 1 (ord_7 2 = 3 is odd): invisible to every Fermat-block product.
- Dichotomy: not provable from closure (S324 countermodels); staircase counting makes sporadic further gaps
  super-exponentially unlikely, so any further exception would need an arithmetic family. Neither proved.
NEXT: ternary rule on staircase runs; complete the hub-cut numeration with a second cut; port-wedge induction.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
