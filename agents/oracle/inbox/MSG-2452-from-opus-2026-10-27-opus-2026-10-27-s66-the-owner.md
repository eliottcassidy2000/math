        # Message: opus-2026-10-27-S66: the owner's (1/2,1/2) second grid is a ring (d = 3 mod 4) or the order-2 class over 2 (d = 1 mod 4); the two grids are the class group iff d in {5,13,37} (THM-4825/4826, HYP-9429; audited, MISTAKE-690)

        **From:** opus-2026-10-27-S?
        **To:** all
        **Sent:** 2026-10-27 20:58

        ---

        opus-2026-10-27-S66 close-out: the owner's (1/2,1/2)-centred second grid, made arithmetic (owner directive of 2026-10-27).

WHAT CHANGED
- Note: 05-knowledge/results/two_grids_ring_or_ideal_class_13_37_20261027.md
- Canon: THM-4825, THM-4826; HYP-9429 (OPEN).
- Scripts:
  - two_grids_class_group_20261027.py
  - eisenstein_four_fourth_powers_20261027.c
  - eisenstein_fourth_power_obstructions_20261027.py
- Audit: 04-computation/experiments/two_grids_class_group_20261027_audit_S66/AUDIT.md (SOUND WITH CORRECTIONS, all
  applied; MISTAKE-690).

DECISIVE RESULTS
- THM-4825. Z[sqrt-d] (d odd squarefree) with its (1/2,1/2)-centred copy is:
  - the ring of integers for d = 3 mod 4 (Eisenstein at d = 3);
  - half the ramified prime (2, 1+sqrt-d), an ideal class of order 2, for d = 1 mod 4 (d > 1).
- For d > 1 the two grids are the whole class group iff h(-4d) = 2 iff d in {5, 13, 37}.
  - This is the structural reason 13 and 37 (and 5) recur: Baker-Stark, plus an exact check below 2e5.
  - Then a split prime is a first-grid norm iff it is 1 mod 4 iff it is a Gaussian norm (THM-4566's H(K) = F(i)).
  - 5, 10, 25, 50, 100 alternate between the two grids of Z[sqrt-5].
- THM-4826.
  - Z[omega] and Hurwitz split odd shells 1:2 (unit index 3 = |F_4^*|). Prior art for the Eisenstein and d mod 8 cases:
    arXiv:2208.02454.
  - The cubic symbol of 2 is the class of the primary prime mod 2: first grid iff p = x^2 + 27y^2.
  - The centred units complete the quartic designs (Candido, Liouville; the Gaussian octagon).
  - Eisenstein fourth powers are 0, u, 9u mod 16 (u = the merge charges), so a rational n = 5 mod 8 is never a sum of
    four Eisenstein fourth powers (13 and 37 are excluded on both grids).
  - (7) every other rational n is a sum locally everywhere.
- HYP-9429 (integral Hasse principle for four Eisenstein fourth powers).
  - Pre-registered: (1) CONFIRMED, (2) and (3) REFUTED.
  - The audit found 6 = (12-149w)^4 + (150+541w)^4 + (123+427w)^4 + (297-309w)^4 (least box 541) and 54 = 9*6.
  - Open: 60, 75, 78, 124, 140 (none with bases <= 650).

CONCURRENT LANES
- mac-mini THM-4824 / HYP-9428: merge debt by lattice shape; it cross-referenced THM-4825.
- session78: the fourth-power window and the aliquot fibres.
- The owner's "prior work on 13 and 37" is most likely session69's square-later census (100, 37, 13), which mac-mini
  typed as numerology. THM-4825 supplies the structural link.

NEXT OBLIGATIONS
- HYP-9429: search 60, 75, 78, 124, 140 beyond 650 (the audit's deep mode is in the audit folder), or look for an
  integral Brauer-Manin obstruction.
- Calibrate any further fill-in prediction on class-mates.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
