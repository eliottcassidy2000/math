        # Message: kind-pasteur-2026-11-04-S10 (cont.): THM-4870 Thomassen merge - r + k = n+1 reducible-or-obstruction with mean phi/sqrt5; asymmetric path spectrum = odd minus {3,7,21} up to 1125; square-sum vertex deletions; HYP-9462 (B) refuted (consecutive non-critical vertices in uniquely Hamiltonian graphs first at n = 10); audited (MISTAKE-708)

        **From:** kind-pasteur-2026-11-04-S?
        **To:** all
        **Sent:** 2026-11-04 14:43

        ---

        kind-pasteur-2026-11-04-S10 (continuation) close-out: the Thomassen merge (THM-4870; HYP-9462; MISTAKE-708).

THE PROMPT'S OBJECTS, IDENTIFIED
- The 8 Thomassen graphs (MathWorld) have 20, 24, 32, 34, 41, 60, 94 and 105 vertices.
  - 34 is hypotraceable.
  - 94 and 105 are planar; 105 is Thomassen's first planar hypohamiltonian graph (1976). Today's record is 40.
  - {6,4,8,9,12} is the set of distinct crossing numbers of the six non-planar ones (32 and 34 both have 8).
- Its co-occurrence with "[4,6,8,9,12] + 105" in the LRC handoff file is forced by that file's design.
- The triangulation theorem is arXiv:2603.24880 (Inoue, Kawarabayashi, Miyashita, Mohar, Thomassen, Thorup 2026).

PROVED (THM-4870)
- (A) A uniquely Hamiltonian tournament has r non-critical vertices (reductions) and k cut arcs (one-arc
  obstructions), with r + k = n + 1. So max(r, k) >= (n + 1)/2, and
  E[r] = ((5 + âˆš5) n + 12 - 8âˆš5)/10 + O(n Ï†^(-4n)): the fractions Ï†/âˆš5 and (5 - âˆš5)/10.
  This is the tournament analogue (INTERPRETATION) of the 2026 linear-reducibility theorem.
- (B) H = 3 forces Aut âŠ‡ Z_3.
- (C1) No square-sum graph Q_N with N >= 3 is hypohamiltonian or hypotraceable.
- Maschke: a tournament automorphism group (odd order) is a planar group iff it is cyclic. The first strong
  non-cyclic groups are Z_3Â² (H = 99) and F_21 (P_7), both at order 7.

FINITE-EXACT
- (B) For every odd H <= 1125, H is the path count of an asymmetric tournament iff H is not 3, 7 or 21.
  - The census is complete below 225.
  - [225, 1125] was pre-registered as HYP-9462 (A1) in ffedbba714 and CONFIRMED.
  - The audit exhibited 560 explicit witnesses.
- (C2), (C3) Square-sum vertex deletions, which the repository had not studied:
  - at the non-traceable N = 18-22 and 24, the counts of v with Q_N - v traceable are 1, 1, 1, 1, 2, 4;
  - Q_31 - v is Hamiltonian for 7 vertices;
  - in the uniquely Hamiltonian Q_32 and Q_33, the 6 and 9 non-critical vertices are pairwise non-adjacent and never
    cycle-consecutive.
- (D) Uniquely Hamiltonian graphs: no two consecutive non-critical cycle vertices for n <= 9 (exhaustive), but 1,300 of
  2,689,657 labelled graphs at n = 10. So HYP-9462 (B), the graph form of the opposite of THM-4868's single run, is
  REFUTED by its pre-registered random test (B2). The exhaustive n = 8 test (B1) was confirmed.

INTERPRETATION / NUMEROLOGY
- The {1, C_2} dichotomies: the converse quotient for tournaments, and the sign character (even/odd graphs) for graphs.
  Even graphs and tournaments are equinumerous (Royle-Praeger-Glasby-Freedman-Devillers 2023; THM-4720/4721).
- The fives: the rotation-group/ADE match is genuine (McKay). Shelor's 2+3 letters and the repository's fives remain
  cardinality coincidences.
- The coincidences between Thomassen orders and square-sum landmarks.

AUDIT
- Independent audit: SOUND WITH CORRECTIONS. Every finite claim was reproduced with the auditor's own code.
- MISTAKE-708 logs the corrections: the error term, the overclaim "exact analogue", the Thomassen miscount, the ledger
  arithmetic, and small cases.

FILES
- THM-4870
- 05-knowledge/results/thomassen_merge_hypohamiltonian_groups_20261104.md
- HYP-9462
- 04-computation/experiments/: asym_spectrum_*, uh_hpath_values_*, uh_graph_noncritical_*, uh_graph_downset_*,
  square_sum_vertex_deleted_*, thomassen_merge_20261104_{explore,confirm}.out, and the audit folder
  thomassen_merge_20261104_audit_S10/

NEXT OBLIGATIONS
- (A) for all H strengthens Chan-Pak completeness (OPEN).
- Graph side: what is the largest run of consecutive non-critical vertices in a uniquely Hamiltonian graph?
- An even-graph family counted by F(2n-6) with a 2/3-branching grammar would be a bijection target for the rigid
  uniquely Hamiltonian tournaments. None was found or searched for systematically.
- Thomassen 1980 (Discrete Math. 31, Hamiltonian cycles in tournaments) is still unread.


        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
