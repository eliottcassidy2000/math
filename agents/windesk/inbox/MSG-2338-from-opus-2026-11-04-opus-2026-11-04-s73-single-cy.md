        # Message: opus-2026-11-04-S73: single-cycle floor (KNOWN) and crevasses as holes of the uniquely Hamiltonian stratum; cycle ceiling through order 10; audited PASS WITH CORRECTIONS (MISTAKE-705)

        **From:** opus-2026-11-04-S?
        **To:** all
        **Sent:** 2026-11-04 11:57

        ---

        opus-2026-11-04-S73 close-out: the single-cycle floor and the crevasses of the Hamiltonian-path spectrum. Owner
inspiration came from three papers: Bannwart (real realization), Grogan (prime pairs) and Ding-Kim (orthogonal matroids).

WHAT CHANGED
- Note: 05-knowledge/results/single_cycle_floor_crevasses_lorentzian_joins_20261104.md
- Canon:
  - THM-4866;
  - HYP-9457: P2, P4, P5 confirmed; P6 refuted; P1 and P3 were sanity checks; (b) CLOSED KNOWN;
  - MISTAKE-705.
- Program: 04-computation/experiments/scf_crevasses_20261104.c (+ census / chains / joins .out).
- Audit: 04-computation/experiments/scf_crevasses_20261104_audit_S73/AUDIT.md. Verdict: PASS WITH CORRECTIONS, all
  applied.

DECISIVE RESULTS
- Single-cycle floor: Busch-extremal strong tournaments are uniquely Hamiltonian. This is KNOWN: interval chains are
  upset tournaments (Brualdi-Li; Poet-Shader; Moon-Yang). The note gives an independent cut proof.
- Through order 10:
  - the cycle ceiling H >= f + 2(hc - 1) holds, with equality only at the floor;
  - the two-cycle floor is 13, 23, 41, 69, 115, 205, always below f(n + 1);
  - the single-cycle band is the bottom part of the uniquely Hamiltonian stratum;
  - 7 and 21 are holes of that stratum: 7 lies between U(4) = {5} and 9, and 21 is the hole of U(6) = {15, 17, 19, 23}.
- One-arc join law:
  - H = x_h y_h + x_R y_R - x_e y_e is KNOWN (session29 ADDENDUM sections 2 and 5);
  - the cycle formula hc = s_a e_b is added here;
  - every uniquely Hamiltonian tournament has a single-vertex one-arc cut. This is KNOWN (Shelor 1977).
- One-arc joins of orders <= 7 cover every odd number in [23, 84705]. This is a statement about the family, not new
  spectrum information. There is a 3-adic bias from triples â‰¡ (0, 0, 0) mod 3.
- For self-converse tournaments, H = core x s^2. This gives parity crevasses: 5 never occurs at odd order, 3 never at
  even order.
- Pfaffian test: null.

NEXT OBLIGATIONS
- HYP-9457 (a): the cycle ceiling for all n.
- HYP-9457 (c): join completeness.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
