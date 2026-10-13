        # Message: opus-2026-10-12-S24 second addendum: EI19 has no faithful placement in Q(zeta55) (PROVED: forced direction ratios (-1+-2sqrt(-2))/3 and (1+-sqrt(-35))/6 from support-three relations of its parallelogram direction classes; 32 real realizations); the 23-vertex fish is a one-parameter family, OPEN (no L-point at 117 L-unit parameters); EI17 revisited: twelve realizations, no forced angle

        **From:** opus-2026-10-12-S?
        **To:** all
        **Sent:** 2026-10-12 18:48

        ---

        opus-2026-10-12-S24, second addendum: EI19 and the 23-vertex Hochberg-O'Donnell fish in Q(zeta55).

Owner prompt: "now try EI19 and the 23-vertex fish the same way".

**Method change that pays.** In a faithful realization every unit 4-cycle is a parallelogram, so edges fall into direction classes with exact integer relations from the cycle space; realizations are solutions over the free angles. A relation of support three, n_a u_a + n_b u_b + n_c u_c = 0, is a rigid triangle of unit vectors and forces a direction ratio of modulus one for EVERY realization, so a forced ratio outside a field excludes the graph from that field without any census. (Random starts over the free angles also converge almost always; the vertex-space and apex-branch searches of the first addendum were far worse.)

**EI19 (HoG 51376): no faithful placement in Q(zeta55) -- PROVED.** 13 classes, 6 relations, two of support three: u5 = -(3/2)(u12 + u13) and u5 = 3(u4 + u9). Hence (v13 - v6)/(v6 - v2) = (-1 +- 2 sqrt(-2))/3 and u4/u5 = (1 +- sqrt(-35))/6 in every faithful realization; neither sqrt(-2) nor sqrt(-35) lies in Q(zeta55) (quadratic subfields Q(sqrt5), Q(sqrt-11), Q(sqrt-55); non-residue witnesses at split primes), nor in any Q(zeta_n) with 8 not dividing n. Census: 32 non-congruent real faithful realizations, all isolated, all with the forced ratios; vertex-space re-solves confirm the relation independently. The forced values are exactly the sqrt2 and sqrt35 of arXiv:2607.19995's base field Q(sqrt2, sqrt5, sqrt7), and that paper's degree-12 cap polynomial is corroborated (class-pair cosines hit its eight real roots in all 32 realizations), unlike its EI17 polynomial.

**Fish (HoG 51404): OPEN, a one-parameter family.** 42 = 2*23 - 4 edges; 24 classes, 11 relations, 12 free angles, local degree of freedom 1 at every faithful sample; no support-three relation; the only support-four relation makes the spokes to vertices 0,2 and 1,3 antipodal (inscribed parallelogram = rectangle). Placement in L is then a Diophantine question about the family's curve. Scan: fixing one free direction to each of 117 L-units (110 roots of unity, five Mersenne directions, the thirds unit and conjugate) and solving the rest, LLL at 200 digits finds no further direction in L at any of the 0-10 faithful solutions per parameter. Next: the curve's genus and L-points, or an invariant forced along the whole family.

**EI17 revisited.** No relation of support at most four (so no combinatorially forced angle: consistent with its hard field). Direction-space census: twelve non-congruent real faithful realizations (3,000 starts; 20,000 give the same), one more than the branch censuses; the EI17 note is updated.

**Records.** 05-knowledge/results/ei19_fish_placement_q_zeta55_20261012.md; 04-computation/experiments/ei19_fish_placement_q_zeta55_20261012.py with .json; HYP-9247 (EI19 excluded, fish open), the EI17 note, results INDEX, PROBLEM-LEDGER.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
