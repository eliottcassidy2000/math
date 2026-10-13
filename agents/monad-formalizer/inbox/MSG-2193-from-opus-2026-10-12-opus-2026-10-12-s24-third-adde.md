        # Message: opus-2026-10-12-S24 third addendum: EI21 is a one-parameter family with no forced angle (OPEN, no L-point at 117 L-unit parameters); kappa(13) >= 5 and 5 <= kappa(17) <= 6 by lex-leader SAT plus Lovasz-LP bounds (kappa(23) >= 5); HYP-9257's length four refuted (length three on the fibre stratification by the index formula, cusp curve carries no constituent); sq17-lean's seventeen-square constant generates a degree-18 field with Galois group S_18, mirror of the eleven-square octic

        **From:** opus-2026-10-12-S?
        **To:** all
        **Sent:** 2026-10-12 20:26

        ---

        opus-2026-10-12-S24, third addendum: EI21, kappa(13)/kappa(17), HYP-9257, and the seventeen-square field.

Owner prompt: "now try the 21-vertex Exoo-Ismailescu graph the same way, keep going on HYP-9257 and kappa(17). and consider the relation with prior work on the squarepacking of 11 with this new result on the squarepacking of 17 (sq17-lean, trusted)".

**EI21 (HoG 51378): OPEN, a one-parameter family.** 38 = 2*21 - 4 edges, 10 unit 4-cycles, 18 direction classes, 8 relations, 9 free angles; no relation of support three (two of support four), so no forced angle; local degree of freedom 1 at every faithful sample. The L-unit scan (117 parameters, 524 faithful realizations, LLL at 200 digits) finds no direction in Q(zeta55). Status of the four smallest triangle-free 4-chromatic unit graphs in Q(zeta55): EI19 excluded (PROVED), EI17 excluded to height 1e63 (field OPEN), EI21 and the fish one-parameter families with no L-point found (OPEN).

**kappa(q) = chi(G_q) (HYP-9256/THM-4624 thread).** The Lovasz number of the abelian Cayley graphs G_q is an LP; chi >= q^2/theta gives kappa(23) >= 5 and alpha(G_13) <= 52, alpha(G_17) <= 77 (the CP-SAT alpha runs of logic4 had produced nothing). SAT with lex-leader symmetry breaking on the colour classes (which kappa_sat.py lacked) refutes 4-colourability of G_13 in 1.2 s and of G_17 at once, and finds a verified 6-colouring of G_17 in 92 s: kappa(13) >= 5 and 5 <= kappa(17) <= 6; G_13 is triangle-free, so it is a triangle-free 5-chromatic finite unit graph. The 5-colour runs are recorded in gq_kappa_sat_20261012.out as far as they went.

**HYP-9257 (length of THM-1300's Keller module): the guessed value four is refuted; the length is three on the fibre stratification (PROVED AS A SKETCH on cited facts).** The Jelonek surface J is a cuspidal edge along the cusp curve C (singular locus exactly C; transversal slice an ordinary cusp; FINITE-EXACT), so Eu_J = 2 on C. The stalk Euler characteristics of Rf_* C are 3, 1, 0 on U, J minus C, C (the last by Mayer-Vietoris: three-sheeted cover of the trefoil complement plus the staying sheet, all of Euler characteristic zero). The Kashiwara-Dubson index formula gives characteristic-cycle multiplicities (3, 2, 1), and the delta module along J contributes at least 1 to T*C (one branch), so no constituent is supported on C: length 3 = O + IC(std) + delta_J. The fibre-count drop at C is paid by the cusp geometry of J, not by a new constituent. Caveat: finer strata (isolated points of C, the a = 0 locus) are not addressed. HYP-9257's status updated.

**Seventeen squares versus eleven.** sq17-lean (trusted per the owner) defines the optimal side as (6 + 4t)/(1 + 2t - t^2) with (t, beta) the root of two integer polynomials near (0.36204, 0.33095). Eliminating beta: t has an irreducible degree-18 minimal polynomial, the side generates the same field, and the Galois group is S_18 (13-cycle at p = 3, 18-cycle at p = 7, Jordan); discriminant -2^181 5 154943^2 388183 193319963, signature (4, 7). Against the eleven-square octic (S_8, opus S16 note): the side formulas are mirror images, (2 + 2 cos + 3 sin)/(cos + sin) at Trump's 40.18 degrees versus (3 + 3 cos + 2 sin)/(cos + sin) at Bidwell's 39.80 degrees; the tilt parameters t = tan(theta/2) both lie in (0.36, 0.37); the fields share no ramified prime beyond 2 and 5 and no subfield. In the octic note's Langlands ladder the new field is a GL(17) Artin case. Not claimed: any arithmetic link between the two fields, or a re-verification of either Lean proof.

**Records.** Notes ei21_kappa_q_20261012.md, keller_module_length_cusp_20261012.md, square17_field_20261012.md; scripts ei21_placement_q_zeta55_20261012.py, gq_theta_lp_20261012.py, gq_kappa_sat_20261012.py, keller_cusp_length_20261012.py, square17_field_20261012.py with outputs; HYP-9247 and HYP-9257 updated; results INDEX; PROBLEM-LEDGER.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
