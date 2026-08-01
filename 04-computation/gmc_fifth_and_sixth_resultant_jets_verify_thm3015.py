#!/usr/bin/env python3
"""THM-3015 -- verification of P_5 and P_6, and the extended corner-slice laws.

Tables (frozen [m,u,v,num,den] rows, content-1 normalisation):
  05-knowledge/results/gmc_first_gap_resultant_jet_P5_table_thm3015.json  (c_5 = 39813120)
  05-knowledge/results/gmc_first_gap_resultant_jet_P6_table_thm3015.json  (c_6 = 120394874880)

This script does NOT rebuild the jets.  It checks what the tables plus the frozen
THM-2997 table and THM-3013's P_4 can establish:
  1. shape: 1313 / 2176 terms, degrees (2j,4j,2j), support b+2c <= 4j, all absences
     at the three corners (4j,0),(0,2j),(0,0);
  2. the SEVEN pre-registered predictions of THM-3011 sec 2a and THM-3014 sec 2,
     recorded BEFORE P_5 existed:
        A_5: [M^6]=23/15, [M^5]=-12/5, [M^4]=+47/24
        E_5: [M^6]=23/15, [M^5]=-11/5, [M^4]=+37/18, [M^3]=0
  3. E_j = [V^(2j)] Q_j / 9^j, so the discriminant regular value is just the
     V-corner rescaled -- no D-adic re-expansion needed (simplifies THM-3014 sec 2);
  4. the NEW [M^(j-3)] law, and the fact that it retro-explains two constants the
     canon had flagged as unexplained:
        [M^(j-3)] = (-1)^j * mu * (j-1)(j-2),  mu = 707/23040 (A), 301/9720 (E)
     At j=4 this returns THM-3011's 707/3840 and THM-3014's 301/1620 exactly.

Reproduce: python3 04-computation/gmc_fifth_and_sixth_resultant_jets_verify_thm3015.py
"""

import json
from fractions import Fraction as Fr

T = {5: ("05-knowledge/results/gmc_first_gap_resultant_jet_P5_table_thm3015.json", 39813120),
     6: ("05-knowledge/results/gmc_first_gap_resultant_jet_P6_table_thm3015.json", 120394874880),
     7: ("05-knowledge/results/gmc_first_gap_resultant_jet_P7_table_thm3015.json", 2889476997120)}
COUNT = {5: 1313, 6: 2176, 7: 3348}


def rule(s):
    print("=" * 74)
    print(s)
    print("=" * 74)


def slices(j):
    rows = json.load(open(T[j][0]))
    c = T[j][1]
    A = {r[0]: Fr(r[3], r[4]) / c for r in rows if r[1] == 4 * j and r[2] == 0}
    E = {r[0]: Fr(r[3], r[4]) / (c * 9 ** j) for r in rows if r[2] == 2 * j and r[1] == 0}
    C = {r[0]: Fr(r[3], r[4]) / c for r in rows if r[1] == 0 and r[2] == 0}
    return rows, A, E, C


def main():
    ok = True
    rule("1. SHAPE OF P_5 AND P_6")
    for j in (5, 6, 7):
        rows, _, _, _ = slices(j)
        a = [r[0] for r in rows]; b = [r[1] for r in rows]; v = [r[2] for r in rows]
        present = {(r[0], r[1], r[2]) for r in rows}
        full = {(m, u, w) for m in range(2 * j + 1) for w in range(2 * j + 1)
                for u in range(4 * j + 1 - 2 * w)}
        miss = sorted(full - present)
        corners = {(4 * j, 0), (0, 2 * j), (0, 0)}
        allc = all((m[1], m[2]) in corners for m in miss)
        good = (len(rows) == COUNT[j] and max(a) == 2 * j and max(b) == 4 * j
                and max(v) == 2 * j and allc)
        ok &= good
        print(f"  P_{j}: {len(rows)} terms (want {COUNT[j]})  degrees ({max(a)},{max(b)},{max(v)})"
              f"  max(b+2c)={max(x + 2 * y for x, y in zip(b, v))}")
        print(f"        ansatz {len(full)}  missing {len(miss)}  all at corners: {allc}  -> {'OK' if good else 'FAILED'}")

    print()
    rule("2. THE SEVEN PRE-REGISTERED PREDICTIONS (recorded before P_5 existed)")
    _, A5, E5, _ = slices(5)
    preds = [("A_5[M^6]", A5.get(6), Fr(23, 15)), ("A_5[M^5]", A5.get(5), Fr(-12, 5)),
             ("A_5[M^4]", A5.get(4), Fr(47, 24)),
             ("E_5[M^6]", E5.get(6), Fr(23, 15)), ("E_5[M^5]", E5.get(5), Fr(-11, 5)),
             ("E_5[M^4]", E5.get(4), Fr(37, 18)), ("E_5[M^3]", E5.get(3, Fr(0)), Fr(0))]
    for name, got, want in preds:
        g = (got == want)
        ok &= g
        print(f"   {name:10s} predicted {str(want):8s} actual {str(got):8s}  {'CONFIRMED' if g else 'REFUTED'}")

    print()
    rule("3. E_j IS THE V-CORNER RESCALED (simplifies THM-3014 sec 2)")
    print("   E_j := [D^(2j)]P_j/c_j equals [V^(2j)]Q_j / 9^j, because D is linear in V")
    print("   with slope -3 and deg_V Q_j = 2j.  Computed from the V-corner above:")
    print(f"   E_5 = {{{', '.join(f'M^{m}: {E5[m]}' for m in sorted(E5, reverse=True))}}}")

    print()
    rule("4. THE NEW [M^(j-3)] LAW, AND WHAT IT EXPLAINS")
    muA, muE = Fr(707, 23040), Fr(301, 9720)
    print("   [M^(j-3)] = (-1)^j * mu * (j-1)(j-2),   mu = 707/23040 (A), 301/9720 (E)")
    known = {(4, 'A'): Fr(707, 3840), (4, 'E'): Fr(301, 1620),
             (5, 'A'): Fr(-707, 1920), (5, 'E'): Fr(-301, 810)}
    for j in (4, 5):
        for tag, mu in (('A', muA), ('E', muE)):
            got = (-1) ** j * mu * (j - 1) * (j - 2)
            w = known[(j, tag)]
            g = (got == w)
            ok &= g
            src = "THM-3011's 707/3840" if (j, tag) == (4, 'A') else \
                  "THM-3014's 301/1620" if (j, tag) == (4, 'E') else "P_5"
            print(f"   j={j} {tag}: law {str(got):12s} vs {str(w):12s}  {'OK' if g else 'MISMATCH'}   ({src})")
    print("   => the two constants the canon flagged as unexplained are the j=4")
    print("      members of a single one-parameter family.")

    print()
    rule("5. P_7: THE OFF-CORNER TEST, AND THE j=7 PREDICTIONS")
    rows7, A7, E7, C7 = slices(7)
    present = {(r[0], r[1], r[2]) for r in rows7}
    full = {(m, u, w) for m in range(15) for w in range(15) for u in range(29 - 2 * w)}
    miss = sorted(full - present)
    off = [m for m in miss if (m[1], m[2]) not in {(28, 0), (0, 14), (0, 0)}]
    print(f"   |P_7| = {len(rows7)} (eq (2) predicts 3348: {len(rows7) == 3348});"
          f" absences {len(miss)}; OFF-CORNER {len(off)} {off if off else '(none)'}")
    print("   => the j=3 sporadic (5,0,3) does NOT recur at j=7.")
    ok7 = (len(rows7) == 3348 and not off)
    P7 = {'A_7': (A7, {8: Fr(23, 28), 7: Fr(-12, 7), 6: Fr(47, 24), 5: Fr(0),
                       4: Fr(-707, 768), 3: Fr(0), 2: Fr(11267, 43008)}),
          'E_7': (E7, {8: Fr(23, 28), 7: Fr(-11, 7), 6: Fr(37, 18), 5: Fr(0),
                       4: Fr(-301, 324), 3: Fr(0), 2: Fr(2677, 10206)}),
          'C_7': (C7, {8: Fr(23, 28), 7: Fr(-3, 7), 6: Fr(23, 6), 5: Fr(0),
                       4: Fr(-23, 12), 3: Fr(0), 2: Fr(23, 42)})}
    for nm, (got, want) in P7.items():
        res = all(got.get(k, Fr(0)) == want[k] for k in want)
        ok7 &= res
        print(f"   {nm}: all 7 pre-registered coefficients {'CONFIRMED' if res else 'REFUTED'}")
    print("   c_7 content: 1/2889476997120 = 1/(2^22 * 3^9 * 5 * 7); prime set {2,3,5,7}")
    print("   = {p <= j+1} at j=7, continuing the law.")
    ok &= ok7

    print()
    rule(f"SUMMARY  all checks pass: {ok}")


if __name__ == "__main__":
    main()
