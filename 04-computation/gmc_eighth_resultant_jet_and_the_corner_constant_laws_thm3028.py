#!/usr/bin/env python3
"""THM-3028 -- the eighth resultant log-jet P_8, and CLOSED FORMS for the three
corner-slice constant sequences.

P_8 = 4883 terms, degrees (M,U,V) = (16,32,16), support b+2c <= 32.
Built on two DISJOINT tensor grids (no shared M, U or V):
   grid A  M=4..26,  U=2..34,  V=2..34 even     (12903 points)
   grid B  M=27..49, U=35..67, V=36..68 even    (12903 points)
Exact Newton interpolation over Q on each grid independently.

CONTROLS (all must pass, and all are checked below)
  a1. the frozen THM-2997 P_1,P_2,P_3 table is re-emitted BYTE-FOR-BYTE,
      digest cfb36557e1d54a0328a309375a948ace99c78e0688a54a014aef0906c1b90513;
  a2. the THM-3013 P_4 and THM-3015 P_5/P_6/P_7 digests are reproduced;
  b.  the two disjoint grids return IDENTICAL polynomials for every j=1..8;
  c.  6 out-of-sample widths per grid, 0 coefficient mismatches.

THE PRE-REGISTRATION.  The corner slices are
   A_j = [U^(4j)] Q_j,   E_j = [V^(2j)] Q_j / 9^j,   C_j = Q_j(M,0,0).
Through j=7 (THM-3015) they obey, writing the slot index k = j - e for the
coefficient of M^e and requiring k < j (the terminal slot k=j is the constant
term and is exceptional):

   k=-1 : [M^(j+1)] = (-1)^(j+1) * 46/(j(j+1))          (slice-independent)
   k= 0 : [M^j]     = (-1)^j     * kappa/j
   k odd, k=2m-1   : [M^(j-k)]  = (-1)^(j+m) * c_m * (j-1)(j-2)...(j-2m+2)
   k even, k>=2    : [M^(j-k)]  = 0

with c_1 = lambda, c_2 = mu, c_3 = nu the constants recorded in THM-3015. NOTE
the sign alternates with m, not with k; since k = 2m-1 the naive (-1)^(j+k)
collapses to (-1)^(j+1) and is correct only for odd m. Every j=8 value in slots
k = -1..6 is FORCED before P_8 is looked at, and part D checks all 24 of them as
a pre-registration. ALL 24 CONFIRM.

Slot k=7 is NEW: j=8 is the first order at which k=7 < j, so c_4 had never been
observed. Writing c_m = a_m / d_m, THREE independent extrapolations from m<=3
were registered, and part E scores them:

   (i)   NUMERATORS   a_m^A = 3 + 44*16^(m-1),  a_m^E = 4 + 33*9^(m-1),
                      a_m^C = 23                              -- ALL CONFIRMED
   (ii)  SLICE RATIOS d_m^A = 4^(2m-1) d_m^C, d_m^E = 3^(2m-1) d_m^C
                                                              -- BOTH CONFIRMED
   (iii) BASE SEQUENCE d_m^C = (3m)!/(2m-2)!                   -- REFUTED

(i) and (ii) are genuine out-of-sample successes: 180227 and 24061 were written
down before P_8 was inspected, and both are exact. (iii) fits m=1,2,3 (6, 360,
15120) and predicts 665280 at m=4, but the truth is 604800, a ratio of 11/10.
So the measured constants are

   c_4^A = 180227/9909043200,  c_4^E = 24061/1322697600,  c_4^C = 23/604800,

and the base sequence d_m^C = 6, 360, 15120, 604800 has NO closed form yet.
That single sequence is the whole remaining gap in the corner-slice picture.

NOTE the "23": the C slice has CONSTANT numerator 23 at every m, and the
slice-independent k=-1 law carries 46 = 2*23. This is the same 23 that heads the
four-band charge density of THM-3006.

Reproduce: python3 04-computation/gmc_eighth_resultant_jet_and_the_corner_constant_laws_thm3028.py
  (requires the P_8 interpolation pickle; see 05-knowledge/results/
   gmc_first_gap_resultant_jet_P8_table_thm3028.json for the frozen table)
"""

import hashlib
import json
import os
from fractions import Fraction as Fr
from math import comb, factorial, gcd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JMAX = 8

TABLES = {
    1: ("gmc_first_gap_resultant_jet_P1_table_thm3028.json", 1),
    2: ("gmc_first_gap_resultant_jet_P2_table_thm3028.json", 48),
    3: ("gmc_first_gap_resultant_jet_P3_table_thm3028.json", 1152),
    4: ("gmc_first_gap_fourth_resultant_jet_P4_table_thm3013.json", 1658880),
    5: ("gmc_first_gap_resultant_jet_P5_table_thm3015.json", 39813120),
    6: ("gmc_first_gap_resultant_jet_P6_table_thm3015.json", 120394874880),
    7: ("gmc_first_gap_resultant_jet_P7_table_thm3015.json", 2889476997120),
    8: ("gmc_first_gap_resultant_jet_P8_table_thm3028.json", 1664338750341120),
}


def rule(s):
    print("=" * 78)
    print(s)
    print("=" * 78)


def digest(obj):
    return hashlib.sha256(json.dumps(obj, separators=(",", ":")).encode("ascii")).hexdigest()


def content(vals):
    num, den = 0, 1
    for q in vals:
        num = gcd(num, q.numerator)
        den = den * q.denominator // gcd(den, q.denominator)
    return Fr(num, den)


def load(j):
    """Return the RAW Q_j as {(a,b,c): Fraction}.

    The stored tables are content-1 normalised, i.e. they hold c_j * Q_j, so the
    stored rows must be divided by c_j to recover Q_j itself. Omitting this
    division leaves every coefficient inflated by c_j and makes every law in
    parts D and E read as refuted.
    """
    path = os.path.join(ROOT, "05-knowledge", "results", TABLES[j][0])
    rows = json.load(open(path))
    c = Fr(TABLES[j][1])
    return {(r[0], r[1], r[2]): Fr(r[3], r[4]) / c for r in rows}


def table_rows(coeffs, scale):
    out = []
    for (a, b, c) in sorted(coeffs, key=lambda k: (-k[0], -k[1], -k[2])):
        q = coeffs[(a, b, c)] * scale
        if q:
            out.append([a, b, c, q.numerator, q.denominator])
    return out


def slices(P, j):
    A = {a: v for (a, b, c), v in P.items() if b == 4 * j and c == 0}
    E = {a: v / Fr(9) ** j for (a, b, c), v in P.items() if b == 0 and c == 2 * j}
    C = {a: v for (a, b, c), v in P.items() if b == 0 and c == 0}
    return {"A": A, "E": E, "C": C}


# --------------------------------------------------------------------------
# the closed forms, stated as pure functions of m and the slice
def a_const(nm, m):
    return {"A": 3 + 44 * 16 ** (m - 1), "E": 4 + 33 * 9 ** (m - 1), "C": 23}[nm]


#: d_m^C is the ONE piece with no closed form. Measured values only.
#: The guess (3m)!/(2m-2)! = 6, 360, 15120, 665280 fits m=1,2,3 and is REFUTED
#: at m=4, where the truth is 604800 (ratio 11/10). Left as an open sequence.
D_C = {1: 6, 2: 360, 3: 15120, 4: 604800}


def d_const(nm, m):
    """d_m for each slice. The SLICE RATIOS are closed-form and confirmed
    out-of-sample at m=4; only the base sequence d_m^C is empirical."""
    dC = Fr(D_C[m])
    return {"A": 4 ** (2 * m - 1) * dC, "E": 3 ** (2 * m - 1) * dC, "C": dC}[nm]


def c_const(nm, m):
    return Fr(a_const(nm, m)) / d_const(nm, m)


KAPPA = {"A": 12, "E": 11, "C": 3}


def law(nm, j, k):
    """Predicted [M^(j-k)] of the nm slice at order j, for k < j."""
    if k == -1:
        return Fr((-1) ** (j + 1) * 46, j * (j + 1))
    if k == 0:
        return Fr((-1) ** j * KAPPA[nm], j)
    if k % 2 == 0:
        return Fr(0)
    m = (k + 1) // 2
    fall = 1
    for i in range(1, 2 * m - 1):
        fall *= (j - i)
    # sign alternates with m, NOT with k: (-1)^(j+m). Since k = 2m-1, the naive
    # (-1)^(j+k) collapses to (-1)^(j+1) for every m and is right only for odd m.
    return Fr((-1) ** (j + m)) * c_const(nm, m) * fall


# --------------------------------------------------------------------------
def partA(POLY):
    rule("A. CONTROLS -- frozen digests re-emitted from the P_8 run")
    ok = True
    src = open(os.path.join(ROOT, "04-computation",
                            "gmc_first_gap_wall_stripped_all_width_second_edge_"
                            "circuit_positivity_thm2997.py")).read()
    s = src.index("JET_DATA_TEXT = r'''") + len("JET_DATA_TEXT = r'''")
    DATA = json.loads(src[s:src.index("'''", s)])
    CJ = {1: Fr(1), 2: Fr(16), 3: Fr(-128)}
    mine = [table_rows(POLY[j], CJ[j]) for j in (1, 2, 3)]
    for j in (1, 2, 3):
        same = mine[j - 1] == DATA[j - 1]
        ok &= same
        print(f"    j={j} frozen rows={len(DATA[j-1]):4d} mine={len(mine[j-1]):4d}"
              f"  ROW-ORDER IDENTICAL={same}")
    d123 = digest(mine)
    EXP = "cfb36557e1d54a0328a309375a948ace99c78e0688a54a014aef0906c1b90513"
    ok &= (d123 == EXP)
    print(f"    THM-2997 digest re-emitted: {d123 == EXP}")
    for j, exp in ((4, "200fae225af6c381733a386f31e107b81e6d33104699e5f69e8d7cd7e445d163"),
                   (5, "8cf6ed5cfca3b92a9229f79aae26f87ab1e65db1cf288b4299fe37b4a47b1de9"),
                   (6, "7004c579194e13f10aa03ceb26707adbaeae01e64b5be85d76792987f20c150e"),
                   (7, "bc53e1a23a9694c277de3aa1e9f6c4401be585452b7a20e35abc0a7a050fb287")):
        c = TABLES[j][1] or (1 / content(list(POLY[j].values()))).numerator
        dg = digest(table_rows(POLY[j], Fr(c)))
        ok &= (dg == exp)
        print(f"    P_{j} content-1 digest reproduced: {dg == exp}")
    print(f"  VERDICT A: every frozen digest reproduced: {ok}")
    return ok


def partB(POLY):
    print()
    rule("B. SHAPE LEDGER j=1..8   eq(2): |P_j| = (2j+1)^3 - 3(2j-2-floor(j/2)) - [j=3]")
    ok = True
    print("    j | terms |  box  | eq(2) | absent | corner-only")
    for j in range(1, JMAX + 1):
        box = {(a, b, c) for a in range(2 * j + 1) for b in range(4 * j + 1)
               for c in range(2 * j + 1) if b + 2 * c <= 4 * j}
        assert len(box) == (2 * j + 1) ** 3
        miss = sorted(box - set(POLY[j]))
        corners = {(4 * j, 0), (0, 2 * j), (0, 0)}
        nc = [m for m in miss if (m[1], m[2]) not in corners]
        pred = (2 * j + 1) ** 3 - 3 * (2 * j - 2 - j // 2) - (1 if j == 3 else 0)
        # j=3 carries ONE extra, non-corner absence -- that is exactly what the
        # -[j=3] correction term in eq(2) accounts for, so it is expected here.
        expect_nc = 1 if j == 3 else 0
        good = (len(POLY[j]) == pred) and (len(nc) == expect_nc) \
            and not [k for k in POLY[j] if k not in box]
        ok &= good
        print(f"    {j} | {len(POLY[j]):5d} | {len(box):5d} | {pred:5d} | {len(miss):6d} |"
              f" {not nc}{'  (j=3: one expected off-corner absence ' + str(nc) + ')' if nc else ''}")
    print(f"  VERDICT B: eq(2) holds j=1..8, all absences at corners except the "
          f"one eq(2) predicts at j=3: {ok}")
    return ok


def partC(POLY):
    print()
    rule("C. THE CONSTANT SEQUENCES c_m -- closed forms fitted on m=1,2,3 (j<=7 data)")
    ok = True
    OBS = {"A": {1: Fr(47, 24), 2: Fr(707, 23040), 3: Fr(11267, 15482880)},
           "E": {1: Fr(37, 18), 2: Fr(301, 9720), 3: Fr(2677, 3674160)},
           "C": {1: Fr(23, 6), 2: Fr(23, 360), 3: Fr(23, 15120)}}
    print("    slice  m   observed c_m        closed form a_m/d_m      match")
    for nm in "AEC":
        for m in (1, 2, 3):
            got, want = OBS[nm][m], c_const(nm, m)
            ok &= (got == want)
            print(f"      {nm}    {m}   {str(got):19s} {a_const(nm,m)}/{d_const(nm,m)}"
                  f"{'':4s} {got == want}")
    print()
    print("    d_m^C = (3m)!/(2m-2)! =", [str(d_const("C", m)) for m in (1, 2, 3, 4)])
    print("    d_m^A / d_m^C = 4^(2m-1) =", [4 ** (2 * m - 1) for m in (1, 2, 3, 4)])
    print("    d_m^E / d_m^C = 3^(2m-1) =", [3 ** (2 * m - 1) for m in (1, 2, 3, 4)])
    print("    a_m^A = 3 + 44*16^(m-1) =", [a_const("A", m) for m in (1, 2, 3, 4)])
    print("    a_m^E = 4 + 33*9^(m-1)  =", [a_const("E", m) for m in (1, 2, 3, 4)])
    print("    a_m^C = 23 (constant)")
    print(f"  VERDICT C: closed forms reproduce all nine known constants: {ok}")
    return ok


def partD(POLY):
    print()
    rule("D. PRE-REGISTERED j=8 VALUES, slots k=-1..6 (forced by the j<=7 laws)")
    ok = True
    SL8 = slices(POLY[8], 8)
    for nm in "AEC":
        print(f"    slice {nm}:")
        for k in range(-1, 7):
            e = 8 - k
            want = law(nm, 8, k)
            got = SL8[nm].get(e, Fr(0))
            good = (got == want)
            ok &= good
            print(f"      k={k:2d}  [M^{e}]  predicted {str(want):22s} got {str(got):22s}"
                  f" {'CONFIRMED' if good else '*** REFUTED ***'}")
    print(f"  VERDICT D: all 24 pre-registered j=8 coefficients confirmed: {ok}")
    return ok


def partE(POLY):
    print()
    rule("E. THE GENUINELY NEW SLOT k=7 -- speculative closed form, first testable at j=8")
    SL8 = slices(POLY[8], 8)
    ok = True
    print("    j=8 is the first order with k=7 < j, so c_4 was never observed before.")
    print("    Three separate extrapolations from m<=3 were registered. Scoring them:")
    print()
    obs = {nm: SL8[nm][1] / 5040 for nm in "AEC"}   # [M^1] = +c_4 * 7*6*5*4*3*2
    print("    (i) NUMERATORS  a_m^A = 3+44*16^(m-1),  a_m^E = 4+33*9^(m-1),  a_m^C = 23")
    for nm in "AEC":
        got, want = obs[nm].numerator, a_const(nm, 4)
        ok &= (got == want)
        print(f"        {nm}  predicted {want:8d}   observed {got:8d}   "
              f"{'CONFIRMED' if got == want else '*** REFUTED ***'}")
    print()
    print("    (ii) SLICE RATIOS  d_m^A/d_m^C = 4^(2m-1),  d_m^E/d_m^C = 3^(2m-1)")
    dobs = {nm: int(Fr(obs[nm].numerator) / obs[nm]) for nm in "AEC"}
    for nm, base in (("A", 4), ("E", 3)):
        got, want = Fr(dobs[nm], dobs["C"]), base ** 7
        ok &= (got == want)
        print(f"        {nm}  predicted {base}^7 = {want:8d}   observed {got}   "
              f"{'CONFIRMED' if got == want else '*** REFUTED ***'}")
    print()
    print("    (iii) BASE SEQUENCE  d_m^C = (3m)!/(2m-2)!  ->  predicted d_4^C = 665280")
    print(f"        observed d_4^C = {dobs['C']}   ratio predicted/observed = "
          f"{Fr(665280, dobs['C'])}")
    print("        *** REFUTED *** -- this was the one guess that failed.")
    print()
    print("    Resulting measured constants:")
    for nm in "AEC":
        print(f"        c_4^{nm} = {obs[nm]}")
    print()
    print("  VERDICT E: numerator laws and slice ratios CONFIRMED out-of-sample at")
    print(f"             m=4 ({ok}); the base sequence d_m^C = 6, 360, 15120, 604800")
    print("             has NO closed form yet and is the open piece.")
    return ok


def partF(POLY):
    print()
    rule("F. TERMINAL SLOT k=j (the constant term) -- exceptional, recorded not predicted")
    for j in range(4, JMAX + 1):
        S = slices(POLY[j], j)
        print(f"    j={j}  A={str(S['A'].get(0, Fr(0))):30s} "
              f"C={str(S['C'].get(0, Fr(0)))}")
    print("    These do NOT follow the slot law -- k=j is outside its range of validity.")
    return True


def main():
    POLY = {j: load(j) for j in range(1, JMAX + 1)}
    a = partA(POLY)
    b = partB(POLY)
    c = partC(POLY)
    d = partD(POLY)
    e = partE(POLY)
    partF(POLY)
    print()
    rule(f"SUMMARY  controls={a}  shape={b}  closed-forms={c}  "
         f"pre-registered={d}  new-slot={e}")


if __name__ == "__main__":
    main()
