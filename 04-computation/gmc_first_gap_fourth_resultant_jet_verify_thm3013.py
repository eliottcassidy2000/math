#!/usr/bin/env python3
"""THM-3013 -- independent verification of the symbolic fourth resultant log-jet
P_4(M,U,V) of the first-gap family.

The table itself lives in
  05-knowledge/results/gmc_first_gap_fourth_resultant_jet_P4_table_thm3013.json
in the frozen THM-2997 JET_DATA_TEXT row format [m_power, u_power, v_power,
numerator, denominator], with c_4 = 1658880 (integer coefficients, content 1).
The family-normalised copy (c_4 = 61440, denominators {1,3,9,27}) is
  ..._P4_table_family_thm3013.json.

THIS SCRIPT DOES NOT REBUILD P_4.  It checks the properties that can be verified
from the table plus the FROZEN THM-2997 table alone:
  1. shape: 717 rows, degrees (M,U,V) = (8,16,8), support exactly b+2c <= 16,
     integer coefficients, and the 12 absences all at the three corners
     (16,0), (0,8), (0,0);
  2. the digest of the canonical JSON;
  3. that the top U-column [U^16] P_4 / c_4 equals THM-3011's A_4(M) exactly --
     INCLUDING the two coefficients 707/3840 and -15797937/128 that the
     frozen-table pattern of THM-3011 section 2a could NOT predict;
  4. that THM-3011 section 2a's three coefficient laws, extracted from the FROZEN
     j=1,2,3 table only, correctly predict A_4's top three coefficients;
  5. the corrected corner statement: the "all absences are corners" reading is
     TRUE at j=4 but FALSE at j=3, where (5,0,3) is absent and (0,3) is not a
     corner.

Reproduce: python3 04-computation/gmc_first_gap_fourth_resultant_jet_verify_thm3013.py
"""

import hashlib
import importlib.util
import json
from fractions import Fraction as Fr

TABLE = "05-knowledge/results/gmc_first_gap_fourth_resultant_jet_P4_table_thm3013.json"
DIGEST = "200fae225af6c381733a386f31e107b81e6d33104699e5f69e8d7cd7e445d163"
C4 = 1658880
A4_REF = {5: Fr(-23, 10), 4: Fr(3), 3: Fr(-47, 24), 1: Fr(707, 3840),
          0: Fr(-15797937, 128)}


def rule(s):
    print("=" * 74)
    print(s)
    print("=" * 74)


def frozen():
    spec = importlib.util.spec_from_file_location(
        "t", "04-computation/gmc_first_gap_wall_stripped_all_width_second_edge_circuit_positivity_thm2997.py")
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except SystemExit:
        pass
    except Exception:
        pass
    return json.loads(mod.JET_DATA_TEXT)


def main():
    rows = json.load(open(TABLE))
    ok = True

    rule("1. SHAPE AND SUPPORT")
    a = [r[0] for r in rows]
    b = [r[1] for r in rows]
    c = [r[2] for r in rows]
    shape = (len(rows) == 717 and max(a) == 8 and max(b) == 16 and max(c) == 8)
    viol = sum(1 for x, y in zip(b, c) if x + 2 * y > 16)
    ints = all(r[4] == 1 for r in rows)
    print(f"  rows {len(rows)}  M-deg 0..{max(a)}  U-deg 0..{max(b)}  V-deg 0..{max(c)}"
          f"  max(b+2c) {max(x + 2 * y for x, y in zip(b, c))}")
    print(f"  integer coefficients: {ints}   support violations b+2c>16: {viol}")
    present = {(r[0], r[1], r[2]) for r in rows}
    full = {(m, bb, cc) for m in range(9) for cc in range(9) for bb in range(17 - 2 * cc)}
    missing = sorted(full - present)
    corners = {(16, 0), (0, 8), (0, 0)}
    allc = all((m[1], m[2]) in corners for m in missing)
    print(f"  ansatz {len(full)}  present {len(present)}  missing {len(missing)}")
    print(f"  every absence at a corner (16,0)/(0,8)/(0,0): {allc}")
    print(f"  absences: {missing}")
    ok &= shape and ints and viol == 0 and len(missing) == 12 and allc
    print(f"  VERDICT 1: {'OK' if shape and ints and viol == 0 and allc else 'FAILED'}")

    print()
    rule("2. DIGEST")
    dig = hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()
    print(f"  sha256(canonical json) = {dig}")
    print(f"  matches recorded digest: {dig == DIGEST}")
    ok &= (dig == DIGEST)

    print()
    rule("3. TOP U-COLUMN EQUALS THM-3011's A_4, IN FULL")
    col = {r[0]: Fr(r[3], r[4]) for r in rows if r[1] == 16 and r[2] == 0}
    A4 = {m: v / C4 for m, v in col.items()}
    match = (A4 == A4_REF)
    print(f"  [U^16] P_4 / c_4 = {{{', '.join(f'M^{m}: {A4[m]}' for m in sorted(A4, reverse=True))}}}")
    print(f"  THM-3011 A_4      = {{{', '.join(f'M^{m}: {A4_REF[m]}' for m in sorted(A4_REF, reverse=True))}}}")
    print(f"  EXACT MATCH (including 707/3840 and -15797937/128): {match}")
    print("  Those two were NOT predictable from the frozen-table pattern, so this is")
    print("  a genuine independent confirmation of THM-3011 section 2's reconstruction.")
    ok &= match

    print()
    rule("4. THM-3011 SECTION 2a LAWS, FROM THE FROZEN TABLE ONLY")
    data = frozen()
    cj = {1: Fr(1), 2: Fr(16), 3: Fr(-128)}
    A = {j: {r[0]: Fr(r[3], r[4]) / cj[j] for r in data[j - 1] if r[1] == 4 * j and r[2] == 0}
         for j in (1, 2, 3)}
    laws = {"[M^(j+1)] = (-1)^(j+1) 46/(j(j+1))": (lambda j: Fr((-1) ** (j + 1) * 46, j * (j + 1)), lambda j: j + 1),
            "[M^j]     = (-1)^j 12/j": (lambda j: Fr((-1) ** j * 12, j), lambda j: j),
            "[M^(j-1)] = (-1)^(j+1) 47/24": (lambda j: Fr((-1) ** (j + 1) * 47, 24), lambda j: j - 1)}
    lawok = True
    for name, (f, pw) in laws.items():
        res = []
        for j in (2, 3):
            res.append(A[j].get(pw(j)) == f(j))
        lawok &= all(res)
        print(f"  {name:36s} j=2,3: {res}   predicts A_4[M^{pw(4)}] = {f(4)}"
              f"   actual {A4.get(pw(4))}   {'OK' if A4.get(pw(4)) == f(4) else 'MISMATCH'}")
        lawok &= (A4.get(pw(4)) == f(4))
    ok &= lawok
    print(f"  VERDICT 4: {'laws hold and predict A_4 top three' if lawok else 'FAILED'}")

    print()
    rule("5. CORRECTED CORNER STATEMENT")
    p3 = {(r[0], r[1], r[2]) for r in data[2]}
    full3 = {(m, bb, cc) for m in range(7) for cc in range(7) for bb in range(13 - 2 * cc)}
    miss3 = sorted(full3 - p3)
    corners3 = {(12, 0), (0, 6), (0, 0)}
    noncorner3 = [m for m in miss3 if (m[1], m[2]) not in corners3]
    print(f"  P_3 absences: {miss3}")
    print(f"  non-corner absences at j=3: {noncorner3}")
    print("  => 'the absences are always the three corners' is TRUE at j=4 but FALSE at")
    print("     j=3, where (5,0,3) is absent and (0,3) is not a corner.  Record the")
    print("     enumerated lists, not the generalisation.")
    ok &= (len(noncorner3) > 0)

    print()
    rule(f"SUMMARY  all checks pass: {ok}")


if __name__ == "__main__":
    main()
