#!/usr/bin/env python3
"""audit_wave2: section 1 and section 3 of monster_levels_square_sum_7_21_20261022.md.

The sporadic orders are rebuilt from the auditor's own recollection of their prime factorisations
(ATLAS), multiplied out, and compared with the decimal constants hard-coded in level_ogg_lists.gp
(read from that file as data). Then |M|/|B|, the owner's number, the valuation classes, the exponent
table and the chain sizes are recomputed from the exponent vectors alone.
Run: python3 -B audit_group_orders.py > audit_group_orders.out
"""
import math
import re
import sys
from pathlib import Path

from sympy import factorint

P = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]
EXP = {  # auditor's recollection (ATLAS)
    "M": {2: 46, 3: 20, 5: 9, 7: 6, 11: 2, 13: 3, 17: 1, 19: 1, 23: 1, 29: 1, 31: 1, 41: 1, 47: 1, 59: 1, 71: 1},
    "B": {2: 41, 3: 13, 5: 6, 7: 2, 11: 1, 13: 1, 17: 1, 19: 1, 23: 1, 31: 1, 47: 1},
    "Fi24'": {2: 21, 3: 16, 5: 2, 7: 3, 11: 1, 13: 1, 17: 1, 23: 1, 29: 1},
    "HN": {2: 14, 3: 6, 5: 6, 7: 1, 11: 1, 19: 1},
    "Fi22": {2: 17, 3: 9, 5: 2, 7: 1, 11: 1, 13: 1},
    "He": {2: 10, 3: 3, 5: 2, 7: 3, 17: 1},
    "HS": {2: 9, 3: 2, 5: 3, 7: 1, 11: 1},
    "M12": {2: 6, 3: 3, 5: 1, 11: 1},
    "L3(3)": {2: 4, 3: 3, 13: 1},
    "Th": {2: 15, 3: 10, 5: 3, 7: 2, 13: 1, 19: 1, 31: 1},
    "Co1": {2: 21, 3: 9, 5: 4, 7: 2, 11: 1, 13: 1, 23: 1},
}
order = {g: math.prod(p ** e for p, e in d.items()) for g, d in EXP.items()}
fails = []

src = Path(sys.argv[1] if len(sys.argv) > 1 else "../level_ogg_lists.gp").read_text()
consts = dict(re.findall(r'\["([^"]+)",\s*(\d+)\]', src))
print("== 1. sporadic orders: auditor's factorisations against the decimals in level_ogg_lists.gp ==")
for g in order:
    if g in consts:
        same = int(consts[g]) == order[g]
        fails += [] if same else [g]
        print(f"  {g:6s} {order[g]}  {'= script' if same else 'DIFFERS from script ' + consts[g]}")
    else:
        print(f"  {g:6s} {order[g]}  (not in the script)")
    if factorint(order[g]) != EXP[g]:
        fails.append("factor " + g)

print("\n== 2. |M|/|B|, the owner's number, the 2A class ==")
r = order["M"] // order["B"]
assert order["M"] % order["B"] == 0
fr = factorint(r)
print("  |M|/|B| =", " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(fr.items())))
num = 2 ** 5 * 3 ** 7 * 5 ** 3 * 7 ** 4 * 11 * 13 ** 2
print("  2^5*3^7*5^3*7^4*11*13^2 =", num, "; |M|/|B| / (29*41*59*71) =", r // (29 * 41 * 59 * 71), "; equal:", r // (29 * 41 * 59 * 71) == num and r % (29 * 41 * 59 * 71) == 0)
fails += [] if num == 39046339332000 and r == num * 29 * 41 * 59 * 71 else ["owner number"]
print("  |2A| = |M|/|2.B| = |M|/(2|B|) =", r // 2, "; so |M|/|B| = 2 |2A| (B is a section of M, not a subgroup: C_M(2A) = 2.B)")
print("\n== 3. valuation classes ==")
cl = {"v_p(M) >= 2": [p for p in P if EXP["M"][p] >= 2],
      "v_p(M) = v_p(B) = 1": [p for p in P if EXP["M"][p] == 1 and EXP["B"].get(p, 0) == 1],
      "v_p(M) = 1, v_p(B) = 0": [p for p in P if EXP["M"][p] == 1 and EXP["B"].get(p, 0) == 0]}
for k, v in cl.items():
    print(f"  {k:24s}: {v}")
fails += [] if list(cl.values()) == [[2, 3, 5, 7, 11, 13], [17, 19, 23, 31, 47], [29, 41, 59, 71]] else ["classes"]
print("  primes of Co1 inside primes of B:", set(EXP["Co1"]) <= set(EXP["B"]), "; primes of Fi22:", sorted(EXP["Fi22"]))
print("\n== 4. exponent table at 2, 3, 5, 7, 11, 13 ==")
sm = [2, 3, 5, 7, 11, 13]
note = {"B": (5, 7, 3, 4, 1, 2), "Fi24'": (25, 4, 7, 3, 1, 2), "Fi22": (29, 11, 7, 5, 1, 2), "Th": (31, 10, 6, 4, 2, 2),
        "HN": (32, 14, 3, 5, 1, 3), "He": (36, 17, 7, 3, 2, 3), "Co1": (25, 11, 5, 4, 1, 2)}
for g, t in note.items():
    v = tuple(EXP["M"][p] - EXP[g].get(p, 0) for p in sm)
    print(f"  |M|/|{g}|: {v} {'= note' if v == t else 'DIFFERS from note ' + str(t)}")
    fails += [] if v == t else ["exp " + g]
v = tuple(EXP["B"].get(p, 0) - EXP["Fi22"].get(p, 0) for p in sm)
print(f"  |B|/|Fi22|: {v} {'= note' if v == (24, 4, 4, 1, 0, 0) else 'DIFFERS'}")
# the 'cycle' reading p -> v_p(|M|/|B|)
m = {p: EXP["M"][p] - EXP["B"].get(p, 0) for p in sm}
print("  p -> v_p(|M|/|B|):", m, "; values form a permutation of", sorted(m.values()))
print("\n== 5. level-chain difference sizes (lists from the note's table, with the primes of the level) ==")
full = {1: set(P), 2: set(EXP["B"]), 3: set(EXP["Fi24'"]), 5: set(EXP["HN"]), 6: set(EXP["Fi22"]), 7: set(EXP["He"]),
        10: set(EXP["HS"]), 21: {2, 3, 5, 7}}
for c in ([1, 2, 6], [1, 2, 10], [1, 3, 6], [1, 7, 21]):
    a, b, d = (full[x] for x in c)
    print(f"  {c}: |L1 minus L2|, |L2 minus L3|, |L3| = {len(a - b)}, {len(b - d)}, {len(d)}; nested: {b <= a and d <= b}")
print("\nfailures:", fails)
if not fails:
    print("ALL AUDIT CHECKS PASSED")
