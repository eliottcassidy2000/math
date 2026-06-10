#!/usr/bin/env python3
# The cubic rung of the additive-relation ladder (HYP-2376; THM-446 frame;
# complements kind-pasteur THM-462/463 + HYP-2370 three-cubes ledger — cited, not duplicated).
# mac-mini-2026-06-10-S1
#
# Placement of CUBES = {1,8,27,...} on the ladder (THM-446: Schur 3-term ⊂ Sidon 4-term ⊂ B_h):
#  * 3-TERM (Schur/cauldron): cubes are SUM-FREE — a^3+b^3=c^3 has no positive solutions
#    (Euler / FLT n=3). Verified in range as sanity. The cubic cauldron never boils.
#  * 4-TERM (Sidon/B_2): FAILS first at 1729 = 1^3+12^3 = 9^3+10^3 (Ramanujan-Hardy;
#    = first C4 of the cubic summand graph = first unit of additive energy above the
#    Sidon floor, in THM-441 language E(S)=||1_S * 1_S||^2).
#  * SIGNED 4-term: x^3+y^3+z^3=w^3 solvable small: 3,4,5,6 (and families).
# Census: taxicab numbers (two-representation sums) up to N^3 bound; additive-energy
# excess growth; signed quadruple count.

import itertools
from collections import Counter

B = 300  # cubes 1..B
cubes = [i**3 for i in range(1, B+1)]
cubeset = set(cubes)

print("=== 3-term: sum-free check (Euler/FLT3 sanity) ===")
viol = 0
for a in range(1, B+1):
    for b in range(a, B+1):
        if a**3 + b**3 in cubeset:
            viol += 1
print(f"   a^3+b^3=c^3 solutions with a,b <= {B}: {viol}   (FLT n=3 => 0)")

print("=== 4-term: taxicab census (pair-sums with >=2 representations) ===")
sums = Counter()
for i in range(len(cubes)):
    for j in range(i, len(cubes)):
        sums[cubes[i]+cubes[j]] += 1
multi = sorted([(s, c) for s, c in sums.items() if c >= 2])
print(f"   pair-sums up to 2*{B}^3: {len(sums)} distinct; with >=2 reps: {len(multi)}")
print(f"   FIRST: {multi[0][0]}  (should be 1729)")
print("   first 8 taxicab numbers:", [s for s, c in multi[:8]])
ta3 = [s for s, c in multi if c >= 3]
print(f"   first with 3 representations (Ta(3)=87539319 if in range): {ta3[0] if ta3 else 'beyond range'}")
# additive energy excess (THM-441): E(S) - (2|S|^2 - |S|) counts ordered quadruple excess
S = cubes
E = 0
for s, c in sums.items():
    # ordered representations: each unordered pair {a,b}, a!=b gives 2; a=b gives 1... reconstruct:
    pass
# direct ordered energy
osums = Counter()
for a in cubes:
    for b in cubes:
        osums[a+b] += 1
E = sum(c*c for c in osums.values())
floor = 2*len(S)**2 - len(S)
print(f"   additive energy E = {E}, Sidon floor = {floor}, EXCESS = {E - floor}")
for Bs in (50, 100, 200, 300):
    cs = [i**3 for i in range(1, Bs+1)]
    oc = Counter()
    for a in cs:
        for b in cs:
            oc[a+b] += 1
    e = sum(c*c for c in oc.values())
    fl = 2*Bs*Bs - Bs
    print(f"   B={Bs}: excess = {e-fl}  (excess/B: {(e-fl)/Bs:.3f})")

print("=== signed 4-term: x^3+y^3+z^3 = w^3, 0<x<=y<=z<w<=120 ===")
qs = []
for w in range(2, 121):
    w3 = w**3
    for x in range(1, w):
        for y in range(x, w):
            r = w3 - x**3 - y**3
            if r < y**3:
                break
            z = round(r ** (1/3))
            for zz in (z-1, z, z+1):
                if zz >= y and zz < w and zz**3 == r:
                    qs.append((x, y, zz, w))
print(f"   solutions: {len(qs)}; first: {qs[:5]}")
print("   (3,4,5,6) present:", (3,4,5,6) in qs)
print()
print("LADDER PLACEMENT: cubes are sum-free (3-term rung: SAFE forever, Euler/FLT3)")
print("but non-Sidon from 1729 (4-term rung: first C4 of the cubic summand graph);")
print("signed relations enter at (3,4,5,6). Cross-refs: THM-446 ladder, THM-441 energy,")
print("kind-pasteur THM-463 (divisor criterion for two-cube reps = the split axis),")
print("HYP-2370 ledger (mod-9 law; 2026 open-list: smallest open k=114; 33/42 fell 2019).")
