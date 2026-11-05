#!/usr/bin/env python3
"""Profile fibers (Session 200, carmichael-fibers) versus strong-class multiplicity (THM-4884 (2)).

Session 200 counts, for the physical source grammar on responses s = (D, E, U, d) starting at t = (1, 0, 2, 1)
with chronological moves

    O(s) = (D, 0, 2U, d),  Z(s) = (D - E + U, d, 2U - d, d),  J(s) = (D - E + U, D - E, D + U, D),

the number pf(n) of distinct numerical profiles s with C1 readout F(s) = D - E + U = n (closing by J gives a
strong tournament with n Hamiltonian paths). THM-4884 (2) counts m(h) = number of strong tournament classes of all
orders with H = h, read from THM-4744's inventories (strong_ears_375.out / strong_ears_625.out, lines "P_d(1)").
This script enumerates the complete physical inventory with U <= CAP (U is strictly increasing along the grammar
since U > d), tabulates pf(n) for odd n <= CAP + 1, and prints the comparison with m(h).

usage: python3 profile_vs_class_multiplicity_20261105.py CAP strong_ears_375.out strong_ears_625.out
"""
import re
import sys

CAP = int(sys.argv[1]) if len(sys.argv) > 1 else 624
EARS = sys.argv[2:]


def moves(s):
    D, E, U, d = s
    yield (D, 0, 2 * U, d)                      # O
    yield (D - E + U, d, 2 * U - d, d)          # Z
    yield (D - E + U, D - E, D + U, D)          # J


t = (1, 0, 2, 1)
seen = {t}
stack = [t]
while stack:
    s = stack.pop()
    for s2 in moves(s):
        if s2[2] <= CAP and s2 not in seen:
            seen.add(s2)
            stack.append(s2)

pf = {}
for (D, E, U, d) in seen:
    n = D - E + U
    pf[n] = pf.get(n, 0) + 1

limit = CAP + 1
print("profiles with U <= %d: %d" % (CAP, len(seen)))
odd = [n for n in range(3, limit + 1, 2)]
missed = [n for n in odd if pf.get(n, 0) == 0]
single = [n for n in odd if pf.get(n, 0) == 1]
print("odd n <= %d with pf(n) = 0 (not a C1 readout of a physical source): %s" % (limit, missed))
print("odd n <= %d with pf(n) = 1: %s" % (limit, single))

m = {}
for path in EARS:
    txt = open(path, encoding="utf-8", errors="replace").read()
    # lines of the form "P_d(1), odd d = A..B: [...]" (375 file) or "P_d(1) for odd d = A..B: [...]" (625 file)
    for a, b, body in re.findall(r"P_d\(1\)(?:,| for) odd d = (\d+)\.\.(\d+): \[([^\]]*)\]", txt):
        vals = [int(x) for x in re.findall(r"-?\d+", body)]
        for i, v in enumerate(vals):
            m[int(a) + 2 * i] = v

if m:
    print("\nh : pf(h) [profiles, C1 readout] : m(h) [strong classes, all orders]  (odd h <= 101)")
    for h in range(3, 102, 2):
        print("%3d : %4d : %4d" % (h, pf.get(h, 0), m.get(h, 0)))
    both = [h for h in odd if h in m]
    print("\nodd h <= %d with m(h) = 1 and pf(h) >= 2: %s" % (limit, [h for h in both if m[h] == 1 and pf.get(h, 0) >= 2]))
    print("odd h <= %d with pf(h) = 1 and m(h) >= 2: %s" % (limit, [h for h in both if pf.get(h, 0) == 1 and m[h] >= 2]))
    print("odd h <= %d with pf(h) = 0 and m(h) >= 1: %s" % (limit, [h for h in both if pf.get(h, 0) == 0 and m[h] >= 1]))
    print("odd h <= %d with pf(h) > m(h): %s" % (limit, [h for h in both if pf.get(h, 0) > m[h]]))
