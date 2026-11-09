"""Auditor's check of HYP-9540 (B)'s construction for the pair (1, k): reverse the single tile (0, 3).
Claim: nu(t) = 1 and nu(t + 1) = k = floor(n/3) for every n >= 4; the line is black for n >= 5."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from staircase_audit import tourn, all_tiles, nu, blue

bad = []
for n in range(4, 21):
    rev = {(0, 3)}
    t, u = tourn(n, rev), tourn(n, all_tiles(n) - rev)
    a, b = nu(t, n), nu(u, n)
    col = "blue" if blue(n, rev) else "black"
    if (a, b) != (1, n // 3) or (col == "blue") != (n == 4):
        bad.append((n, a, b, col))
print(f"single reversed tile (0,3), n = 4..20: nu-pair (1, floor(n/3)), black for n >= 5: failures {bad}")
