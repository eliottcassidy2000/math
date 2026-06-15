#!/usr/bin/env python3
"""Lean driver for ocf_two_dimensions_monad.run with smaller sample/cap so the
pure-python phase-2 (H_dp + cycle enum) finishes. Rank of level-sums stabilizes
quickly, so modest sampling suffices to confirm dim(H)=floor(n/3)."""
from ocf_two_dimensions_monad import run
import sys

plan = {9: (70000, 40), 10: (250000, 12), 11: (700000, 6), 12: (1500000, 4)}
res = []
for n in [int(a) for a in sys.argv[1:]]:
    s, c = plan[n]
    res.append((n,) + run(n, s, c))
print("\nSUMMARY  n  carrier-dim  H-func-dim  q(n)-3")
for n, rc, rl, g in res:
    print(f"  {n}  {rc}  {rl}  {g}")
