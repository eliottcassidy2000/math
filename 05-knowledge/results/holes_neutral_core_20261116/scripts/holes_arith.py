#!/usr/bin/env python3
"""holes_arith.py: the borrowed-response pencils through the holes and the 7p lane against the census of counts to order 9
(counts.txt of the toothpick packet)."""
S = set(int(l.split()[1]) for l in open('counts.txt'))
print("counts <= 60 realised at orders <= 9:", sorted(c for c in S if c <= 60))
print("pencil 1 + 2^k*3 (base count 1 with the borrowed response 3):", [(k, 1 + 2**k * 3, (1 + 2**k * 3) in S) for k in range(1, 7)])
print("pencil 3 + 2^k*9 (base count 3 with the borrowed response 9):", [(k, 3 + 2**k * 9, (3 + 2**k * 9) in S) for k in range(1, 7)])
print("7p lane p + 6p for p = 1, 3, 5, 7, 9, 11, 13:", [(p, 7 * p, (7 * p) in S) for p in (1, 3, 5, 7, 9, 11, 13)])
