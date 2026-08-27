#!/usr/bin/env python3
"""Emit every strong fixed-gauge two-reversal presentation at one order."""

import sys


def arcs(n):
    return tuple((i, j) for i in range(n) for j in range(i + 1, n))


def shell(n):
    z = n - 1
    rows = set()
    for edge in arcs(n):
        if edge not in ((0, z), (0, 1), (z - 1, z)):
            rows.add(tuple(sorted(((0, z), edge))))
    for b in range(1, z):
        for c in range(1, b + 1):
            rows.add(tuple(sorted(((0, b), (c, z)))))
    return tuple(sorted(rows))


def encode(n, reversals):
    reverse = frozenset(reversals)
    return "".join("0" if edge in reverse else "1" for edge in arcs(n))


if len(sys.argv) != 2:
    raise SystemExit("usage: tournament_two_reversal_words_thm4239.py ORDER")
n = int(sys.argv[1])
if n < 3:
    raise SystemExit("ORDER must be at least 3")
rows = shell(n)
if len(rows) != (n - 1) ** 2 - 3:
    raise RuntimeError("presentation count mismatch")
for reversals in rows:
    print(encode(n, reversals))
