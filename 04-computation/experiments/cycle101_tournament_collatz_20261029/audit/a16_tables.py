#!/usr/bin/env python3
"""Compare the audit's own reduced denominators and U-traces of every tournament cycle at q = 7, 11, 13 with the
session's out/circulant_collatz.out (THM-4849 FINITE-EXACT dictionary tables)."""
import re, os
from fractions import Fraction as Fr
from collections import Counter
here = os.path.dirname(os.path.abspath(__file__))
def C_of(w):
    k = len(w); s = 0; t = 0
    for j in range(k):
        t += 3 ** (k - 1 - j) * 2 ** s; s += w[j]
    return t
def canon(S, q): return min(tuple(sorted((x + t) % q for x in S)) for t in range(q))
mine = {}
for q in (7, 11, 13):
    k = (q - 1) // 2
    seen = {}
    for mask in range(2 ** k):
        S = frozenset((x if (mask >> (x - 1)) & 1 == 0 else q - x) for x in range(1, k + 1))
        seen.setdefault(canon(S, q), S)
    for neck in seen:
        s = sorted(neck); w = tuple(((s[(i + 1) % k] - s[i]) % q) or q for i in range(k))
        x = Fr(C_of(w), 2 ** q - 3 ** k)
        trace = sum(Fr(C_of(w[i:] + w[:i]), 2 ** q - 3 ** k) for i in range(k))
        mine[(q, neck)] = (x.denominator, trace)
theirs = {}; q = None
for line in open(os.path.join(here, "..", "out", "circulant_collatz.out")):
    m = re.match(r"q=(\d+):", line)
    if m: q = int(m.group(1)); continue
    m = re.search(r"class \(([\d, ]+)\).*reduced denominator (\d+) .*\| sum ([\d/]+) \|", line)
    if m and q in (7, 11, 13):
        neck = tuple(int(t) for t in m.group(1).replace(" ", "").split(",") if t != "")
        theirs[(q, neck)] = (int(m.group(2)), Fr(m.group(3)))
print("tournament cycles compared:", len(mine), "; session entries:", len(theirs),
      "; reduced denominators and U-traces agree:", all(mine[key] == theirs.get(key) for key in mine))
for q in (11, 13):
    print(f"q={q}: reduced denominators {dict(Counter(v[0] for kk, v in mine.items() if kk[0] == q))}")
