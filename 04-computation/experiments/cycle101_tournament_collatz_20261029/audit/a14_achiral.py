#!/usr/bin/env python3
"""THM-4849 (B) 'two base points iff reversal-symmetric': does it hold for all shape-(k,q) cycles, or only for
tournament cycles?  Count achiral necklaces with 0 base points."""
from itertools import combinations
def canon(S, q): return min(tuple(sorted((x + t) % q for x in S)) for t in range(q))
def is_cs(S, q): return 0 not in S and all((x in S) != ((-x) % q in S) for x in range(1, q))
for q in (5, 7, 11, 13, 17, 19):
    k = (q - 1) // 2
    necks = {canon(c, q) for c in combinations(range(q), k)}
    stats = {}
    for N in necks:
        nb = sum(1 for t in range(q) if is_cs({(x + t) % q for x in N}, q))
        ach = canon({(-x) % q for x in N}, q) == N
        stats[(nb, ach)] = stats.get((nb, ach), 0) + 1
    ex = None
    for N in sorted(necks):
        if canon({(-x) % q for x in N}, q) == N and not any(is_cs({(x + t) % q for x in N}, q) for t in range(q)):
            s = sorted(N); ex = tuple(((s[(i + 1) % k] - s[i]) % q) for i in range(k)); break
    print(f"q={q}: (base points, achiral) -> necklaces {dict(sorted(stats.items()))}; first achiral word with no base point: {ex}")
