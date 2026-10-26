# Cross-check against the vitali-response atlas table: total strong presentations and the largest k = 1 target fibre.
import sys
from collections import Counter
from rooted import tournaments, canon_rooted, ham_paths, partner, is_strong, relabel
for n in (3, 4, 5, 6, 7):
    tot = 0; mx = 0
    for a in tournaments(n):
        if not is_strong(a): continue
        P = ham_paths(a)
        pres = Counter()
        for p in P:
            perm = [0] * n
            for i, v in enumerate(p): perm[v] = i
            pres[relabel(a, perm)] += 1
        aut = len(P) // len(pres); tot += len(pres)
        F = Counter(canon_rooted(partner(a, p), p[0], p[-1]) for p in P)
        mx = max(mx, max(c // aut for c in F.values()))
    print(f"n = {n}: strong presentations {tot}, largest joint target fibre {mx}")
