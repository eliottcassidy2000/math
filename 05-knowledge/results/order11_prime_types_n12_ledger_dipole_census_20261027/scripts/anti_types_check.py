#!/usr/bin/env python3
"""Which cycle types can an anti-automorphism of an n-tournament have?  For each partition of n
into cycle lengths, build psi and propagate u->v  =>  psi(v)->psi(u) around every pair orbit;
a contradiction (the orbit returns to the same pair reversed) rules the type out.  Prints the
admissible types and checks the law: all cycles of length = 2 (mod 4) plus at most one fixed point."""
import sys
def partitions(n, m=None):
    if m is None: m = n
    if n == 0: yield []; return
    for k in range(min(n, m), 0, -1):
        for rest in partitions(n-k, k): yield [k] + rest
def consistent(n, lens):
    psi = list(range(n)); start = 0
    for L in lens:
        for i in range(L): psi[start+i] = start + (i+1) % L
        start += L
    seen = set()
    for u in range(n):
        for v in range(u+1, n):
            if (u, v) in seen: continue
            x, y = u, v
            while True:
                k = (min(x,y), max(x,y))
                if k in seen:
                    if (x, y) != (u, v): return False
                    break
                seen.add(k); x, y = psi[y], psi[x]
    return True
for n in range(2, int(sys.argv[1]) + 1):
    ok = [p for p in partitions(n) if consistent(n, p)]
    law = [p for p in partitions(n) if sum(1 for L in p if L == 1) <= 1 and all(L == 1 or L % 4 == 2 for L in p)]
    print(n, 'admissible:', ok, 'law holds:', ok == law)
