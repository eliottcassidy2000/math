#!/usr/bin/env python3
"""Do the two Gersonides carries (unit (b+1,b)~(b,b+2), Catalan (b+3,*,b)~(b,*,b+6)) derive the primitive collisions at -1?  For each primitive pair, BFS from the padded short side.  Usage: python3 carry_derivation.py E CAP"""
# For each primitive pair (s, t) of canonical increasing representations with K(s) = K(t), test whether the padded s
# (leading 1's to the length of t) reaches t by unit/Catalan carries (both directions, any positions), exponents <= CAP.
import itertools, sys
from collections import defaultdict, deque
E, CAP = int(sys.argv[1]), int(sys.argv[2])
def K_of(s):
    K = -1
    for e in s: K = 3 * K + (1 << e)
    return K
def nbrs(B):
    m = len(B); out = []
    for j in range(m - 1):
        x, y = B[j], B[j + 1]
        if x == y + 1 and y + 2 <= CAP: out.append(B[:j] + (y, y + 2) + B[j + 2:])
        if y == x + 2 and x >= 0: out.append(B[:j] + (x + 1, x) + B[j + 2:])
    for j in range(m - 2):
        x, y = B[j], B[j + 2]
        if x == y + 3 and y + 6 <= CAP: out.append(B[:j] + (y,) + B[j + 1:j + 2] + (y + 6,) + B[j + 3:])
        if y == x + 6: out.append(B[:j] + (x + 3,) + B[j + 1:j + 2] + (x,) + B[j + 3:])
    return out
reps = defaultdict(list)
for n in range(E + 1):
    for s in itertools.combinations(range(2, E + 1), n):
        reps[K_of(s)].append(s)
ok = bad = 0; badlist = []
for k, rs in reps.items():
    if k == -1: continue
    for s, t in itertools.combinations(rs, 2):
        if s[-1] == t[-1]: continue
        if len(s) > len(t): s, t = t, s
        src = (1,) * (len(t) - len(s)) + s
        seen = {src}; dq = deque([src]); found = False
        while dq and not found:
            B = dq.popleft()
            for C in nbrs(B):
                if C not in seen:
                    if C == t: found = True; break
                    seen.add(C); dq.append(C)
            if len(seen) > 200000: break
        if found: ok += 1
        else:
            bad += 1
            if len(badlist) < 6: badlist.append((k, s, t, len(seen)))
print(f"E <= {E}, carries within exponent {CAP}: {ok} primitive pairs derivable from the padded short side, {bad} not")
for b in badlist: print("   not derivable:", b)
