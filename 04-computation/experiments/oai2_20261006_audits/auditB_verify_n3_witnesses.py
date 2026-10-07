#!/usr/bin/env python3
"""Independent verifier (audit B) for the reader's fully gap-determined witnesses E on [t]^3
(finv_witness_n3_t{4,5,6}.txt). Checks:
  (0) every d in E is lex-positive with |d_i| < t;
  (1) triangle-free on [t]^3: no d1, d2 in E with d1 + d2 in E and the triple realizable in [t]^3
      (realizable iff per coordinate max(|a|,|b|,|a+b|) <= t-1, checked here directly by points);
  (2) every binary subgrid of [t]^3 contains a pair whose gap is in E.
(2) is done hierarchically: a binary subgrid = (a0 < a1, X, Y) with X, Y binary subgrids of [t]^2;
it is hit iff X or Y is internally hit (gap (0, *, *) in E) or some x in X, y in Y has (a1-a0, y-x) in E.
Only g = a1 - a0 matters (gap-determined), and the number of (a0, a1) with a given g is t - g.
"""
import itertools, sys
import numpy as np

path = sys.argv[1]
t = int(sys.argv[2])
E = set(tuple(map(int, l.split())) for l in open(path) if l.strip())
print(f"file {path.split('/')[-1]}: |E| = {len(E)}")

def lexpos(d):
    for x in d:
        if x:
            return x > 0
    return False

assert all(len(d) == 3 and lexpos(d) and all(abs(x) < t for x in d) for d in E), "bad gap vector"

# (1) triangle-free, by explicit points
pts = list(itertools.product(range(t), repeat=3))
tri = 0
for x in pts:
    for d1 in E:
        y = tuple(x[i] + d1[i] for i in range(3))
        if not all(0 <= c < t for c in y):
            continue
        for d2 in E:
            z = tuple(y[i] + d2[i] for i in range(3))
            if all(0 <= c < t for c in z) and tuple(d1[i] + d2[i] for i in range(3)) in E:
                tri += 1
print("triangles:", tri)

# (2) subgrids
P2 = [(b, c) for b in range(t) for c in range(t)]
pidx = {p: i for i, p in enumerate(P2)}
pairs = list(itertools.combinations(range(t), 2))
subs2 = []
for b0, b1 in pairs:
    for c in pairs:
        for c2 in pairs:
            subs2.append(((b0, c[0]), (b0, c[1]), (b1, c2[0]), (b1, c2[1])))
def internally_hit(X):
    for p, q in itertools.combinations(X, 2):   # p <lex q by construction order
        if (0, q[0] - p[0], q[1] - p[1]) in E:
            return True
    return False
ih = np.array([internally_hit(X) for X in subs2])
masks = np.array([sum(1 << pidx[p] for p in X) for X in subs2], dtype=object)
# use python ints (t^2 <= 36 bits) -> int64 is fine
masks = np.array([int(m) for m in masks], dtype=np.int64)
free = np.where(~ih)[0]
print(f"2-level subgrids: {len(subs2)}, internally independent: {len(free)}")
missed_total = 0
total = 0
for g in range(1, t):
    nrow = t - g
    for xi in free:
        X = subs2[xi]
        A = 0
        for y in P2:
            for x in X:
                if (g, y[0] - x[0], y[1] - x[1]) in E:
                    A |= 1 << pidx[y]
                    break
        missed = int(np.count_nonzero((masks[free] & A) == 0))
        missed_total += missed * nrow
    total += nrow * len(subs2) ** 2
print(f"binary subgrids of [t]^3: {total} (= C(t,2)^7 = {len(pairs)**7}); independent (missed): {missed_total}")
print("VERDICT:", "WITNESS OK" if tri == 0 and missed_total == 0 else "WITNESS FAILS")
