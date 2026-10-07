"""Brute force (no SAT) the bare (sign, v_2) gap-determined algebra F2 on [t]^2, t = 3, 4: all 2^12 class sets."""
import itertools
import auditB_n2_games as G
def v2(x):
    x = abs(x); v = 0
    while x % 2 == 0: x //= 2; v += 1
    return v
feat = lambda d: ('=',) if d == 0 else ((1 if d > 0 else -1), v2(d))
def key(x, y):
    if y < x: x, y = y, x
    return (feat(y[0] - x[0]), feat(y[1] - x[1]))
for t in (3, 4):
    P = G.points(t)
    classes = sorted({key(x, y) for x, y in itertools.combinations(P, 2)}, key=repr)
    ci = {c: i for i, c in enumerate(classes)}
    tri = {frozenset((ci[key(x, y)], ci[key(x, z)], ci[key(y, z)])) for x, y, z in itertools.combinations(P, 3)}
    hit = {frozenset(ci[key(a, b)] for a, b in itertools.combinations(g, 2)) for g in G.subgrids(t)}
    trim = [sum(1 << i for i in s) for s in tri]; hitm = [sum(1 << i for i in s) for s in hit]
    win = [m for m in range(1 << len(classes)) if all((m & c) != c for c in trim) and all(m & h for h in hitm)]
    print(f"F2 n=2 t={t}: {len(classes)} classes, winning class sets: {len(win)}")
