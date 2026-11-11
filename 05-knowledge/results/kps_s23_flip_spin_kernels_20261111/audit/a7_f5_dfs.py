"""a7_f5_dfs.py (AUDIT) - an independent re-implementation of the SEARCH DESIGN of cd_f5_orbits.py (own equations,
own elimination, own orbit representatives under the full GL(4,2)): next to five odd planes, all consistent choices
of k further planes of defect exactly 2 with a prescribed defect-2 violation vector each.  The counts per level are
compared with runs/cd_f5_orbits.txt; the conclusion (no 4 such planes) is also obtained by complete enumeration of
all patterns in a7_enum_five_reps.out.

Also: the 21 defect-2 violation vectors of a plane are in bijection with the 21 pairs of lines of the plane."""
import os, sys, itertools, collections
from pg32 import *

HERE = os.path.dirname(os.path.abspath(__file__))
out = []


def P(s):
    print(s)
    out.append(s)
    sys.stdout.flush()


# realizable violation vectors of one plane: image of the affine map restricted to the plane's 14 rows
def realizable(p):
    rows = plane_rows(p)
    base = 0
    for j, (m, c) in enumerate(rows):
        base |= c << j
    used = 0
    for m, c in rows:
        used |= m
    cols = []
    for v in range(NV):
        if used >> v & 1:
            col = 0
            for j, (m, c) in enumerate(rows):
                if m >> v & 1:
                    col |= 1 << j
            cols.append(col)
    span = {0}
    for col in independent(cols):
        span |= {s ^ col for s in span}
    return sorted(base ^ s for s in span)


def dfct(w):
    return bin((w | (w >> 1)) & 0x1555).count("1")


R = [realizable(p) for p in range(15)]
hist = collections.Counter(dfct(w) for w in R[0])
P("one plane: %d realizable violation vectors; by defect %s" % (len(R[0]), sorted(hist.items())))
w2 = [[w for w in R[p] if dfct(w) == 2] for p in range(15)]
ok21 = all(len(x) == 21 for x in w2)
pairs_ok = True
for p in range(15):
    seen = set()
    for w in w2[p]:
        qs = tuple(q for q in range(7) if (w >> (2 * q)) & 3)
        seen.add(qs)
    pairs_ok &= (len(seen) == 21)
P("every plane has exactly 21 defect-2 vectors: %s; they realise each of the 21 pairs of 4-sets (= pairs of lines) exactly once: %s" % (ok21, pairs_ok))


def rows_with(p, w):
    return [(m, c ^ (w >> j & 1)) for j, (m, c) in enumerate(plane_rows(p))]


reps = [(31, 315), (143, 1680), (158, 840)]          # from a7_pg32.out (orbits under the full GL(4,2))
overall = 0
for mask, size in reps:
    e = Elim()
    for i in range(15):
        if mask >> i & 1:
            assert e.add_rows(plane_rows(i))
    others = [i for i in range(15) if not mask >> i & 1]
    levels = collections.Counter()
    best = [0]

    def dfs(start, el, depth):
        for idx in range(start, len(others)):
            i = others[idx]
            for w in w2[i]:
                e2 = el.copy()
                if e2.add_rows(rows_with(i, w)):
                    levels[depth + 1] += 1
                    best[0] = max(best[0], depth + 1)
                    if depth + 1 < 5:
                        dfs(idx + 1, e2, depth + 1)
    dfs(0, e, 0)
    overall = max(overall, best[0])
    P("orbit of %d sets (representative mask %d, rank %d): consistent additions of k planes of defect exactly 2: %s; largest k = %d"
      % (size, mask, e.rank(), sorted(levels.items()), best[0]))
P("largest number of planes of defect exactly 2 next to five odd planes: %d" % overall)
P("author (runs/cd_f5_orbits.txt): 315: [(1, 98), (2, 1405), (3, 1816), (4, 0)]; 1680: [(1, 84), (2, 1038), (3, 3456), (4, 0)]; 840: [(1, 66), (2, 819), (3, 1759), (4, 0)]")
open(os.path.join(HERE, "a7_f5_dfs.out"), "w").write("\n".join(out) + "\n")
