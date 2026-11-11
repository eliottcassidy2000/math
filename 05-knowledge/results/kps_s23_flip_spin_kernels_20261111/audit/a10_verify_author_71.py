"""a10_verify_author_71.py (AUDIT) - independent direct recount of the certificate data/pg42_71_odd_planes.txt that
appeared in the author's working tree AFTER the audited commit bdb2f8edd8 (file dated 10:48; not part of the
commit under audit).  Convention of that file: pairs u < v in lexicographic order, bit k = 1 iff v -> u.
No linear algebra; cyclic triangles are counted directly."""
import itertools, os
HERE = os.path.dirname(os.path.abspath(__file__))
x = None
for l in open(os.path.join(HERE, "..", "data", "pg42_71_odd_planes.txt")):
    if l.startswith("x ="):
        x = int(l.split("=")[1])
PTS = list(range(1, 32))
arc = {}
k = 0
for u in PTS:
    for v in PTS:
        if u < v:
            b = (x >> k) & 1
            arc[(u, v)] = not b
            arc[(v, u)] = bool(b)
            k += 1
assert x >> k == 0
def cyc(u, v, w):
    return (arc[(u, v)] and arc[(v, w)] and arc[(w, u)]) or (arc[(v, u)] and arc[(w, v)] and arc[(u, w)])
planes = sorted({tuple(sorted({a, b, c, a ^ b, a ^ c, b ^ c, a ^ b ^ c})) for a in PTS for b in PTS for c in PTS if a < b < c and c != a ^ b})
lines = sorted({tuple(sorted((a, b, a ^ b))) for a in PTS for b in PTS if a < b})
cl = sum(1 for l in lines if cyc(*l))
dist = {}
ev = 0
oddset = []
for pl in planes:
    bad = 0
    for q in itertools.combinations(pl, 4):
        if q[0] ^ q[1] ^ q[2] ^ q[3]:
            continue
        kk = sum(1 for tr in itertools.combinations(q, 3) if cyc(*tr))
        bad += (kk != 1)
    dist[bad] = dist.get(bad, 0) + 1
    ev += bad
    if bad == 0:
        oddset.append(pl)
L = {1, 2, 3}
by = {}
for pl in oddset:
    by[len(set(pl) & L)] = by.get(len(set(pl) & L), 0) + 1
res = "author's certificate (post-commit): cyclic lines %d/155; planes by number of even 4-sets %s; odd planes %d; even zero-sum 4-sets %d; odd planes by number of points on the line {1,2,3}: %s" % (cl, sorted(dist.items()), len(oddset), ev, sorted(by.items()))
print(res)
open(os.path.join(HERE, "a10_verify_author_71.out"), "w").write(res + "\n")
