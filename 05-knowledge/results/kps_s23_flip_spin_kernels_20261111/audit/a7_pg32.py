"""a7_pg32.py (AUDIT, independent) - item 7: Theorem 3.4 on PG(3,2), own linear algebra over F_2.

 0. validation of the linear model against DIRECT counts of cyclic triangles on random tournaments;
 1. the 210 parity equations: rank, consistency; kernel of the parity map = span of the 30 moves;
 2. every subset of the 15 planes: can all its planes be odd at once?  counts by size; the 8-sets and their types;
 3. type A and type B: violation patterns of the other seven planes (problem file for a7_enum.c);
 4. 'all lines cyclic': consistency with each type, dimensions; a plane with all lines cyclic has defect 0 or 6;
 5. all feasible 5-sets of planes with their orbits under the FULL group GL(4,2) (20160 matrices), and one
    problem per 5-set for a7_enum.c (the enumeration of ALL violation patterns next to five odd planes)."""
import random, sys, os, itertools, collections
from pg32 import *

HERE = os.path.dirname(os.path.abspath(__file__))
out = []


def P(s):
    print(s)
    out.append(s)
    sys.stdout.flush()


rng = random.Random(7)
# ---- 0. model validation
bad = 0
dist3 = collections.Counter()
for trial in range(300):
    t = rng.getrandbits(NV)
    T = arcs_from_bits(t)
    w = violation_vector(t)
    k = 0
    for p in range(15):
        for q in PLANE_QUADS[p]:
            v1 = w >> k & 1
            v2 = w >> (k + 1) & 1
            k += 2
            c = c3_in(T, q)
            dist3[c] += 1
            if (c % 2 == 1) != (v1 == 0 and v2 == 0):
                bad += 1
P("0. linear model vs direct counts of cyclic triangles: 300 random tournaments x 105 zero-sum 4-sets, mismatches %d; cyclic triangles per 4-set: %s"
  % (bad, sorted(dist3.items())))

# ---- 1. rank and consistency
E = Elim()
ok = E.add_rows(ROWS)
E2 = Elim()
E2.add_rows([(m, 0) for m, c in ROWS])
P("1. 'all 105 zero-sum 4-sets odd': %d equations in %d unknowns, rank %d, consistent: %s  (kernel of the parity map: dimension %d)"
  % (len(ROWS), NV, E2.rank(), ok, NV - E2.rank()))
# two independent equations per 4-set?
ind2 = all(Elim().add_rows([(ROWS[2 * i][0], 0), (ROWS[2 * i + 1][0], 0)]) and len(independent([ROWS[2 * i][0], ROWS[2 * i + 1][0]])) == 2 for i in range(105))
P("   the two equations of each 4-set are independent: %s" % ind2)
# moves
moves = []
for u in PTS:
    g = 0
    for v in PTS:
        if v != u:
            g |= 1 << PIDX[(min(u, v), max(u, v))]
    moves.append(g)
for f in PTS:
    g = 0
    for v in PTS:
        w_ = v ^ f
        if w_ > v:
            g |= 1 << PIDX[(v, w_)]
    moves.append(g)
inker = all(lin_image(g) == 0 for g in moves)
P("   15 switchings + 15 line coboundaries: all in the kernel of the linear part: %s; rank of the 30 moves: %d (kernel dimension %d)"
  % (inker, len(independent(moves)), NV - E2.rank()))
# rescaling e_u -> -e_u = switching at u + coboundary at u reverses all three arcs of every line through u
u = 5
resc = moves[u - 1] ^ moves[15 + u - 1]
through = [l for l in LINES if u in l]
arcs_rev = {PAIRS[i] for i in range(NV) if resc >> i & 1}
expect = set()
for l in through:
    for a, b in itertools.combinations(l, 2):
        expect.add((min(a, b), max(a, b)))
P("   rescaling e_u -> -e_u (switching + coboundary at u) reverses exactly the three arcs of each of the 7 lines through u: %s (%d arcs)"
  % (arcs_rev == expect, len(arcs_rev)))

# ---- 2. subsets of planes
feasible = {}          # mask -> rank


def dfs(start, el, mask):
    feasible[mask] = el.rank()
    for i in range(start, 15):
        e2 = el.copy()
        if e2.add_rows(plane_rows(i)):
            dfs(i + 1, e2, mask | 1 << i)


dfs(0, Elim(), 0)
bysize = collections.Counter(bin(m).count("1") for m in feasible)
P("2. subsets of the 15 planes that can be odd at once, by size: %s ; largest %d" % (sorted(bysize.items()), max(bysize)))
# monotonicity check on a sample: every subset of a feasible set is feasible; an infeasible set stays infeasible
eight = sorted(m for m in feasible if bin(m).count("1") == 8)
typeA, typeB, other = [], [], []
for m in eight:
    sel = [PLANES[i] for i in range(15) if m >> i & 1]
    missing = [p for p in PTS if all(p not in s for s in sel)]
    if len(missing) == 1:
        # the planes off the point 'missing'
        assert sorted(sel) == sorted(pl for pl in PLANES if missing[0] not in pl)
        typeA.append(m)
        continue
    okB = False
    for i0 in range(15):
        if not m >> i0 & 1:
            continue
        rest = [PLANES[i] for i in range(15) if m >> i & 1 and i != i0]
        common = [p for p in PTS if all(p in s for s in rest)]
        if len(common) == 1 and common[0] not in PLANES[i0]:
            assert sorted(rest) == sorted(pl for pl in PLANES if common[0] in pl)
            okB = True
    (typeB if okB else other).append(m)
P("   sets of 8: %d = %d of type A (all planes off a point) + %d of type B (a plane and the 7 planes through a point off it) + %d other"
  % (len(eight), len(typeA), len(typeB), len(other)))
P("   ranks (solution dimension = 105 - rank): type A %s, type B %s" % (sorted({feasible[m] for m in typeA}), sorted({feasible[m] for m in typeB})))
# brute-force cross-check of the DFS on random subsets (fresh elimination each time)
chk = 0
for trial in range(400):
    m = rng.getrandbits(15)
    e = Elim()
    okm = True
    for i in range(15):
        if m >> i & 1:
            okm &= e.add_rows(plane_rows(i))
    if okm != (m in feasible):
        chk += 1
P("   cross-check of the subset search on 400 random subsets with a fresh elimination: disagreements %d" % chk)


# ---- 3. problems for the enumerator
def problem(mask):
    e = Elim()
    for i in range(15):
        if mask >> i & 1:
            assert e.add_rows(plane_rows(i))
    t0, basis = e.solve()
    w0 = violation_vector(t0)
    # sanity: the selected planes have zero violation
    for i in range(15):
        if mask >> i & 1:
            assert (w0 >> (14 * i)) & 0x3FFF == 0
    imgs = independent([lin_image(x) for x in basis])
    imgs = [v for v in imgs if v]
    return t0, basis, w0, imgs


probs = []
for m in typeA:
    t0, basis, w0, imgs = problem(m)
    probs.append(("A_%d" % m, w0, imgs))
for m in typeB:
    t0, basis, w0, imgs = problem(m)
    probs.append(("B_%d" % m, w0, imgs))
P("3. type A: %d sets, dimension of the pattern space %s; type B: %d sets, dimension %s"
  % (len(typeA), sorted({len(p[2]) for p in probs if p[0][0] == "A"}), len(typeB), sorted({len(p[2]) for p in probs if p[0][0] == "B"})))
with open(os.path.join(HERE, "a7_problems_AB.txt"), "w") as fh:
    for label, w0, imgs in probs:
        fh.write("%s %d %x %s\n" % (label, len(imgs), w0, " ".join("%x" % v for v in imgs)))

# ---- 4. all lines cyclic
line_rows = []
for (a, b, c) in LINES:
    # cyclic iff [a->b] = [b->c] = [c->a]
    m1, c1 = arc_term(a, b)
    m2, c2 = arc_term(b, c)
    m3, c3 = arc_term(c, a)
    line_rows.append((m1 ^ m2, c1 ^ c2))
    line_rows.append((m2 ^ m3, c2 ^ c3))
for name, sets in (("A", typeA), ("B", typeB)):
    dims = collections.Counter()
    for m in sets:
        e = Elim()
        for i in range(15):
            if m >> i & 1:
                assert e.add_rows(plane_rows(i))
        okl = e.add_rows(line_rows)
        dims[(okl, NV - e.rank() if okl else None)] += 1
    P("4. type %s + all 35 lines cyclic: (consistent, solution dimension) -> number of sets: %s" % (name, sorted(dims.items(), key=str)))
el = Elim()
okl = el.add_rows(line_rows)
P("   all 35 lines cyclic alone: consistent %s, solution dimension %d (= 2^35 orientations: %s)" % (okl, NV - el.rank(), NV - el.rank() == 35))
# with all lines cyclic every plane has defect 0 or 6: check on random all-cyclic tournaments + number of odd planes
t0, basis = el.solve()
hist = collections.Counter()
maxodd = 0
for trial in range(20000):
    t = t0
    for x in basis:
        if rng.random() < 0.5:
            t ^= x
    w = violation_vector(t)
    ds = tuple(sorted(defect_of_plane(w, p) for p in range(15)))
    assert set(ds) <= {0, 6}
    hist[ds.count(0)] += 1
P("   20000 random tournaments with all lines cyclic: every plane has defect 0 or 6; number of odd planes -> count: %s" % sorted(hist.items()))

# ---- 5. feasible 5-sets, orbits under GL(4,2)
G = gl42()
P("5. |GL(4,2)| enumerated: %d" % len(G))
plane_index = {tuple(sorted(pl)): i for i, pl in enumerate(PLANES)}
perms = set()
for imgs in G:
    perm = tuple(plane_index[tuple(sorted(apply_lin(imgs, v) for v in pl))] for pl in PLANES)
    perms.add(perm)
perms = sorted(perms)
P("   distinct permutations of the 15 planes induced: %d" % len(perms))
five = sorted(m for m in feasible if bin(m).count("1") == 5)
fset = set(five)
orbit_of = {}
orbits = []
for m in five:
    if m in orbit_of:
        continue
    idxs = [i for i in range(15) if m >> i & 1]
    orb = set()
    for perm in perms:
        mm = 0
        for i in idxs:
            mm |= 1 << perm[i]
        orb.add(mm)
    assert orb <= fset
    for mm in orb:
        orbit_of[mm] = len(orbits)
    orbits.append((m, len(orb), feasible[m]))
P("   feasible sets of 5 planes: %d, in %d orbits under GL(4,2): (representative mask, size, rank) %s" % (len(five), len(orbits), orbits))
# same for all sizes (orbit counts), as extra information
for size in (6, 7, 8):
    sets = sorted(m for m in feasible if bin(m).count("1") == size)
    seen = set()
    norb = []
    for m in sets:
        if m in seen:
            continue
        idxs = [i for i in range(15) if m >> i & 1]
        orb = set()
        for perm in perms:
            mm = 0
            for i in idxs:
                mm |= 1 << perm[i]
            orb.add(mm)
        seen |= orb
        norb.append((len(orb), feasible[m]))
    P("   feasible sets of %d planes: %d, orbits (size, rank): %s" % (size, len(sets), norb))
with open(os.path.join(HERE, "a7_problems_five.txt"), "w") as fh:
    dimc = collections.Counter()
    for m in five:
        t0, basis, w0, imgs = problem(m)
        dimc[(orbit_of[m], len(imgs))] += 1
        fh.write("F%d_%d %d %x %s\n" % (orbit_of[m], m, len(imgs), w0, " ".join("%x" % v for v in imgs)))
P("   problems written for all %d feasible 5-sets; (orbit, dimension of the pattern space) -> count: %s" % (len(five), sorted(dimc.items())))
with open(os.path.join(HERE, "a7_problems_five_reps.txt"), "w") as fh:
    for (m, size, rk) in orbits:
        t0, basis, w0, imgs = problem(m)
        fh.write("R%d_%d %d %x %s\n" % (orbit_of[m], m, len(imgs), w0, " ".join("%x" % v for v in imgs)))
# the empty set: the whole pattern space (dimension 80) cannot be enumerated; record w0 for the annealer
open(os.path.join(HERE, "a7_pg32.out"), "w").write("\n".join(out) + "\n")
