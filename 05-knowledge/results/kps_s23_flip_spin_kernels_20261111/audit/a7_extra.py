"""a7_extra.py (AUDIT) - further checks of Theorem 3.4 and of the surrounding statements, own code.

 a. type B: per-plane defects of the other seven planes (all 120 sets);
 b. type A: the 64 patterns with 21 even 4-sets: per-plane defects; an explicit tournament, recounted directly;
 c. the sedenion sign tournament (own Cayley-Dickson doubling): odd planes, type, even 4-sets;
 d. the minimiser printed in runs/cd_minimum.txt (105-bit word): direct recount;
 e. is 21 attained only with 8 odd planes?  (from a7_enum_five_reps.out);
 f. a second, different proof of 'k2 <= 3 next to five odd planes' is in a7_enum_five_reps.out (full enumeration)."""
import os, sys, itertools, collections, random
from pg32 import *
from aud_common import ham_c

HERE = os.path.dirname(os.path.abspath(__file__))
out = []


def P(s):
    print(s)
    out.append(s)
    sys.stdout.flush()


def feasible_mask(mask):
    e = Elim()
    ok = True
    for i in range(15):
        if mask >> i & 1:
            ok &= e.add_rows(plane_rows(i))
    return ok, e


# all 8-subsets that are feasible (recomputed here)
eights = []
for comb in itertools.combinations(range(15), 8):
    m = sum(1 << i for i in comb)
    ok, e = feasible_mask(m)
    if ok:
        eights.append((m, e))
P("feasible 8-sets recomputed by direct enumeration of all C(15,8) = 6435 subsets: %d" % len(eights))
typeA, typeB = [], []
for m, e in eights:
    sel = [PLANES[i] for i in range(15) if m >> i & 1]
    missing = [p for p in PTS if all(p not in s for s in sel)]
    (typeA if missing else typeB).append((m, e))
# a. type B
pat = collections.Counter()
for m, e in typeB:
    t0, basis = e.solve()
    w = violation_vector(t0)
    # the pattern does not depend on the solution: check on random solutions
    for trial in range(3):
        t = t0
        for x in basis:
            if random.random() < 0.5:
                t ^= x
        assert violation_vector(t) == w
    pat[tuple(sorted(defect_of_plane(w, p) for p in range(15)))] += 1
P("a. type B (120 sets): per-plane defects (sorted) -> number of sets: %s" % sorted(pat.items()))
# in type B, which 4-set stays odd in each of the seven other planes?
m, e = typeB[0]
t0, basis = e.solve()
w = violation_vector(t0)
sel = [i for i in range(15) if m >> i & 1]
rest = [i for i in range(15) if not m >> i & 1]
P("   example: odd planes %s; other planes %s; per other plane the number of odd 4-sets: %s"
  % (sel, rest, [7 - defect_of_plane(w, p) for p in rest]))
# b. type A, total 21
m, e = typeA[0]
t0, basis = e.solve()
w0 = violation_vector(t0)
imgs = []
src = []
el = {}
for x in basis:
    v = lin_image(x)
    t = v
    while t:
        hb = t.bit_length() - 1
        if hb in el:
            t ^= el[hb]
        else:
            el[hb] = t
            imgs.append(v)
            src.append(x)
            break
assert len(imgs) == 19
cnt21 = 0
shapes = collections.Counter()
witness = None
w = w0
tcur = t0
for g in range(1 << 19):
    if g:
        bit = (g & -g).bit_length() - 1
        w ^= imgs[bit]
        tcur ^= src[bit]
    tot = sum(defect_of_plane(w, p) for p in range(15))
    if tot == 21:
        cnt21 += 1
        shapes[tuple(sorted(defect_of_plane(w, p) for p in range(15)))] += 1
        if witness is None:
            witness = tcur
P("b. type A (one set, 2^19 patterns): patterns with 21 even 4-sets: %d; per-plane defects: %s" % (cnt21, sorted(shapes.items())))
T = arcs_from_bits(witness)
ev = even_quads_direct(T)
cyc = sum(1 for l in LINES if c3_in(T, l) == 1)
P("   explicit tournament with 21 even 4-sets (direct recount of cyclic triangles: %d even 4-sets); cyclic lines %d/35; scores %s; H = %d"
  % (len(ev), cyc, sorted(len(s) for s in T), ham_c([T])[0]))
P("   arc word of this audit witness (bit i = pair i in lexicographic order of (u < v), 1 iff u -> v): %d" % witness)

# c. sedenions by own doubling
def conj(x):
    if len(x) == 1:
        return x[:]
    h = len(x) // 2
    return conj(x[:h]) + [-t for t in x[h:]]


def mul(x, y):
    if len(x) == 1:
        return [x[0] * y[0]]
    h = len(x) // 2
    a, b, c, d = x[:h], x[h:], y[:h], y[h:]
    first = [p - q for p, q in zip(mul(a, c), mul(conj(d), b))]
    second = [p + q for p, q in zip(mul(d, a), mul(b, conj(c)))]
    return first + second


g = [[0] * 16 for _ in range(16)]
for i in range(16):
    for j in range(16):
        x = [0] * 16
        y = [0] * 16
        x[i] = 1
        y[j] = 1
        pr = mul(x, y)
        k = [t for t in range(16) if pr[t]]
        assert k == [i ^ j]
        g[i][j] = pr[i ^ j]
t = 0
for (u, v), i in PIDX.items():
    if g[u][v] == 1:
        t |= 1 << i
w = violation_vector(t)
T = arcs_from_bits(t)
defs = [defect_of_plane(w, p) for p in range(15)]
oddp = [p for p in range(15) if defs[p] == 0]
sel = [PLANES[p] for p in oddp]
missing = [p for p in PTS if all(p not in s for s in sel)]
desc = "A" if missing else "?"
if not missing:
    for p0 in oddp:
        others = [PLANES[p] for p in oddp if p != p0]
        common = [q for q in PTS if all(q in s for s in others)]
        if len(common) == 1 and common[0] not in PLANES[p0]:
            desc = "B: plane with normal %d (points %s) and the seven planes through the point %d" % (p0 + 1, PLANES[p0], common[0])
P("c. sedenion sign tournament: per-plane defects %s; odd planes %d; type %s; even 4-sets (direct) %d; cyclic lines %d; H = %d"
  % (sorted(defs), len(oddp), desc, len(even_quads_direct(T)), sum(1 for l in LINES if c3_in(T, l) == 1), ham_c([T])[0]))

# d. the packet's minimiser
x = 2211813282662872962504903186636          # bit = 1 iff v -> u for the pair u < v (lexicographic order)
t = 0
for i in range(NV):
    if not (x >> i & 1):
        t |= 1 << i
T = arcs_from_bits(t)
ev = even_quads_direct(T)
w = violation_vector(t)
P("d. the minimiser of runs/cd_minimum.txt: even 4-sets by direct recount %d; per-plane defects %s; cyclic lines %d/35; scores %s; H = %d"
  % (len(ev), sorted(defect_of_plane(w, p) for p in range(15)), sum(1 for l in LINES if c3_in(T, l) == 1), sorted(len(s) for s in T), ham_c([T])[0]))
open(os.path.join(HERE, "a7_extra.out"), "w").write("\n".join(out) + "\n")
