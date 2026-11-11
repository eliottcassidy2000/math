"""cd_min.py - the least number of even zero-sum 4-sets of a tournament on PG(3,2), exactly (S23).

Vocabulary: a zero-sum 4-set Q = {a,b,c,d} (a^b^c^d = 0, points of PG(3,2)) is EVEN if it contains an even number
of cyclic triangles (it is transitive or strong), ODD otherwise.  Q lies in one plane; a plane is ODD if its seven
4-sets (line complements) are odd.  The defect of a plane = number of its even 4-sets.

Facts used (all recomputed here):
 F1. Fano plane: the defect is never 1 (distribution over all 2^21 tournaments, via the image of the affine map).
 F2. PG(3,2): at most 8 planes are odd at once; the 135 sets of 8 are 15 of type A (planes off a point) and 120 of
     type B (a plane and the seven planes through a point off it).
 F3. Type B: the remaining seven planes have defect 6 each, whatever the solution (42 even 4-sets).
 F4. Type A: all violation vectors of the remaining seven planes enumerated; least total defect.
 F5. At most 7 odd planes: 8 or more planes have defect >= 2; total <= 20 would need 4 planes of defect exactly 2;
     depth-first search over sets of planes of defect exactly 2.
Also: all-cyclic Fano orientations (128) by defect; a tournament with the least number, described."""
import itertools, sys
from cd_more import F2System, planes, quads, c3_of
from cd_signs import ham_paths


def popc(x):
    return bin(x).count("1")


class Elim:
    """incremental Gaussian elimination over F_2 on (mask, const) rows; copyable."""
    def __init__(self, piv=None):
        self.piv = dict(piv) if piv else {}

    def copy(self):
        return Elim(self.piv)

    def add(self, mask, const):
        while mask:
            hb = mask.bit_length() - 1
            if hb in self.piv:
                pm, pc = self.piv[hb]
                mask ^= pm
                const ^= pc
            else:
                self.piv[hb] = (mask, const)
                return True
        return const == 0

    def add_all(self, rows):
        ok = True
        for m, c in rows:
            ok &= self.add(m, c)
        return ok


def plane_data(S, P, Q4):
    """per plane: list of 7 quads, 14 rows (mask, const) in the order (q0 pairing 1, q0 pairing 2, q1 ...)."""
    out = []
    for pl in P:
        s = set(pl)
        qs = [q for q in Q4 if s >= set(q)]
        rows = []
        for q in qs:
            rows += S.quad_rows(q)
        out.append((qs, rows))
    return out


def realizable(rows, nv):
    """all violation vectors (14 bits) of one plane: viol_j = parity(x & mask_j) ^ const_j."""
    base = 0
    for j, (m, c) in enumerate(rows):
        base |= c << j
    cols = []
    used = 0
    for m, c in rows:
        used |= m
    v = 0
    while used >> v:
        if used >> v & 1:
            col = 0
            for j, (m, c) in enumerate(rows):
                if m >> v & 1:
                    col |= 1 << j
            cols.append(col)
        v += 1
    span = {0}
    for col in cols:
        if col not in span:
            span |= {s ^ col for s in span}
    return sorted(base ^ s for s in span)


def defect(w):
    return sum(1 for i in range(7) if (w >> (2 * i)) & 3)


if __name__ == "__main__":
    S = F2System(4)
    nv = len(S.var)
    P = planes(4)
    Q4 = quads(4)
    PD = plane_data(S, P, Q4)
    R = [realizable(rows, nv) for (qs, rows) in PD]
    dist = {}
    for w in R[0]:
        dist[defect(w)] = dist.get(defect(w), 0) + 1
    print("F1. one plane: realizable violation vectors %d, by defect: %s" % (len(R[0]), sorted(dist.items())))
    assert all(len(r) == 2048 for r in R)

    def rows_with(i, w):
        return [(m, c ^ (w >> j & 1)) for j, (m, c) in enumerate(PD[i][1])]
    # F2
    odd_sets = []
    for mask in range(1 << 15):
        if popc(mask) < 8:
            continue
        e = Elim()
        ok = True
        for i in range(15):
            if mask >> i & 1:
                ok &= e.add_all(rows_with(i, 0))
        if ok:
            odd_sets.append(mask)
    print("F2. sets of at least 8 planes that can be odd at once: %d, sizes %s" % (len(odd_sets), sorted(set(popc(m) for m in odd_sets))))
    V = range(1, 16)
    typeA = [m for m in odd_sets if any(all(p not in P[i] for i in range(15) if m >> i & 1) for p in V)]
    typeB = [m for m in odd_sets if m not in typeA]
    print("    type A (the planes off a point): %d ; type B: %d" % (len(typeA), len(typeB)))
    # F3 / F4: enumerate the violation vectors of the other planes over the whole solution space
    def enumerate_rest(mask):
        e = Elim()
        for i in range(15):
            if mask >> i & 1:
                assert e.add_all(rows_with(i, 0))
        # particular solution and homogeneous basis
        piv = e.piv
        free = [v for v in range(nv) if v not in piv]

        def back(x):
            for hb in sorted(piv):
                m, c = piv[hb]
                val = c ^ (popc(x & (m ^ (1 << hb))) & 1)
                x = (x | (1 << hb)) if val else (x & ~(1 << hb))
            return x
        x0 = back(0)
        hom = []
        for fv in free:
            # homogeneous: set free var fv = 1, others 0, constants 0
            x = 1 << fv
            for hb in sorted(piv):
                m, c = piv[hb]
                val = popc(x & (m ^ (1 << hb))) & 1
                x = (x | (1 << hb)) if val else (x & ~(1 << hb))
            hom.append(x)
        rest = [i for i in range(15) if not (mask >> i & 1)]

        def viol(x, homog):
            out = 0
            k = 0
            for i in rest:
                for (m, c) in PD[i][1]:
                    out |= ((popc(x & m) & 1) ^ (0 if homog else c)) << k
                    k += 1
            return out
        w0 = viol(x0, False)
        span_basis = []
        seen = Elim()
        for hx in hom:
            hv = viol(hx, True)
            if hv and seen.add(hv, 0):
                # keep reduced independent vector
                span_basis.append(hv)
        # independent set via elimination on the images
        el = {}
        indep = []
        for hx in hom:
            hv = viol(hx, True)
            t = hv
            while t:
                hb = t.bit_length() - 1
                if hb in el:
                    t ^= el[hb]
                else:
                    el[hb] = t
                    indep.append(hv)
                    break
        nq = 7 * len(rest)
        best = None
        hist = {}
        shapes21 = {}
        w = w0
        total = 1 << len(indep)
        # Gray code enumeration
        lowmask = sum(1 << (2 * i) for i in range(nq))
        for g in range(total):
            if g:
                bit = (g & -g).bit_length() - 1
                w ^= indep[bit]
            d = popc((w | (w >> 1)) & lowmask)
            hist[d] = hist.get(d, 0) + 1
            if d == 21:
                shape = tuple(sorted(popc(((w >> (14 * pi)) | (w >> (14 * pi + 1))) & 0x1555) for pi in range(len(rest))))
                shapes21[shape] = shapes21.get(shape, 0) + 1
            if best is None or d < best[0]:
                best = (d, w, g)
        if shapes21:
            print("    patterns with total defect 21, by the defects of the other planes:", sorted(shapes21.items()))
        return len(indep), hist, best, x0, hom, indep, rest
    kB, histB, bestB, *_ = enumerate_rest(typeB[0])
    print("F3. type B: violation patterns of the other 7 planes: 2^%d, total defect distribution %s" % (kB, sorted(histB.items())))
    kA, histA, bestA, x0, hom, indep, rest = enumerate_rest(typeA[0])
    print("F4. type A: violation patterns of the other 7 planes: 2^%d, total defect distribution %s" % (kA, sorted(histA.items())))
    # F5: sets of planes of defect exactly 2 (and the search for 4 of them)
    w2 = [[w for w in R[i] if defect(w) == 2] for i in range(15)]
    assert all(len(x) == 21 for x in w2)
    MAXD = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    sys.stdout.flush()
    counts = {1: 0, 2: 0, 3: 0, 4: 0}

    def dfs(start, e, depth):
        for i in range(start, 15):
            for w in w2[i]:
                e2 = e.copy()
                if e2.add_all(rows_with(i, w)):
                    counts[depth + 1] += 1
                    if depth + 1 < MAXD:
                        dfs(i + 1, e2, depth + 1)
    dfs(0, Elim(), 0)
    print("F5 (preliminary; the search next to five odd planes is cd_f5_orbits.py). consistent choices of k planes of defect exactly 2, with the pair of lines, searched to depth %d: %s" % (MAXD, [(k, counts[k]) for k in range(1, MAXD + 1)]))
    # all-cyclic Fano orientations by defect
    S3 = F2System(3)
    lines3 = sorted({tuple(sorted((a, b, a ^ b))) for a in range(1, 8) for b in range(1, 8) if a != b})
    Q3 = quads(3)
    cyc_hist = {}
    for bits in range(128):
        adj = {}
        for k, (a, b, c) in enumerate(lines3):
            tri = (a, b, c) if bits >> k & 1 else (a, c, b)
            for u, v in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                adj[(u, v)] = True
        arc = lambda u, v: (u, v) in adj
        d = sum(1 for q in Q3 if c3_of(arc, q) % 2 == 0)
        A = [[1 if (i != j and arc(i + 1, j + 1)) else 0 for j in range(7)] for i in range(7)]
        key = (d, ham_paths(A))
        cyc_hist[key] = cyc_hist.get(key, 0) + 1
    print("all-cyclic orientations of the Fano plane by (defect, H):", sorted(cyc_hist.items()))
    # describe a minimiser of type A
    d, w, g = bestA
    x = x0
    # reconstruct x from the Gray index: need the hom vectors that produced 'indep' in order
    el = {}
    chosen = []
    def violh(xh):
        out = 0
        k = 0
        for i in rest:
            for (m, c) in PD[i][1]:
                out |= (popc(xh & m) & 1) << k
                k += 1
        return out
    for hx in hom:
        t = violh(hx)
        while t:
            hb = t.bit_length() - 1
            if hb in el:
                t ^= el[hb]
            else:
                el[hb] = t
                chosen.append(hx)
                break
    gray = g ^ (g >> 1)
    for b in range(len(chosen)):
        if gray >> b & 1:
            x ^= chosen[b]

    def arc(u, v):
        if u < v:
            return not (x >> S.var[(u, v)] & 1)
        return bool(x >> S.var[(v, u)] & 1)
    ev = [q for q in Q4 if c3_of(arc, q) % 2 == 0]
    lines = sorted({tuple(sorted((a, b, a ^ b))) for a in V for b in V if a != b})
    per_plane = sorted(sum(1 for q in ev if set(q) <= set(pl)) for pl in P)
    scores = sorted(sum(1 for j in V if j != i and arc(i, j)) for i in V)
    A = [[1 if (i != j and arc(i + 1, j + 1)) else 0 for j in range(15)] for i in range(15)]
    print("a minimiser: even 4-sets %d, per plane %s, cyclic lines %d/35, scores %s, H = %d" % (len(ev), per_plane, sum(1 for L in lines if c3_of(arc, L) == 1), scores, ham_paths(A)))
    print("   arcs u->v as a 105-bit word (pairs u<v in lexicographic order, bit 1 iff v -> u):", x)
