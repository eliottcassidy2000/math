"""pg32.py (AUDIT, independent) - model of sign tournaments on PG(3,2) used by the a7_* programs.

Points: 1..15 (nonzero vectors of Z_2^4).  Arc variable t[(u, v)], u < v: 1 iff u -> v (bit index = position of the
pair in lexicographic order).  Zero-sum 4-set Q = (a < b < c < d), a^b^c^d = 0; pairings pi1 = ab|cd, pi2 = ac|bd.
f_pi(t) = number of arcs from the first side to the second, mod 2 (an affine function of t).
Violation bits v1 = f_pi1 + 1, v2 = f_pi2 + 1;  Q is ODD (exactly one cyclic triangle) iff v1 = v2 = 0.
The 210 violation bits are ordered plane by plane: bit 14 p + 2 q + j for the q-th 4-set of plane p, pairing j.
Nothing is imported from the author's scripts."""
import itertools

PTS = list(range(1, 16))
PAIRS = [(u, v) for u in PTS for v in PTS if u < v]
PIDX = {p: i for i, p in enumerate(PAIRS)}
NV = len(PAIRS)                                   # 105


def dot(l, v):
    return bin(l & v).count("1") & 1


PLANES = [tuple(v for v in PTS if dot(l, v) == 0) for l in PTS]          # plane l = kernel of the functional l
assert all(len(p) == 7 for p in PLANES)
QUADS = sorted({tuple(sorted((a, b, c, a ^ b ^ c))) for a in PTS for b in PTS for c in PTS
                if a < b < c and (a ^ b ^ c) not in (0, a, b, c)})
assert len(QUADS) == 105
LINES = sorted({tuple(sorted((a, b, a ^ b))) for a in PTS for b in PTS if a < b})
assert len(LINES) == 35
PLANE_QUADS = [[q for q in QUADS if set(q) <= set(pl)] for pl in PLANES]
assert all(len(x) == 7 for x in PLANE_QUADS)
assert sorted(q for x in PLANE_QUADS for q in x) == QUADS                  # each 4-set in exactly one plane


def arc_term(u, v):
    """affine form of [u -> v] in the variables t: (mask, const)."""
    if u < v:
        return (1 << PIDX[(u, v)], 0)
    return (1 << PIDX[(v, u)], 1)


def pairing_form(side1, side2):
    m, c = 0, 0
    for u in side1:
        for v in side2:
            mm, cc = arc_term(u, v)
            m ^= mm
            c ^= cc
    return m, c


def quad_rows(q):
    """two rows (mask, const) with violation bit = parity(t & mask) ^ const  (= f_pi + 1)."""
    a, b, c, d = q
    m1, c1 = pairing_form((a, b), (c, d))
    m2, c2 = pairing_form((a, c), (b, d))
    return [(m1, c1 ^ 1), (m2, c2 ^ 1)]


ROWS = []            # 210 rows in the order plane, quad, pairing
for p in range(15):
    for q in PLANE_QUADS[p]:
        ROWS += quad_rows(q)
assert len(ROWS) == 210


def parity(x):
    return bin(x).count("1") & 1


def violation_vector(t):
    w = 0
    for k, (m, c) in enumerate(ROWS):
        if parity(t & m) ^ c:
            w |= 1 << k
    return w


def arcs_from_bits(t):
    """out-sets on vertices 0..14 (vertex i = point i + 1)."""
    out = [set() for _ in range(15)]
    for (u, v), i in PIDX.items():
        if t >> i & 1:
            out[u - 1].add(v - 1)
        else:
            out[v - 1].add(u - 1)
    return out


def c3_in(out, S):
    """number of cyclic triangles inside the point set S (points 1..15) of the tournament 'out' (0-based)."""
    c = 0
    for u, v, w in itertools.combinations(S, 3):
        a, b, d = u - 1, v - 1, w - 1
        if (b in out[a] and d in out[b] and a in out[d]) or (a in out[b] and b in out[d] and d in out[a]):
            c += 1
    return c


def even_quads_direct(out):
    """list of even zero-sum 4-sets by DIRECT count of cyclic triangles."""
    return [q for q in QUADS if c3_in(out, q) % 2 == 0]


class Elim:
    """incremental Gaussian elimination over F_2 for rows (mask, const): parity(t & mask) = const."""
    __slots__ = ("piv",)

    def __init__(self, piv=None):
        self.piv = dict(piv) if piv else {}

    def copy(self):
        return Elim(self.piv)

    def add(self, m, c):
        piv = self.piv
        while m:
            hb = m.bit_length() - 1
            pr = piv.get(hb)
            if pr is None:
                piv[hb] = (m, c)
                return True
            m ^= pr[0]
            c ^= pr[1]
        return c == 0

    def add_rows(self, rows):
        ok = True
        for m, c in rows:
            if not self.add(m, c):
                ok = False
        return ok

    def rank(self):
        return len(self.piv)

    def solve(self):
        """one solution (free variables 0) and a basis of the homogeneous solutions."""
        piv = self.piv
        # reduce to reduced echelon form
        order = sorted(piv)                       # increasing pivot bit
        red = {}
        for hb in order:
            m, c = piv[hb]
            # eliminate lower pivots from m
            for lb in order:
                if lb >= hb:
                    break
                if m >> lb & 1:
                    m ^= red[lb][0]
                    c ^= red[lb][1]
            red[hb] = (m, c)
        t0 = 0
        for hb, (m, c) in red.items():
            if c:
                t0 |= 1 << hb
        free = [v for v in range(NV) if v not in red]
        basis = []
        for fv in free:
            x = 1 << fv
            for hb, (m, c) in red.items():
                if m >> fv & 1:
                    x |= 1 << hb
            basis.append(x)
        return t0, basis


def plane_rows(p):
    """the 14 equations 'violation = 0' of plane p."""
    return ROWS[14 * p:14 * p + 14]


def lin_image(x):
    """linear part of the violation map applied to x."""
    w = 0
    for k, (m, c) in enumerate(ROWS):
        if parity(x & m):
            w |= 1 << k
    return w


def independent(vectors):
    """an independent subfamily spanning the same space (list of ints)."""
    el = {}
    res = []
    for v in vectors:
        t = v
        while t:
            hb = t.bit_length() - 1
            if hb in el:
                t ^= el[hb]
            else:
                el[hb] = t
                res.append(v)
                break
    return res


def defect_of_plane(w, p):
    f = (w >> (14 * p)) & 0x3FFF
    return bin((f | (f >> 1)) & 0x1555).count("1")


def gl42():
    """all invertible 4x4 matrices over F_2 as tuples of images of the basis vectors 1, 2, 4, 8."""
    res = []
    for imgs in itertools.product(range(1, 16), repeat=4):
        span = {0}
        ok = True
        for im in imgs:
            if im in span:
                ok = False
                break
            span |= {s ^ im for s in span}
        if ok:
            res.append(imgs)
    return res


def apply_lin(imgs, v):
    o = 0
    for b in range(4):
        if v >> b & 1:
            o ^= imgs[b]
    return o
