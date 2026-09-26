"""procgen_fencelim_20260926_geom.py -- exact planar engine for Friedman fence configurations.

Session collatz-procgen-20260922, lane "fencelim" (2026-09-26).

A configuration is a list of fences ((x0,y0),(x1,y1)), each of length 1.  Coordinates may be
Fractions (exact), QS numbers a+b*sqrt(D) with rational a,b (exact signs), or mpmath mpf
(high precision; signs decided with a tolerance, used only for the numerically optimised n=6
record).  analyze() checks the rules, builds the plane graph, and returns junction and face
data: e_j, t_j, sectors with the classification theta<pi / =pi / >pi, convex-corner counts
kappa_f, holes, outer-face data, exact field areas and perimeters.
"""
from fractions import Fraction as Fr
from functools import cmp_to_key
import math

try:
    import mpmath
except ImportError:  # pragma: no cover
    mpmath = None


# ------------------------------------------------------------------ Q(sqrt D) with exact sign
class QS:
    __slots__ = ('a', 'b', 'D')

    def __init__(self, a, b=0, D=3):
        self.a = Fr(a)
        self.b = Fr(b)
        self.D = D

    def _c(self, o):
        if isinstance(o, QS):
            assert o.D == self.D
            return o
        return QS(o, 0, self.D)

    def __add__(self, o):
        o = self._c(o)
        return QS(self.a + o.a, self.b + o.b, self.D)
    __radd__ = __add__

    def __sub__(self, o):
        o = self._c(o)
        return QS(self.a - o.a, self.b - o.b, self.D)

    def __rsub__(self, o):
        return self._c(o) - self

    def __neg__(self):
        return QS(-self.a, -self.b, self.D)

    def __mul__(self, o):
        o = self._c(o)
        return QS(self.a * o.a + self.D * self.b * o.b, self.a * o.b + self.b * o.a, self.D)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._c(o)
        n = o.a * o.a - self.D * o.b * o.b
        assert n != 0
        conj = QS(o.a / n, -o.b / n, self.D)
        return self * conj

    def __rtruediv__(self, o):
        return self._c(o) / self

    def sign(self):
        a, b = self.a, self.b
        if b == 0:
            return (a > 0) - (a < 0)
        if a == 0:
            return (b > 0) - (b < 0)
        if a > 0 and b > 0:
            return 1
        if a < 0 and b < 0:
            return -1
        # opposite signs: compare a^2 with D b^2
        c = a * a - self.D * b * b
        s = (c > 0) - (c < 0)
        return s if a > 0 else -s

    def __eq__(self, o):
        return (self - o).sign() == 0

    def __lt__(self, o):
        return (self - o).sign() < 0

    def __le__(self, o):
        return (self - o).sign() <= 0

    def __gt__(self, o):
        return (self - o).sign() > 0

    def __ge__(self, o):
        return (self - o).sign() >= 0

    def __hash__(self):
        return hash((self.a, self.b, self.D))

    def __float__(self):
        return float(self.a) + float(self.b) * math.sqrt(self.D)

    def mp(self):
        return mpmath.mpf(self.a.numerator) / self.a.denominator + \
            mpmath.mpf(self.b.numerator) / self.b.denominator * mpmath.sqrt(self.D)

    def __repr__(self):
        return f'({self.a}+{self.b}*sqrt{self.D})'


MP_TOL = None  # set by callers using mpf (e.g. mpmath.mpf(10)**-30)


def sgn(x):
    if isinstance(x, QS):
        return x.sign()
    if mpmath is not None and isinstance(x, mpmath.mpf):
        if abs(x) <= MP_TOL:
            return 0
        return 1 if x > 0 else -1
    return (x > 0) - (x < 0)


def tomp(x):
    if isinstance(x, QS):
        return x.mp()
    if isinstance(x, Fr):
        return mpmath.mpf(x.numerator) / x.denominator
    return mpmath.mpf(x)


def key(p):
    """hashable key of a point (exact types: itself; mpf: rounded)."""
    if mpmath is not None and isinstance(p[0], mpmath.mpf):
        return (mpmath.nstr(p[0], 22), mpmath.nstr(p[1], 22))
    return (p[0], p[1])


def sub(p, q):
    return (p[0] - q[0], p[1] - q[1])


def cross(u, v):
    return u[0] * v[1] - u[1] * v[0]


def dot(u, v):
    return u[0] * v[0] + u[1] * v[1]


def orient(a, b, c):
    return cross(sub(b, a), sub(c, a))


def param_on(p, a, b):
    """if p lies on closed segment ab return t with p=a+t(b-a), else None."""
    if sgn(orient(a, b, p)) != 0:
        return None
    d = sub(b, a)
    if sgn(d[0]) != 0:
        t = (p[0] - a[0]) / d[0]
    else:
        t = (p[1] - a[1]) / d[1]
    if sgn(t) < 0 or sgn(t - 1) > 0:
        return None
    return t


def seg_meet(f, g):
    """classify the intersection of closed segments f, g: ('none',), ('point', P), ('overlap',)."""
    a, b = f
    c, d = g
    d1, d2 = sgn(orient(c, d, a)), sgn(orient(c, d, b))
    d3, d4 = sgn(orient(a, b, c)), sgn(orient(a, b, d))
    if d1 == 0 and d2 == 0:
        pts = []
        for p in (a, b):
            if param_on(p, c, d) is not None:
                pts.append(p)
        for p in (c, d):
            if param_on(p, a, b) is not None:
                pts.append(p)
        ks = {}
        for p in pts:
            ks[key(p)] = p
        if not ks:
            return ('none',)
        if len(ks) == 1:
            return ('point', list(ks.values())[0])
        return ('overlap',)
    if d1 * d2 > 0 or d3 * d4 > 0:
        return ('none',)
    den = cross(sub(b, a), sub(d, c))
    t = cross(sub(c, a), sub(d, c)) / den
    return ('point', (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))


def half(v):
    """0 for directions in [0,pi), 1 for [pi,2pi)."""
    sx, sy = sgn(v[0]), sgn(v[1])
    return 0 if (sy > 0 or (sy == 0 and sx > 0)) else 1


def cmp_dir(u, v):
    hu, hv = half(u), half(v)
    if hu != hv:
        return hu - hv
    c = sgn(cross(u, v))
    return -c  # u before v if v is counterclockwise of u


def sector_class(u, v):
    """angle of the sector from ray u counterclockwise to ray v: -1 (<pi), 0 (=pi), +1 (>pi)."""
    c = sgn(cross(u, v))
    if c > 0:
        return -1
    if c < 0:
        return 1
    if sgn(dot(u, v)) < 0:
        return 0
    return 2  # same direction: full turn (degree-1 vertex) -- cannot occur in valid configs


def sector_angle(u, v):
    """numeric angle (mpmath) of the sector from u ccw to v, in (0, 2pi]."""
    a = mpmath.atan2(tomp(cross(u, v)), tomp(dot(u, v)))
    if a <= 0:
        a += 2 * mpmath.pi
    return a


def length2(f):
    d = sub(f[1], f[0])
    return dot(d, d)


def analyze(fences, check_unit=True):
    """Check the fence rules and compute the plane graph data.  Returns a dict; raises
    ValueError on a rule violation (crossing, overlap, loose end, non-unit fence)."""
    n = len(fences)
    if check_unit:
        for i, f in enumerate(fences):
            if sgn(length2(f) - 1) != 0:
                raise ValueError(f'fence {i} not of unit length')
    # points on each fence: its ends plus meeting points
    onf = [dict() for _ in range(n)]  # key -> (t, point)
    for i, f in enumerate(fences):
        onf[i][key(f[0])] = (0, f[0])
        onf[i][key(f[1])] = (1, f[1])
    supported = [[False, False] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            m = seg_meet(fences[i], fences[j])
            if m[0] == 'none':
                continue
            if m[0] == 'overlap':
                raise ValueError(f'fences {i},{j} overlap')
            P = m[1]
            kp = key(P)
            ends_i = [e for e in (0, 1) if key(fences[i][e]) == kp]
            ends_j = [e for e in (0, 1) if key(fences[j][e]) == kp]
            if not ends_i and not ends_j:
                raise ValueError(f'fences {i},{j} cross')
            for e in ends_i:
                supported[i][e] = True
            for e in ends_j:
                supported[j][e] = True
            if kp not in onf[i]:
                onf[i][kp] = (param_on(P, *fences[i]), P)
            if kp not in onf[j]:
                onf[j][kp] = (param_on(P, *fences[j]), P)
    for i in range(n):
        for e in (0, 1):
            if not supported[i][e]:
                raise ValueError(f'end {e} of fence {i} is loose')
    # vertices, edges
    pts = {}
    ends_at = {}
    through_at = {}
    edges = []  # (ku, kv, fence)
    for i in range(n):
        lst = sorted(onf[i].items(), key=cmp_to_key(lambda A, B: sgn(A[1][0] - B[1][0])))
        for k, (t, P) in lst:
            pts[k] = P
        ends_at.setdefault(lst[0][0], []).append(i)
        ends_at.setdefault(lst[-1][0], []).append(i)
        for k, _ in lst[1:-1]:
            through_at.setdefault(k, []).append(i)
        for a_, b_ in zip(lst[:-1], lst[1:]):
            edges.append((a_[0], b_[0], i))
    for k, fl in through_at.items():
        if len(fl) > 1:
            raise ValueError('two fences through one point')
    # rays
    rays = {k: [] for k in pts}
    for (ku, kv, i) in edges:
        rays[ku].append(kv)
        rays[kv].append(ku)
    order = {}
    for k, lst in rays.items():
        P = pts[k]
        srt = sorted(lst, key=cmp_to_key(lambda A, B: cmp_dir(sub(pts[A], P), sub(pts[B], P))))
        order[k] = srt
    # junction data
    junctions = {}
    for k in pts:
        P = pts[k]
        rs = order[k]
        d = len(rs)
        secs = []
        for idx in range(d):
            u = sub(pts[rs[idx]], P)
            v = sub(pts[rs[(idx + 1) % d]], P)
            secs.append(sector_class(u, v))
        e = len(ends_at.get(k, []))
        t = 1 if k in through_at else 0
        kk = sum(1 for s in secs if s == -1)
        r = d - kk
        junctions[k] = dict(e=e, t=t, d=d, k=kk, r=r, alpha=2 - r, secs=secs)
    # faces by half-edge traversal
    used = set()
    cycles = []
    for (ku, kv, i) in edges:
        for (a_, b_) in ((ku, kv), (kv, ku)):
            if (a_, b_) in used:
                continue
            cyc = []
            cur = (a_, b_)
            while cur not in used:
                used.add(cur)
                u_, v_ = cur
                rs = order[v_]
                idx = rs.index(u_)
                w_ = rs[(idx - 1) % len(rs)]
                P = pts[v_]
                cls = sector_class(sub(pts[w_], P), sub(pts[u_], P))
                cyc.append((u_, v_, cls, w_))
                cur = (v_, w_)
            cycles.append(cyc)
    info = []
    for cyc in cycles:
        vs = [pts[c[0]] for c in cyc]
        A2 = 0
        per = mpmath.mpf(0)
        for idx in range(len(vs)):
            p, q = vs[idx], vs[(idx + 1) % len(vs)]
            A2 = A2 + cross(p, q)
            per += mpmath.sqrt(tomp(dot(sub(q, p), sub(q, p))))
        kappa = sum(1 for c in cyc if c[2] == -1)
        info.append(dict(cyc=cyc, A2=A2, per=per, kappa=kappa, pts=vs))
    pos = [c for c in info if sgn(c['A2']) > 0]
    neg = [c for c in info if sgn(c['A2']) < 0]
    if any(sgn(c['A2']) == 0 for c in info):
        raise ValueError('degenerate boundary cycle')

    def inside(pt, poly):
        # strict point-in-polygon by ray casting with exact arithmetic (pt not on boundary)
        cnt = 0
        m = len(poly)
        for idx in range(m):
            a_, b_ = poly[idx], poly[(idx + 1) % m]
            if (sgn(a_[1] - pt[1]) > 0) != (sgn(b_[1] - pt[1]) > 0):
                # x-coordinate of crossing compared with pt: sign of orient
                o = sgn(orient(a_, b_, pt))
                if sgn(b_[1] - a_[1]) > 0:
                    if o > 0:
                        cnt += 1
                else:
                    if o < 0:
                        cnt += 1
        return cnt % 2 == 1

    # connected components (union-find over edges); each component has exactly one negative cycle
    par = {k: k for k in pts}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for (ku, kv, i) in edges:
        ru, rv = find(ku), find(kv)
        if ru != rv:
            par[ru] = rv
    for c in info:
        c['comp'] = find(c['cyc'][0][0])
    ncomp = len({find(k) for k in pts})
    if len(neg) != ncomp:
        raise ValueError('negative cycles != components')
    faces = [dict(outer=c, holes=[]) for c in pos]
    outer_face = []
    for h in neg:
        pt = h['pts'][0]
        best = None
        for fi, f in enumerate(faces):
            if f['outer']['comp'] == h['comp']:
                continue
            if inside(pt, f['outer']['pts']):
                if best is None or sgn(f['outer']['A2'] - faces[best]['outer']['A2']) < 0:
                    best = fi
        if best is None:
            outer_face.append(h)
        else:
            faces[best]['holes'].append(h)
    for f in faces:
        f['A2'] = f['outer']['A2'] + sum((h['A2'] for h in f['holes']), 0)
        f['kappa'] = f['outer']['kappa'] + sum(h['kappa'] for h in f['holes'])
        f['h'] = len(f['holes'])
        f['per'] = f['outer']['per'] + sum(h['per'] for h in f['holes'])
    kappa_o = sum(h['kappa'] for h in outer_face)
    c_o = len(outer_face)
    P_o = sum((h['per'] for h in outer_face), mpmath.mpf(0))
    return dict(n=n, junctions=junctions, faces=faces, kappa_o=kappa_o, c_o=c_o, P_o=P_o,
                outer=outer_face, pts=pts, order=order, edges=edges)


def junction_type(J):
    e, t, r, d = J['e'], J['t'], J['r'], J['d']
    if (e, t, r) == (1, 1, 1):
        return 'T'
    if (e, t, r) == (2, 1, 0):
        return 'X'
    if (e, t, r) == (3, 0, 0):
        return 'Y'
    if (e, t, r) == (2, 0, 1):
        return 'L'
    return f'e{e}t{t}r{r}'


def corner_report(R):
    """Exact integer identities of the audit for an analysed configuration."""
    J = R['junctions']
    n = R['n']
    C = len(R['faces'])
    H = sum(f['h'] for f in R['faces'])
    c_o = R['c_o']
    kappa_o = R['kappa_o']
    s_k_minus_a = sum(j['k'] - j['alpha'] for j in J.values())
    s_3a_minus_k = sum(3 * j['alpha'] - j['k'] for j in J.values())
    s_e = sum(j['e'] for j in J.values())
    lhs3 = sum(f['kappa'] - 3 for f in R['faces'])
    ok1 = all(3 * j['alpha'] - j['k'] <= j['e'] for j in J.values())
    ok_e = (s_e == 2 * n)
    ok2 = (2 * C == s_k_minus_a + 2 * H + 2 * c_o)
    ok3 = (2 * lhs3 == s_3a_minus_k - 2 * kappa_o - 6 * H - 6 * c_o)
    ok3b = lhs3 <= n - 3
    types = {}
    for j in J.values():
        tname = junction_type(j)
        types[tname] = types.get(tname, 0) + 1
    tight = all((3 * j['alpha'] - j['k'] == j['e']) == (junction_type(j) in 'TXYL') for j in J.values())
    return dict(C=C, H=H, c_o=c_o, kappa_o=kappa_o, lhs3=lhs3, rhs3=Fr(s_3a_minus_k, 2) - kappa_o - 3 * H - 3 * c_o,
                ok1=ok1, ok_e=ok_e, ok2=ok2, ok3=ok3, ok3b=ok3b, types=types, tight_types_exact=tight)


def angle_identity_check(R):
    """numeric (mpmath) check of the face angle identity for every face and for the outer face."""
    pi = mpmath.pi
    pts = R['pts']
    worst = mpmath.mpf(0)

    def corner_val(c):
        u_, v_, cls, w_ = c
        P = pts[v_]
        th = sector_angle(sub(pts[w_], P), sub(pts[u_], P))
        return th if cls == -1 else th - pi

    for f in R['faces']:
        s = sum(corner_val(c) for cyc in [f['outer']] + f['holes'] for c in cyc['cyc'])
        target = pi * (f['kappa'] - 2 + 2 * f['h'])
        worst = max(worst, abs(s - target))
    s = sum(corner_val(c) for cyc in R['outer'] for c in cyc['cyc'])
    worst = max(worst, abs(s - pi * (R['kappa_o'] + 2 * R['c_o'])))
    return worst
