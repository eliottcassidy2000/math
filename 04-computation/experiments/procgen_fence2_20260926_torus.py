"""procgen_fence2_20260926_torus.py -- exact fence engine for periodic (torus) and finite (plane)
configurations, with the per-fence corner/side typing of lane "fence2".

Session collatz-procgen-20260922, lane "fence2" (2026-09-26).

A configuration is a list of fences ((x0,y0),(x1,y1)) with exact coordinates (Fractions, or the
QS numbers a+b*sqrt(D) of procgen_fencelim_20260926_geom), optionally with a lattice (v1, v2): then
the fences are representatives of the periodic configuration  F + Z v1 + Z v2  (the torus).

analyze(fences, lattice=None) checks the rules (unit lengths, no crossing, no overlap, no loose end,
no two fences through a point) and returns the plane/torus graph with, for every face walk, its
corners typed O / I / E / R (reflex) / S0 (straight, through) / S2 (straight join), its sides with
e = o(start) + iota(end), exact squared side lengths, and the fence-side pieces.

Typing (walk with the face on its left; at a corner, iota = [incoming fence ends here],
o = [outgoing fence ends here]):  O = (o,iota) = (1,0), I = (0,1), E = (1,1) convex, R = reflex.
"""
from fractions import Fraction as Fr
from functools import cmp_to_key
import math

import procgen_fencelim_20260926_geom as G
from procgen_fencelim_20260926_geom import sgn, sub, cross, dot, seg_meet, cmp_dir, sector_class


def _floor(x):
    if isinstance(x, Fr):
        return math.floor(x)
    f = math.floor(float(x))
    # exact correction (x is exact: Fraction or QS)
    while sgn(x - f) < 0:
        f -= 1
    while sgn(x - (f + 1)) >= 0:
        f += 1
    return f


class Lattice:
    def __init__(self, v1, v2):
        self.v1, self.v2 = v1, v2
        self.det = cross(v1, v2)
        if sgn(self.det) == 0:
            raise ValueError('degenerate lattice')

    def coords(self, p):
        a = cross(p, self.v2) / self.det
        b = cross(self.v1, p) / self.det
        return a, b

    def canon(self, p):
        a, b = self.coords(p)
        fa, fb = _floor(a), _floor(b)
        q = (p[0] - fa * self.v1[0] - fb * self.v2[0], p[1] - fa * self.v1[1] - fb * self.v2[1])
        return q, (fa, fb)

    def shift(self, p, off):
        return (p[0] + off[0] * self.v1[0] + off[1] * self.v2[0], p[1] + off[0] * self.v1[1] + off[1] * self.v2[1])


def _shift_f(f, L, off):
    return (L.shift(f[0], off), L.shift(f[1], off))


def _bbox(f):
    xs = sorted((f[0][0], f[1][0]), key=float)
    ys = sorted((f[0][1], f[1][1]), key=float)
    return float(xs[0]), float(xs[1]), float(ys[0]), float(ys[1])


def analyze(fences, lattice=None, R=2, check_unit=True):
    n = len(fences)
    L = Lattice(*lattice) if lattice is not None else None
    if check_unit:
        for i, f in enumerate(fences):
            d = sub(f[1], f[0])
            if sgn(dot(d, d) - 1) != 0:
                raise ValueError(f'fence {i} not of unit length')
    offs = [(0, 0)] if L is None else [(a, b) for a in range(-R, R + 1) for b in range(-R, R + 1)]
    inst = []  # (j, off, segment)
    for j, f in enumerate(fences):
        for off in offs:
            inst.append((j, off, f if L is None else _shift_f(f, L, off)))
    boxes = [_bbox(s) for (_, _, s) in inst]
    onf = [dict() for _ in range(n)]   # key -> (t, point, list of (j,off,endidx) ending there)
    for i, f in enumerate(fences):
        onf[i][G.key(f[0])] = [Fr(0), f[0]]
        onf[i][G.key(f[1])] = [Fr(1), f[1]]
    supported = [[False, False] for _ in range(n)]
    for i, f in enumerate(fences):
        bx = _bbox(f)
        for idx, (j, off, g) in enumerate(inst):
            if j == i and off == (0, 0):
                continue
            by = boxes[idx]
            if by[0] > bx[1] + 1e-9 or by[1] < bx[0] - 1e-9 or by[2] > bx[3] + 1e-9 or by[3] < bx[2] - 1e-9:
                continue
            m = seg_meet(f, g)
            if m[0] == 'none':
                continue
            if m[0] == 'overlap':
                raise ValueError(f'fences {i},{j}{off} overlap')
            P = m[1]
            kp = G.key(P)
            ends_i = [e for e in (0, 1) if G.key(f[e]) == kp]
            ends_g = [e for e in (0, 1) if G.key(g[e]) == kp]
            if not ends_i and not ends_g:
                raise ValueError(f'fences {i},{j}{off} cross')
            for e in ends_i:
                supported[i][e] = True
            if kp not in onf[i]:
                onf[i][kp] = [G.param_on(P, *f), P]
    for i in range(n):
        for e in (0, 1):
            if not supported[i][e]:
                raise ValueError(f'end {e} of fence {i} is loose')
    # canonical vertices and edges
    def ck(p):
        if L is None:
            return G.key(p), (0, 0)
        q, off = L.canon(p)
        return G.key(q), off
    vpos = {}          # canonical key -> canonical point
    fence_pts = []     # per fence: sorted list of points
    edges = []         # (i, k, p_k, p_k+1)
    for i in range(n):
        lst = sorted(onf[i].values(), key=cmp_to_key(lambda A, B: sgn(A[0] - B[0])))
        pts = [p for (_, p) in lst]
        fence_pts.append(pts)
        for k in range(len(pts) - 1):
            edges.append((i, k))
        for p in pts:
            kk, off = ck(p)
            if kk not in vpos:
                vpos[kk] = p if off == (0, 0) else L.shift(p, (-off[0], -off[1]))
    # rays: half-edge h = (i, k, s): s=+1 leaves p_k towards p_k+1, s=-1 leaves p_k+1 towards p_k
    rays = {kk: [] for kk in vpos}
    ends_at = {kk: 0 for kk in vpos}
    thr_at = {kk: set() for kk in vpos}
    for i in range(n):
        pts = fence_pts[i]
        m = len(pts) - 1
        ends_at[ck(pts[0])[0]] += 1
        ends_at[ck(pts[-1])[0]] += 1
        for k in range(1, m):
            thr_at[ck(pts[k])[0]].add(i)
        for k in range(m):
            a, b = pts[k], pts[k + 1]
            rays[ck(a)[0]].append((sub(b, a), (i, k, 1)))
            rays[ck(b)[0]].append((sub(a, b), (i, k, -1)))
    for kk, s in thr_at.items():
        if len(s) > 1:
            raise ValueError('two fences through one point')
    order = {}
    for kk, lst in rays.items():
        srt = sorted(lst, key=cmp_to_key(lambda A, B: cmp_dir(A[0], B[0])))
        order[kk] = srt
    pos_in = {}
    for kk, lst in order.items():
        for idx, (d, h) in enumerate(lst):
            pos_in[h] = (kk, idx)

    def h_start(h):
        i, k, s = h
        return fence_pts[i][k] if s == 1 else fence_pts[i][k + 1]

    def h_end(h):
        i, k, s = h
        return fence_pts[i][k + 1] if s == 1 else fence_pts[i][k]

    def h_rev(h):
        return (h[0], h[1], -h[2])

    def o_leave(h):   # does fence of h end at the start of h
        i, k, s = h
        m = len(fence_pts[i]) - 1
        return int(k == 0) if s == 1 else int(k + 1 == m)

    def iota_arrive(h):  # does fence of h end at the end of h
        i, k, s = h
        m = len(fence_pts[i]) - 1
        return int(k + 1 == m) if s == 1 else int(k == 0)

    used = set()
    walks = []
    for (i, k) in edges:
        for s in (1, -1):
            h0 = (i, k, s)
            if h0 in used:
                continue
            walk = []
            h = h0
            cur = h_start(h)
            path = [cur]
            while h not in used:
                used.add(h)
                d = sub(h_end(h), h_start(h))
                nxt_pt = (cur[0] + d[0], cur[1] + d[1])
                kk, idx = pos_in[h_rev(h)]
                lst = order[kk]
                hn = lst[(idx - 1) % len(lst)][1]
                back = lst[idx][0]
                outd = lst[(idx - 1) % len(lst)][0]
                cls = sector_class(outd, back)
                walk.append(dict(h_in=h, h_out=hn, cls=cls, iota=iota_arrive(h), o=o_leave(hn),
                                 vkey=kk, pt=nxt_pt, back=back, outd=outd))
                cur = nxt_pt
                path.append(cur)
                h = hn
            closed = (sgn(path[-1][0] - path[0][0]) == 0 and sgn(path[-1][1] - path[0][1]) == 0)
            A2 = 0
            for a, b in zip(path[:-1], path[1:]):
                A2 = A2 + cross(a, b)
            walks.append(dict(corners=walk, closed=closed, A2=A2, start=path[0]))
    # junction data
    junctions = {}
    for kk, lst in order.items():
        d = len(lst)
        secs = [sector_class(lst[idx][0], lst[(idx + 1) % d][0]) for idx in range(d)]
        junctions[kk] = dict(e=ends_at[kk], t=len(thr_at[kk]), d=d, secs=secs)
    return dict(n=n, lattice=L, fences=fences, fence_pts=fence_pts, walks=walks, junctions=junctions,
                order=order, vpos=vpos, h_start=h_start, h_end=h_end)


# ------------------------------------------------------------------------------------ typing
def corner_type(c):
    if c['cls'] == 1:
        return 'R'
    if c['cls'] == 0:
        return 'S2' if (c['o'], c['iota']) == (1, 1) else ('S0' if (c['o'], c['iota']) == (0, 0) else 'S?')
    return {(1, 0): 'O', (0, 1): 'I', (1, 1): 'E', (0, 0): 'C00'}[(c['o'], c['iota'])]


def walk_sides(w):
    """sides of a walk between consecutive non-straight corners.  Returns list of dicts with
    start/end corner index, e, j (straight joins inside), squared length (exact), vector."""
    cs = w['corners']
    m = len(cs)
    ns = [idx for idx in range(m) if cs[idx]['cls'] != 0]
    out = []
    if not ns:
        return out
    for a_i, a in enumerate(ns):
        b = ns[(a_i + 1) % len(ns)]
        # corners strictly between a and b (cyclically) are straight
        j = 0
        idx = (a + 1) % m
        while idx != b:
            if corner_type(cs[idx]) == 'S2':
                j += 1
            idx = (idx + 1) % m
        pa, pb = cs[a]['pt'], cs[b]['pt']
        # unrolled positions: walk positions are consecutive; recompute vector by summing edges
        vec = (0, 0)
        idx = a
        while True:
            h = cs[(idx + 1) % m]['h_in']
            idx = (idx + 1) % m
            # edge vector of h
            vec = (vec[0] + (cs[idx]['pt'][0] - cs[(idx - 1) % m]['pt'][0]), vec[1] + (cs[idx]['pt'][1] - cs[(idx - 1) % m]['pt'][1]))
            if idx == b:
                break
        e = cs[a]['o'] + cs[b]['iota']
        out.append(dict(a=a, b=b, e=e, j=j, len2=dot(vec, vec), vec=vec,
                        ta=corner_type(cs[a]), tb=corner_type(cs[b])))
    return out
