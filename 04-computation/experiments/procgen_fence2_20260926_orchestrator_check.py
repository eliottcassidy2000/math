#!/usr/bin/env python3
"""Orchestrator audit of lane `fence2` (per-fence accounting for Friedman's fences). Independent code: the
orchestrator's own exact half-edge face computation (from procgen_fencelim_20260926_orchestrator_check.py),
extended to track which fence carries each edge; the lane's scripts were not read.

For every boundary walk of every field (outer walks and hole walks), with corners = non-straight vertices and
sides = maximal straight runs between corners:
  (P1) every convex corner has an adjacent side whose fence ends there; at every reflex corner both do;
  (P2 identity) #whole sides - #through sides = #convex corners where both adjacent fences end + #reflex corners,
       where e(side) = number of its two corners at which the fence carrying it ends, whole: e = 2, through: e = 0;
  (length) a whole side without an interior straight join has length exactly 1;
  on 9 exact rational configurations (T, X, Y, L junctions, a hole, two components, straight joins).
"""
import math
from fractions import Fraction as Fr


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg, flush=True)


# ------------------------------------------------------------------ 2. exact faces
def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def sub(p, q):
    return (p[0] - q[0], p[1] - q[1])


def half(v):          # 0 for angle in [0, pi), 1 for [pi, 2pi)
    return 0 if (v[1] > 0 or (v[1] == 0 and v[0] > 0)) else 1


def ang_key(v):
    return half(v)


def sort_ccw(vs):
    import functools

    def cmp(a, b):
        ha, hb = half(a), half(b)
        if ha != hb:
            return ha - hb
        c = cross(a, b)
        return -1 if c > 0 else (1 if c < 0 else 0)
    return sorted(vs, key=functools.cmp_to_key(cmp))


def on_segment_interior(p, a, b):
    if cross(sub(b, a), sub(p, a)) != 0:
        return False
    dot = (p[0] - a[0]) * (b[0] - a[0]) + (p[1] - a[1]) * (b[1] - a[1])
    L2 = (b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2
    return 0 < dot < L2


def analyse(segs):
    n = len(segs)
    for (a, b) in segs:
        assert (b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2 == 1, "not unit"
    ends = [p for s in segs for p in s]
    junc = sorted(set(ends))
    # rules: every end on another fence; pairwise intersections only at an end of one of them
    for i, (a, b) in enumerate(segs):
        for p in (a, b):
            assert any(j != i and (p in segs[j] or on_segment_interior(p, *segs[j])) for j in range(n)), ("dangling", p)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = segs[i]
            c, d = segs[j]
            d1, d2 = cross(sub(b, a), sub(c, a)), cross(sub(b, a), sub(d, a))
            d3, d4 = cross(sub(d, c), sub(a, c)), cross(sub(d, c), sub(b, c))
            if d1 * d2 < 0 and d3 * d4 < 0:
                raise AssertionError("crossing")
            if d1 == 0 and d2 == 0:                   # collinear: may only touch at ends
                assert not (on_segment_interior(c, a, b) or on_segment_interior(d, a, b) or on_segment_interior(a, c, d) or on_segment_interior(b, c, d) or {a, b} == {c, d})
    # split fences at junction points
    edges = set()
    for (a, b) in segs:
        pts = [a, b] + [p for p in junc if on_segment_interior(p, a, b)]
        pts.sort(key=lambda p: (p[0] - a[0]) * (b[0] - a[0]) + (p[1] - a[1]) * (b[1] - a[1]))
        for u, v in zip(pts, pts[1:]):
            edges.add((u, v))
            edges.add((v, u))
    out = {}
    for (u, v) in edges:
        out.setdefault(u, []).append(v)
    order = {u: sort_ccw([sub(v, u) for v in vs]) for u, vs in out.items()}
    nbr = {u: [(u[0] + d[0], u[1] + d[1]) for d in order[u]] for u in out}
    # junction statistics
    stats = {}
    Sk_alpha = 0
    S3a_k = 0
    V0 = 0
    for u in out:
        ds = order[u]
        dcount = len(ds)
        e = sum(1 for s in segs for p in s if p == u)
        t = 1 if any(on_segment_interior(u, *s) for s in segs) else 0
        assert dcount == e + 2 * t
        r = 0
        for i in range(dcount):
            a, b = ds[i], ds[(i + 1) % dcount]
            if dcount == 1:
                r += 1
                continue
            c = cross(a, b)
            if c < 0 or (c == 0 and (a[0] * b[0] + a[1] * b[1]) < 0) or (dcount == 2 and c == 0):
                r += 1
            elif c == 0 and a == b:
                r += 1
        # sector from a to b CCW is >= pi iff cross(a,b) <= 0 (with the 2-direction straight case)
        k = dcount - r
        alpha = 2 - r
        Sk_alpha += k - alpha
        S3a_k += 3 * alpha - k
        if t == 0:
            V0 += 1
        typ = ("T" if (e, t, r) == (1, 1, 1) else "X" if (e, t, r) == (2, 1, 0) else "Y" if (e, t, r) == (3, 0, 0) else "L" if (e, t, r) == (2, 0, 1) else f"e{e}t{t}r{r}")
        stats[typ] = stats.get(typ, 0) + 1
    # faces: next half-edge = previous (clockwise) neighbour of the twin at v
    used = set()
    walks = []
    for h in edges:
        if h in used:
            continue
        walk = []
        cur = h
        while cur not in used:
            used.add(cur)
            walk.append(cur)
            u, v = cur
            lst = nbr[v]
            i = lst.index(u)
            w = lst[(i - 1) % len(lst)]
            cur = (v, w)
        walks.append(walk)
    info = []
    for walk in walks:
        area2 = sum(cross(u, v) for (u, v) in walk)
        kappa = 0
        for (u, v), (v2, w) in zip(walk, walk[1:] + walk[:1]):
            a, b = sub(w, v), sub(u, v)                     # face sector from (v->w) CCW to (v->u)
            c = cross(a, b)
            if c > 0:
                kappa += 1
        info.append((Fr(area2, 2), kappa, walk))
    pos = [x for x in info if x[0] > 0]
    neg = [x for x in info if x[0] < 0]

    def inside(p, walk):
        cnt = 0
        for (u, v) in walk:
            if (u[1] > p[1]) != (v[1] > p[1]):
                xint = u[0] + (p[1] - u[1]) * (v[0] - u[0]) / (v[1] - u[1])
                if xint > p[0]:
                    cnt += 1
        return cnt % 2 == 1
    C = len(pos)
    kap = {id(x): x[1] for x in pos}
    H = 0
    c_o = 0
    kappa_o = 0
    for (ar, kk, walk) in neg:
        p = walk[0][0]
        cands = [x for x in pos if not any(p == u for (u, v) in x[2]) and inside(p, x[2])]
        if cands:
            host = min(cands, key=lambda x: x[0])
            kap[id(host)] += kk
            H += 1
        else:
            c_o += 1
            kappa_o += kk
    corner_sum = sum(kap[id(x)] - 3 for x in pos)
    # components (union-find over edges)
    par = {u: u for u in out}

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for (u, v) in edges:
        par[f(u)] = f(v)
    comps = len({f(u) for u in out})
    A = sum(x[0] for x in pos)
    return dict(n=n, C=C, H=H, c_o=c_o, kappa_o=kappa_o, corner_sum=corner_sum, Sk_alpha=Sk_alpha, S3a_k=S3a_k,
                V0=V0, comps=comps, A=A, stats=stats)


F0, F1 = Fr(0), Fr(1)


def P(x, y):
    return (Fr(x), Fr(y))


def seg(a, b):
    return (P(*a), P(*b))


def grid(w, h):
    s = []
    for j in range(h + 1):
        for i in range(w):
            s.append(seg((i, j), (i + 1, j)))
    for i in range(w + 1):
        for j in range(h):
            s.append(seg((i, j), (i, j + 1)))
    return s


configs = {
    "unit square": grid(1, 1),
    "square + mid chord": grid(1, 1) + [seg((Fr(1, 2), 0), (Fr(1, 2), 1))],
    "n=5 record: square + chord (0,3/5)-(4/5,0)": grid(1, 1) + [seg((0, Fr(3, 5)), (Fr(4, 5), 0))],
    "grid 1x2": grid(2, 1),
    "grid 2x2": grid(2, 2),
    "Y junction + parallelogram (n=7)": [seg((0, 0), (1, 0)), seg((0, 0), (0, 1)), seg((0, 0), (Fr(-3, 5), Fr(-4, 5))),
                                          seg((0, 1), (1, 1)), seg((1, 0), (1, 1)),
                                          seg((Fr(-3, 5), Fr(-4, 5)), (Fr(2, 5), Fr(-4, 5))), seg((Fr(2, 5), Fr(-4, 5)), (1, 0))],
    "X junction: 2x1 rectangle, middle vertical fence with two stems at its midpoint": [seg((0, 0), (1, 0)), seg((1, 0), (2, 0)), seg((0, 1), (1, 1)), seg((1, 1), (2, 1)),
                                                           seg((0, 0), (0, 1)), seg((2, 0), (2, 1)), seg((1, 0), (1, 1)),
                                                           seg((0, Fr(1, 2)), (1, Fr(1, 2))), seg((1, Fr(1, 2)), (2, Fr(1, 2)))],
    "nested: unit square inside a 3x3 frame (hole)": [seg((i, 0), (i + 1, 0)) for i in range(3)] + [seg((i, 3), (i + 1, 3)) for i in range(3)] +
                                                     [seg((0, j), (0, j + 1)) for j in range(3)] + [seg((3, j), (3, j + 1)) for j in range(3)] + grid(1, 1)[:0] +
                                                     [seg((1, 1), (2, 1)), seg((1, 2), (2, 2)), seg((1, 1), (1, 2)), seg((2, 1), (2, 2))],
    "two disjoint squares": grid(1, 1) + [seg((3, 0), (4, 0)), seg((3, 1), (4, 1)), seg((3, 0), (3, 1)), seg((4, 0), (4, 1))],
}


def walks_with_fences(segs):
    """half-edge faces as in analyse(), keeping the fence index of every edge."""
    n = len(segs)
    ends = [p for s in segs for p in s]
    junc = sorted(set(ends))
    edges = {}
    for i, (a, b) in enumerate(segs):
        pts = [a, b] + [p for p in junc if on_segment_interior(p, a, b)]
        pts.sort(key=lambda p: (p[0] - a[0]) * (b[0] - a[0]) + (p[1] - a[1]) * (b[1] - a[1]))
        for u, v in zip(pts, pts[1:]):
            edges[(u, v)] = i
            edges[(v, u)] = i
    out = {}
    for (u, v) in edges:
        out.setdefault(u, []).append(v)
    order = {u: sort_ccw([sub(v, u) for v in vs]) for u, vs in out.items()}
    nbr = {u: [(u[0] + d[0], u[1] + d[1]) for d in order[u]] for u in out}
    used = set()
    walks = []
    for h in edges:
        if h in used:
            continue
        walk = []
        cur = h
        while cur not in used:
            used.add(cur)
            walk.append(cur)
            u, v = cur
            lst = nbr[v]
            i = lst.index(u)
            w = lst[(i - 1) % len(lst)]
            cur = (v, w)
        walks.append(walk)
    return edges, walks


def audit_walk(walk, edges, segs):
    """returns (ok_P1, lhs, rhs, whole_len_ok) for one boundary walk (face on the left)."""
    m = len(walk)
    kind = []                       # at the vertex between walk[i] and walk[i+1]
    for i in range(m):
        (u, v), (v2, w) = walk[i], walk[(i + 1) % m]
        c = cross(sub(v, u), sub(w, v))
        if c > 0:
            kind.append("convex")
        elif c < 0:
            kind.append("reflex")
        else:
            d = (v[0] - u[0]) * (w[0] - v[0]) + (v[1] - u[1]) * (w[1] - v[1])
            kind.append("straight" if d > 0 else "turnaround")
    corners = [i for i in range(m) if kind[i] != "straight"]
    if not corners:
        return True, 0, 0, True
    ends_at = lambda edge, pt: pt in segs[edges[edge]]
    ok = True
    both = 0
    refl = 0
    for i in corners:
        e_in = walk[i]
        e_out = walk[(i + 1) % m]
        v = e_in[1]
        a = ends_at(e_in, v)
        b = ends_at(e_out, v)
        if kind[i] == "convex":
            ok &= (a or b)
            both += (a and b)
        elif kind[i] == "reflex":
            ok &= (a and b)
            refl += 1
    # sides: from corner i (outgoing edge walk[i+1]) to the next corner j (incoming edge walk[j])
    whole = through = 0
    len_ok = True
    K = len(corners)
    for t in range(K):
        i, j = corners[t], corners[(t + 1) % K]
        first = walk[(i + 1) % m]
        last = walk[j]
        e = int(ends_at(first, first[0])) + int(ends_at(last, last[1]))
        # interior straight joins: a fence change inside the side
        idx = []
        k = (i + 1) % m
        while True:
            idx.append(edges[walk[k]])
            if k == j:
                break
            k = (k + 1) % m
        joins = sum(1 for x, y in zip(idx, idx[1:]) if x != y)
        if e == 2:
            whole += 1
            if joins == 0:
                p, q = first[0], last[1]
                len_ok &= ((q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2 == 1)
        elif e == 0:
            through += 1
    return ok, whole - through, both + refl, len_ok



tot = 0
for name, s in configs.items():
    r = analyse(s)                      # validates the configuration (rules, exact faces)
    edges, walks = walks_with_fences(s)
    for walk in walks:
        area2 = sum(cross(u, v) for (u, v) in walk)
        ok, lhs, rhs, lok = audit_walk(walk, edges, s)
        if area2 > 0 or True:           # every walk: fields (area > 0), holes and outer walks (area < 0)
            assert ok, (name, "P1")
            assert lhs == rhs, (name, lhs, rhs)
            assert lok, (name, "whole side length")
            tot += 1
check(True, f"P1, the P2 walk identity (#whole - #through = #convex both-end corners + #reflex corners) and whole-side length 1 hold on all {tot} boundary walks of {len(configs)} exact configurations")
