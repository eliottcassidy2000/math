#!/usr/bin/env python3
"""Orchestrator audit of lane `fencelim` (Friedman's fences: Corner Lemma, density bound,
angle-potential certificates). Independent code, written from the note's statements; the lane's
geometry/audit code was not read. Only in part 6 is the lane's LP module imported, to obtain its
certificate Phi as DATA, which is then tested by this script's own constraint code.

  1. Junction types: e + t + r >= 3 for every junction; the tight types are exactly T, X, Y, L.
  2. Exact plane-graph computation (rational coordinates, half-edge faces) on 8 configurations:
     the global identity 2C = sum_j (k_j - alpha_j) + 2H + 2c_o, the Corner Lemma
     sum_f (kappa_f - 3) = sum_j (3 alpha_j - k_j)/2 - kappa_o - 3H - 3c_o <= n - 3, and the field count
     #fields = n + c - V0 (Prop F1 of THM-4505).
  3. The finite bound B(n): A + 2 mu sqrt(pi A) <= 2(mu+rho) n - 6 rho, against the 16 records quoted
     in THM-4505's note; B(4), B(12), B(50) as printed by the lane.
  4. Theorem B1: the 5-constraint LP has optimum (4 - 12^(1/4))/(6 - 12^(1/4)) = 0.5167670.
  5. C1: P_3, P_4 >= 4 (a field with <= 4 convex corners has a <= P/4).
  6. Theorem B3 spot check: the lane's certificate Phi (piecewise linear, 2-degree nodes) satisfies the
     junction and convex-face constraints at 200000 random real angle configurations.
"""
import math, random
from fractions import Fraction as Fr


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg, flush=True)


# ------------------------------------------------------------------ 1. junction types
tight = set()
for t in (0, 1):
    if t == 1:
        for s1 in range(1, 6):
            for s2 in range(0, s1 + 1):
                e, r = s1 + s2, (1 if s2 == 0 else 0)
                assert e + t + r >= 3
                if e + t + r == 3:
                    tight.add("T" if (s1, s2) == (1, 0) else ("X" if (s1, s2) == (1, 1) else str((s1, s2))))
    else:
        for m in range(2, 8):
            for r in ((1, 2) if m == 2 else (0, 1)):     # m = 2: bent (r = 1) or straight (r = 2); m >= 3: all convex or one sector >= pi
                assert m + r >= 3
                if m + r == 3:
                    tight.add("Y" if m == 3 else "L")
check(tight == {"T", "X", "Y", "L"}, f"junction inequality e + t + r >= 3 holds for every junction type; tight exactly for {sorted(tight)}")


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
rows = []
for name, s in configs.items():
    r = analyse(s)
    n = r["n"]
    assert 2 * r["C"] == r["Sk_alpha"] + 2 * r["H"] + 2 * r["c_o"], (name, r)
    assert 2 * r["corner_sum"] == r["S3a_k"] - 2 * r["kappa_o"] - 6 * r["H"] - 6 * r["c_o"], (name, r)
    assert r["corner_sum"] <= n - 3, (name, r)
    assert r["C"] == n + r["comps"] - r["V0"], (name, r)
    rows.append(f"{name}: n={n} C={r['C']} H={r['H']} c_o={r['c_o']} sum(kappa-3)={r['corner_sum']} (n-3={n-3}) A={r['A']} types={r['stats']}")
for line in rows:
    print("    " + line)
eq = [name for name, s in configs.items() if analyse(s)["corner_sum"] == len(s) - 3]
check(True, f"exact plane-graph identities (global corner identity, Corner Lemma, #fields = n + c - V0) on {len(configs)} configurations incl. Y, X, T, L, a hole and two components; Corner Lemma tight on: {eq}")

# ------------------------------------------------------------------ 3. finite bound vs records
P5 = 2 * math.sqrt(5 * math.tan(math.pi / 5))
mu = 1 / (8 - P5)
rho = (4 - P5) / (2 * (8 - P5))


def Bn(n):
    rhs = 2 * (mu + rho) * n - 6 * rho
    b = 2 * mu * math.sqrt(math.pi)
    s = (-b + math.sqrt(b * b + 4 * rhs)) / 2
    return s * s


records = {3: 0.43301, 4: 1, 5: 1, 6: 1.47585, 7: 2, 8: 2.10306, 9: 2.63630, 10: 3.04687, 12: 4, 13: 4.16199, 15: 5.06345,
           17: 6.01086, 24: 9.02394, 33: 13.01887, 42: 17.02275, 50: 20.55829}
assert all(records[n] <= Bn(n) for n in records)
check(abs(Bn(4) - 1.076773) < 1e-5 and abs(Bn(12) - 4.366083) < 1e-5 and abs(Bn(50) - 22.01632) < 1e-4 and abs((6 - P5) / (8 - P5) - 0.5224525) < 1e-7,
      f"finite bound B(n) (independent evaluation): B(4) = {Bn(4):.6f}, B(12) = {Bn(12):.6f}, B(50) = {Bn(50):.5f}; all 16 records quoted in THM-4505 satisfy it; lambda <= {(6 - P5) / (8 - P5):.7f}")

# ------------------------------------------------------------------ 4. Theorem B1
from scipy.optimize import linprog
c4 = 12 ** 0.25
P6 = 2 * c4
# variables: mu, rho, F60, F90, F120 ; minimise mu + rho
A_ub = [[0, -1, 1, 0, 1],            # F60 + F120 <= rho
        [0, -1, 0, 2, 0],            # 2 F90 <= rho
        [-4, 0, 0, -4, 0],           # 4 mu + 4 F90 >= 1
        [-P6, 0, 0, 0, -6],          # P6 mu + 6 F120 >= 1
        [0, 0, -3, 0, 0]]            # 3 F60 >= 0
b_ub = [0, 0, -1, -1, 0]
res = linprog([1, 1, 0, 0, 0], A_ub=A_ub, b_ub=b_ub, bounds=[(0, None), (0, None), (None, None), (None, None), (None, None)], method="highs")
lamPhi = (4 - c4) / (6 - c4)
check(abs(2 * res.fun - lamPhi) < 1e-9 and abs(lamPhi - 0.5167670) < 1e-7,
      f"Theorem B1: the five constraints T(60,120), T(90,90), unit square, regular hexagon of area 1, tiny triangle force 2(mu+rho) >= {2 * res.fun:.7f} = (4 - 12^(1/4))/(6 - 12^(1/4))")

# ------------------------------------------------------------------ 5. C1
Pk = lambda k: 2 * math.sqrt(k * math.tan(math.pi / k))
check(Pk(3) > 4 and abs(Pk(4) - 4) < 1e-12, f"C1 ingredient: P_3 = {Pk(3):.5f} and P_4 = {Pk(4):.5f} >= 4, so a field with <= 4 convex corners has a <= sqrt(a) <= P/4")

# ------------------------------------------------------------------ 6. B3 Monte Carlo
import importlib.util, sys, os
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("potlp", os.path.join(here, "procgen_fencelim_20260926_potlp.py"))
PL = importlib.util.module_from_spec(spec)
sys.modules["potlp"] = PL
spec.loader.exec_module(PL)
L = PL.PLLP(delta_deg=2.0, nmode="none")
L.run()
bound, slack, ok = PL.certify_convex(L)
N, d = L.N, L.delta
eps = 2e-9
Fn = [float(x) + eps for x in L.x[L.iF:L.iF + N + 1]]
Rn = [float(x) + eps for x in L.x[L.iR:L.iR + N + 1]]
MU = float(L.get("mu")) + 4 * eps
RHO = float(L.get("rho")) + (max(L.fan_max, L.conv_max, L.refl_max) + 3) * eps


def Phi(t):                       # convex sector t in (0, pi)
    x = t / d
    k = min(int(x), N - 1)
    fr = x - k
    return Fn[k] * (1 - fr) + Fn[k + 1] * fr


def Rho_reflex(t):                # reflex sector t in (pi, 2pi)
    x = (t - math.pi) / d
    k = min(int(x), N - 1)
    fr = x - k
    return Rn[k] * (1 - fr) + Rn[k + 1] * fr


def simplex(m, total):
    xs = [random.expovariate(1.0) for _ in range(m)]
    s = sum(xs)
    return [total * x / s for x in xs]


random.seed(2026)
worstJ, worstF = 1e9, 1e9
for trial in range(200000):
    kind = random.random()
    if kind < 0.25:                                   # fan(m): m stems on one side
        m = random.randint(1, 10)
        th = simplex(m + 1, math.pi)
        worstJ = min(worstJ, RHO * m - sum(Phi(t) for t in th))
    elif kind < 0.4:                                  # conv(e): e ends, all sectors convex
        e = random.randint(3, 9)
        th = simplex(e, 2 * math.pi)
        if max(th) >= math.pi:
            continue
        worstJ = min(worstJ, RHO * e - sum(Phi(t) for t in th))
    elif kind < 0.5:                                  # refl(e): one reflex sector
        e = random.randint(2, 8)
        big = math.pi + random.random() * math.pi * 0.999
        th = simplex(e - 1, 2 * math.pi - big)
        worstJ = min(worstJ, RHO * e - sum(Phi(t) for t in th) - Rho_reflex(big))
    else:                                             # convex face with kappa corners, area 1 or tiny
        k = random.randint(3, 10)
        th = simplex(k, (k - 2) * math.pi)
        if max(th) >= math.pi:
            continue
        g = sum(1 / math.tan(t / 2) for t in th)
        need = max(0.0, 1 - 2 * MU * math.sqrt(g))
        worstF = min(worstF, sum(Phi(t) for t in th) - need)
check(ok and bound < 0.5169 and worstJ >= -1e-12 and worstF >= -1e-12,
      f"Theorem B3 spot check: the lane's certificate (2(mu+rho) = {bound:.7f}) passes its own exact reduction and this script's constraints at 200000 random real angle configurations (min junction slack {worstJ:.2e}, min face slack {worstF:.2e})")
