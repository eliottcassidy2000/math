"""procgen_fencelim_20260926_typed.py -- Task C helpers: typed pentagons, unit sides, Cairo geometry,
unit-side mean field.  Session collatz-procgen-20260922, lane "fencelim" (2026-09-26).

Corner types (tight junctions): at a Y junction every corner is end-end (EE); at a T (or X) junction
the two corners beside the stem are end-through (one side is the stem, which ends there, the other
side is the through fence).  A face side whose fence ends at both of its corners is a whole fence,
so it has length exactly 1 (at least 1 in general: a concatenation of whole fences).
"""
import itertools
import math
import numpy as np
from scipy.optimize import minimize, brentq, linprog


def polygon_from_angles(angles_deg, lengths):
    """vertices of the polygon with the given interior angles (in order) and side lengths
    (side i joins corner i to corner i+1); returns (closure_error, area)."""
    k = len(angles_deg)
    d = 0.0
    P = [np.zeros(2)]
    for i in range(k - 1):
        P.append(P[-1] + lengths[i] * np.array([math.cos(math.radians(d)), math.sin(math.radians(d))]))
        d += 180 - angles_deg[i + 1]
    last = P[-1] + lengths[k - 1] * np.array([math.cos(math.radians(d)), math.sin(math.radians(d))])
    P = np.array(P)
    A = 0.5 * np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1])
    return np.linalg.norm(last), A


def solve_two_free(angles_deg, unit):
    """side lengths with unit[i] -> 1 and exactly two free sides fixed by closure."""
    k = len(angles_deg)
    d = 0.0
    dirs = []
    for i in range(k):
        dirs.append(d)
        d += 180 - angles_deg[(i + 1) % k]
    u = np.array([[math.cos(math.radians(t)), math.sin(math.radians(t))] for t in dirs])
    free = [i for i in range(k) if not unit[i]]
    if len(free) != 2:
        return None
    fixed = sum(u[i] for i in range(k) if unit[i])
    M = np.array([u[free[0]], u[free[1]]]).T
    if abs(np.linalg.det(M)) < 1e-12:
        return None
    x = np.linalg.solve(M, -fixed)
    L = np.ones(k)
    L[free[0]], L[free[1]] = x
    P = [np.zeros(2)]
    for i in range(k - 1):
        P.append(P[-1] + L[i] * u[i])
    P = np.array(P)
    A = 0.5 * np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1])
    return L, A


def typed_patterns():
    """YYTYT / YYYTT pentagons: for each T corner choose which adjacent side is the stem; the unit
    sides are those whose fence ends at both corners.  Keep patterns with exactly 3 unit sides."""
    out = []
    for order in (['Y', 'Y', 'T', 'Y', 'T'], ['Y', 'Y', 'Y', 'T', 'T']):
        Tidx = [i for i, c in enumerate(order) if c == 'T']
        for stems in itertools.product((0, 1), repeat=2):
            def isE(c, as_start):
                if order[c] == 'Y':
                    return True
                st = stems[Tidx.index(c)]
                return (st == 1) if as_start else (st == 0)
            unit = [isE(i, True) and isE((i + 1) % 5, False) for i in range(5)]
            if sum(unit) == 3:
                out.append((order, stems, unit))
    return out


def best_balanced_typed_pentagon(starts=120, seed=0):
    """max a/P over typed pentagons (3 unit sides forced by the corner types) whose Y angles sum to
    360 and T angles to 180 (a Cairo-like face whose own corners close up its junctions), a <= 1."""
    rng = np.random.default_rng(seed)
    results = []
    for order, stems, unit in typed_patterns():
        Yi = [i for i, c in enumerate(order) if c == 'Y']
        Ti = [i for i, c in enumerate(order) if c == 'T']

        def neg(v):
            ang = [0.0] * 5
            ang[Yi[0]], ang[Yi[1]], ang[Yi[2]] = v[0], v[1], 360 - v[0] - v[1]
            ang[Ti[0]], ang[Ti[1]] = v[2], 180 - v[2]
            if min(ang) <= 0.5 or max(ang) >= 179.5:
                return 1.0
            r = solve_two_free(ang, unit)
            if r is None:
                return 1.0
            L, A = r
            if np.any(L <= 1e-9) or np.any(L > 1 + 1e-12) or A <= 0:
                return 1.0
            if A > 1:
                return 1.0 + (A - 1)
            return -A / np.sum(L)
        best = 0.0
        for _ in range(starts):
            v0 = [rng.uniform(60, 170), rng.uniform(60, 170), rng.uniform(20, 160)]
            r = minimize(neg, v0, method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-13, 'maxiter': 4000})
            best = max(best, -r.fun)
        results.append((''.join(order), stems, best))
    return results


def cyclic_area(sides):
    s = np.array(sides, float)
    M = s.max()
    f = lambda R: np.sum(2 * np.arcsin(np.minimum(1, s / (2 * R)))) - 2 * np.pi
    if f(M / 2 * (1 + 1e-15)) >= 0:
        R = brentq(f, M / 2 * (1 + 1e-15), 1e6)
        ang = 2 * np.arcsin(s / (2 * R))
        return 0.5 * R * R * np.sum(np.sin(ang))
    g = lambda R: np.sum(2 * np.arcsin(np.minimum(1, s / (2 * R)))) - 4 * np.arcsin(min(1, M / (2 * R)))
    R = brentq(g, M / 2 * (1 + 1e-15), 1e6)
    ang = 2 * np.arcsin(s / (2 * R))
    i = int(np.argmax(s))
    return 0.5 * R * R * (np.sum(np.sin(ang)) - 2 * np.sin(ang[i]))


def min_perimeter_unit_sides(k, u):
    """least perimeter of a k-gon of area 1 with u sides of length 1 (the other k-u sides equal):
    by the cyclic-polygon theorem (max area for given sides is the cyclic one; CITED) this is the
    root w of cyclic_area([1]*u + [w/(k-u)]*(k-u)) = 1."""
    lo = 1e-3 if u >= 2 else 1.0 + 1e-9     # the polygon must close: longest side < sum of the others
    w = brentq(lambda w: cyclic_area([1.0] * u + [w / (k - u)] * (k - u)) - 1, lo, (k - u) * 0.999)
    return u + w


def unit_side_mean_field(P3, P2):
    """LP over face types (per fence): P(3,2) pentagons [cost 2 fences, perimeter allowance 4, P=P3],
    P(2,3) pentagons [cost 1.75, allowance 3.5, P=P2, corner excess +1], unit squares [cost 2,
    allowance 4, P=4, corner excess -4]; maximise area (all fields area 1)."""
    # variables x3, x2, x4 >= 0; cost: 2x3 + 1.75x2 + 2x4 = 1; excess: x2 - 4x4 <= 0;
    # perimeter: (P3-4) x3 + (P2-3.5) x2 + 0 x4 <= 0
    res = linprog([-1, -1, -1], A_ub=[[0, 1, -4], [P3 - 4, P2 - 3.5, 0]], b_ub=[0, 0],
                  A_eq=[[2, 1.75, 2]], b_eq=[1], bounds=[(0, None)] * 3, method='highs')
    return -res.fun, res.x


def cairo_data():
    """Cairo tiling with 4-valent vertices at (0,0),(1,1) (period 2) and 3-valent A=(a,a-1),
    a=(3+sqrt3)/6: angles 90/120, pentagon area 1."""
    a = (3 + math.sqrt(3)) / 6
    A = np.array([a, a - 1])
    Ap = np.array([2 - a, 1 - a])
    B = np.array([1 - a, a])            # rotation of A by 90 deg about the origin
    xy = np.linalg.norm(A)
    yy = np.linalg.norm(A - Ap)
    u1, u2, u3 = -A, np.array([1, -1]) - A, Ap - A

    def ang(p, q):
        return math.degrees(math.acos(np.dot(p, q) / np.linalg.norm(p) / np.linalg.norm(q)))
    return dict(a=a, xy=xy, yy=yy, angles_Y=(ang(u1, u2), ang(u2, u3), ang(u3, u1)),
                angle_X=ang(A, B), through=2 * xy, perimeter=4 * xy + yy)
