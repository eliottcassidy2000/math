"""procgen_smallgraph_20260926_t2.py -- T2: Friedman's Maximum Fence Area problem read against
the repository (planar graphs, Euler, isoperimetry, grids = Sundaram's sieve, polyomino thresholds).

Session collatz-procgen-20260922, lane "smallgraph" (2026-09-26).  Called by the runner.
The record table below was transcribed from https://erich-friedman.github.io/packing/fence/
(fetched 2026-09-26 with a generic user agent; page sha256 recorded in the note).  The
records are BEST KNOWN values (lower bounds for the true maxima), not proven optima.
"""
from fractions import Fraction as Fr
from math import sqrt, pi, isqrt, ceil
import math

from procgen_smallgraph_20260926_lib import check

A_REC = {3: 0.43301, 4: 1.0, 5: 1.0, 6: 1.47585, 7: 2.0, 8: 2.10306, 9: 2.63630, 10: 3.04687,
         11: 3.53721, 12: 4.0, 13: 4.16199, 14: 4.74494, 15: 5.06345, 16: 5.53706, 17: 6.01086,
         18: 6.24791, 19: 6.87498, 20: 7.14732, 21: 7.72489, 22: 8.07384, 23: 8.53463, 24: 9.02394,
         25: 9.30077, 26: 9.89324, 27: 10.16136, 28: 10.73270, 29: 11.08218, 30: 11.53492,
         31: 12.03575, 32: 12.43064, 33: 13.01887, 34: 13.28913, 35: 13.90296, 36: 14.18217,
         37: 14.76529, 38: 15.08152, 39: 15.53725, 40: 16.04172, 41: 16.44841, 42: 17.02275,
         43: 17.34095, 44: 17.91952, 45: 18.17998, 46: 18.82057, 47: 19.15056, 48: 19.75072,
         49: 20.07383, 50: 20.55829}
TRIVIAL = [3, 4, 5, 7, 12]   # marked "Trivial." on the page
TRUNC = 1e-5                 # the page prints values truncated to 5 decimals ("+")


def upper_bound(n):
    """A(n) <= ((sqrt(1 + 4n/sqrt(pi)) - 1)/2)^2  (isoperimetry per field + outer boundary)."""
    return ((sqrt(1 + 4 * n / sqrt(pi)) - 1) / 2) ** 2


# ---------------------------------------------------------------------------------------
# exact planar arrangement engine for configurations with rational coordinates
# ---------------------------------------------------------------------------------------
def orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def on_segment(p, a, b):
    if orient(a, b, p) != 0:
        return False
    return min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])


def seg_intersection(a, b, c, d):
    """Intersection of closed segments ab, cd: returns ('none'), ('point', P), ('overlap')."""
    d1, d2, d3, d4 = orient(c, d, a), orient(c, d, b), orient(a, b, c), orient(a, b, d)
    if d1 == 0 and d2 == 0:  # collinear
        pts = [p for p in (a, b) if on_segment(p, c, d)] + [p for p in (c, d) if on_segment(p, a, b)]
        pts = sorted(set(pts))
        if not pts:
            return ('none',)
        if len(pts) == 1:
            return ('point', pts[0])
        return ('overlap',)
    if (d1 > 0 and d2 > 0) or (d1 < 0 and d2 < 0) or (d3 > 0 and d4 > 0) or (d3 < 0 and d4 < 0):
        return ('none',)
    # proper or touching intersection: solve
    den = (b[0] - a[0]) * (d[1] - c[1]) - (b[1] - a[1]) * (d[0] - c[0])
    t = Fr((c[0] - a[0]) * (d[1] - c[1]) - (c[1] - a[1]) * (d[0] - c[0])) / den
    P = (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))
    return ('point', P)


def analyse(fences):
    """fences: list of ((x1,y1),(x2,y2)) with Fraction coordinates.
    Returns dict with rule checks, V, E, V0, T, components, bounded face areas (exact)."""
    n = len(fences)
    fences = [(tuple(map(Fr, a)), tuple(map(Fr, b))) for a, b in fences]
    res = {'n': n}
    res['unit'] = all((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2 == 1 for a, b in fences)
    cross = False
    pts_on = [set([a, b]) for a, b in fences]
    for i in range(n):
        for j in range(i + 1, n):
            r = seg_intersection(*fences[i], *fences[j])
            if r[0] == 'overlap':
                cross = True
            elif r[0] == 'point':
                P = r[1]
                interior_i = P not in (fences[i][0], fences[i][1])
                interior_j = P not in (fences[j][0], fences[j][1])
                if interior_i and interior_j:
                    cross = True
                pts_on[i].add(P)
                pts_on[j].add(P)
    res['no_cross'] = not cross
    ends_ok = True
    for i, (a, b) in enumerate(fences):
        for e in (a, b):
            if not any(on_segment(e, *fences[j]) for j in range(n) if j != i):
                ends_ok = False
    res['ends_on_fences'] = ends_ok
    # vertices and pieces
    V = set()
    for s in pts_on:
        V |= s
    V = sorted(V)
    vid = {p: k for k, p in enumerate(V)}
    pieces = set()
    for i, (a, b) in enumerate(fences):
        d = (b[0] - a[0], b[1] - a[1])
        pl = sorted(pts_on[i], key=lambda p: (p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1])
        for p, q in zip(pl, pl[1:]):
            pieces.add(tuple(sorted((vid[p], vid[q]))))
    pieces = sorted(pieces)
    # through-vertices: points interior to some fence
    through = set()
    for i, (a, b) in enumerate(fences):
        for p in pts_on[i]:
            if p != a and p != b:
                through.add(p)
    res['V'] = len(V)
    res['E'] = len(pieces)
    res['T'] = len(through)
    res['V0'] = len(V) - len(through)
    # components
    parent = list(range(len(V)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for u, v in pieces:
        parent[find(u)] = find(v)
    res['components'] = len({find(x) for x in range(len(V))})
    # faces via half-edges
    out = {k: [] for k in range(len(V))}
    for u, v in pieces:
        out[u].append(v)
        out[v].append(u)
    for u in out:
        pu = V[u]
        out[u].sort(key=lambda w: math.atan2(float(V[w][1] - pu[1]), float(V[w][0] - pu[0])))
    visited = set()
    faces = []
    for u in out:
        for v in out[u]:
            if (u, v) in visited:
                continue
            face = []
            a, b = u, v
            while (a, b) not in visited:
                visited.add((a, b))
                face.append(a)
                lst = out[b]
                k = lst.index(a)
                c = lst[(k - 1) % len(lst)]   # next edge: clockwise from the reverse edge
                a, b = b, c
            area2 = Fr(0)
            for k in range(len(face)):
                p, q = V[face[k]], V[face[(k + 1) % len(face)]]
                area2 += p[0] * q[1] - q[0] * p[1]
            faces.append(area2 / 2)
    bounded = sorted(a for a in faces if a > 0)
    res['faces_total'] = len(faces)
    res['bounded_areas'] = bounded
    res['area'] = sum(bounded)
    return res


def grid_config(i, j):
    F = []
    for y in range(j + 1):
        for x in range(i):
            F.append(((x, y), (x + 1, y)))
    for x in range(i + 1):
        for y in range(j):
            F.append(((x, y), (x, y + 1)))
    return F


def polyomino_cells(k):
    w = isqrt(k - 1) + 1 if k > 0 else 0      # ceil(sqrt k)
    cells = [(c % w, c // w) for c in range(k)]
    return cells, w, -(-k // w)


def polyomino_config(k):
    cells, w, h = polyomino_cells(k)
    edges = set()
    for (x, y) in cells:
        edges |= {((x, y), (x + 1, y)), ((x, y + 1), (x + 1, y + 1)), ((x, y), (x, y + 1)), ((x + 1, y), (x + 1, y + 1))}
    return sorted(edges)


# ---------------------------------------------------------------------------------------
def t2():
    print('== T2 fences: records, bounds, grids, polyominoes, Euler ==')
    ns = sorted(A_REC)
    check('T2.0 transcribed table', len(ns) == 48 and ns[0] == 3 and ns[-1] == 50 and
          abs(A_REC[3] - sqrt(3) / 4) < TRUNC, 'A(3..50) best-known values; A(3)=sqrt(3)/4')
    # --- upper bound
    viol = [n for n in ns if A_REC[n] > upper_bound(n)]
    ratios = {n: A_REC[n] / upper_bound(n) for n in ns}
    check('T2.A1 isoperimetric upper bound U(n)=((sqrt(1+4n/sqrt(pi))-1)/2)^2 respected by every record', not viol,
          'A_rec/U: ' + ' '.join(f'{n}:{ratios[n]:.3f}' for n in (3, 4, 6, 10, 12, 15, 20, 30, 40, 50)))
    check('T2.A2 U(n)/n -> 1/sqrt(pi)', abs(upper_bound(10 ** 8) / 10 ** 8 - 1 / sqrt(pi)) < 1e-3,
          f'1/sqrt(pi)={1 / sqrt(pi):.5f}; U(50)={upper_bound(50):.3f} vs record {A_REC[50]}; records A/n rise to {A_REC[50] / 50:.4f} at n=50')
    check('T2.A3 Euler bound A <= #fields <= n-2 respected', all(A_REC[n] <= n - 2 + 1e-12 for n in ns), 'weak but exact')
    # --- grids = Sundaram
    grid = {}
    for n in range(1, 51):
        best = None
        for i in range(1, n + 1):
            for j in range(i, n + 1):
                if i + j + 2 * i * j == n:
                    best = max(best or 0, i * j)
        grid[n] = best
    comp = {n for n in range(1, 51) if any((2 * n + 1) % d == 0 for d in range(3, isqrt(2 * n + 1) + 1, 2))}
    check('T2.B1 an i x j unit grid uses i+j+2ij fences, so a full grid exists iff 2n+1 is composite (Sundaram)',
          {n for n in grid if grid[n]} == comp,
          'grid n <= 50: ' + ' '.join(f'{n}({grid[n]})' for n in sorted(comp)))
    eq = [n for n in ns if grid.get(n) and abs(A_REC[n] - grid[n]) < TRUNC]
    check('T2.B2 records equal the best grid exactly at n in {4,7,12}', eq == [4, 7, 12],
          'elsewhere with composite 2n+1 the record beats the grid: ' +
          ' '.join(f'{n}:+{A_REC[n] - grid[n]:.3f}' for n in ns if grid.get(n) and n not in (4, 7, 12)))
    check('T2.B3 first non-grid record n=6 has 2n+1=13 prime; trivial values are 4,7,12 (grids) and 3,5',
          6 not in comp and set(TRIVIAL) - {3, 5} == {4, 7, 12}, 'first "found" record at n=6')
    # --- polyominoes (Harary-Harborth perimeter; construction checked)
    ok = True
    for k in range(1, 100001):
        w = isqrt(k - 1) + 1
        h = -(-k // w)
        c2 = isqrt(4 * k - 1) + 1          # ceil(2 sqrt k) for k>=1 (4k is never ... handled below)
        if isqrt(4 * k) ** 2 == 4 * k:
            c2 = isqrt(4 * k)
        if w + h != c2:
            ok = False
    check('T2.C1 quasi-square polyomino: ceil(sqrt k) + ceil(k/ceil(sqrt k)) = ceil(2 sqrt k)', ok,
          'k <= 1e5; so a k-omino with perimeter 2*ceil(2 sqrt k) uses 2k+ceil(2 sqrt k) unit fences')
    for k in (1, 2, 3, 5, 7, 10, 13, 20, 30):
        r = analyse(polyomino_config(k))
        c2 = ceil(2 * sqrt(k) - 1e-12)
        check(f'T2.C2 polyomino k={k} is a legal configuration', r['unit'] and r['no_cross'] and r['ends_on_fences']
              and r['n'] == 2 * k + c2 and r['area'] == k and max(r['bounded_areas']) == 1,
              f"fences={r['n']}=2k+ceil(2sqrt k), fields={len(r['bounded_areas'])} of area 1")
    Npoly = {k: 2 * k + ceil(2 * sqrt(k) - 1e-12) for k in range(1, 21)}
    Nrec = {k: min(n for n in ns if A_REC[n] >= k - 1e-12) for k in range(1, 21)}
    diff = {k: Npoly[k] - Nrec[k] for k in Npoly}
    check('T2.C3 least n with record >= k equals 2k+ceil(2 sqrt k) for k<=20 except k=13,17 (one fence fewer)',
          {k for k in diff if diff[k] != 0} == {13, 17} and diff[13] == 1 and diff[17] == 1 and all(d >= 0 for d in diff.values()),
          ' '.join(f'{k}:{Nrec[k]}' for k in Nrec) + ' (poly: 34 and 43 at k=13,17; records 33 and 42)')
    lower = {n: max([k for k in range(0, 60) if k == 0 or 2 * k + ceil(2 * sqrt(k) - 1e-12) <= n]) for n in ns}
    check('T2.C4 every record >= polyomino bound max{k: 2k+ceil(2 sqrt k) <= n}', all(A_REC[n] >= lower[n] for n in ns),
          'records exceed it except at n=4,7,12 (equal): ' + ' '.join(f'{n}:{A_REC[n] - lower[n]:.2f}' for n in (6, 10, 20, 30, 40, 50)))
    # --- Euler identity on exact configurations
    configs = {f'grid{i}x{j}': grid_config(i, j) for i in range(1, 5) for j in range(i, 5)}
    for k in (1, 3, 6, 11, 17, 26):
        configs[f'poly{k}'] = polyomino_config(k)
    configs['square+chord (n=5)'] = grid_config(1, 1) + [((0, Fr(3, 5)), (Fr(4, 5), 0))]
    configs['square+chord+chord'] = grid_config(1, 1) + [((0, Fr(3, 5)), (Fr(4, 5), 0)), ((Fr(1, 5), 1), (1, Fr(2, 5)))]
    control = grid_config(1, 1) + [((0, Fr(3, 5)), (Fr(4, 5), 0)), ((1, Fr(3, 5)), (Fr(1, 5), 0))]
    ok = True
    rows = []
    for name, F in configs.items():
        r = analyse(F)
        fb = len(r['bounded_areas'])
        good = r['unit'] and r['no_cross'] and fb == r['n'] + r['components'] - r['V0'] and r['E'] == r['n'] + r['T']
        ok = ok and good
        rows.append(f"{name}: n={r['n']} V0={r['V0']} T={r['T']} fields={fb}")
    check('T2.D1 Euler identity #fields = n + c - V0 and #pieces = n + T on exact configurations', ok, '; '.join(rows))
    r = analyse(configs['square+chord (n=5)'])
    check('T2.D2 the n=5 record configuration (unit square + chord (0,3/5)-(4/5,0))',
          r['ends_on_fences'] and r['bounded_areas'] == [Fr(6, 25), Fr(19, 25)] and r['area'] == 1,
          'fields 6/25 and 19/25; total 1 (exact)')
    r2 = analyse(control)
    check('T2.D3 control: two crossing chords are rejected', not r2['no_cross'], 'crossing detected')
    # --- superadditivity of the records
    bad = [(m, k) for m in ns for k in ns if m <= k and m + k in A_REC and A_REC[m + k] < A_REC[m] + A_REC[k] - TRUNC]
    check('T2.E records are superadditive (else a disjoint union would beat a record)', not bad,
          'no pair m+k<=50 violates A(m+k) >= A(m)+A(k)')
    # --- empirical growth
    xs = [n for n in ns if n >= 12]
    import numpy as np
    M = np.array([[sqrt(n), 1.0] for n in xs])
    y = np.array([A_REC[n] - n / 2 for n in xs])
    coef, *_ = np.linalg.lstsq(M, y, rcond=None)
    resid = float(max(abs(M @ coef - y)))
    check('T2.F EMPIRICAL fit A_rec(n) ~ n/2 - a sqrt(n) + b on 12<=n<=50', resid < 0.2,
          f'a={-coef[0]:.4f}, b={coef[1]:.4f}, max residual {resid:.3f} (polyomino: a=1/sqrt2=0.7071)')
    check('T2.G integer thresholds of the records (for the numerology table)',
          Nrec[5] == 15 and Nrec[4] == 12 and Nrec[2] == 7 and Nrec[1] == 4,
          'area >= 5 first at n=15 (same number as the square-sum threshold: NUMEROLOGY, different mechanisms)')
