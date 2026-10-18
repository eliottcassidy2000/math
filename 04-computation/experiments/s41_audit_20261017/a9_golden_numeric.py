#!/usr/bin/env python3
"""S41 audit A9 (independent referee code): the golden source from first principles in floating point.

1. Points: ports P (8, literal), spindle M (7), A = P + M (ports first), C = M u P[:5]; u = (w+3-j)/(w+3+j),
   d = 1 + zeta5, Z = d + u, N = |Z|^2, s = i sqrt(4/N - 1), v = Z(-1+s)/2, in the certificate embedding
   (zeta15 = e^(2 pi i/15), j = +i sqrt11, 1 < w < 2). Checks |v| = |Z+v| = 1 and 37/15 < N < 59/20.
2. All 14.6M pairs of the 5400 points of A + uC + vC at unit distance (tolerance 1e-9), compared with the receipt's
   edge list (same address order) and with the rule axis + mixed (dq, q, q) (the classification G1 applied to e = 1).
   Counts 146 / 15 / 56 of the base recipe.
3. Float census of every difference dA + u d1 + v d2 + v^2 d3 (dA in D(A) u {0}, d1, d2 in D(C) u {0}, d3 in D(C)):
   all candidates with ||x|^2 - 1| < 1e-7, their classification, and the distance of the nearest non-unit to 1.
   A true unit has |x|^2 = 1 exactly, so a float census with error << 1e-7 cannot miss one.
Usage: python a9_golden_numeric.py
"""
import os, sys, json, time
sys.dont_write_bytecode = True
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(os.environ.get('TEMP', '/tmp'), 's41_audit_numba_cache'))
import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
REC = os.path.join(os.path.dirname(HERE), 'root_positive_transport_20261017g', 'golden_shared.json')


def uniq(lst, tol=1e-9):
    out = []
    for x in lst:
        if all(abs(x - y) > tol for y in out):
            out.append(x)
    return out


def diffset(X, tol=1e-9):
    D = sorted((a - b for a in X for b in X if abs(a - b) > tol), key=lambda x: (x.real, x.imag))
    out = []
    for x in D:
        dup = False
        for y in reversed(out):
            if x.real - y.real > tol:
                break
            if abs(x - y) < tol:
                dup = True
                break
        if not dup:
            out.append(x)
    return np.array(out)


@njit(cache=True)
def census(a, b, c, e, tol, near):
    hits = []
    best = 1e9
    for i3 in range(e.shape[0]):
        for i2 in range(c.shape[0]):
            for i1 in range(b.shape[0]):
                pr = e[i3].real + c[i2].real + b[i1].real
                pi = e[i3].imag + c[i2].imag + b[i1].imag
                for iA in range(a.shape[0]):
                    xr = pr + a[iA].real
                    xi = pi + a[iA].imag
                    val = xr * xr + xi * xi - 1.0
                    av = abs(val)
                    if av < tol:
                        hits.append((iA, i1, i2, i3, av))
                    elif av < best:
                        best = av
    return hits, best


def main():
    t0 = time.time()
    z = np.exp(2j * np.pi / 15)
    al = z ** 3
    j = 1j * np.sqrt(11)
    roots = np.roots([1, 2, -4, -4])
    wv = float([r.real for r in roots if abs(r.imag) < 1e-12 and 1 < r.real < 2][0])
    phi = (1 + al + al ** 4).real
    r5 = 2 * phi - 1
    P = [0, 1, 1 + al, 1 + al + al ** 2, -al ** 4, (3 - j * r5) / 8, 3 / 4, (3 + j * r5) / 8]
    P = [complex(p) for p in P]
    eta = 1 + z ** 5
    rho = (5 + j) / 6
    M = [0, 1, eta, 1 + eta, rho, rho * eta, rho * (1 + eta)]
    M = [complex(m) for m in M]
    A = uniq(P + [p + m for p in P for m in M])
    C = uniq(M + P[:5])
    u = (wv + 3 - j) / (wv + 3 + j)
    d = 1 + al
    Z = d + u
    N = abs(Z) ** 2
    s = 1j * np.sqrt(4 / N - 1)
    v = Z * (-1 + s) / 2
    res = dict(w=wv, N=N, absv=abs(v), absZv=abs(Z + v), sizeA=len(A), sizeC=len(C))
    print('1. |A| = %d, |C| = %d; w = %.6f; N = |Z|^2 = %.6f in (37/15, 59/20): %s; |v| = %.15f, |Z+v| = %.15f; |d|^2 = %.6f'
          % (len(A), len(C), wv, N, 37 / 15 < N < 59 / 20, abs(v), abs(Z + v), abs(d) ** 2))
    # 2. all unit-distance pairs
    R = json.load(open(REC))
    addr = [tuple(x) for x in R['source_addresses']]
    pts = np.array([A[a] + u * C[b] + v * C[c] for (a, b, c) in addr])
    n = len(pts)
    found = set()
    gap = 1e9
    for i in range(n - 1):
        dd = np.abs(np.abs(pts[i + 1:] - pts[i]) - 1)
        k = np.flatnonzero(dd < 1e-9)
        for kk in k:
            found.add((i, i + 1 + int(kk)))
        m = dd[dd >= 1e-9]
        if m.size:
            gap = min(gap, float(m.min()))
    rec_edges = set(tuple(sorted(e)) for e in R['edges'])
    # rule-based edges
    def unit_pairs(X):
        return [(a, b) for a in range(len(X)) for b in range(a + 1, len(X)) if abs(abs(X[a] - X[b]) - 1) < 1e-9]
    EA = unit_pairs(A)
    EC = unit_pairs(C)
    mixed = set()
    ordered = [(b, b2) for b, b2 in EC] + [(b2, b) for b, b2 in EC]
    for (b, b2) in ordered:
        q = C[b2] - C[b]
        for (c, c2) in ordered:
            if abs((C[c2] - C[c]) - q) > 1e-9:
                continue
            for a in range(len(A)):
                for a2 in range(len(A)):
                    if abs(A[a2] - A[a] - d * q) < 1e-9:
                        mixed.add(tuple(sorted(((a, b, c), (a2, b2, c2)))))
    idx = {x: i for i, x in enumerate(addr)}
    rule = set()
    nC = len(C)
    for a, a2 in EA:
        for b in range(nC):
            for c in range(nC):
                rule.add(tuple(sorted((idx[(a, b, c)], idx[(a2, b, c)]))))
    for b, b2 in EC:
        for a in range(len(A)):
            for c in range(nC):
                rule.add(tuple(sorted((idx[(a, b, c)], idx[(a, b2, c)]))))
                rule.add(tuple(sorted((idx[(a, c, b)], idx[(a, c, b2)]))))
    for x, y in mixed:
        rule.add(tuple(sorted((idx[x], idx[y]))))
    res.update(numeric_unit_pairs=len(found), receipt_edges=len(rec_edges), equal_to_receipt=(found == rec_edges),
               rule_edges=len(rule), rule_equals_numeric=(rule == found), EA=len(EA), EC=len(EC), mixed=len(mixed),
               nearest_nonunit_gap=gap)
    print('2. numerically unit pairs among the 5400 points: %d (receipt %d, identical %s); rule axis + mixed (dq, q, q): %d '
          '(identical %s); A edges %d, C edges %d, mixed %d; closest non-unit distance | |x-y| - 1 | = %.2e [%.0f s]'
          % (len(found), len(rec_edges), found == rec_edges, len(rule), rule == found, len(EA), len(EC), len(mixed), gap,
             time.time() - t0), flush=True)
    # 3. v^2 census
    DA = diffset(A)
    DC = diffset(C)
    a_ = np.concatenate([[0], DA]).astype(np.complex128)
    b_ = (u * np.concatenate([[0], DC])).astype(np.complex128)
    c_ = (v * np.concatenate([[0], DC])).astype(np.complex128)
    e_ = (v * v * DC).astype(np.complex128)
    print('3. |D(A)| = %d, |D(C)| = %d; candidates %d' % (len(DA), len(DC), len(a_) * len(b_) * len(c_) * len(e_)), flush=True)
    hits, best = census(a_, b_, c_, e_, 1e-6, 0.0)
    pure = [h for h in hits if h[0] == 0 and h[1] == 0 and h[2] == 0]
    other = sorted(h[4] for h in hits if not (h[0] == 0 and h[1] == 0 and h[2] == 0))
    unitDC = int(np.sum(np.abs(np.abs(DC) - 1) < 1e-9))
    hist = {t: sum(1 for h in hits if h[4] < t) for t in (1e-12, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6)}
    res.update(DA=len(DA), DC=len(DC), census_hits_1e6=len(hits), pure_axis=len(pure), pure_axis_max_dev=max(h[4] for h in pure),
               unit_differences_C=unitDC, hist=hist, smallest_non_axis_deviation=(other[0] if other else None),
               nearest_outside_window=best)
    print('   float census: cumulative counts of ||x|^2 - 1| below 1e-12, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6: %s' % list(hist.values()))
    print('   pure v^2-axis candidates %d (largest deviation %.1e); unit differences of C: %d; smallest deviation among all '
          'other candidates %.3e [%.0f s]' % (len(pure), max(h[4] for h in pure), unitDC, other[0] if other else -1,
                                             time.time() - t0))
    json.dump(res, open(os.path.join(HERE, 'a9_golden_numeric.json'), 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
