#!/usr/bin/env python3
"""S44 audit: the (5,6)/(5,7)/(6,7) 'stability' in (m = D - 2N, s - N) coordinates.
(1) Forced part: per-weight totals.  dim ker pi_D = dim (C[x]/J)_D - dim L_D = chi_(2,2N+1)(D) - chi_(p,q)(D).  From the
    Rocha-Caridi/BGG expansion of the vacuum characters, for m = D - 2N < min(2(p+q), 4N+1):
        dim ker pi_(2N+m) = p(m-3) - p(m-2p+1) - p(m-2q+1)        (p(n) = partition number, p(n<0) = 0, p < q)
    -- N-independent, and depends on q only through 2q-1.  Checked against every audited model.
(2) Unforced part: bidegree (s - N) and transfer-order refinement; compared across models.
"""
import json, glob


def pnum(n, cache={}):
    if n < 0:
        return 0
    if n not in cache:
        a = [1] + [0] * n
        for k in range(1, n + 1):
            for j in range(k, n + 1):
                a[j] += a[j - k]
        cache[n] = a[n]
    return cache[n]


def load(p, q, D):
    fn = 'audit_orders_%d_%d_D%d_P2147483647.json' % (p, q, D)
    d = json.load(open(fn))
    N = (p - 1) * (q - 1) // 2
    rows = {}
    tot = {}
    for x in d['rows']:
        if not x['ker']:
            continue
        m = x['D'] - 2 * N
        rows[(m, x['s'] - N)] = (x['ker'], tuple(sorted((int(k), v) for k, v in x['orders'].items())))
        tot[m] = tot.get(m, 0) + x['ker']
    return N, rows, tot


def main():
    models = [(3, 4, 30), (3, 5, 30), (3, 7, 28), (3, 8, 30), (3, 10, 30), (4, 5, 28), (4, 7, 30), (5, 6, 32), (5, 7, 36), (6, 7, 43)]
    allok = True
    data = {}
    for (p, q, D) in models:
        N, rows, tot = load(p, q, D)
        data[(p, q)] = (N, rows, tot, D)
        bad = []
        for m in range(0, D - 2 * N + 1):
            if m >= min(2 * (p + q), 4 * N + 1):
                break
            pred = pnum(m - 3) - pnum(m - 2 * p + 1) - pnum(m - 2 * q + 1)
            if tot.get(m, 0) != pred:
                bad.append((m, tot.get(m, 0), pred))
        allok &= not bad
        print('(%d,%d) N=%d m<=%d: totals = p(m-3) - p(m-2p+1) - p(m-2q+1) at every m: %s %s' % (p, q, N, D - 2 * N, not bad, bad[:5]))
    print('FORCED TOTALS FORMULA HOLDS FOR ALL:', allok)

    def cmp(a, b, mmax):
        ra, rb = data[a][1], data[b][1]
        keys = sorted(k for k in set(ra) | set(rb) if k[0] <= mmax)
        diff = [(k, ra.get(k), rb.get(k)) for k in keys if ra.get(k) != rb.get(k)]
        return diff
    for a, b in [((5, 6), (5, 7)), ((5, 6), (6, 7)), ((5, 7), (6, 7)), ((3, 7), (3, 8)), ((3, 8), (3, 10)), ((4, 5), (4, 7))]:
        mmax = min(data[a][3] - 2 * data[a][0], data[b][3] - 2 * data[b][0])
        diff = cmp(a, b, mmax)
        first = diff[0][0] if diff else None
        print('%s vs %s (m <= %d): bidegree+order agreement up to first difference at (m, s-N) = %s; first diffs: %s'
              % (a, b, mmax, first, diff[:3]))


if __name__ == '__main__':
    main()


def universal_free_region():
    """Compare every pair of models on m < min(2p-1, 2p'-1, 2q-1, 2q'-1) (the common Verma-free range)."""
    models = [(3, 4, 30), (3, 5, 30), (3, 7, 28), (3, 8, 30), (3, 10, 30), (4, 5, 28), (4, 7, 30), (5, 6, 32), (5, 7, 36), (6, 7, 43)]
    data = {}
    for (p, q, D) in models:
        N, rows, tot = load(p, q, D)
        data[(p, q)] = (N, rows, D)
    ok = True
    for i, (p, q, D) in enumerate(models):
        for (p2, q2, D2) in models[i + 1:]:
            lim = min(2 * p - 1, 2 * p2 - 1, 2 * q - 1, 2 * q2 - 1, D - 2 * data[(p, q)][0] + 1, D2 - 2 * data[(p2, q2)][0] + 1)
            ra, rb = data[(p, q)][1], data[(p2, q2)][1]
            keys = [k for k in set(ra) | set(rb) if k[0] < lim]
            diff = [k for k in keys if ra.get(k) != rb.get(k)]
            if diff:
                ok = False
                print('  DIFF in common free range', (p, q), (p2, q2), 'm <', lim, sorted(diff)[:3])
    print('refined (bidegree, order) data identical for every pair of models on the common Verma-free range m < min(2p-1, 2q-1): %s' % ok)


universal_free_region()
