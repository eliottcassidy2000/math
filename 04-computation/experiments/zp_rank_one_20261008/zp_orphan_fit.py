"""Repunit orphans: exponent fits with intervals, and a valid refutation rule (session opus-2026-10-08-S22; added after the audit).
Note: 05-knowledge/results/zp_rank_one_coalescence_20261008.md (section 4); HYP-9245 (b).
Reads the committed data repunit_orbits_p{7,11,13}_K8000.txt of zp_repunit_lines.py (columns K n n_mult entry cycle_min;
n = steps from R_K = p^K - 1 to the first value on a cycle, n_mult = multiplication steps among them, entry = that value).

  (O1) FINITE-EXACT: the pair (height n - K, entry) determines n_mult and the cycle for every K <= 8000, so partner classes
       (equal height and entry, i.e. C^(n(K') + K - K')(R_K) = C^(n(K'))(R_K'): an exit R_K ~> R_K' completed by cycle entry)
       are the classes of zp_repunit_lines.py's key.  Orphans = least members of their class.
  (O2) NUMERICAL: orphan counts per dyadic window with Poisson errors, times sqrt(K); local exponents with standard errors.
  (O3) NUMERICAL: Poisson maximum-likelihood fit of the orphan rate c K^(-alpha) on K in [800, 8000], treating orphan events
       as independent, with profile-likelihood 95% and 99% intervals.
  (O4) NUMERICAL (simulation under the null of independent Bernoulli(c K^(-1/2)) orphans):
       * the original HYP-9245 (b) rule (refute if both local exponents of the last two dyadic windows ending at K = 30000
         fall outside [0.35, 0.65]) rejects a true K^(-1/2) law about 40% of the time: invalid;
       * the replacement rule (refute if the 99% profile-likelihood interval for alpha on K in [10^4, 10^5] excludes 1/2)
         rejects a true law at about the nominal 1%, and has high power against alpha = 1 and moderate power against 1/4.
Prints ALL CHECKS PASSED.  Runtime about 10 seconds."""
import os, math, time
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = []
EDGES = [100, 200, 400, 800, 1600, 3200, 6400, 8001]


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def load(p):
    rows = {}
    with open(os.path.join(HERE, f'repunit_orbits_p{p}_K8000.txt')) as f:
        for line in f:
            K, n, nm, x, c = map(int, line.split())
            rows[K] = (n, nm, x, c)
    return rows


def profile_fit(Ks, hit, grid):
    """Poisson rate c K^(-a): with c profiled out, ll(a) = n ln(n/Z(a)) - a S - n, Z(a) = sum K^(-a), S = sum_hits ln K."""
    lk = np.log(Ks)
    Z = np.array([np.exp(-a * lk).sum() for a in grid])
    n = hit.sum(); S = lk[hit].sum()
    ll = n * np.log(n / Z) - grid * S - n
    return ll, n, Z


def interval(ll, grid, drop):
    ok = grid[ll >= ll.max() - drop]
    return float(ok.min()), float(ok.max())


if __name__ == '__main__':
    t0 = time.time()
    print('(O1) partner classes')
    data = {p: load(p) for p in (7, 11, 13)}
    ok1 = True
    bottoms = {}
    for p, rows in data.items():
        full = defaultdict(set); cls = defaultdict(list)
        for K, (n, nm, x, c) in rows.items():
            full[(n - K, x)].add((nm, c)); cls[(n - K, x)].append(K)
        ok1 &= all(len(v) == 1 for v in full.values())
        bottoms[p] = {min(v) for v in cls.values()}
        print(f'     p = {p}: {len(cls)} classes; (height, entry) determines (n_mult, cycle): {all(len(v) == 1 for v in full.values())}')
    check(ok1, 'for p = 7, 11, 13 and K <= 8000 the height n - K and the cycle entry determine n_mult and the cycle (FINITE-EXACT)')

    print('(O2) orphan counts per window (fraction +- Poisson s.e.; x sqrt(K_mid)) and local exponents')
    for p in data:
        pts = []; line = []
        for lo, hi in zip(EDGES, EDGES[1:]):
            cnt = sum(1 for K in range(lo, hi) if K in bottoms[p])
            fr = cnt / (hi - lo); km = math.sqrt(lo * hi)
            se = fr / math.sqrt(cnt) if cnt else float('nan')
            pts.append((km, fr, cnt))
            line.append(f'[{lo},{hi}): {cnt} ({fr * math.sqrt(km):.3f} +- {se * math.sqrt(km):.3f})')
        ex = []
        for (k1, o1, c1), (k2, o2, c2) in zip(pts, pts[1:]):
            if c1 and c2:
                ex.append(f'{math.log(o1 / o2) / math.log(k2 / k1):.2f} +- {math.sqrt(1 / c1 + 1 / c2) / math.log(k2 / k1):.2f}')
        print(f'     p = {p}: ' + '; '.join(line))
        print(f'            local exponents: {ex}')

    print('(O3) Poisson ML fit of c K^(-alpha) on K in [800, 8000]')
    grid = np.round(np.arange(0, 2.0001, 0.001), 3)
    Ks = np.arange(800, 8001, dtype=float)
    fits = {}
    for p in data:
        hit = np.array([int(K) in bottoms[p] for K in Ks])
        ll, n, Z = profile_fit(Ks, hit, grid)
        a = float(grid[ll.argmax()]); c = n / Z[ll.argmax()]
        i95 = interval(ll, grid, 1.92); i99 = interval(ll, grid, 3.32)
        fits[p] = (a, i95, i99)
        print(f'     p = {p}: {int(n)} orphans; alpha = {a:.2f}, 95% [{i95[0]:.2f}, {i95[1]:.2f}], 99% [{i99[0]:.2f}, {i99[1]:.2f}]; c = {c:.3f}')
    check(all(f[1][0] < 0.5 < f[1][1] and f[1][1] - f[1][0] > 0.5 for f in fits.values()),
          'every 95% interval contains 1/2 and is wider than 0.5: K^(-1/2) is consistent with all three bases and resolved for none (NUMERICAL)')

    print('(O4) refutation rules under the null (independent Bernoulli(c K^(-1/2)) orphans)')
    rng = np.random.default_rng(20261008)
    trials = 4000
    old_rates = []
    for c in (0.56, 0.65):
        edges = [3750, 7500, 15000, 30001]
        rej = 0
        for _ in range(trials):
            fr = []
            for lo, hi in zip(edges, edges[1:]):
                k = np.arange(lo, hi, dtype=float)
                fr.append((rng.random(k.size) < c / np.sqrt(k)).sum() / (hi - lo))
            mids = [math.sqrt(lo * hi) for lo, hi in zip(edges, edges[1:])]
            ex = [math.log(fr[i] / fr[i + 1]) / math.log(mids[i + 1] / mids[i]) if fr[i] and fr[i + 1] else float('nan')
                  for i in range(2)]
            rej += all(not (0.35 <= x <= 0.65) for x in ex)
        old_rates.append(rej / trials)
        print(f'     original rule, c = {c}: a true K^(-1/2) law is declared refuted in {rej / trials:.3f} of {trials} trials')
    Kw = np.arange(10 ** 4, 10 ** 5 + 1, dtype=float)
    lk = np.log(Kw)
    g2 = np.round(np.arange(0, 1.6001, 0.0025), 4)
    Z = np.array([np.exp(-a * lk).sum() for a in g2])
    kmid = math.sqrt(1e9)
    new = {}
    for label, alpha in (('null alpha = 1/2', 0.5), ('alpha = 1', 1.0), ('alpha = 1/4', 0.25)):
        c = 0.6 * kmid ** (alpha - 0.5)           # same rate as 0.6 K^(-1/2) at the window's geometric middle
        rej = 0
        for _ in range(trials // 2):
            hit = rng.random(Kw.size) < c * Kw ** (-alpha)
            n = hit.sum(); S = lk[hit].sum()
            ll = n * np.log(n / Z) - g2 * S - n
            lo_, hi_ = interval(ll, g2, 3.32)
            rej += not (lo_ <= 0.5 <= hi_)
        new[label] = rej / (trials // 2)
        print(f'     replacement rule (99% interval on [10^4, 10^5] excludes 1/2), {label}: rejected in {rej / (trials // 2):.3f}')
    check(min(old_rates) > 0.3 and new['null alpha = 1/2'] < 0.025 and new['alpha = 1'] > 0.9,
          'the original rule rejects a true law in more than 30% of trials (invalid); the replacement rejects it in under 2.5% '
          '(nominal 1%) and rejects alpha = 1 in over 90% (NUMERICAL)')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
