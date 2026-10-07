"""The 2-adic Haar model of Mersenne switching and the debt walk (session opus-2026-10-06-S19).
Note: 05-knowledge/results/collatz_cycles_tubes_debt_walk_openai_20261006.md, section 1.
Usage:  python3 <this> [N_A M_A DMAX SEED_A]      defaults 1200 12000 41 2026 (about 2 minutes on 12 cores)
        the large run of the note:  python3 <this> 3000 20000 61 7   (about 15 minutes; output in <this>_large.out)
  A. Haar model of THM-4556 (iv): odd a <-> X = 3^(a-1) Haar on 1 + 8Z_2; post-run starts p = 2X - 1, q_D = 2X 3^(-D) - 1;
     switch at lag D iff U^i(p) = U^(i+D)(q_D).  Monte Carlo of the non-switch share q(T) and the lag-1 non-merge share q1(T),
     least-squares exponent on T >= 400 with a bootstrap interval, sqrt(T) q1(T) with a bootstrap interval, least-total lags.
  B. Exact debt bookkeeping for lag-1 pairs x0 = 3 2^v y0 + 1 (x_t against y_(t+1)), cf. mac-mini's HYP-9214 model:
     relation x = 2^L y + Delta, debt D = 3 Delta + 1 - 2^L, transition b = v_2(2^(L+a) y' + D);
     D = 0 iff the next relation is the identity (L' = 0, Delta' = 0); a merge with D != 0 needs the value coincidence
     y' = -Delta'/(2^(L') - 1) (Haar measure 0).  Kicked law (integer D, i.e. L >= 0; v_2(D) < L + a; L' > w'):
     D'_odd = U(D_odd) - 2^(L' - w'), w' = v_2(3 D_odd + 1).
  C. The L-walk: increments by L-range; the two exponent streams (debt exponent w vs partner exponent a) in generic steps.
  D. Actual exponents a <= 2001: sigma(M_a) census vs the Haar prediction 1 - mean q(T_post(a)).
Exact claims print ok/FAIL; Monte Carlo numbers are NUMERICAL (seeded)."""
import sys, random, math, time
from fractions import Fraction
from multiprocessing import Pool
from collections import Counter

FAIL = []


def check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg)
    if not cond:
        FAIL.append(msg)


def v2(n):
    return (n & -n).bit_length() - 1


def U(x):
    y = 3 * x + 1
    return y >> v2(y)


# ------------------------------------------------------------------ A. Haar model
_args = [int(a) for a in sys.argv[1:5]]
N_A, M_A, DMAX, SEED_A = (_args + [1200, 12000, 41, 2026][len(_args):])[:4]
MARGIN = 64
TGRID = [27, 50, 100, 200, 400, 800, 1600, 3200, 6400, 9600, 12800, 19000]


def orbit(x, prec):
    vals, precs, tots = [x], [prec], [0]
    tot = 0
    while True:
        y = 3 * x + 1
        if y & ((1 << prec) - 1) == 0:
            break
        e = v2(y)
        prec -= e
        if prec < MARGIN:
            break
        x = (y >> e) & ((1 << prec) - 1)
        tot += e
        vals.append(x); precs.append(prec); tots.append(tot)
    return vals, precs, tots


def first_merge(P, Q, D):
    pv, pp, pt = P
    qv, qp, qt = Q
    for i in range(len(pv)):
        j = i + D
        if j >= len(qv):
            return None
        m = min(pp[i], qp[j])
        if ((pv[i] - qv[j]) & ((1 << m) - 1)) == 0:
            return pt[i]
    return None


def haar_sample(s):
    rng = random.Random(SEED_A * 1000003 + s)
    mod = 1 << M_A
    X = (rng.getrandbits(M_A) & ~7) | 1
    P = orbit((2 * X - 1) & (mod - 1), M_A)
    inv3 = pow(3, -1, mod)
    res = {}
    for D in range(1, DMAX + 1, 2):
        Q = orbit((2 * X * pow(inv3, D, mod) - 1) & (mod - 1), M_A)
        res[D] = first_merge(P, Q, D)
    return res, P[2][-1]


# ------------------------------------------------------------------ B/C. lag-1 debt pairs
BITS_B, N_B, SEED_B = 12000, 1000, 4242


def debt_pair(s, exact_steps=60):
    rng = random.Random(SEED_B * 7919 + s)
    y = rng.getrandbits(BITS_B) | 1
    v = 2
    while rng.random() < 0.5:
        v += 1
    x = 3 * (1 << v) * y + 1
    t3 = 3 * y + 1
    a0 = v2(t3)
    yy = t3 >> a0
    L = v + a0
    stats = {'trans': True, 'kick_tested': 0, 'kick_ok': True, 'merge_ok': None}
    incs, Ls, streams, ab = [], [L], [], []
    merged = None
    maxsteps = int(BITS_B / 2.5)
    prev = None
    for t in range(maxsteps):
        if x == yy:
            merged = t
            Lp, Dp, Deltap = prev
            stats['merge_ok'] = (Dp == 0 and Lp % 2 == 0 and Lp != 0 and Deltap == (Fraction(2) ** Lp - 1) / 3)
            break
        Delta = Fraction(x) - Fraction(2) ** L * yy
        D = 3 * Delta + 1 - Fraction(2) ** L
        prev = (L, D, Delta)
        bx, ay = v2(3 * x + 1), v2(3 * yy + 1)
        ny = (3 * yy + 1) >> ay
        if D != 0 and D.denominator == 1 and L >= 8:
            w = v2(abs(D.numerator))
            if w < L + ay:
                streams.append((ay, w, bx))           # generic step: b = w
        if t < exact_steps:
            Fr = Fraction(2) ** (L + ay) * ny + D       # transition law: b = v_2(2^(L+a) y' + D)
            vb = v2(Fr.numerator) - v2(Fr.denominator)
            if vb != bx:
                stats['trans'] = False
            if D != 0 and D.denominator == 1:          # kicked law (integer debt)
                vD = v2(abs(D.numerator))
                Lnew = L + ay - bx
                if vD < L + ay:
                    Dodd = D.numerator >> vD if D.numerator > 0 else -((-D.numerator) >> vD)
                    w = v2(3 * Dodd + 1)
                    if Lnew > w:
                        stats['kick_tested'] += 1
                        Dnext = 3 * (Fraction((3 * x + 1) >> bx) - Fraction(2) ** Lnew * ny) + 1 - Fraction(2) ** Lnew
                        vn = v2(Dnext.numerator) - v2(Dnext.denominator)
                        if not ((vn == w) and (Dnext / Fraction(2) ** vn == U(Dodd) - 2 ** (Lnew - w))):
                            stats['kick_ok'] = False
        x = (3 * x + 1) >> bx
        yy = ny
        L = L + ay - bx
        incs.append((Ls[-1], ay - bx))
        ab.append((ay, bx))
        Ls.append(L)
    return merged, stats, incs, Ls, streams, ab


# ------------------------------------------------------------------ D. actual exponents
def sigma_and_post_total(a):
    """odd-step count of M_a = 2^a - 1 to 1, and the exponent total of its post-run orbit."""
    x = 2 * 3 ** (a - 1) - 1         # post-run start, reached after a - 1 odd steps
    s, tot = a - 1, 0
    while x != 1:
        y = 3 * x + 1
        e = v2(y)
        x = y >> e
        s += 1
        tot += e
    return s, tot


def ls_fit(points):
    n = len(points)
    mx = sum(p[0] for p in points) / n; my = sum(p[1] for p in points) / n
    a = -sum((p[0] - mx) * (p[1] - my) for p in points) / sum((p[0] - mx) ** 2 for p in points)
    return a, math.exp(my + a * mx)


def shares(least, least1, grid):
    n = len(least)
    q = [sum(1 for v in least if v is None or v > T) / n for T in grid]
    q1 = [sum(1 for v in least1 if v is None or v > T) / n for T in grid]
    return q, q1


if __name__ == '__main__':
    t0 = time.time()
    print('A. Haar model of Mersenne switching (NUMERICAL, seed %d, N = %d, %d-bit 2-adic precision, odd D <= %d)' % (SEED_A, N_A, M_A, DMAX))
    with Pool(12) as pool:
        out = pool.map(haar_sample, range(N_A), chunksize=4)
        pairs = pool.map(debt_pair, range(N_B), chunksize=8)
        sig = pool.map(sigma_and_post_total, range(2, 2002), chunksize=16)
    reach = min(o[1] for o in out)
    grid = [T for T in TGRID if T <= reach]
    least = [min([v for v in r.values() if v is not None], default=None) for r, _ in out]
    least1 = [r[1] for r, _ in out]
    q, q1 = shares(least, least1, grid)
    for T, a, b in zip(grid, q, q1):
        print(f'     T = {T:5d}: non-switch q(T) = {a:.4f}   lag-1 non-merge q1(T) = {b:.4f}   sqrt(T) q1 = {math.sqrt(T) * b:6.2f}')
    def median_T(qs):
        for (T1, v1), (T2, v2) in zip(zip(grid, qs), list(zip(grid, qs))[1:]):
            if v1 >= 0.5 >= v2:
                return math.exp(math.log(T1) + (math.log(T2) - math.log(T1)) * (v1 - 0.5) / (v1 - v2))
        return None
    print(f'     smallest reachable total over samples: {reach}; median merge total (log-interpolated): any lag {median_T(q):.0f}, lag 1 {median_T(q1):.0f}')
    s27 = 1 - q[0]
    sd = math.sqrt(s27 * (1 - s27) / N_A)
    check(abs(s27 - 0.1556) < 3 * sd + 0.005, f'Haar switch share at T = 27 is {s27:.4f} +- {sd:.4f}; certified exact share at K = 27 is 0.1556 (consistent)')
    fitT = [T for T in grid if T >= 400]
    alpha, Cq = ls_fit([(math.log(T), math.log(v)) for T, v in zip(grid, q) if T >= 400])
    tailT = [T for T in grid if T >= 3200]
    rng = random.Random(1)
    boots, boots_c1, boots_tail = [], [], []
    for _ in range(400):
        idx = [rng.randrange(N_A) for _ in range(N_A)]
        bq, bq1 = shares([least[i] for i in idx], [least1[i] for i in idx], grid)
        if min(v for T, v in zip(grid, bq) if T >= 400) > 0:
            boots.append(ls_fit([(math.log(T), math.log(v)) for T, v in zip(grid, bq) if T >= 400])[0])
        boots_c1.append(math.sqrt(grid[-1]) * bq1[-1])
        if len(tailT) >= 2:
            boots_tail.append(ls_fit([(math.log(T), math.log(v)) for T, v in zip(grid, bq1) if T >= 3200])[0])
    boots.sort(); boots_c1.sort(); boots_tail.sort()
    ci = lambda L: (L[int(0.025 * len(L))], L[int(0.975 * len(L)) - 1])
    print(f'     LS fit of q(T) on T in [{fitT[0]}, {fitT[-1]}]: q ~ {Cq:.2f} T^-{alpha:.3f}, bootstrap 95% [{ci(boots)[0]:.3f}, {ci(boots)[1]:.3f}]')
    print(f'     sqrt(T) q1(T) at T = {grid[-1]}: {math.sqrt(grid[-1]) * q1[-1]:.2f}, bootstrap 95% [{ci(boots_c1)[0]:.2f}, {ci(boots_c1)[1]:.2f}]')
    if boots_tail:
        a1t, _ = ls_fit([(math.log(T), math.log(v)) for T, v in zip(grid, q1) if T >= 3200])
        print(f'     lag-1 tail exponent on T >= 3200: {a1t:.3f}, bootstrap 95% [{ci(boots_tail)[0]:.3f}, {ci(boots_tail)[1]:.3f}]')
    a1all, _ = ls_fit([(math.log(T), math.log(v)) for T, v in zip(grid, q1) if T >= 400])
    print(f'     lag-1 exponent on T >= 400: {a1all:.3f}')
    lagc = Counter()
    for r, _ in out:
        best = None
        for D, v in r.items():
            if v is not None and (best is None or v < best[1]):
                best = (D, v)
        if best:
            lagc[best[0]] += 1
    tot = sum(lagc.values())
    print('     least-total lag among switching samples: ' + ', '.join(f'D={D}: {c / tot:.3f}' for D, c in sorted(lagc.items())[:6]))

    print('B. Exact debt bookkeeping for lag-1 pairs x0 = 3 2^v y0 + 1')
    check(all(p[1]['trans'] for p in pairs), f'transition law b_t = v_2(2^(L+a) y\' + D) exact on the first 60 steps of all {N_B} pairs')
    kt = sum(p[1]['kick_tested'] for p in pairs)
    check(all(p[1]['kick_ok'] for p in pairs) and kt > 1000,
          f'kicked-Collatz law D\'_odd = U(D_odd) - 2^(L\'-w\'), v_2(D\') = w\' = v_2(3 D_odd + 1), in all {kt} generic integer-debt steps tested')
    mg = [p for p in pairs if p[0] is not None]
    check(all(p[1]['merge_ok'] for p in mg),
          f'all {len(mg)} sampled merges pass through D = 0 one step earlier (L even, nonzero, Delta = (2^L - 1)/3: a sibling pair); '
          'value-coincidence merges with D != 0 have Haar measure 0')

    print('C. The L-walk and the two exponent streams (NUMERICAL)')
    by = {}
    for p in pairs:
        for Lp, d in p[2]:
            key = '<= 0' if Lp <= 0 else ('1-6' if Lp <= 6 else ('7-15' if Lp <= 15 else '>= 16'))
            by.setdefault(key, []).append(d)
    for k in ['<= 0', '1-6', '7-15', '>= 16']:
        d = by[k]; m = sum(d) / len(d); var = sum((x - m) ** 2 for x in d) / len(d)
        print(f'     L {k:5s}: {len(d):8d} steps, mean increment {m:+.4f} (s.e. {math.sqrt(var / len(d)):.4f}), variance {var:.3f}')
    allinc = [d for v in by.values() for d in v]
    m = sum(allinc) / len(allinc); var = sum((x - m) ** 2 for x in allinc) / len(allinc)
    check(abs(m) < 0.02 and abs(var - 4) < 0.1, f'all increments: mean {m:+.4f}, variance {var:.3f} (Var of a difference of two Geom(1/2) = 4)')
    st = [s for p in pairs for s in p[4]]
    nst = len(st)
    W = Counter(s[1] for s in st); Aa = Counter(s[0] for s in st)
    print(f'     generic integer-debt steps with L >= 8: {nst}; P(w = k), P(a = k), 2^-k for k = 1..5:')
    for k in range(1, 6):
        print(f'        k = {k}: {W[k] / nst:.4f}  {Aa[k] / nst:.4f}  {2 ** -k:.4f}')
    J = Counter((min(s[0], 4), min(s[1], 4)) for s in st)
    mA = Counter(min(s[0], 4) for s in st); mW = Counter(min(s[1], 4) for s in st)
    chi = sum((J[(i, j)] - mA[i] * mW[j] / nst) ** 2 / (mA[i] * mW[j] / nst) for i in mA for j in mW)
    ws = [s[1] for s in st]
    mw = sum(ws) / len(ws); vw = sum((x - mw) ** 2 for x in ws) / len(ws)
    ac = sum((ws[i] - mw) * (ws[i + 1] - mw) for i in range(len(ws) - 1)) / (len(ws) - 1) / vw
    print(f'     independence of a and w: chi^2 = {chi:.1f} on 9 dof; mean w = {mw:.4f}, var w = {vw:.3f}, lag-1 autocorrelation of w = {ac:+.4f}')
    print('     (these are consistency checks of exact Haar facts: in generic steps w = b is the x-orbit\'s own exponent, and x_0 = 3 2^v y_0 + 1')
    print('      is Haar on a coset, so its exponent stream is i.i.d. Geom(1/2); same-step independence of a and w is also exact)')
    check(all(s[2] == s[1] for s in st), 'b = w (= v_2(D)) in every generic step, as Theorem D(c) says')
    # the open question: long-lag joint law of the two exponent streams. L_t - L_0 = (sum of y-exponents) - (sum of x-exponents)
    print('     block variances Var(L_(t+K) - L_t)/K over disjoint blocks starting at L >= 40 (i.i.d. increments give 4 for every K);')
    print('     in brackets: blocks required to stay at L >= 16 throughout (a conditioning that biases the variance down for large K):')
    bv = []
    for K in (1, 4, 16, 64, 256):
        diffs, diffs_c = [], []
        for p in pairs:
            Ls = p[3]
            for t0b in range(0, len(Ls) - K, K):
                if Ls[t0b] >= 40:
                    diffs.append(Ls[t0b + K] - Ls[t0b])
                if Ls[t0b] >= 16 and min(Ls[t0b:t0b + K + 1]) >= 16:
                    diffs_c.append(Ls[t0b + K] - Ls[t0b])
        if len(diffs) > 30:
            mdf = sum(diffs) / len(diffs)
            vK = sum((d - mdf) ** 2 for d in diffs) / len(diffs) / K
            se = vK * math.sqrt(2 / len(diffs))
            mdc = sum(diffs_c) / len(diffs_c)
            vc = sum((d - mdc) ** 2 for d in diffs_c) / len(diffs_c) / K
            bv.append((K, vK, se, len(diffs)))
            print(f'        K = {K:3d}: {vK:.3f} (+- {se:.3f}, {len(diffs)} blocks)   [{vc:.3f}, {len(diffs_c)} blocks]')
    # lagged cross-covariances of the two streams inside the L >= 16 regime: Cov(a_s, b_(s+k)) and Cov(b_s, a_(s+k))
    def cov(u, v):
        mu, mv = sum(u) / len(u), sum(v) / len(v)
        return sum((x - mu) * (y - mv) for x, y in zip(u, v)) / len(u)
    lagc_out = []
    for k in (1, 2, 4, 8, 16, 32, 64):
        A1, B1, A2, B2 = [], [], [], []
        for p in pairs:
            Ls, ab = p[3], p[5]
            for s in range(len(ab) - k):
                if Ls[s] >= 16 and Ls[s + k] >= 16:
                    A1.append(ab[s][0]); B1.append(ab[s + k][1])
                    A2.append(ab[s + k][0]); B2.append(ab[s][1])
        if len(A1) > 100:
            lagc_out.append((k, cov(A1, B1), cov(B2, A2), 2 / math.sqrt(len(A1))))
    print('     cross-covariances in the L >= 16 regime, Cov(a_s, b_(s+k)) / Cov(b_s, a_(s+k)) (s.e. ~ 2/sqrt(n)):')
    print('        ' + '; '.join(f'k={k}: {c1:+.4f}/{c2:+.4f} (s.e. {se:.4f})' for k, c1, c2, se in lagc_out))
    check(all(abs(v - 4) < 4 * se + 0.15 for K, v, se, n in bv), 'block variances consistent with 4 for K up to 256 (NUMERICAL)')
    print('     L one step before merge:', Counter(p[3][p[0] - 1] for p in mg).most_common(8))
    for T in [100, 300, 1000, 3000, 4800]:
        nm = sum(1 for p in pairs if p[0] is None or p[0] > T) / N_B
        print(f'     P(no merge by step {T:5d}) = {nm:.4f}   sqrt(T) P = {math.sqrt(T) * nm:5.2f}')

    print('D. Actual exponents: sigma(M_a) census for a <= 2001 vs the Haar prediction')
    sig_of = {a: s for a, (s, _) in zip(range(2, 2002), sig)}
    post = {a: tot for a, (_, tot) in zip(range(2, 2002), sig)}
    seen = {0}                        # sigma(M_1) = sigma(1) = 0 (the convention of THM-4556 (vi) counts a >= 1)
    levels, new_level = [], {}
    for a in range(2, 2002):
        new_level[a] = sig_of[a] not in seen
        seen.add(sig_of[a])
        levels.append(len(seen))
    odd = [a for a in range(1001, 2002, 2)]
    share = sum(1 for a in odd if not new_level[a]) / len(odd)
    lag1 = sum(1 for a in odd if sig_of[a] == sig_of[a - 1]) / len(odd)
    ratio = sum(post[a] / a for a in odd) / len(odd)
    pred = 1 - sum(min(1.0, Cq * post[a] ** -alpha) for a in odd) / len(odd)
    c1h = math.sqrt(grid[-1]) * q1[-1]
    pred1 = 1 - sum(min(1.0, c1h / math.sqrt(post[a])) for a in odd) / len(odd)
    print(f'     odd a in [1001, 2001]: share with a smaller sigma-partner {share:.3f} (Haar prediction {pred:.3f}); '
          f'share with sigma(M_(a-1)) {lag1:.3f} (Haar prediction {pred1:.3f} with sqrt(T) q1 = {c1h:.2f})')
    print(f'     post-run template total / a = {ratio:.3f};  #distinct sigma values for a <= 100, 400, 1000, 2001: '
          f'{levels[98]}, {levels[398]}, {levels[998]}, {levels[-1]}')
    predlev = {A: 2 + sum(min(1.0, Cq * post[a] ** -alpha) for a in range(3, A + 1, 2)) for A in (100, 400, 1000, 2001)}
    print('     Haar-predicted level counts:', {A: round(v, 1) for A, v in predlev.items()})
    obs = {100: levels[98], 400: levels[398], 1000: levels[998], 2001: levels[-1]}
    def logslope(d):
        pts = [(math.log(A), math.log(v)) for A, v in d.items()]
        n_ = len(pts); mx_ = sum(u for u, _ in pts) / n_; my_ = sum(w for _, w in pts) / n_
        return sum((u - mx_) * (w - my_) for u, w in pts) / sum((u - mx_) ** 2 for u, _ in pts)
    print(f'     log-log slope of level counts over A = 100..2001: observed {logslope(obs):.3f}, predicted {logslope(predlev):.3f} (fit of this run)')
    check(abs(share - pred) < 0.03, 'Haar model predicts the switching share of actual exponents within 0.03 (NUMERICAL agreement)')
    check(levels[98] == 23 and levels[398] == 37, 'sigma(M_a) takes 23 and 37 values for a <= 100, 400 (THM-4556 (vi), HYP-9213)')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
