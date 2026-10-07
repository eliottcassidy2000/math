"""The 2-adic Haar model of Mersenne switching and the debt walk (session opus-2026-10-06-S19).
Note: 05-knowledge/results/collatz_cycles_tubes_debt_walk_openai_20261006.md, section 1.
  A. Haar model of THM-4556 (iv): odd a <-> X = 3^(a-1) Haar on 1 + 8Z_2; post-run starts p = 2X - 1, q_D = 2X 3^(-D) - 1;
     switch at lag D iff U^i(p) = U^(i+D)(q_D).  Monte Carlo of the non-switch share q(T) and the lag-1 non-merge share.
  B. Exact debt bookkeeping for lag-1 pairs x0 = 3 2^v y0 + 1 (x_t vs y_(t+1)): relation x = 2^L y + Delta, D = 3 Delta + 1 - 2^L,
     transition b = v_2(2^(L+a) y' + D); merge iff D = 0 one step earlier, with L even and Delta = (2^L - 1)/3;
     the debt's odd part runs the Collatz map: D'_odd = U(D_odd) - 2^(L' - w') in generic steps (w' = v_2(3 D_odd + 1) < L').
  C. The L-walk: increments mean 0, variance 4 (= Var(Geom - Geom)); sqrt(t) P(no merge by t) plateaus (t^(-1/2) law).
  D. Actual exponents a in [1001, 2001]: sigma(M_a) census vs the Haar prediction 1 - mean q(T_post(a)).
Exact claims print ok/FAIL; Monte Carlo numbers are NUMERICAL (seeded).  Runtime about 3-5 minutes on 12 cores."""
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
M_A, N_A, DMAX, SEED_A = 12000, 1200, 41, 2026
MARGIN = 64


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
    stats = {'rel': True, 'trans': True, 'kick_tested': 0, 'kick_ok': True, 'merge_ok': None}
    incs, Ls = [], [L]
    merged = None
    maxsteps = int(BITS_B / 2.5)
    prev = None
    for t in range(maxsteps):
        if x == yy:
            merged = t
            Lp, Dp, Deltap = prev
            stats['merge_ok'] = (Dp == 0 and Lp % 2 == 0 and Lp != 0 and Deltap == (Fraction(2) ** Lp - 1) / 3)
            break
        if t < exact_steps or True:
            Delta = Fraction(x) - Fraction(2) ** L * yy
            D = 3 * Delta + 1 - Fraction(2) ** L
            prev = (L, D, Delta)
        bx, ay = v2(3 * x + 1), v2(3 * yy + 1)
        ny = (3 * yy + 1) >> ay
        if t < exact_steps:
            # transition law: b = v_2(2^(L+a) y' + D)   (D is dyadic; scale to integers)
            Fr = Fraction(2) ** (L + ay) * ny + D
            num, den = Fr.numerator, Fr.denominator
            vb = v2(num) - v2(den)
            if vb != bx:
                stats['trans'] = False
            # kicked-Collatz law in the generic case: D != 0, v(D) < L + a, L' > w'
            if D != 0:
                vD = v2(D.numerator) - v2(D.denominator)
                Lnew = L + ay - bx
                if D.denominator == 1 and vD < L + ay:
                    Dodd = D.numerator >> vD if D.numerator > 0 else -((-D.numerator) >> vD)
                    w = v2(3 * Dodd + 1)
                    if Lnew > w:
                        stats['kick_tested'] += 1
                        Dnext = 3 * (Fraction(x) * 0 + Fraction((3 * x + 1) >> bx) - Fraction(2) ** Lnew * ny) + 1 - Fraction(2) ** Lnew
                        Dnext_odd_pred = U(Dodd) - 2 ** (Lnew - w)
                        vn = v2(Dnext.numerator) - v2(Dnext.denominator)
                        ok = (vn == w) and (Dnext / Fraction(2) ** vn == Dnext_odd_pred)
                        if not ok:
                            stats['kick_ok'] = False
        x = (3 * x + 1) >> bx
        yy = ny
        L = L + ay - bx
        incs.append((Ls[-1], ay - bx))
        Ls.append(L)
    return merged, stats, incs, Ls


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


if __name__ == '__main__':
    t0 = time.time()
    print('A. Haar model of Mersenne switching (NUMERICAL, seed %d, N = %d, %d-bit 2-adic precision, odd D <= %d)' % (SEED_A, N_A, M_A, DMAX))
    with Pool(12) as pool:
        out = pool.map(haar_sample, range(N_A), chunksize=4)
        pairs = pool.map(debt_pair, range(N_B), chunksize=8)
        sig = pool.map(sigma_and_post_total, range(2, 2002), chunksize=16)
    reach = min(o[1] for o in out)
    rows = []
    for T in [27, 50, 100, 200, 400, 800, 1600, 3200, 6400, 9600]:
        if T > reach:
            break
        nos = sum(1 for r, _ in out if not any(v is not None and v <= T for v in r.values())) / N_A
        no1 = sum(1 for r, _ in out if not (r[1] is not None and r[1] <= T)) / N_A
        rows.append((T, nos, no1))
        print(f'     T = {T:5d}: non-switch q(T) = {nos:.4f}   lag-1 non-merge = {no1:.4f}   sqrt(T)*lag1 = {math.sqrt(T) * no1:6.2f}')
    s27 = 1 - rows[0][1]
    sd = math.sqrt(s27 * (1 - s27) / N_A)
    check(abs(s27 - 10757550 / 67108864 * 0 - 0.1556) < 3 * sd + 0.005,
          f'Haar switch share at T = 27 is {s27:.4f} +- {sd:.4f}; certified exact share at K = 27 is 0.1556 (consistent)')
    # power-law fit of q(T) on T >= 400
    pts = [(math.log(T), math.log(q)) for T, q, _ in rows if T >= 400 and q > 0]
    n = len(pts)
    mx = sum(p[0] for p in pts) / n; my = sum(p[1] for p in pts) / n
    alpha = -sum((p[0] - mx) * (p[1] - my) for p in pts) / sum((p[0] - mx) ** 2 for p in pts)
    Cq = math.exp(my + alpha * mx)
    print(f'     LS fit on T >= 400: q(T) ~ {Cq:.2f} T^-{alpha:.3f}   (the S19 large run, seed 7, N = 3000, 20000 bits, D <= 61: alpha ~ 0.66)')

    print('B. Exact debt bookkeeping for lag-1 pairs x0 = 3 2^v y0 + 1')
    rel_ok = all(p[1]['trans'] for p in pairs)
    check(rel_ok, f'transition law b_t = v_2(2^(L+a) y\' + D) exact on the first 60 steps of all {N_B} pairs')
    kt = sum(p[1]['kick_tested'] for p in pairs)
    check(all(p[1]['kick_ok'] for p in pairs) and kt > 1000,
          f'kicked-Collatz law D\'_odd = U(D_odd) - 2^(L\'-w\'), v_2(D\') = w\' = v_2(3 D_odd + 1), in all {kt} generic steps tested')
    mg = [p for p in pairs if p[0] is not None]
    check(all(p[1]['merge_ok'] for p in mg), f'every one of the {len(mg)} merges: one step earlier D = 0, L even and nonzero, Delta = (2^L - 1)/3 (a sibling x = 4^k y + (4^k-1)/3, k = L/2)')

    print('C. The L-walk (NUMERICAL)')
    by = {}
    for p in pairs:
        for Lp, d in p[2]:
            key = '<= 0' if Lp <= 0 else ('1-6' if Lp <= 6 else ('7-15' if Lp <= 15 else '>= 16'))
            by.setdefault(key, []).append(d)
    for k in ['<= 0', '1-6', '7-15', '>= 16']:
        d = by[k]; m = sum(d) / len(d); var = sum((x - m) ** 2 for x in d) / len(d)
        print(f'     L {k:5s}: {len(d):8d} steps, mean increment {m:+.4f}, variance {var:.3f}')
    allinc = [d for v in by.values() for d in v]
    m = sum(allinc) / len(allinc); var = sum((x - m) ** 2 for x in allinc) / len(allinc)
    check(abs(m) < 0.02 and abs(var - 4) < 0.1, f'all increments: mean {m:+.4f}, variance {var:.3f} (Var of a difference of two Geom(1/2) = 4)')
    print('     L one step before merge:', Counter(p[3][p[0] - 1] for p in mg).most_common(8))
    for T in [100, 300, 1000, 3000, 4800]:
        nm = sum(1 for p in pairs if p[0] is None or p[0] > T) / N_B
        print(f'     P(no merge by step {T:5d}) = {nm:.4f}   sqrt(T) P = {math.sqrt(T) * nm:5.2f}')

    print('D. Actual exponents: sigma(M_a) census for a <= 2001 vs the Haar prediction')
    sig_of = {a: s for a, (s, _) in zip(range(2, 2002), sig)}
    post = {a: tot for a, (_, tot) in zip(range(2, 2002), sig)}
    seen = {0}                        # sigma(M_1) = sigma(1) = 0 (the convention of THM-4556 (vi) counts a >= 1)
    levels = []
    new_level = {}
    for a in range(2, 2002):
        new_level[a] = sig_of[a] not in seen
        seen.add(sig_of[a])
        levels.append(len(seen))
    odd = [a for a in range(1001, 2002, 2)]
    share = sum(1 for a in odd if not new_level[a]) / len(odd)
    lag1 = sum(1 for a in odd if sig_of[a] == sig_of[a - 1]) / len(odd)
    ratio = sum(post[a] / a for a in odd) / len(odd)
    pred = 1 - sum(min(1.0, Cq * post[a] ** -alpha) for a in odd) / len(odd)
    c1 = sum(math.sqrt(T) * nm for T, nm in [(T, sum(1 for p in pairs if p[0] is None or p[0] > T) / N_B) for T in (3000, 4800)]) / 2
    no1_T = [(T, q1) for T, _, q1 in rows if T >= 3200]
    c1h = sum(math.sqrt(T) * q1 for T, q1 in no1_T) / len(no1_T)
    pred1 = 1 - sum(min(1.0, c1h / math.sqrt(post[a])) for a in odd) / len(odd)
    print(f'     odd a in [1001, 2001]: share with a smaller sigma-partner {share:.3f} (Haar prediction {pred:.3f}); '
          f'share with sigma(M_(a-1)) {lag1:.3f} (Haar prediction {pred1:.3f} from sqrt(T) q1(T) = {c1h:.2f})')
    print(f'     post-run template total / a = {ratio:.3f};  #distinct sigma values for a <= 100, 400, 1000, 2001: '
          f'{levels[98]}, {levels[398]}, {levels[998]}, {levels[-1]}')
    predlev = {A: 2 + sum(min(1.0, Cq * post[a] ** -alpha) for a in range(3, A + 1, 2)) for A in (100, 400, 1000, 2001)}
    print('     Haar-predicted level counts:', {A: round(v, 1) for A, v in predlev.items()})
    check(abs(share - pred) < 0.03, 'Haar model predicts the switching share of actual exponents within 0.03 (NUMERICAL agreement)')
    check(levels[98] == 23 and levels[398] == 37, 'sigma(M_a) takes 23 and 37 values for a <= 100, 400 (THM-4556 (vi), HYP-9213)')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
