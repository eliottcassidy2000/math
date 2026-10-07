"""How the two coupled orbits' exponents correlate (session opus-2026-10-07-S20).
Note: 05-knowledge/results/two_readers_exponent_correlation_20261007.md.
Pair: y Haar odd 2-adic (a random BITS-bit odd integer), v >= 2 with P(v = k) = 2^-(k-1), x = 3 2^v y + 1.
A_s = e(y_s), B_k = e(x_k), e(z) = v2(3z+1), S_y(s) = A_0 + ... + A_(s-1), S_x(k) = B_0 + ... + B_(k-1).
For a pair (s, k): lam = v + S_y(s) - S_x(k) (offset of y_s's reading window against x_k's in x's bit frame),
j = k + 1 - s, E = (3 x_k + 1) - 2^lam 3^j (3 y_s + 1) (2-adically), mu = v2(E), saturation depth M = mu - max(0, lam).
Exact statements checked:
  (K)  kernel: for z Haar odd and z' = 2^l 3^j z + Delta, Cov(e(z), e(z')) = (2 - 6 2^-m) 1{m >= 1}, m = v2(3 Delta + 1 - 2^l 3^j) - l
       (m = infinity, kernel 2, when that constant vanishes); exact enumeration over z mod 2^22.
  (R)  one source, two readers: the leader's next exponent is Geom(1/2); the lagger's is predicted from the pinned debt.
  (G)  every pair (s, k): the windows overlap iff M >= 1; reading identity  B_k = min(lam + A_s, mu)  (lam >= 0) or
       A_s = -lam + min(B_k, mu)  (lam < 0) away from ties.  Then the law of M on overlaps by time lag d = k - s, and the
       joint law of (fresh exponent, re-read exponent) on overlaps (independent iff M is Geom: the memoryless-competition lemma).
  (W)  window-overlap identity Cov(A_s, B_k) = E[(A_s - 2)(B_k - 2) 1{overlap}] (NUMERICAL agreement; exact by proof).
  (S)  same-step kernel Cov(A_(t+1), B_t | L_t-bin) = E[kappa(m_t)] (NUMERICAL agreement); saturation law P(m_t >= 1 | L_t).
  (C)  cross-covariances Cov(A_(t+1), B_(t+d) | L_t in bin), conditioned on the past only.
Prints ALL CHECKS PASSED.  Runtime about 4-8 minutes on 12 cores."""
import sys, random, math, time, bisect
from multiprocessing import Pool
from collections import Counter
import numpy as np

FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m)
    if not c:
        FAIL.append(m)


def v2(n):
    return (n & -n).bit_length() - 1


def e(z):
    return v2(3 * z + 1)


INF = 10 ** 9
W = 96                        # 2-adic working precision beyond the lag; M >= W is reported as infinity
BITS, NS, TMAX, SEED = 3000, 12000, 600, 2007
DMAX = 40
CLASSES = ['d=-1 merged', 'd=-1 pre-merge', 'd=0', 'd=-2,+1', '2<=|d+1/2|<=5', '6<=|d+1/2|<=12', '13<=|d+1/2|<=40']


def kappa(m):
    if m >= INF:
        return 2.0
    return 2 - 6 * 2.0 ** -m if m >= 1 else 0.0


# ------------------------------------------------------------------ (K) exact kernel by enumeration
def kernel_exact(l, j, Delta, N=22):
    sa = sb = sab = cnt = 0
    for z in range(1, 1 << N, 2):
        zp = (3 ** j) * z * (1 << l) + Delta
        a, b = e(z), e(zp)
        if a >= N - 4 or b >= N + l - 2:
            continue          # unresolved at this precision (Haar mass <= 2^-(N-5))
        sa += a; sb += b; sab += a * b; cnt += 1
    cov = sab / cnt - (sa / cnt) * (sb / cnt)
    E = 3 * Delta + 1 - (1 << l) * 3 ** j
    m = INF if E == 0 else v2(E) - l
    return l, j, Delta, m, cov, kappa(m)


KCASES = [(0, 0, 2), (0, 0, 4), (0, 0, 8), (0, 0, 16), (0, 0, -2), (0, 0, -4), (0, 1, 2), (0, 1, 10), (0, 2, 10), (0, 3, 8),
          (1, 0, 1), (2, 0, 5), (2, 0, -3), (2, 1, 3), (3, 1, 7), (1, 1, 3), (2, 2, -9), (0, 0, 0), (4, 0, 5), (6, 0, 21)]


# ------------------------------------------------------------------ comparison constant
def mu_val(X, Y, lam, j):
    """v2 of E = (3X+1) - 2^lam 3^j (3Y+1), computed 2-adically (3^j a unit); INF if >= max(0, lam) + W."""
    a3 = max(0, -j); b3 = max(0, j)
    if lam >= 0:
        nb = lam + W; mod = (1 << nb) - 1
        val = (pow(3, a3, 1 << nb) * ((3 * X + 1) & mod) - ((pow(3, b3, 1 << nb) * ((3 * Y + 1) & mod)) << lam)) & mod
        return v2(val) if val else INF
    L_ = -lam; nb = L_ + W; mod = (1 << nb) - 1
    val = (((pow(3, a3, 1 << nb) * ((3 * X + 1) & mod)) << L_) - pow(3, b3, 1 << nb) * ((3 * Y + 1) & mod)) & mod
    return (v2(val) - L_) if val else INF


def dclass(d, merged):
    if d == -1:
        return 0 if merged else 1
    if d == 0:
        return 2
    if d in (-2, 1):
        return 3
    h = abs(d + 0.5)
    if h <= 5:
        return 4
    if h <= 12:
        return 5
    return 6


# ------------------------------------------------------------------ trajectories
def run(s):
    rng = random.Random(SEED * 1000003 + s)
    y = rng.getrandbits(BITS) | 1
    v = 2
    while rng.random() < 0.5:
        v += 1
    x = 3 * (1 << v) * y + 1
    A, B, Y, X = [], [], [], []
    yy, xx = y, x
    for t in range(TMAX + DMAX + 2):
        a = e(yy); b = e(xx)
        A.append(a); B.append(b); Y.append(yy); X.append(xx)
        yy = (3 * yy + 1) >> a; xx = (3 * xx + 1) >> b
    Sy = [0]; Sx = [0]
    for a in A:
        Sy.append(Sy[-1] + a)
    for b in B:
        Sx.append(Sx[-1] + b)
    tau = next((t for t in range(TMAX) if X[t] == Y[t + 1]), None)
    # debt bookkeeping (up to the merge): x_t = 2^L y_(t+1) + Delta_t
    rec = []
    pred_ok = True
    for t in range(TMAX if tau is None else tau + 1):
        L = v + Sy[t + 1] - Sx[t]
        m = None
        if L >= 0:
            Delta = X[t] - (Y[t + 1] << L)
            D = 3 * Delta + 1 - (1 << L)
            m = (v2(D) - L) if D != 0 else INF
            vd = v2(D) if D != 0 else INF
            s_ = L + A[t + 1]
            if vd < s_:
                pred_ok &= (B[t] == vd)
            elif vd > s_:
                pred_ok &= (B[t] == s_)
        rec.append((t, L, m))
    # (G) all overlapping pairs (s, k) with 10 <= s < TMAX - DMAX, |k - s| <= DMAX, plus the adjacent non-overlapping k
    nc = len(CLASSES)
    g_n = np.zeros(nc); g_prod = np.zeros(nc); g_kap = np.zeros(nc); g_p2 = np.zeros(nc)
    g_M = np.zeros((nc, 9)); g_joint = np.zeros((nc, 6, 6))
    ov_by_d = np.zeros(2 * DMAX + 1); prod_by_d = np.zeros(2 * DMAX + 1); kap_by_d = np.zeros(2 * DMAX + 1)
    dev2_by_d = np.zeros(2 * DMAX + 1)
    # alignment split of pre-merge overlaps: 0 = exact alignment (lam = 0), 1 = offset overlap (lam != 0); long gaps separately
    a_n = np.zeros((2, 2)); a_kap = np.zeros((2, 2)); a_prod = np.zeros((2, 2)); a_info = np.zeros((2, 2)); a_absk = np.zeros((2, 2))
    a_M = np.zeros((2, 2, 9))
    iff_ok = read_ok = bsc_ok = True
    nonov_checked = 0
    for s in range(10, TMAX - DMAX):
        lo = Sy[s] + v; hi = Sy[s + 1] + v                      # y_s reads (lo, hi] in x's frame
        kf = bisect.bisect_right(Sx, lo) - 1
        kl = bisect.bisect_left(Sx, hi) - 1
        for k in range(max(kf - 1, 0), kl + 2):
            d = k - s
            if abs(d) > DMAX:
                continue
            ov = Sx[k] < hi and Sx[k + 1] > lo
            lam = v + Sy[s] - Sx[k]
            j = k + 1 - s
            mu = mu_val(X[k], Y[s], lam, j)
            M = INF if mu >= INF else mu - max(0, lam)
            iff_ok &= (ov == (M >= 1))
            if lam >= 0:
                if mu != lam + A[s]:
                    read_ok &= (B[k] == min(lam + A[s], mu))
            else:
                if mu != B[k]:
                    read_ok &= (A[s] == -lam + min(B[k], mu))
            if not ov:
                nonov_checked += 1
                continue
            merged = tau is not None and k >= tau and s >= tau + 1
            c = dclass(d, merged)
            prod = (A[s] - 2) * (B[k] - 2)
            kap = kappa(M)
            mbin = 8 if M >= INF else min(M, 7)
            g_n[c] += 1; g_prod[c] += prod; g_kap[c] += kap; g_p2[c] += 0.0 if M >= INF else 2.0 ** -M
            g_M[c, mbin] += 1
            F, R = (A[s], B[k] - lam) if lam >= 0 else (B[k], A[s] + lam)     # fresh (later window) and re-read exponent
            g_joint[c, min(F, 5), min(R, 5)] += 1
            bsc_ok &= ((R == 1) == ((F == 1) != (M == 1)))                     # first re-read bit = dictator XOR 1{M = 1}
            ov_by_d[d + DMAX] += 1; prod_by_d[d + DMAX] += prod; kap_by_d[d + DMAX] += kap; dev2_by_d[d + DMAX] += (prod - kap) ** 2
            if not (tau is not None and k >= tau - 1 and s >= tau):          # pre-merge (excludes the merging step and after)
                ai = 0 if lam == 0 else 1
                gi = 1 if abs(d + 0.5) >= 6 else 0
                a_n[ai, gi] += 1; a_kap[ai, gi] += kap; a_prod[ai, gi] += prod; a_absk[ai, gi] += abs(kap)
                a_info[ai, gi] += 2.0 if M >= INF else 2 * (1 - 2.0 ** -M)
                a_M[ai, gi, mbin] += 1
    # full covariance numerator at each lag d (all pairs, overlapping or not), same s range
    a = np.array(A, dtype=np.float64) - 2.0; b = np.array(B, dtype=np.float64) - 2.0
    srange = np.arange(10, TMAX - DMAX)
    full_by_d = np.array([np.sum(a[srange] * b[srange + d]) for d in range(-DMAX, DMAX + 1)])
    ns = len(srange)
    gstats = (g_n, g_prod, g_kap, g_p2, g_M, g_joint, ov_by_d, prod_by_d, kap_by_d, full_by_d, ns, iff_ok, read_ok, nonov_checked,
              dev2_by_d, (a_n, a_kap, a_prod, a_info, a_absk, a_M), bsc_ok)
    return v, np.array(A[:TMAX + 2], dtype=np.int16), np.array(B[:TMAX + 2], dtype=np.int16), rec, tau, pred_ok, gstats


if __name__ == '__main__':
    t0 = time.time()
    with Pool(12) as pool:
        kres = pool.starmap(kernel_exact, KCASES)
        print('(K) universal comparison kernel, exact enumeration over z mod 2^22')
        for l, j, Delta, m, cov, kap in kres:
            print('     (l, j, Delta) = (%d, %d, %d): m = %s, Cov = %+.4f, kernel = %+.4f' % (l, j, Delta, 'inf' if m >= INF else m, cov, kap))
        check(all(abs(cov - kap) < 3e-3 for l, j, Delta, m, cov, kap in kres),
              'Cov(e(z), e(2^l 3^j z + Delta)) = (2 - 6 2^-m) 1{m >= 1} (m = inf -> 2) for all listed maps (truncation error < 3e-3)')
        res = pool.map(run, range(NS), chunksize=20)
    print(f'     ({NS} pairs, {BITS}-bit y, {TMAX} steps; merged within {TMAX} steps: {sum(1 for r in res if r[4] is not None) / NS:.4f})')

    print('(R) one source, two readers')
    check(all(r[5] for r in res), 'whenever y leads (L_t >= 0) and the lagger does not tie, B_t equals its prediction from the pinned debt (all steps, all pairs)')
    lead = Counter(); lead2 = Counter()
    for v, A, B, rec, tau, _, _ in res:
        for (t, L, m) in rec:
            if t < 3 or (tau is not None and t >= tau):
                continue
            if L >= 0:
                lead[min(int(A[t + 1]), 6)] += 1
            else:
                lead2[min(int(B[t]), 6)] += 1
    n1, n2 = sum(lead.values()), sum(lead2.values())
    print('     leader y (L>=0): P(A=k) k=1..5: ' + ', '.join(f'{lead[k] / n1:.4f}' for k in range(1, 6)) +
          ';  leader x (L<0): P(B=k): ' + ', '.join(f'{lead2[k] / n2:.4f}' for k in range(1, 6)))
    check(all(abs(lead[k] / n1 - 2.0 ** -k) < 0.004 for k in range(1, 6)) and all(abs(lead2[k] / n2 - 2.0 ** -k) < 0.006 for k in range(1, 6)),
          'the leader\'s exponent is Geom(1/2) in both regimes (NUMERICAL, within 0.004/0.006)')

    print('(G) all pairs (s, k), |k - s| <= %d, 10 <= s < %d' % (DMAX, TMAX - DMAX))
    G = [r[6] for r in res]
    check(all(g[11] for g in G), 'pointwise: the reading windows overlap iff the saturation depth M >= 1 (every overlapping pair and its non-overlapping neighbours: %d overlaps, %d neighbours)'
          % (int(sum(g[0].sum() for g in G)), sum(g[13] for g in G)))
    check(all(g[12] for g in G), 'pointwise reading identity away from ties: B_k = min(lam + A_s, mu) (lam >= 0), A_s = -lam + min(B_k, mu) (lam < 0)')
    check(all(g[16] for g in G), 'pointwise BSC identity on every overlap: 1{R = 1} = 1{F = 1} XOR 1{M = 1} (F fresh exponent, R re-read exponent minus offset)')
    g_n = sum(g[0] for g in G); g_prod = sum(g[1] for g in G); g_kap = sum(g[2] for g in G); g_p2 = sum(g[3] for g in G)
    g_M = sum(g[4] for g in G); g_joint = sum(g[5] for g in G)
    ov_by_d = sum(g[6] for g in G); prod_by_d = sum(g[7] for g in G); kap_by_d = sum(g[8] for g in G); full_by_d = sum(g[9] for g in G)
    ns_tot = sum(g[10] for g in G); dev2_by_d = sum(g[14] for g in G)
    a_n, a_kap, a_prod, a_info, a_absk, a_M = [sum(g[15][i] for g in G) for i in range(6)]
    print('     law of M on overlapping pairs by time-lag class (Geom: 0.5, 0.25, 0.125, 0.0625, P(7 <= M < inf) = 1/64; E[2^-M] = 1/3; E[kappa] = 0):')
    long_ok = True
    for c, name in enumerate(CLASSES):
        n = g_n[c]
        if n == 0:
            continue
        pm = g_M[c] / n
        J = g_joint[c] / n
        PF = J.sum(axis=1); PR = J.sum(axis=0)
        dep = max(abs(J[i, k] - PF[i] * PR[k]) for i in range(1, 5) for k in range(1, 5))
        al = pm[1]
        h2 = 0.0 if al in (0.0, 1.0) else -(al * math.log2(al) + (1 - al) * math.log2(1 - al))
        print(f'     {name:17s} n = {int(n):9d}: P(M=1..4) = ' + ', '.join(f'{pm[m]:.4f}' for m in range(1, 5)) +
              f', P(7<=M<inf) = {pm[7]:.4f}, P(M=inf) = {pm[8]:.4f}; E[2^-M] = {g_p2[c] / n:.4f}; E[kappa] = {g_kap[c] / n:+.4f},'
              f' E[(A-2)(B-2)] = {g_prod[c] / n:+.4f}; max|P(F,R) - P(F)P(R)| = {dep:.4f}; first-bit BSC 1 - h2(P(M=1)) = {1 - h2:.5f} bits')
        if c >= 5:
            se = 1 / math.sqrt(n)
            long_ok &= abs(pm[1] - 0.5) < 4 * se and abs(g_p2[c] / n - 1 / 3) < 4 * se and pm[8] == 0.0
    check(long_ok, 'at time lags |d + 1/2| >= 6 the saturation depth on overlaps is Geom(1/2) within 4 s.e. and never infinite (NUMERICAL)')
    print('     pre-merge overlaps by alignment type (lam = 0: exact alignment, the case of THM-4564; lam != 0: offset overlap):')
    for gi, gname in ((0, '|d+1/2| < 6'), (1, '|d+1/2| >= 6')):
        for ai, aname in ((0, 'lam = 0 '), (1, 'lam != 0')):
            n = a_n[ai, gi]
            pm = a_M[ai, gi] / n
            print(f'       {gname:12s} {aname}: n = {int(n):8d} (share {n / a_n[:, gi].sum():.3f}); P(M=1..3) = ' + ', '.join(f'{pm[m]:.4f}' for m in range(1, 4)) +
                  f'; E[kappa] = {a_kap[ai, gi] / n:+.4f}, E[(A-2)(B-2)] = {a_prod[ai, gi] / n:+.4f}; E|kappa| = {a_absk[ai, gi] / n:.3f},'
                  f' E[I(F;R | past)] = {a_info[ai, gi] / n:.3f} bits')
    check(a_n[1].sum() > a_n[0].sum() and all(a_absk[1, gi] / a_n[1, gi] > 0.5 for gi in (0, 1)),
          'offset overlaps (lam != 0) outnumber exact alignments and carry conditional coupling E|kappa(M)| > 0.5 (exact kernel; counts NUMERICAL)')
    print('     covariance at lag d (pooled over s): full Cov(A_s, B_(s+d)), its overlap part, the kernel prediction, P(overlap):')
    rows = []
    for d in list(range(-8, 8)) + [-12, 11, -20, 19, -30, 29, -40]:
        i = d + DMAX
        rows.append((d, full_by_d[i] / ns_tot, prod_by_d[i] / ns_tot, kap_by_d[i] / ns_tot, ov_by_d[i] / ns_tot, math.sqrt(dev2_by_d[i]) / ns_tot))
    for d, cf, co, ck, po, sd in rows:
        print(f'       d = {d:+3d}: Cov = {cf:+.4f}, overlap part = {co:+.4f}, E[kappa(M) 1_O] = {ck:+.4f} (s.e. of difference {sd:.4f}), P(O) = {po:.4f}')
    se = 2 / math.sqrt(ns_tot)
    check(all(abs(cf - co) < 4 * se for d, cf, co, ck, po, sd in rows), f'window-overlap identity at every listed lag (agreement within 4 s.e. = {4 * se:.4f}; exact by proof)')
    check(all(abs(co - ck) < 4 * sd + 1e-4 for d, cf, co, ck, po, sd in rows), 'overlap part = E[kappa(M) 1_O] at every listed lag (within 4 s.e. of the pairwise difference; exact by proof)')

    print('(W) window-overlap identity over t in [8, 200] (independent tally, includes s < 10)')
    wl = []
    for d in (-6, -3, -1, 0, 1, 6):
        s_all = s_ov = 0.0; n = n_ov = 0
        for v, A, B, rec, tau, _, _ in res:
            Sy = np.concatenate([[0], np.cumsum(A)]); Sx = np.concatenate([[0], np.cumsum(B)])
            for t in range(8, 200):
                k = t + d
                prod = (int(A[t]) - 2) * (int(B[k]) - 2)
                ov = (Sx[k] < Sy[t + 1] + v) and (Sx[k + 1] > Sy[t] + v)
                s_all += prod; n += 1
                if ov:
                    s_ov += prod; n_ov += 1
        wl.append((d, s_all / n, s_ov / n, n_ov / n, 2 / math.sqrt(n)))
        print(f'     d = {d:+d}: Cov = {s_all / n:+.4f}, E[(A-2)(B-2)1_O] = {s_ov / n:+.4f}, P(overlap) = {n_ov / n:.4f} (s.e. {2 / math.sqrt(n):.4f})')
    check(all(abs(c - o) < 4 * se for d, c, o, p, se in wl), 'Cov and its overlap part agree within 4 s.e. at every tested lag')

    print('(S) same-step kernel and saturation law (strictly pre-merge steps t < tau - 1 with L_t >= 0, t >= 3; tau = merge time)')
    byL = {}
    mergeL = Counter()
    for v, A, B, rec, tau, _, _ in res:
        for (t, L, m) in rec:
            if tau is not None and t == tau - 1:
                mergeL[L] += 1                     # the merging step: D_t = 0, m_t = infinity
            if m is None or t < 3 or (tau is not None and t >= tau - 1):
                continue
            byL.setdefault(L, []).append((m, int(A[t + 1]), int(B[t])))
    print('     merging steps (D_t = 0) by L_t: ' + ', '.join(f'L={l}: {mergeL[l]}' for l in sorted(mergeL)[:8]))
    check(all(l != 0 and l % 2 == 0 for l in mergeL), 'every merging step has L_t even and nonzero: the relation one step before the merge is a sibling relation, x_t = 4^k y_(t+1) + (4^k - 1)/3 (L_t = 2k) or y_(t+1) = 4^k x_t + (4^k - 1)/3 (L_t = -2k)')
    srows = []
    for lo, hi in ((0, 0), (1, 4), (5, 12), (13, 10 ** 6)):
        R = [r for L in byL if lo <= L <= hi for r in byL[L]]
        pred = sum(kappa(m) for m, _, _ in R) / len(R)
        ma = sum(a for _, a, _ in R) / len(R); mb = sum(b for _, _, b in R) / len(R)
        cov = sum(a * b for _, a, b in R) / len(R) - ma * mb
        sat = [m for m, _, _ in R if m >= 1]
        pm = [sum(1 for m in sat if m == k) / max(len(sat), 1) for k in (1, 2, 3, 4)]
        srows.append((lo, hi, len(R), pred, cov))
        al = pm[0]
        ck = 1 + (al * math.log2(al) + (1 - al) * math.log2(1 - al)) if 0 < al < 1 else 1.0
        print(f'     L in [{lo},{hi}]: n = {len(R)}, P(m>=1) = {len(sat) / len(R):.4f}, P(m=1..4 | m>=1) = ' + ', '.join(f'{p:.3f}' for p in pm) +
              f'; predicted E[kappa(m)] = {pred:+.4f}, measured Cov(A_(t+1), B_t) = {cov:+.4f}; first-bit BSC 1 - h2(P(m=1 | m>=1)) = {ck:.4f} bits')
    check(all(abs(p - c) < 4 * 2 / math.sqrt(n) for lo, hi, n, p, c in srows), 'same-step covariance = E[kappa(m_t)] in every L-bin (within 4 s.e.)')
    sat = []
    for l in range(0, 13):
        R = byL.get(l, [])
        if len(R) < 2000:
            continue
        sat.append((l, sum(1 for m, _, _ in R if m >= 1) / len(R)))
    print('     P(m >= 1 | L = l): ' + ', '.join(f'l={l}: {p:.4f}' for l, p in sat))
    check(all(p <= 2.0 ** -l * 2 + 0.01 for l, p in sat if l >= 1), 'saturation P(m >= 1 | L = l) is of order 2^-l (at most 2^(1-l) + 0.01)')

    print('(C) Cov(A_(t+1), B_(t+d) | L_t in bin, not merged by t), d = -8..8, past conditioning only, t in [10, %d)' % (TMAX - 10))
    bins = [(-10 ** 6, -9), (-8, -1), (0, 0), (1, 8), (9, 10 ** 6)]
    for lo, hi in bins:
        rowc = []; nn = 0
        for d in range(-8, 9):
            s_ab = s_a = s_b = 0.0; n = 0
            for v, A, B, rec, tau, _, _ in res:
                Sy = np.concatenate([[0], np.cumsum(A)]); Sx = np.concatenate([[0], np.cumsum(B)])
                t = np.arange(10, TMAX - 10 if tau is None else min(TMAX - 10, tau))
                L = v + Sy[t + 1] - Sx[t]
                sel = t[(L >= lo) & (L <= hi)]
                aa = A[sel + 1].astype(np.float64); bb = B[sel + d].astype(np.float64)
                s_ab += float(np.sum(aa * bb)); s_a += float(np.sum(aa)); s_b += float(np.sum(bb)); n += len(sel)
            rowc.append(s_ab / n - (s_a / n) * (s_b / n)); nn = n
        print(f'     L in [{lo},{hi}] (n = {nn}): ' + ' '.join(f'{c:+.3f}' for c in rowc))

    print('(P) the near-synchronous coupling resolved by single L_t = l (not merged by t): Cov(A_(t+1), B_(t+d)), d = -2..8')
    far = []
    for l in list(range(-4, 0)) + list(range(1, 13)):
        acc = {d: [0.0, 0.0, 0.0, 0] for d in range(-2, 9)}
        for v, A, B, rec, tau, _, _ in res:
            Sy = np.concatenate([[0], np.cumsum(A)]); Sx = np.concatenate([[0], np.cumsum(B)])
            t = np.arange(10, TMAX - 10 if tau is None else min(TMAX - 10, tau))
            L = v + Sy[t + 1] - Sx[t]
            sel = t[L == l]
            for d in acc:
                aa = A[sel + 1].astype(np.float64); bb = B[sel + d].astype(np.float64)
                acc[d][0] += float(np.sum(aa * bb)); acc[d][1] += float(np.sum(aa)); acc[d][2] += float(np.sum(bb)); acc[d][3] += len(sel)
        n = acc[0][3]
        row = [acc[d][0] / n - (acc[d][1] / n) * (acc[d][2] / n) for d in acc]
        se = 2 / math.sqrt(n)
        print(f'     L = {l:+3d} (n = {n:6d}, s.e. {se:.4f}): ' + ' '.join(f'{c:+.3f}' for c in row) + f'   max |Cov|/s.e. = {max(abs(c) for c in row) / se:.1f}')
        if l >= 9:
            far.append(max(abs(c) for c in row) / se)
    check(max(far) < 4.5, 'for L_t = 9..12 no cross-covariance at lags -2..8 exceeds 4.5 s.e. (NUMERICAL)')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
