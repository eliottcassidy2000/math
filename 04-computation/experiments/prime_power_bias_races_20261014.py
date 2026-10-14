#!/usr/bin/env python3
"""Prime powers, not primes: the biases behind Chebyshev's race, Skewes' crossing and Polya's conjecture, and the
unbiased Mertens function, computed and compared.

Session opus-2026-10-14-S32.  Note: 05-knowledge/results/prime_power_bias_collatz_collector_20261014.md

Part A (segmented sieve to N_A = 10^10, every prime visited in order):
  * exact sign-change regions of pi(x;4,3) - pi(x;4,1) and of pi(x;3,2) - pi(x;3,1) (the races);
  * pi(26861; 4, 1) and pi(26861; 4, 3);
  * pi(x) < li(x) at every prime p <= N_A (li anchored by a convergent series each segment, Simpson + Kahan inside);
  * samples on a log grid (step 1/1000 in log x) of the normalised races
        E_4(x) = (log x/sqrt x)(pi(x;4,3) - pi(x;4,1)),   E_3(x) = (log x/sqrt x)(pi(x;3,2) - pi(x;3,1)),
        S(x)   = (log x/sqrt x)(pi(x) - li(x)),
    and of the prime-power-corrected versions, in which every prime power p^k counts with weight 1/k:
        Pi(x; chi) = sum_{p^k <= x} chi(p^k)/k,   Pi(x) = sum_{p^k <= x} 1/k   (Riemann's prime-power count).
    The corrected races are what the explicit formula sees; the difference is the prime squares (and cubes ...).
Part B (segmented Omega-parity / squarefree sieve to N_B = 10^9):
  * Liouville L(x) = sum lambda(n) and Mertens M(x) = sum mu(n);
  * the first n >= 2 with L(n) > 0 (Polya; Tanaka 1980 found 906150257), the extent and maximum of that first region;
  * max |M(x)|/sqrt x on [10^4, N_B];
  * log-grid samples of L(x)/sqrt x and M(x)/sqrt x (their log-averages: 1/zeta(1/2) = -0.6848 vs 0 expected).

Reproduce: python3 04-computation/experiments/prime_power_bias_races_20261014.py   (numba; several minutes)
"""
import math, time, sys, json
import numpy as np
import numba as nb
import mpmath as mp

N_A = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 10
N_B = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10 ** 9
PI_TABLE = {10 ** 6: 78498, 10 ** 7: 664579, 10 ** 8: 5761455, 10 ** 9: 50847534, 10 ** 10: 455052511}
SEG = 1 << 24
GRID_STEP = 1e-3          # log-grid spacing
FAIL = []


def check(label, ok):
    if not ok:
        FAIL.append(label)
    return ok


def small_primes(n):
    s = np.ones(n + 1, dtype=np.bool_)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


@nb.njit(cache=True)
def li_series(x):
    """li(x) for x > 1 by Ramanujan's series (float64)."""
    g = 0.5772156649015329
    lx = math.log(x)
    tot = 0.0
    fact = 1.0
    inner = 0.0
    term_pow = 1.0
    for n in range(1, 400):
        fact *= n
        term_pow *= lx
        if (n - 1) % 2 == 0:
            inner += 1.0 / (2 * ((n - 1) // 2) + 1)
        sign = 1.0 if (n - 1) % 2 == 0 else -1.0
        t = sign * term_pow / (fact * 2.0 ** (n - 1)) * inner
        tot += t
        if n > 10 and abs(t) < 1e-18 * abs(tot):
            break
    return g + math.log(lx) + math.sqrt(x) * tot


@nb.njit(cache=True)
def sweep_races(N, seg, base_primes, grid_x, out_grid, regions4, regions3, info):
    """Visit every prime <= N in order.  out_grid[i] = (pi, pi41, pi43, pi31, pi32, li) at grid_x[i] (counts of primes <= x)."""
    n_grid = grid_x.shape[0]
    gi = 0
    while gi < n_grid and grid_x[gi] < 2.0:
        gi += 1
    pi = 0
    p41 = 0
    p43 = 0
    p31 = 0
    p32 = 0
    # race states: sign of (pi(x;4,1) - pi(x;4,3)); regions where > 0 recorded as (start, end, max lead)
    lead4 = False
    lead3 = False
    nr4 = 0
    nr3 = 0
    max_lead4 = 0
    max_lead3 = 0
    # Skewes check
    li_prev_x = 2.0
    li_val = li_series(2.0)
    comp = 0.0                      # Kahan compensation
    max_pi_minus_li = -1e300
    argmax_pml = 0.0
    sieve = np.ones(seg, dtype=np.uint8)
    lo = 0
    while lo <= N:
        hi = min(lo + seg - 1, N)
        L = hi - lo + 1
        for i in range(L):
            sieve[i] = 1
        if lo == 0:
            sieve[0] = 0
            if L > 1:
                sieve[1] = 0
        for bi in range(base_primes.shape[0]):
            p = base_primes[bi]
            if p * p > hi:
                break
            start = (lo + p - 1) // p * p
            if start < p * p:
                start = p * p
            for m in range(start, hi + 1, p):
                sieve[m - lo] = 0
        # re-anchor li at the segment start (exact series), unless at the very beginning
        if lo > 2:
            li_val = li_series(float(li_prev_x))
            comp = 0.0
        for i in range(L):
            if sieve[i] == 0:
                continue
            p = lo + i
            # grid points strictly below p see the counts before p
            while gi < n_grid and grid_x[gi] < p:
                # li at grid point by Simpson from li_prev_x
                a = li_prev_x
                b = grid_x[gi]
                if b > a:
                    mid = 0.5 * (a + b)
                    inc = (b - a) / 6.0 * (1.0 / math.log(a) + 4.0 / math.log(mid) + 1.0 / math.log(b))
                else:
                    inc = 0.0
                out_grid[gi, 0] = pi
                out_grid[gi, 1] = p41
                out_grid[gi, 2] = p43
                out_grid[gi, 3] = p31
                out_grid[gi, 4] = p32
                out_grid[gi, 5] = li_val + inc
                gi += 1
            # advance li to p
            a = li_prev_x
            b = float(p)
            if b > a:
                mid = 0.5 * (a + b)
                inc = (b - a) / 6.0 * (1.0 / math.log(a) + 4.0 / math.log(mid) + 1.0 / math.log(b))
                y = inc - comp
                t = li_val + y
                comp = (t - li_val) - y
                li_val = t
                li_prev_x = b
            pi += 1
            r4 = p % 4
            if r4 == 1:
                p41 += 1
            elif r4 == 3:
                p43 += 1
            r3 = p % 3
            if r3 == 1:
                p31 += 1
            elif r3 == 2:
                p32 += 1
            d = pi - li_val
            if p >= 3 and d > max_pi_minus_li:
                max_pi_minus_li = d
                argmax_pml = p
            # race mod 4: lead when pi(x;4,1) > pi(x;4,3)
            if p41 > p43:
                if not lead4:
                    lead4 = True
                    if nr4 < regions4.shape[0]:
                        regions4[nr4, 0] = p
                        regions4[nr4, 2] = p41 - p43
                        regions4[nr4, 3] = p41
                        regions4[nr4, 4] = p43
                    max_lead4 = p41 - p43
                else:
                    if p41 - p43 > max_lead4:
                        max_lead4 = p41 - p43
                        if nr4 < regions4.shape[0]:
                            regions4[nr4, 2] = max_lead4
            else:
                if lead4:
                    lead4 = False
                    if nr4 < regions4.shape[0]:
                        regions4[nr4, 1] = p        # first prime at which the lead is lost
                    nr4 += 1
            if p31 > p32:
                if not lead3:
                    lead3 = True
                    if nr3 < regions3.shape[0]:
                        regions3[nr3, 0] = p
                        regions3[nr3, 2] = p31 - p32
                    max_lead3 = p31 - p32
                else:
                    if p31 - p32 > max_lead3:
                        max_lead3 = p31 - p32
                        if nr3 < regions3.shape[0]:
                            regions3[nr3, 2] = max_lead3
            else:
                if lead3:
                    lead3 = False
                    if nr3 < regions3.shape[0]:
                        regions3[nr3, 1] = p
                    nr3 += 1
        lo = hi + 1
    while gi < n_grid:
        out_grid[gi, 0] = pi
        out_grid[gi, 1] = p41
        out_grid[gi, 2] = p43
        out_grid[gi, 3] = p31
        out_grid[gi, 4] = p32
        out_grid[gi, 5] = li_series(grid_x[gi])
        gi += 1
    info[0] = pi
    info[1] = p41
    info[2] = p43
    info[3] = p31
    info[4] = p32
    info[5] = nr4 + (1 if lead4 else 0)
    info[6] = nr3 + (1 if lead3 else 0)
    info[7] = max_pi_minus_li
    info[8] = argmax_pml


@nb.njit(cache=True)
def sweep_liouville_mertens(N, seg, base_primes, grid_x, out_grid, info, first_region, runs):
    """lambda and mu by sieving Omega parity and squarefreeness; a leftover prime > sqrt(N) is detected with log sums."""
    n_grid = grid_x.shape[0]
    gi = 0
    L = 0
    M = 0
    first_pos = 0          # first n >= 2 with L(n) > 0
    region_end = 0         # first n after first_pos with L(n) <= 0
    region_max = 0
    region_argmax = 0
    in_region = False
    nruns = 0
    max_absM = 0.0
    argmax_absM = 0
    max_L_norm = -1e300
    min_L_norm = 1e300
    parity = np.zeros(seg, dtype=np.uint8)
    sqf = np.ones(seg, dtype=np.uint8)
    logsum = np.zeros(seg, dtype=np.float32)
    lo = 1
    while lo <= N:
        hi = min(lo + seg - 1, N)
        Lseg = hi - lo + 1
        for i in range(Lseg):
            parity[i] = 0
            sqf[i] = 1
            logsum[i] = 0.0
        for bi in range(base_primes.shape[0]):
            p = base_primes[bi]
            if p * p > hi and p > hi:
                break
            if p > hi:
                break
            lp = math.log(p)
            pk = p
            k = 1
            while pk <= hi:
                start = (lo + pk - 1) // pk * pk
                for m in range(start, hi + 1, pk):
                    j = m - lo
                    parity[j] ^= 1
                    logsum[j] += lp
                    if k >= 2:
                        sqf[j] = 0
                if pk > hi // p:
                    break
                pk *= p
                k += 1
        for i in range(Lseg):
            n = lo + i
            par = parity[i]
            if n > 1 and math.log(n) - logsum[i] > 0.5:
                par ^= 1           # one prime factor above sqrt(N)
            lam = 1 if par == 0 else -1
            L += lam
            if sqf[i] == 1:
                M += lam
            if n >= 2:
                if L > 0:
                    if first_pos == 0:
                        first_pos = n
                    if not in_region:
                        in_region = True
                        if nruns < runs.shape[0]:
                            runs[nruns, 0] = n
                            runs[nruns, 2] = L
                            runs[nruns, 3] = n
                    if nruns < runs.shape[0] and L > runs[nruns, 2]:
                        runs[nruns, 2] = L
                        runs[nruns, 3] = n
                    if L > region_max:
                        region_max = L
                        region_argmax = n
                elif in_region:
                    in_region = False
                    if nruns < runs.shape[0]:
                        runs[nruns, 1] = n
                    nruns += 1
                    if region_end == 0:
                        region_end = n
            if n >= 10000:
                a = abs(M) / math.sqrt(n)
                if a > max_absM:
                    max_absM = a
                    argmax_absM = n
                ln = L / math.sqrt(n)
                if ln > max_L_norm:
                    max_L_norm = ln
                if ln < min_L_norm:
                    min_L_norm = ln
            while gi < n_grid and grid_x[gi] < n + 1:
                if grid_x[gi] >= n:
                    out_grid[gi, 0] = L
                    out_grid[gi, 1] = M
                gi += 1
        lo = hi + 1
    info[0] = L
    info[1] = M
    info[2] = max_absM
    info[3] = argmax_absM
    info[4] = max_L_norm
    info[5] = min_L_norm
    first_region[0] = first_pos
    first_region[1] = region_end
    first_region[2] = region_max
    first_region[3] = region_argmax
    first_region[4] = nruns


def chi4(p):
    return 0 if p % 2 == 0 else (1 if p % 4 == 1 else -1)


def chi3(p):
    return 0 if p % 3 == 0 else (1 if p % 3 == 1 else -1)


def power_corrections(xs, smallp):
    """For each x: sum_{k>=2} (1/k) sum_{p^k <= x} chi(p)^k for chi = chi_{-4}, chi_{-3}, and the trivial character (vectorised)."""
    pr = np.asarray(smallp, dtype=np.int64)
    c4 = np.array([chi4(int(q)) for q in pr], dtype=np.int64)
    c3 = np.array([chi3(int(q)) for q in pr], dtype=np.int64)
    out = np.zeros((len(xs), 3))
    xs = np.asarray(xs, dtype=np.float64)
    kmax = int(math.log2(xs.max())) + 1
    for k in range(2, kmax + 1):
        r = np.floor(xs ** (1.0 / k) + 1e-9)
        r = np.where((r + 1) ** k <= xs, r + 1, r)
        r = np.where(r ** k > xs, r - 1, r)
        idx = np.searchsorted(pr, r, side='right')          # number of primes <= r
        pre4 = np.concatenate([[0], np.cumsum(c4 ** k)])
        pre3 = np.concatenate([[0], np.cumsum(c3 ** k)])
        out[:, 0] += pre4[idx] / k
        out[:, 1] += pre3[idx] / k
        out[:, 2] += idx / k
    return out


def main():
    t0 = time.time()
    print('== Part A: races and Skewes to N_A = %.0e' % N_A)
    base = small_primes(int(N_A ** 0.5) + 10)
    ys = np.arange(math.log(1000.0), math.log(float(N_A)) - 1e-9, GRID_STEP)
    grid_x = np.exp(ys)
    out_grid = np.zeros((len(grid_x), 6))
    regions4 = np.zeros((100000, 5), dtype=np.int64)
    regions3 = np.zeros((100000, 3), dtype=np.int64)
    info = np.zeros(9)
    sweep_races(N_A, SEG, base, grid_x, out_grid, regions4, regions3, info)
    print('   sweep done in %.1f s' % (time.time() - t0))
    pi_N, p41, p43, p31, p32 = (int(v) for v in info[:5])
    print(f'   pi(N_A) = {pi_N} (table: {PI_TABLE.get(N_A)}: {check("pi(N_A)", PI_TABLE.get(N_A, pi_N) == pi_N)});  pi(N;4,1) = {p41}, pi(N;4,3) = {p43};  pi(N;3,1) = {p31}, pi(N;3,2) = {p32}')
    n4 = int(info[5])
    print(f'   maximal runs of primes with pi(x;4,1) > pi(x;4,3), x <= N_A: {n4}')
    R = regions4[:n4]
    for r in R[:3]:
        print(f'     run from prime {r[0]} to (first prime with the lead lost) {r[1]};  max lead {r[2]};  at the entry pi(x;4,1) = {r[3]}, pi(x;4,3) = {r[4]}')
    # clusters: runs closer than 1% of x belong to one excursion
    clusters = []
    for r in R:
        if clusters and r[0] < clusters[-1][1] * 1.01:
            c = clusters[-1]
            clusters[-1] = (c[0], int(max(c[1], r[1])), int(max(c[2], r[2])), c[3] + 1)
        else:
            clusters.append((int(r[0]), int(r[1]), int(r[2]), 1))
    print(f'   excursions (runs merged when closer than 1% of x): {len(clusters)}')
    for c in clusters:
        print(f'     excursion {c[0]} .. {c[1]}:  {c[3]} runs, max lead {c[2]}')
    check('26861', n4 >= 1 and regions4[0, 0] == 26861)
    print(f'   regions with pi(x;3,1) > pi(x;3,2), x <= N_A: {int(info[6])}  (Bays-Hudson: the first is at 608981813029 > N_A)')
    check('mod3 none', int(info[6]) == 0 or N_A > 608981813029)
    print(f'   max over primes 3 <= p <= N_A of pi(p) - li(p): {info[7]:.4f} at p = {int(info[8])}  (negative: pi(x) < li(x) throughout; Skewes)')
    check('skewes', info[7] < 0)
    # exact check of the smallest primes with mpmath (the float sweep's margin there is small)
    mp.mp.dps = 30
    worst = min(float(mp.li(p)) - (i + 1) for i, p in enumerate(small_primes(10 ** 5)))
    print(f'   exact (mpmath) min over primes p <= 10^5 of li(p) - pi(p): {worst:.6f}  (> 0)', check('skewes small', worst > 0))

    # normalised races on the log grid
    xs = grid_x
    lx = np.log(xs)
    sq = np.sqrt(xs)
    pi_g, p41_g, p43_g, p31_g, p32_g, li_g = (out_grid[:, i] for i in range(6))
    corr = power_corrections(xs, base)          # (c4, c3, c1) = prime-power terms k >= 2
    D4 = p43_g - p41_g                           # = -pi(x; chi_-4)
    D3 = p32_g - p31_g                           # = -pi(x; chi_-3)
    D4c = D4 - corr[:, 0]                        # = -Pi(x; chi_-4): prime powers counted with weight 1/k
    D3c = D3 - corr[:, 1]
    Sk = pi_g - li_g
    Skc = pi_g + corr[:, 2] - li_g               # Pi(x) - li(x)
    norm = lx / sq
    series = {'E4 = race mod 4 (primes only)': D4 * norm, 'E4 corrected (prime powers, weight 1/k)': D4c * norm,
              'E3 = race mod 3 (primes only)': D3 * norm, 'E3 corrected': D3c * norm,
              'Skewes S = pi - li': Sk * norm, 'Skewes corrected Pi - li': Skc * norm}
    print('\n   log-grid statistics on [10^3, 10^10] (step 1e-3 in log x; mean = logarithmic average):')
    stats = {}
    for name, v in series.items():
        pos = float(np.mean(v > 0))
        sc = int(np.sum(np.sign(v[1:]) * np.sign(v[:-1]) < 0))
        stats[name] = (float(np.mean(v)), float(np.std(v)), pos, sc)
        print(f'     {name:45s} mean {np.mean(v):+.4f}  sd {np.std(v):.4f}  share > 0: {pos:.5f}  grid sign changes: {sc}')
    # the square term alone, normalised: (log x/sqrt x) * (1/2) * #odd prime squares <= x  -> 1
    print('     (prime-square term alone: (log x/sqrt x) * corr_4(x) at x = 10^10 =', round(float(corr[-1, 0] * norm[-1]), 4), '; at 10^6 =',
          round(float(np.interp(math.log(1e6), ys, corr[:, 0] * norm)), 4), ')')
    sq4 = corr[:, 0] * norm
    sq1 = corr[:, 2] * norm
    print(f'     log-mean of the prime-power term alone: mod 4 {np.mean(sq4):+.4f}, mod 3 {np.mean(corr[:, 1] * norm):+.4f}, Skewes {np.mean(sq1):+.4f}'
          '  (= uncorrected mean - corrected mean, by construction; > 1 at finite x: (1/2)pi(sqrt x) log x/sqrt x = 1 + 2/log x + ...)')
    for nm, nc in (('E4 = race mod 4 (primes only)', 'E4 corrected (prime powers, weight 1/k)'), ('E3 = race mod 3 (primes only)', 'E3 corrected'),
                   ('Skewes S = pi - li', 'Skewes corrected Pi - li')):
        mu_u, _, pos_u, sc_u = stats[nm]
        mu_c, _, pos_c, sc_c = stats[nc]
        check('bias ' + nm, abs(mu_u) > 0.8 and abs(mu_c) < 0.25 and 0.3 < pos_c < 0.7 and sc_c > 20 * max(sc_u, 1))
    # 26861 details
    pr = small_primes(30000)
    c41 = sum(1 for p in pr if p <= 26861 and p % 4 == 1)
    c43 = sum(1 for p in pr if p <= 26861 and p % 4 == 3)
    sq_odd = sum(1 for p in pr if p > 2 and p * p <= 26861)
    print(f'\n   x = 26861: pi(x;4,1) = {c41}, pi(x;4,3) = {c43};  odd prime squares <= x: {sq_odd} (largest 163^2 = 26569);  sqrt(x)/log x = {math.sqrt(26861) / math.log(26861):.2f}')
    check('1473', c41 == 1473 and c43 == 1472)
    ties, a41, a43 = [], 0, 0
    for p in pr:
        if p > 26861:
            break
        a41 += p % 4 == 1
        a43 += p % 4 == 3
        if a41 == a43 and p > 2:
            ties.append(int(p))
    print(f'   primes p < 26861 at which pi(p;4,1) = pi(p;4,3) (ties; class 3 leads at every other prime >= 3): {[t for t in ties if t < 26861]}')

    print('\n== Part B: Liouville and Mertens to N_B = %.0e' % N_B)
    t1 = time.time()
    baseB = small_primes(int(N_B ** 0.5) + 10)
    ysB = np.arange(math.log(1e4), math.log(float(N_B)) - 1e-9, GRID_STEP)
    gxB = np.ceil(np.exp(ysB))
    ogB = np.zeros((len(gxB), 2))
    infoB = np.zeros(6)
    fr = np.zeros(5, dtype=np.int64)
    runsB = np.zeros((100000, 4), dtype=np.int64)
    sweep_liouville_mertens(N_B, SEG, baseB, gxB, ogB, infoB, fr, runsB)
    print('   sweep done in %.1f s' % (time.time() - t1))
    print(f'   L(N_B) = {int(infoB[0])}, M(N_B) = {int(infoB[1])}')
    print(f'   first n >= 2 with L(n) > 0: {fr[0]} (Tanaka 1980: 906150257)', check('polya', fr[0] == 906150257 or N_B < 906150257),
          f';  the first run ends at n = {fr[1]} (first n with L(n) <= 0 again);  max L over n <= N_B: {fr[2]} at n = {fr[3]};  runs with L > 0: {fr[4]}')
    clB = []
    for r in runsB[:min(int(fr[4]), runsB.shape[0])]:
        if clB and r[0] < clB[-1][1] * 1.01:
            c = clB[-1]
            clB[-1] = (c[0], int(r[1]), max(c[2], int(r[2])), c[3] if c[2] >= r[2] else int(r[3]), c[4] + 1)
        else:
            clB.append((int(r[0]), int(r[1]), int(r[2]), int(r[3]), 1))
    for c in clB:
        print(f'     excursion of L > 0: {c[0]} .. {c[1]} ({c[4]} runs), max L = {c[2]} at n = {c[3]}')
    print(f'   max |M(x)|/sqrt x on [10^4, 10^9]: {infoB[2]:.6f} at x = {int(infoB[3])};  L(x)/sqrt x ranges over [{infoB[5]:.4f}, {infoB[4]:.4f}]')
    Ln = ogB[:, 0] / np.sqrt(gxB)
    Mn = ogB[:, 1] / np.sqrt(gxB)
    print(f'   log-averages on [10^4, 10^9]: L(x)/sqrt x {np.mean(Ln):+.4f} (1/zeta(1/2) = {float(1 / mp.zeta(0.5)):+.4f});  M(x)/sqrt x {np.mean(Mn):+.4f} (0 expected)')
    print(f'   share of the log-grid with L(x) > 0: {np.mean(ogB[:, 0] > 0):.6f};  with |M(x)| > sqrt x: {np.mean(np.abs(Mn) > 1):.6f}')
    check('M bounded', infoB[2] < 1)
    json.dump({'stats': stats, 'excursions4': clusters, 'runs4': int(n4), 'polya': fr.tolist(), 'polya_excursions': clB, 'mertens_max': [float(infoB[2]), int(infoB[3])]},
              open('05-knowledge/results/prime_power_bias_races_20261014.json', 'w'), indent=1)
    print('\nFAILED CHECKS:', FAIL if FAIL else 'none')
    print('ALL CHECKS PASSED' if not FAIL else 'SOME CHECKS FAILED')
    print('total time %.1f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
