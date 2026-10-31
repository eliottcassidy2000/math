#!/usr/bin/env python3
"""Independent audit of THM-4853 (PRIMEGAME step counts).  Own code; nothing imported from the session's scripts.

Run from the repository root:
    python3 04-computation/experiments/groups_of_five_primegame_20261031/audit/audit_primegame.py

Checks
  A1  literal FRACTRAN on true integers (n % den == 0), both lists, first 25 primes; phase entries (the configuration
      right after 1/17) equal 2^N 7^(D(N)-1); the only pure powers of 2 are 2^prime; state-free configurations that
      are NOT phase entries exist (HYP-9448 wording); the owner's orbit reproduces OEIS A007542 (b-file 0..8102, hash).
  A2  per-trial block costs (non-divisor 6N+2q+2, divisor 4N+2N/d+2, end configurations), the state-free segment
      (L_G: N steps, L_O: N+m steps) and the reload (2N-1 steps), on the literal machine, 2 <= N <= 60, 1 <= d < N.
  A3  an accelerated counter-machine simulator (generic 2-cycle / self-loop jumps, own design), validated against A1,
      run to every prime < 1000 for both lists.
  A4  the closed forms cost_G, cost_O (fast O(spf) evaluation), against A3, the 34 OEIS data terms, and the OEIS
      b-files (A007546: 5000 terms, A. P. Heinz; A007547: 1000 terms, B. Dobbelaere), and the difference law.
  A5  the constant c = sum_p p^-2 prod_{q<p}(1-1/q) with a rigorous tail bound; the growth ratios; comparison with
      OEIS A032742's formula (A. Eldar, 2022: sum_{k<=n} A032742(k) ~ 0.165049... n^2).
"""
import hashlib
import math
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FAILS = []


def check(name, ok, info=''):
    print(f"  [{'OK' if ok else 'FAIL'}] {name}" + (f"  -- {info}" if info else ''))
    if not ok:
        FAILS.append(name)


HEAD = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17), (11, 13), (13, 11)]
LG = HEAD + [(15, 14), (15, 2), (55, 1)]   # the 15/14 list
LO = HEAD + [(15, 2), (1, 7), (55, 1)]     # the owner's list (Conway 1987 Theorem 1 per session126)
STATE_PRODUCT = 11 * 13 * 17 * 19 * 23 * 29

# OEIS JSON API data fields, retrieved 2026-10-31 (all listed terms)
OEIS_A007546 = [19, 69, 280, 707, 2363, 3876, 8068, 11319, 19201, 36866, 45551, 75224, 101112, 117831, 152025,
                215384, 293375, 327020, 428553, 507519, 555694, 700063, 808331, 989526, 1273490, 1434366, 1530213,
                1710923, 1818254, 2019962, 2833089, 3104685, 3546320, 3720785]
OEIS_A007547 = [19, 69, 281, 710, 2375, 3893, 8102, 11361, 19268, 36981, 45680, 75417, 101354, 118093, 152344,
                215797, 293897, 327571, 429229, 508284, 556494, 701008, 809381, 990746, 1274952, 1435957, 1531854,
                1712701, 1820085, 2021938, 2835628, 3107393, 3549288, 3723821]
# sha256 of the A007542 b-file (n = 0..8102) after keeping only lines "n a(n)" joined by '\n'
A007542_BFILE_NORMALIZED_SHA256 = 'a7359b7a0e7c068b656a16c2e1a047ab1e23fdb9ff32281af18186948609aab2'


def sieve_spf(n):
    spf = list(range(n + 1))
    for i in range(2, int(n ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


SPF = sieve_spf(60000)


def D(n):
    return 1 if n == 1 else n // SPF[n]


def primes_below(n):
    return [p for p in range(2, n) if SPF[p] == p]


# ------------------------------------------------------------------ A1 literal integer FRACTRAN
def literal(prog, nprimes, keep_orbit=0):
    """returns emissions [(p, step)], phase entries [(step, value)] (value right after 1/17),
    the first 8 state-free configurations, and optionally the orbit prefix."""
    n, step = 2, 0
    em, entries, statefree, orbit = [], [], [], [2]
    i17 = prog.index((1, 17))
    while len(em) < nprimes or len(orbit) < keep_orbit:
        for i, (a, b) in enumerate(prog):
            if n % b == 0:
                n = n // b * a
                break
        step += 1
        if len(orbit) < keep_orbit:
            orbit.append(n)
        if n & (n - 1) == 0 and len(em) < nprimes:
            em.append((n.bit_length() - 1, step))
        if i == i17:
            entries.append((step, n))
        if len(statefree) < 8 and math.gcd(n, STATE_PRODUCT) == 1:
            statefree.append((step, n))
    return em, entries, statefree, orbit


def A1():
    print('=' * 100)
    print('A1  literal integer FRACTRAN (both lists, first 25 primes)')
    print('=' * 100)
    out = {}
    for name, prog in (('L_G (15/14 list)', LG), ('L_O (owner list)', LO)):
        t0 = time.time()
        em, entries, sf, _ = literal(prog, 25)
        out[name] = (em, entries)
        P = primes_below(100)
        check(f'{name}: pure powers of 2 after step 0 are exactly 2^p, p the first 25 primes, in order',
              [p for p, _ in em] == P[:25], f'{time.time() - t0:.1f}s, last step {em[-1][1]}')
        ok = all(v == 2 ** N * 7 ** (D(N) - 1) for (N, (_, v)) in zip(range(2, 10 ** 9), entries))
        check(f'{name}: the k-th configuration produced by 1/17 is 2^N 7^(D(N)-1), N = k+1 (all {len(entries)} '
              f'phase entries up to the last emission)', ok)
        print(f'     first state-free configurations (step, value): {sf}')
        check(f'{name}: a state-free configuration that is not a phase entry occurs (step 1: 15 = 3*5)',
              sf[0] == (1, 15))
    check('L_G steps = OEIS A007546 (first 25 terms)', [s for _, s in out['L_G (15/14 list)'][0]] == OEIS_A007546[:25])
    check('L_O steps = OEIS A007547 (first 25 terms)', [s for _, s in out['L_O (owner list)'][0]] == OEIS_A007547[:25])
    # A007542 orbit (owner list) vs the OEIS b-file 0..8102
    _, _, _, orb = literal(LO, 1, keep_orbit=8103)
    lines = '\n'.join(f'{i} {v}' for i, v in enumerate(orb[:8103]))
    h = hashlib.sha256(lines.encode()).hexdigest()
    check('owner-list orbit a(0..8102) reproduces the OEIS A007542 b-file exactly (normalized sha256)',
          h == A007542_BFILE_NORMALIZED_SHA256, f'a(130..133) = {orb[130:134]}')
    _, _, _, orbg = literal(LG, 1, keep_orbit=140)
    first_diff = next(i for i in range(140) if orb[i] != orbg[i])
    print(f'     first index where the two orbits differ: {first_diff} (owner {orb[first_diff]}, 15/14 list '
          f'{orbg[first_diff]}); OEIS A007542 says Guy\'s version differs "starting from n=132"')
    return out


# ------------------------------------------------------------------ the counter machine (own derivation from the list)
REG = [2, 3, 5, 7]
STATES = [11, 13, 17, 19, 23, 29]


def factor_over(x, primes):
    e = []
    for p in primes:
        k = 0
        while x % p == 0:
            x //= p
            k += 1
        e.append(k)
    assert x == 1
    return e


def machine(prog):
    """state (0 = none) -> ordered list of (need-mask over registers as tuple of indices, delta, next state, index)."""
    rules = []
    for idx, (a, b) in enumerate(prog):
        ea, eb = factor_over(a, REG + STATES), factor_over(b, REG + STATES)
        sin = [STATES[i] for i in range(6) if eb[4 + i]]
        sout = [STATES[i] for i in range(6) if ea[4 + i]]
        assert len(sin) <= 1 and len(sout) <= 1 and max(eb[:4] + ea[:4]) <= 1
        need = tuple(i for i in range(4) if eb[i])
        delta = tuple(ea[i] - eb[i] for i in range(4))
        rules.append((sin[0] if sin else 0, need, delta, sout[0] if sout else 0, idx))
    tab = {}
    for s in [0] + STATES:
        # FRACTRAN priority: a rule can fire in state s iff its denominator's state part is s or empty
        tab[s] = [(need, delta, nxt, idx) for (sin, need, delta, nxt, idx) in rules if sin in (s, 0)]
    return tab


def first(tab, s, r):
    for k, (need, delta, nxt, idx) in enumerate(tab[s]):
        if all(r[j] >= 1 for j in need):
            if s != 0 and nxt != 0 and idx >= 11:
                raise RuntimeError('state-free rule fired with a state present')
            return k
    raise RuntimeError('no rule')


def run_trial(prog, N, d):
    """literal counter machine from state 13, (0,0,N,d); stops at state 13 with r7 = d-1 (non-divisor) or at the first
    state-free configuration (divisor). returns (steps, final state, final registers)."""
    tab = machine(prog)
    s, r, steps = 13, [0, 0, N, d], 0
    while True:
        k = first(tab, s, r)
        need, delta, nxt, idx = tab[s][k]
        r = [r[j] + delta[j] for j in range(4)]
        s = nxt
        steps += 1
        if s == 0:
            return steps, s, r
        if s == 13 and r[3] == d - 1 and r[0] == 0 and r[1] == 0:
            return steps, s, r
        if steps > 10 ** 6:
            raise RuntimeError('runaway')


def run_segment(prog, N, m):
    """from state-free 2^(N-1) 7^m: count steps to state 11 (segment) and then to state 13 (reload)."""
    tab = machine(prog)
    s, r, steps = 0, [N - 1, 0, 0, m], 0
    seg = None
    while True:
        k = first(tab, s, r)
        need, delta, nxt, idx = tab[s][k]
        r = [r[j] + delta[j] for j in range(4)]
        s = nxt
        steps += 1
        if s == 11 and seg is None:
            seg = (steps, tuple(r))
        if s == 13:
            return seg, (steps - seg[0], tuple(r))


def A2():
    print('=' * 100)
    print('A2  block costs on the literal counter machine (2 <= N <= 60, 1 <= d < N), segment and reload')
    print('=' * 100)
    bad = 0
    pairs = 0
    for prog in (LG, LO):
        for N in range(2, 61):
            for d in range(1, N):
                pairs += 1
                steps, s, r = run_trial(prog, N, d)
                q = N // d
                if N % d:
                    ok = (s == 13 and r == [0, 0, N, d - 1] and steps == 6 * N + 2 * q + 2)
                else:
                    ok = (s == 0 and r == [N, 0, 0, d - 1] and steps == 4 * N + 2 * q + 2)
                bad += not ok
    check('trial from state 13, (0,0,N,d): non-divisor -> (0,0,N,d-1) in 6N+2floor(N/d)+2 steps; '
          'divisor -> state-free 2^N 7^(d-1) in 4N+2N/d+2 steps (both lists)', bad == 0, f'{pairs} (list,N,d) cases')
    bad = 0
    for N in range(2, 61):
        m = D(N - 1) - 1
        (sg, rg), (rlg, rrg) = run_segment(LG, N, m)
        (so, ro), (rlo, rro) = run_segment(LO, N, m)
        ok = (sg == N and so == N + m and rg == ro == (0, N - 1, N, 0) and rlg == rlo == 2 * N - 1
              and rrg == rro == (0, 0, N, N - 1))
        bad += not ok
    check('state-free segment from 2^(N-1) 7^(D(N-1)-1): L_G N steps, L_O N+m steps, both to 3^(N-1)5^N 11; '
          'reload 2N-1 steps to state 13 with (0,0,N,N-1)', bad == 0, 'N = 2..60')
    # any m <= N-1 works for L_G (needs r7 <= r2); show the segment formula is about m = D(N-1)-1 only through m
    (sg, _), _ = run_segment(LG, 7, 6)
    (so, _), _ = run_segment(LO, 7, 6)
    print(f'     (the segment counts depend on m only: N=7, m=6 gives L_G {sg} = N, L_O {so} = N+m)')


# ------------------------------------------------------------------ A3 accelerated simulator (own design)
def accelerated(prog, pmax, record_entries=False):
    tab = machine(prog)
    s, r, steps = 0, [1, 0, 0, 0], 0
    em, entries = [], []
    while True:
        k = first(tab, s, r)
        need, delta, nxt, idx = tab[s][k]
        jumped = False
        if nxt == s:
            # self-loop: largest K with rule k first-applicable for K consecutive applications
            K = None
            for j in need:
                if delta[j] < 0:
                    kk = (r[j] - 1) // (-delta[j]) + 1
                    K = kk if K is None else min(K, kk)
            for need2, delta2, _, _ in tab[s][:k]:
                zero = [j for j in need2 if r[j] == 0]
                if not any(delta[j] <= 0 for j in zero):
                    K = 1
            if K is None:
                raise RuntimeError('non-terminating self-loop')
            if K > 1:
                # a pure power of 2 strictly inside the run would need r3=r5=r7=0 at some m < K
                ms = set()
                okm = True
                for j in (1, 2, 3):
                    if delta[j] == 0:
                        if r[j] != 0:
                            okm = False
                    else:
                        if (-r[j]) % delta[j] == 0 and (-r[j]) // delta[j] >= 1:
                            ms.add((-r[j]) // delta[j])
                        else:
                            okm = False
                if okm and s == 0 and len(ms) <= 1 and ms and min(ms) < K:
                    raise RuntimeError('pure power inside an accelerated self-loop (not expected)')
                r = [r[j] + (K - 1) * delta[j] for j in range(4)]
                steps += K - 1
                # the last application is done as an ordinary step below
        else:
            mid = [r[j] + delta[j] for j in range(4)]
            k2 = first(tab, nxt, mid)
            need2, delta2, nxt2, idx2 = tab[nxt][k2]
            if nxt2 == s and s != 0 and nxt != 0:
                Dv = [delta[j] + delta2[j] for j in range(4)]
                K = None
                for j in need:          # guard of rule k at r + m Dv, m = 0..K-1
                    if Dv[j] < 0:
                        kk = (r[j] - 1) // (-Dv[j]) + 1
                        K = kk if K is None else min(K, kk)
                for j in need2:         # guard of rule k2 at mid + m Dv
                    if Dv[j] < 0:
                        kk = (mid[j] - 1) // (-Dv[j]) + 1
                        K = kk if K is None else min(K, kk)
                for base, st, kk0 in ((r, s, k), (mid, nxt, k2)):
                    for needE, _, _, _ in tab[st][:kk0]:
                        zero = [j for j in needE if base[j] == 0]
                        if any(Dv[j] <= 0 for j in zero):
                            continue
                        # earlier rule becomes applicable at the first m >= 1 where all its registers are >= 1
                        mE = max((0 if base[j] >= 1 else -(-(1 - base[j]) // Dv[j])) for j in needE)
                        K = mE if K is None else min(K, mE)
                if K is None:
                    raise RuntimeError('non-terminating 2-cycle')
                if K > 1:
                    r = [r[j] + (K - 1) * Dv[j] for j in range(4)]
                    steps += 2 * (K - 1)
                    jumped = True
        # ordinary single step (re-evaluate the first applicable rule after any jump)
        k = first(tab, s, r)
        need, delta, nxt, idx = tab[s][k]
        r = [r[j] + delta[j] for j in range(4)]
        s = nxt
        steps += 1
        if record_entries and idx == 8:
            entries.append((steps, tuple(r)))
        if s == 0 and r[0] >= 1 and r[1] == r[2] == r[3] == 0:
            em.append((r[0], steps))
            if r[0] >= pmax:
                return em, entries


def A3(lit):
    print('=' * 100)
    print('A3  accelerated counter-machine simulator: validation against A1, then every prime < 1000')
    print('=' * 100)
    res = {}
    for name, prog in (('L_G (15/14 list)', LG), ('L_O (owner list)', LO)):
        em, entries = accelerated(prog, 97, record_entries=True)
        lem, lentries = lit[name]
        conv = [(st, tuple(factor_over(v, REG + STATES)[:4])) for st, v in lentries]
        check(f'{name}: accelerated = literal at every phase entry (step and configuration) and every emission, '
              f'first 25 primes', em == lem and entries == conv[:len(entries)] and len(entries) == len(conv),
              f'{len(entries)} phase entries')
        t0 = time.time()
        em1000, _ = accelerated(prog, 997)
        res[name] = em1000
        check(f'{name}: emits exactly the 168 primes < 1000 in order', [p for p, _ in em1000] == primes_below(1000),
              f'{time.time() - t0:.1f}s; step at 2^997 = {em1000[-1][1]}')
    return res


# ------------------------------------------------------------------ A4 closed forms
def closed_forms(Xmax):
    """returns dict p -> (s_G, s_O) for all primes p <= Xmax, with cost_G evaluated in O(spf) or via the divisor
    summatory function for primes."""
    tau = [0] * (Xmax + 1)
    for dd in range(1, Xmax + 1):
        for mlt in range(dd, Xmax + 1, dd):
            tau[mlt] += 1
    S = 0
    sG = sO = 0
    out = {}
    for N in range(1, Xmax + 1):
        S += tau[N]                      # S = sum_{d<=N} floor(N/d)
        if N == 1:
            continue
        p = SPF[N]
        a = N // p
        if p == N:
            F = S - N - 1                # sum_{d=2}^{N-1} floor(N/d)
        else:
            F = sum(N // k for k in range(1, p)) - (p - 1) * a - 1
        cG = 7 * N + 1 + 2 * p + (N - 1 - a) * (6 * N + 2) + 2 * F
        sG += cG
        sO += cG + D(N - 1) - 1
        if p == N:
            out[N] = (sG, sO)
    return out


def cost_G_brute(N):
    a = D(N)
    return 7 * N + 1 + 2 * SPF[N] + (N - 1 - a) * (6 * N + 2) + 2 * sum(N // d for d in range(a + 1, N))


def read_bfile(fn):
    vals = {}
    with open(fn, encoding='utf-8', errors='replace') as f:
        for line in f:
            m = re.fullmatch(r'(\d+) (\d+)', line.strip())
            if m:
                vals[int(m.group(1))] = int(m.group(2))
    return [vals[i] for i in range(1, len(vals) + 1)]


def A4(sim):
    print('=' * 100)
    print('A4  closed forms vs simulation, OEIS data fields and OEIS b-files; difference law')
    print('=' * 100)
    # fast evaluation = literal formula
    acc_fast = closed_forms(2000)
    accG = 0
    okf = True
    for N in range(2, 2001):
        accG += cost_G_brute(N)
        if SPF[N] == N and acc_fast[N][0] != accG:
            okf = False
    check('fast closed-form evaluator = the formula of THM-4853 (a) summed term by term, N <= 2000', okf)
    X = 48611
    t0 = time.time()
    cf = closed_forms(X)
    P = sorted(cf)
    print(f'     closed forms evaluated for all {len(P)} primes <= {X} in {time.time() - t0:.1f}s')
    g1000 = [s for _, s in sim['L_G (15/14 list)']]
    o1000 = [s for _, s in sim['L_O (owner list)']]
    check('closed forms = accelerated simulator, all 168 primes < 1000, both lists',
          [cf[p][0] for p in P[:168]] == g1000 and [cf[p][1] for p in P[:168]] == o1000,
          f'p=997: s_G {g1000[-1]}, s_O {o1000[-1]}, diff {o1000[-1] - g1000[-1]}')
    check('closed form s_G = OEIS A007546 data (34 terms)', [cf[p][0] for p in P[:34]] == OEIS_A007546)
    check('closed form s_O = OEIS A007547 data (34 terms)', [cf[p][1] for p in P[:34]] == OEIS_A007547)
    b46 = read_bfile(os.path.join(HERE, 'b007546.txt'))
    b47 = read_bfile(os.path.join(HERE, 'b007547.txt'))
    check('b-files start with the JSON data fields', b46[:34] == OEIS_A007546 and b47[:34] == OEIS_A007547,
          f'b007546 has {len(b46)} terms, b007547 has {len(b47)} terms')
    mism46 = [k + 1 for k in range(len(b46)) if cf[P[k]][0] != b46[k]]
    mism47 = [k + 1 for k in range(len(b47)) if cf[P[k]][1] != b47[k]]
    check(f'closed form s_G = A007546 b-file, all {len(b46)} terms (p_k up to {P[len(b46) - 1]})', not mism46,
          f'mismatches {mism46[:5]}')
    check(f'closed form s_O = A007547 b-file, all {len(b47)} terms (p_k up to {P[len(b47) - 1]})', not mism47,
          f'mismatches {mism47[:5]}')
    law = []
    acc = 0
    n = 2
    for p in P[:len(b47)]:
        while n <= p - 1:
            acc += D(n) - 1
            n += 1
        law.append(acc)
    check('difference law A007547(k) - A007546(k) = sum_{n=2}^{p_k-1}(D(n)-1) on the b-files, k <= 1000',
          [b47[k] - b46[k] for k in range(len(b47))] == law, f'k=1000: {b47[-1] - b46[999]}')
    return cf, P


# ------------------------------------------------------------------ A5 constants and growth
def A5(cf, P):
    print('=' * 100)
    print('A5  the constant c and the growth claims')
    print('=' * 100)
    LIM = 10 ** 7
    sv = bytearray([1]) * (LIM + 1)
    sv[0:2] = b'\x00\x00'
    for i in range(2, int(LIM ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = bytearray(len(range(i * i, LIM + 1, i)))
    c = 0.0
    comp = 0.0
    prod = 1.0
    for p in range(2, LIM + 1):
        if sv[p]:
            term = prod / (p * p)
            y = term - comp          # Kahan summation
            t = c + y
            comp = (t - c) - y
            c = t
            prod *= 1.0 - 1.0 / p
    tail = prod / LIM               # sum_{p>LIM} p^-2 prod_{q<p}(1-1/q) < prod_{q<=LIM}(1-1/q) * sum_{n>LIM} n^-2
    print(f'     c (primes <= 1e7) = {c:.10f}; tail bound < {tail:.2e}; so c = {c:.8f}..{c + tail:.8f}')
    check('c = 0.330098... (THM-4853)', abs(c - 0.330098) < 5e-7)
    check('c/2 = 0.165049... = the constant of OEIS A032742 (Eldar 2022: sum_{k<=n} a(k) ~ c n^2)',
          abs(c / 2 - 0.165049) < 5e-7, f'c/2 = {c / 2:.9f}')
    m = 10 ** 6
    mean = sum(1.0 / SPF_big(k) for k in range(2, m + 1)) / m
    print(f'     mean of 1/spf(n), n <= 1e6: {mean:.6f} (c = {c:.6f})')
    for p in (101, 499, 997, 2999, 7919, 48611):
        sG, sO = cf[p]
        law = sO - sG
        print(f'     p = {p:>5}: s_G / (2(1-c)p^3) = {sG / (2 * (1 - c) * p ** 3):.5f};  '
              f'(s_O - s_G) / ((c/2)p^2) = {law / (c / 2 * p ** 2):.5f};  share (s_O-s_G)/s_O = {law / sO:.3e}')
    sG, sO = cf[2999]
    check('ratio at p = 2999 is 1.0015 (THM-4853)', round(sG / (2 * (1 - c) * 2999 ** 3), 4) == 1.0015,
          f's_G(2999) = {sG}')
    sG, sO = cf[997]
    check('law(997)/((c/2) 997^2) = 0.9932 and share 1.2e-4 (THM-4853)',
          round((sO - sG) / (c / 2 * 997 ** 2), 4) == 0.9932 and round((sO - sG) / sO, 5) == 0.00012,
          f'law(997) = {sO - sG}')


_SPF_BIG = None


def SPF_big(k):
    global _SPF_BIG
    if _SPF_BIG is None:
        _SPF_BIG = sieve_spf(10 ** 6)
    return _SPF_BIG[k]


def main():
    t0 = time.time()
    lit = A1()
    A2()
    sim = A3(lit)
    cf, P = A4(sim)
    A5(cf, P)
    print()
    print(f'total {time.time() - t0:.0f}s')
    if FAILS:
        print('FAILED:', FAILS)
        return 1
    print('ALL AUDIT CHECKS PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
