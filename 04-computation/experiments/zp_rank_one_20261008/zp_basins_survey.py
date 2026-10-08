"""Cycles, basins and basin boundaries of the base-p Collatz maps (session opus-2026-10-08-S22).
Note: 05-knowledge/results/zp_rank_one_coalescence_20261008.md (section 3).

C_p(x) = x/p if p | x, else ((p+1)x + p - (x mod p))/p (= ceil((p+1)x/p)), on positive integers.
  (S1) for every prime p <= 31: the cycles reached by starts n <= 10^6 and their basin shares.
  (S2) for p = 3, 11, 17, 23: basin indicator b(n) for n <= 10^7; share s of the second basin and boundary density
       #{n < N : b(n) != b(n+1)}/N at N = 10^3..10^7.  THM-4610 5(e) proves the boundary density tends to 0; here it is
       compared with 2s(1-s), its value if n and n+1 landed independently: the ratio q_eff = boundary/(2s(1-s)) is the
       share of consecutive pairs still unmerged when they reach the cycles (NUMERICAL; compared with the Haar merge tail
       at the typical orbit length in zp_tails.py).
  (S3) for p = 3: the fraction of n <= 10^6 whose orbit merges with that of n+1 at equal time before either enters a
       cycle (density-one merging, THM-4610 5(a), at finite size).
Prints ALL CHECKS PASSED.  Runtime about 2 minutes."""
import sys, os, math, time

FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def C(p, x):
    i = x % p
    return x // p if i == 0 else ((p + 1) * x + p - i) // p


def survey(p, N):
    term = {}; cycmin = {}; share = {}
    lim = 50 * N
    for n in range(1, N + 1):
        if n in term:
            share[term[n]] = share.get(term[n], 0) + 1
            continue
        path = []; pos = {}; x = n
        while x not in term and x not in pos:
            pos[x] = len(path); path.append(x); x = C(p, x)
        if x in term:
            m = term[x]
        else:
            c = path[pos[x]:]; m = min(c); cycmin[m] = len(c)
        for z in path:
            if z <= lim:
                term[z] = m
        share[m] = share.get(m, 0) + 1
    return {m: (cycmin.get(m), share[m] / N) for m in sorted(share)}


def basin_array(p, N):
    term = bytearray(N + 2)            # 1 = trivial cycle, 2 = other cycle(s)
    cyc_cache = {}

    def cyc_min(x):
        seen = set()
        while x not in seen:
            seen.add(x); x = C(p, x)
        c = [x]; z = C(p, x)
        while z != x:
            c.append(z); z = C(p, z)
        return min(c)
    for n in range(1, N + 1):
        if term[n]:
            continue
        path = []; x = n
        while True:
            if x <= N and term[x]:
                t = term[x]; break
            if x <= 10 ** 4:
                if x not in cyc_cache:
                    cyc_cache[x] = cyc_min(x)
                t = 1 if cyc_cache[x] == 1 else 2
                break
            path.append(x); x = C(p, x)
        for z in path:
            if z <= N:
                term[z] = t
        term[n] = t
    return term


if __name__ == '__main__':
    t0 = time.time()
    print('(S1) cycles and basin shares, starts <= 10^6')
    table = {}
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        table[p] = survey(p, 10 ** 6)
        print(f'     p = {p}: ' + ', '.join(f'min {m} (length {L}): {s:.5f}' for m, (L, s) in table[p].items()), flush=True)
    single = [p for p in table if list(table[p]) == [1]]
    double = [p for p in table if len(table[p]) == 2]
    check(single == [5, 7, 13, 19] and double == [3, 11, 17, 23, 29, 31] and table[11][642][0] == 57,
          'single trivial cycle for p = 5, 7, 13, 19; a second cycle with a positive basin share for p = 3, 11, 17, 23, 29, 31 '
          '(p = 11: minimum 642, length 57)')
    print('(S2) basin boundaries, n <= 10^7')
    ok2 = True; qeffs = {}
    for p in (3, 11, 17, 23):
        b = basin_array(p, 10 ** 7)
        rows = []
        for j in range(3, 8):
            M = 10 ** j
            sec = b[1:M + 1].count(2) / M
            bnd = sum(1 for n in range(1, M) if b[n] != b[n + 1]) / M
            rows.append((M, round(sec, 5), round(bnd, 5), round(bnd * math.sqrt(math.log(M)), 4)))
        s_ = rows[-1][1]; q_eff = rows[-1][2] / (2 * s_ * (1 - s_))
        print(f'     p = {p}: (N, second-basin share, boundary density, boundary x sqrt(ln N)) = {rows}; '
              f'q_eff = boundary/(2s(1-s)) at 10^7 = {q_eff:.3f}', flush=True)
        qeffs[p] = q_eff
        ok2 &= abs(rows[-1][1] - rows[-2][1]) < 0.005 and rows[-1][1] > 0.01 and 0 < q_eff < 1
    check(ok2 and qeffs[3] < qeffs[11] < qeffs[17] and qeffs[3] < 0.7 and min(qeffs[11], qeffs[17], qeffs[23]) > 0.8,
          'second-basin shares are stable (10^6 vs 10^7) and positive; the boundary density at 10^7 is a fraction q_eff < 1 of '
          '2s(1-s), about 0.58 for p = 3 and above 0.8 for p = 11, 17, 23 (pairs mostly unmerged at these sizes) (NUMERICAL)')
    print('(S3) consecutive pairs merging at equal time before a cycle, p = 3')
    p = 3; N = 10 ** 6
    cyc_elems = set()
    for x0 in (1, 7):
        z = x0
        while True:
            cyc_elems.add(z); z = C(p, z)
            if z == x0:
                break
    merged = 0; diff_basin = 0
    b3 = basin_array(3, N + 1)
    for n in range(1, N + 1):
        a, b2 = n, n + 1
        for _ in range(2000):
            if a in cyc_elems or b2 in cyc_elems:
                break
            if a == b2:
                merged += 1; break
            a, b2 = C(p, a), C(p, b2)
        diff_basin += b3[n] != b3[n + 1]
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import zp_chain_checks as Z
    import random
    lam3 = (1 / 3) * math.log(1 / 3) + (2 / 3) * math.log(4 / 3)
    T6 = round(math.log(10 ** 6) / abs(lam3))
    rng = random.Random(6); alive = 0; paths = 40000
    for _ in range(paths):
        k, A = 0, 1
        for t in range(T6):
            if k == 0 and A == 0:
                break
            k, A = Z.istep(3, k, A, rng.randrange(3))
        alive += not (k == 0 and A == 0)
    q6 = alive / paths
    print(f'     n <= 10^6: merged at equal time above the cycles {merged / N:.4f}; different basins {diff_basin / N:.4f}; '
          f'Haar prediction 1 - q_3(T) at the typical orbit length T = ln(10^6)/|Lambda_3| = {T6}: {1 - q6:.4f} ({paths} paths)')
    check(0 < merged / N < 1 - diff_basin / N and abs(merged / N - (1 - q6)) < 0.08,
          'the share of consecutive pairs merging above the cycles is within 0.08 of the Haar prediction at the typical orbit length '
          '(finite-size face of THM-4610 5(a)); pairs in different basins never merge (NUMERICAL)')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
