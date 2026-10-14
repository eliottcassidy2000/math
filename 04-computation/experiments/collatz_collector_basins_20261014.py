#!/usr/bin/env python3
"""The Collatz collector: residues mod p of px+-1 trajectories, the 2:1 visit bias mod 3, the balanced backward tree,
the three basins of 3x-1 (a Collatz race), the integer-cycle census for both signs, and the carry budget.

Session opus-2026-10-14-S32.  Note: 05-knowledge/results/prime_power_bias_collatz_collector_20261014.md

Maps (shortcut form), s = +1 or -1:   T_s(n) = n/2 (n even),  (p n + s)/2 (n odd).   Collatz: p = 3, s = +1;  "3x-1": s = -1.
Sections:
 1. Exact law of T_s^k(n) mod p for n uniform mod p 2^k (brute force) against the closed form
        P_k(r) = sum_{j<k} 2^-(j+1) [r = s 2^-(j+1)] + 2^-k / p,
    limit: mass 2^(d-m)/(2^d - 1) on s 2^-m (m = 1..d, d = ord_p 2), zero elsewhere.
 2. The Syracuse map n -> (3n+1)/2^v (odd to odd): law mod 3 is (1/3, 2/3) on (1, 2).
 3. The 3-adic visit law: stationary law of T_+ mod 3^m (m <= 6) from the fair-coin chain, exact rationals.
 4. Trajectory visit counts mod 3 for n <= 10^7 under T_+ and T_-.
 5. The backward tree of T_+ from 1: layer counts by residue mod 3 to depth 50 (balanced).
 6. The three basins of T_- on 1..N_C (N_C = 10^9): every n reaches {1}, {5,7,10} or the 11-cycle through 17; counts,
    densities, and every N at which the order of the three counts changes (the race).
 7. Integer cycles of T_+ on Z for every word of length L <= 24 (both signs; T_- cycles are their negatives).
 8. The carry budget: for l odd steps and L total steps, a cycle's least odd element n_min satisfies
        n_min <= 1/(3 (e^(|delta|/l) - 1))  (above the line, s = +1),   n_min <= 1/(3 (1 - e^(-|delta|/l)))  (below, s = -1),
    delta = L log 2 - l log 3;  table for l <= 12 at the two nearest L, with the cycles found.

Reproduce: python3 04-computation/experiments/collatz_collector_basins_20261014.py [N_C]
"""
import math, sys, time
from fractions import Fraction as Fr
import numpy as np
import numba as nb
import sympy as sp

N_C = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 9
FAIL = []


def check(label, ok):
    if not ok:
        FAIL.append(label)
    return ok


@nb.njit(cache=True)
def law_bruteforce(p, s, k):
    """counts[r] = #{0 <= n < p 2^k : T_s^k(n) = r mod p}"""
    M = p * (1 << k)
    counts = np.zeros(p, dtype=np.int64)
    for n0 in range(M):
        x = n0
        for step in range(k):
            x = x % (p << (k - step))          # T^k(n) mod p only depends on x mod p 2^(k-step): no overflow
            if x % 2 == 0:
                x //= 2
            else:
                x = (p * x + s) // 2
        counts[x % p] += 1
    return counts


def law_formula(p, s, k):
    inv2 = pow(2, -1, p)
    cnt = [1] * p                      # the no-odd-step term: 2^-k/p * p 2^k = 1 per class
    for j in range(k):
        r = (s * pow(inv2, j + 1, p)) % p
        cnt[r] += p * 2 ** (k - j - 1)
    return cnt


def limit_law(p, s):
    d = sp.n_order(2, p)
    inv2 = pow(2, -1, p)
    return {(s * pow(inv2, m, p)) % p: Fr(2 ** (d - m), 2 ** d - 1) for m in range(1, d + 1)}, d


@nb.njit(cache=True)
def syracuse_law(k):
    """law of Syr(n) mod 3 over odd n mod 3 2^(k+1) with v_2(3n+1) <= k (exact part)"""
    M = 3 * (1 << (k + 1))
    counts = np.zeros(3, dtype=np.int64)
    tail = 0
    for n in range(1, M, 2):
        y = 3 * n + 1
        v = 0
        while y % 2 == 0 and v <= k:
            y //= 2
            v += 1
        if v > k:
            tail += 1
            continue
        counts[y % 3] += 1
    return counts, tail


@nb.njit(cache=True)
def visit_counts(Nmax, s, cyc_mark):
    """residue-mod-3 counts of all values T_s^k(n), k >= 1, until the orbit enters a marked cycle element, n = 1..Nmax"""
    counts = np.zeros(3, dtype=np.int64)
    first_odd_seen = np.zeros(3, dtype=np.int64)    # counts restricted to values after the first odd step
    for n in range(1, Nmax + 1):
        x = n
        seen_odd = False
        steps = 0
        while True:
            if x <= 200 and cyc_mark[x] == 1 and steps > 0:
                break
            if x % 2 == 0:
                x //= 2
            else:
                x = (3 * x + s) // 2
                seen_odd = True
            steps += 1
            counts[x % 3] += 1
            if seen_odd:
                first_odd_seen[x % 3] += 1
            if steps > 100000:
                break
    return counts, first_odd_seen


@nb.njit(cache=True)
def basins(N, memo, race_n, race_order, checkpoints, cp_counts):
    """memo[n] in {1,2,3}: the T_- cycle reached from n.  Cycles: {1}; {5,7,10}; {17,25,37,55,82,41,61,91,136,68,34}."""
    c = np.zeros(4, dtype=np.int64)
    nrace = 0
    prev_key = -1
    cpi = 0
    maxval = 0
    for n in range(1, N + 1):
        if memo[n] == 0:
            x = n
            while x >= n and memo[x if x <= 136 else 0] == 0 or (x >= n and x > 136):
                if x % 2 == 0:
                    x //= 2
                else:
                    x = (3 * x - 1) // 2
                if x > maxval:
                    maxval = x
                if x < n or (x <= 136 and memo[x] != 0):
                    break
            memo[n] = memo[x]
        b = memo[n]
        c[b] += 1
        # order key of the three counts (descending, ties broken by id)
        a1, a2, a3 = c[1], c[2], c[3]
        key = 0
        # rank encoding: which id is first, second
        if a1 >= a2 and a1 >= a3:
            first = 1
        elif a2 >= a1 and a2 >= a3:
            first = 2
        else:
            first = 3
        if first == 1:
            second = 2 if a2 >= a3 else 3
        elif first == 2:
            second = 1 if a1 >= a3 else 3
        else:
            second = 1 if a1 >= a2 else 2
        key = first * 10 + second
        if key != prev_key:
            if nrace < race_n.shape[0]:
                race_n[nrace] = n
                race_order[nrace] = key
            nrace += 1
            prev_key = key
        while cpi < checkpoints.shape[0] and checkpoints[cpi] == n:
            cp_counts[cpi, 0] = c[1]
            cp_counts[cpi, 1] = c[2]
            cp_counts[cpi, 2] = c[3]
            cpi += 1
    return nrace, maxval


@nb.njit(cache=True)
def cycle_census(Lmax, out_vals, out_L):
    """all integer x0 = c_w/(2^L - 3^l) over binary words of length L <= Lmax"""
    nout = 0
    for L in range(1, Lmax + 1):
        twoL = np.int64(1) << L
        for w in range(1 << L):
            l = 0
            c = np.int64(0)
            # c = sum_k w_k 3^(ones after k) 2^k: scan from the last letter backwards
            p3 = np.int64(1)
            for k in range(L - 1, -1, -1):
                if (w >> k) & 1:
                    c += p3 * (np.int64(1) << k)
                    p3 *= 3
                    l += 1
            D = twoL - p3
            if D != 0 and c % D == 0:
                if nout < out_vals.shape[0]:
                    out_vals[nout] = c // D
                    out_L[nout] = L
                nout += 1
    return nout


def main():
    t0 = time.time()
    print('== 1. Exact law of T_s^k(n) mod p (n uniform mod p 2^k) against the closed form')
    for p in (3, 5, 7, 11, 13, 17, 23, 31):
        for s in (1, -1):
            ok = True
            for k in (1, 2, 3, 5, 8, 12):
                if p * 2 ** k > 6 * 10 ** 6:
                    continue
                bf = law_bruteforce(p, s, k)
                ok &= list(bf) == law_formula(p, s, k)
            lim, d = limit_law(p, s)
            check(f'law {p} {s}', ok)
            supp = sorted(lim)
            qr = sorted({x * x % p for x in range(1, p)})
            tag = ''
            if supp == qr:
                tag = '  [support = the quadratic residues]'
            elif sorted((-x) % p for x in qr) == supp:
                tag = '  [support = the non-residues]'
            elif len(supp) == p - 1:
                tag = '  [support = all units: 2 is a primitive root]'
            print(f'   p = {p:2d}, s = {s:+d}: brute force = formula for k <= 12: {ok};  ord_p 2 = {d};  limit law ' +
                  ', '.join(f'{r}: {lim[r]}' for r in sorted(lim, key=lambda r: -lim[r])) + tag)
    print('   p = 3 reading: T_+ visits class 2 with frequency 2/3 and class 1 with 1/3; T_- the reverse; class 0 never after the first odd step.')

    print('\n== 2. The Syracuse map n -> (3n+1)/2^v on odd n: law mod 3')
    for k in (4, 8, 12, 16):
        cnts, tail = syracuse_law(k)
        tot = int(cnts.sum())
        print(f'   v <= {k:2d}: classes (0, 1, 2): {list(cnts)} of {tot} (+{tail} with v > {k});  share of class 2: {Fr(int(cnts[2]), tot)} = {cnts[2] / tot:.6f}')
    check('syracuse', abs(syracuse_law(16)[0][2] / syracuse_law(16)[0].sum() - 2 / 3) < 1e-4)

    print('\n== 3. The 3-adic visit law of T_+: stationary law of the fair-coin chain mod 3^m')
    for m in (1, 2):
        q = 3 ** m
        inv2 = pow(2, -1, q)
        P = sp.zeros(q, q)
        for r in range(q):
            P[r, (3 * r + 1) * inv2 % q] += sp.Rational(1, 2)
            P[r, r * inv2 % q] += sp.Rational(1, 2)
        A = (P.T - sp.eye(q)).row_insert(q, sp.ones(1, q))
        b = sp.zeros(q, 1).row_insert(q, sp.Matrix([1]))
        sol = A.solve_least_squares(b)
        vals = [sol[i] for i in range(q)]
        print(f'   mod {q} (exact): ' + ', '.join(f'{i}: {vals[i]}' for i in range(q) if vals[i] != 0))
    Hprev = 0.0
    for m in range(1, 10):
        q = 3 ** m
        inv2 = pow(2, -1, q)
        rr = np.arange(q)
        odd_t = (3 * rr + 1) * inv2 % q
        even_t = rr * inv2 % q
        v = np.full(q, 1.0 / q)
        for _ in range(5000):
            w = np.zeros(q)
            np.add.at(w, odd_t, 0.5 * v)
            np.add.at(w, even_t, 0.5 * v)
            if np.abs(w - v).sum() < 1e-15:
                v = w
                break
            v = w
        nz = v[v > 0]
        H = float(-(nz * np.log(nz)).sum() / math.log(3))
        q1 = q // 3
        tot = v[:q1] + v[q1:2 * q1] + v[2 * q1:]
        ok = tot > 0
        cond = np.stack([v[:q1], v[q1:2 * q1], v[2 * q1:]])[:, ok] / tot[ok]
        l1 = np.abs(cond - 1 / 3).sum(axis=0)
        print(f'   mod 3^{m}: support {int(np.count_nonzero(v > 1e-14))} of {q};  entropy of the new digit {H - Hprev:.4f} (log 3 = 1);  '
              f'top digit given the lower ones: mass-weighted L1 deviation {float((tot[ok] * l1).sum()):.4f}, worst {float(l1.max()):.4f}')
        Hprev = H
    print('   reading: biased at the coarse 3-adic scale (1/3 : 2/3 mod 3); finer digits approach uniformity on average '
          '(conditional entropy rises toward 1, the mean deviation decays slowly), not at every residue (the worst stays near 0.41)')

    print('\n== 4. Visit counts mod 3 along actual trajectories, n <= 10^7')
    cyc_plus = np.zeros(201, dtype=np.int64)
    cyc_plus[1] = cyc_plus[2] = 1
    cyc_minus = np.zeros(201, dtype=np.int64)
    for v in (1, 5, 7, 10, 17, 25, 37, 55, 82, 41, 61, 91, 136, 68, 34):
        cyc_minus[v] = 1
    for s, cm, name in ((1, cyc_plus, 'T_+ (3x+1), until 1 or 2'), (-1, cyc_minus, 'T_- (3x-1), until a cycle element')):
        cnt, after = visit_counts(10 ** 7, s, cm)
        tot, tot2 = cnt.sum(), after.sum()
        print(f'   {name}: all visits (0,1,2) shares {cnt[0] / tot:.5f}, {cnt[1] / tot:.5f}, {cnt[2] / tot:.5f};  after the first odd step: '
              f'{after[0] / tot2:.6f}, {after[1] / tot2:.5f}, {after[2] / tot2:.5f}  ({tot2} visits)')
        check(f'visits {s}', after[0] == 0 and abs(after[2 if s == 1 else 1] / tot2 - 2 / 3) < 0.02)

    print('\n== 5. Backward tree of T_+ from 1 (preimages 2n, and (2n-1)/3 when n = 2 mod 3): layer counts by class mod 3')
    layer = [1]
    seen = {1}
    rows = []
    for depth in range(1, 51):
        nxt = []
        for n in layer:
            for m in ([2 * n] + ([(2 * n - 1) // 3] if n % 3 == 2 else [])):
                if m not in seen:
                    seen.add(m)
                    nxt.append(m)
        layer = nxt
        c = [sum(1 for v in layer if v % 3 == r) for r in range(3)]
        rows.append((depth, len(layer), c))
    for depth, sz, c in rows:
        if depth in (5, 10, 20, 30, 40, 45, 50):
            print(f'   depth {depth:2d}: {sz:8d} nodes, classes (0,1,2) = {c}, shares {c[0] / sz:.4f} {c[1] / sz:.4f} {c[2] / sz:.4f};  (c1 - c2)/size = {(c[1] - c[2]) / sz:+.5f}')
    g = [rows[i][1] / rows[i - 1][1] for i in range(40, 50)]
    print(f'   layer growth ratio at depth 41..50: {min(g):.4f} .. {max(g):.4f} (4/3 = 1.3333)')
    sh = rows[-1][2]
    check('tree balance', all(abs(x / rows[-1][1] - 1 / 3) < 0.03 for x in sh))

    print(f'\n== 6. The three basins of T_- (3x-1) on 1..N_C = {N_C:.0e}')
    t1 = time.time()
    memo = np.zeros(N_C + 1, dtype=np.int8)
    for v in (1,):
        memo[v] = 1
    for v in (5, 7, 10):
        memo[v] = 2
    for v in (17, 25, 37, 55, 82, 41, 61, 91, 136, 68, 34):
        memo[v] = 3
    cps = []
    for e in range(1, 10):
        for mlt in (1, 2, 5):
            if mlt * 10 ** e <= N_C:
                cps.append(mlt * 10 ** e)
    cps = np.array(sorted(set(cps)), dtype=np.int64)
    cpc = np.zeros((len(cps), 3), dtype=np.int64)
    race_n = np.zeros(2000000, dtype=np.int64)
    race_o = np.zeros(2000000, dtype=np.int64)
    nrace, maxval = basins(N_C, memo, race_n, race_o, cps, cpc)
    print(f'   done in {time.time() - t1:.1f} s; every n <= N_C reached one of the three cycles (no overflow; largest value met {maxval})')
    check('basins all', int(np.count_nonzero(memo[1:] == 0)) == 0)
    for N, (a, b, c) in zip(cps, cpc):
        if N >= 100:
            print(f'   N = {N:>11d}:  basin {{1}} {a / N:.6f},  basin {{5,7,10}} {b / N:.6f},  basin {{17,...}} {c / N:.6f}')
    names = {1: '{1}', 2: '{5,7,10}', 3: '{17,...}'}
    nr = min(nrace, race_n.shape[0])
    print(f'   changes in the order of the three counts (leader, second) as N grows: {nrace} in total')
    lead_changes = [i for i in range(1, nr) if race_o[i] // 10 != race_o[i - 1] // 10]
    second_changes = [i for i in range(1, nr) if race_o[i] % 10 != race_o[i - 1] % 10]
    print(f'   leader changes: {len(lead_changes)}, the last at N = {race_n[lead_changes[-1]] if lead_changes else None};  second-place changes: {len(second_changes)}, the last at N = {race_n[second_changes[-1]] if second_changes else None}')
    for i in range(nr):
        if i < 6 or race_n[i] >= 10000 or i in lead_changes[-3:]:
            key = race_o[i]
            print(f'     from N = {race_n[i]}: leader {names[key // 10]}, second {names[key % 10]}')
    del memo

    print('\n== 7. Integer cycles of T_+ on Z from every word of length L <= 24')
    vals = np.zeros(200000, dtype=np.int64)
    Ls = np.zeros(200000, dtype=np.int64)
    nout = cycle_census(24, vals, Ls)
    pts = sorted(set(int(v) for v in vals[:nout]))
    # group into cycles
    def Tplus(n):
        return n // 2 if n % 2 == 0 else (3 * n + 1) // 2
    cycles = []
    done = set()
    for x in pts:
        if x in done:
            continue
        cyc = [x]
        y = Tplus(x)
        while y != x:
            cyc.append(y)
            y = Tplus(y)
        done.update(cyc)
        cycles.append(cyc)
    print(f'   {nout} (word, integer point) solutions; distinct integer cycle points {pts}')
    for cyc in cycles:
        print(f'   cycle {cyc}: period {len(cyc)}, odd steps {sum(1 for v in cyc if v % 2)}')
    check('census', sorted(max(c) for c in cycles) == [-17, -5, -1, 0, 2])
    print('   T_- (3x-1) on the positive integers has exactly the negatives: {1}, {5,7,10}, {17,...,136}; T_+ on the positives only {1,2}.')

    print('\n== 8. The carry budget: which (l, L) can carry a cycle, and how small its least odd element must be')
    print('   l   L  side   2^L - 3^l        |delta|/l      budget for n_min   cycles found')
    found = {(len(c), sum(1 for v in c if v % 2)): c for c in cycles if c != [0]}
    for l in range(1, 13):
        for L in (math.floor(l * math.log2(3)), math.ceil(l * math.log2(3))):
            D = 2 ** L - 3 ** l
            delta = L * math.log(2) - l * math.log(3)
            t = abs(delta) / l
            if D > 0:
                bud = 1 / (3 * (math.exp(t) - 1))
                side = 'above (s=+1)'
            else:
                bud = 1 / (3 * (1 - math.exp(-t)))
                side = 'below (s=-1)'
            hit = found.get((L, l))
            print(f'   {l:2d} {L:3d}  {side:13s} {D:>14d}   {t:.6f}      {bud:12.2f}       {hit if hit else "-"}')
    print('\nFAILED CHECKS:', FAIL if FAIL else 'none')
    print('ALL CHECKS PASSED' if not FAIL else 'SOME CHECKS FAILED')
    print('total time %.1f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
