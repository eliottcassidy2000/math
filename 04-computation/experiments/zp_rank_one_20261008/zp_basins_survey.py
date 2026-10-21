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
  (S4) the sign-conjugate side: cycles of C_p on the negative integers, i.e. of F_p(x) = -C_p(-x) = x/p or
       floor((p+1)x/p) on the positive integers, from the starts -10^5..-1, p <= 31.  F_p is conjugate to C_p on Z_p by
       x -> -x, so it satisfies everything THM-4610 proves, and it fixes 1, ..., p-1 (the 3x - 1 barrier, KNOWN).
  (S5) the 1-cycle equation: a cycle with one run of divisions, read from the start of its multiplication run x, satisfies
       x (p^(k1+k2) - (p+1)^k2) = c_w with c_w = sum_j (p+1)^(k2-1-j) p^j (p - i_j) over the run's digits i_j in {1..p-1}
       (checked for every extra 1-cycle of (S1) and (S4)).  For p = 2 (the Terras map) the run's digits are all 1, so each
       shape has the single candidate c = 3^k2 - 2^k2; for k2 <= 2000 and every k1 with |D| <= c it is divisible only for the
       trivial cycle (1, 1) on the positive side and for (0, 1), (1, 2) (the cycles -1, -5) on the negative side, as
       Steiner's theorem says for the positive side.  For p >= 3 each shape has (p-1)^k2 candidates.
       HEURISTIC count: if each candidate is divisible with probability 1/|D|, the expected number of 1-cycles with at least
       one division is sum (p-1)^k2/|D| over shapes with |D| <= (p-1)((p+1)^k2 - p^k2) (k1 <= 11, k2 < 60p); for p = 2 the
       same sum over primitive shapes with 1 < |D| <= 3^k2 - 2^k2, k2 <= 2000, split at k2 = 5.
  The cycles found in (S1) and (S4) are compared with Carnielli's Table 1 (AMEN 15 (2015) 207-215; cycles computed by
  K. R. Matthews for d <= 150, |x| <= 50000), and each cycle's shape is printed: k1 division steps, k2 multiplication
  steps, lambda = (1/p)^k1 ((p+1)/p)^k2 < 1 (the near-balance at which cycles sit; cf. Tremblay arXiv:1907.13302).
Prints ALL CHECKS PASSED.  Runtime about 1 minute."""
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


CARNIELLI_POS = {3: [(7, 9)], 5: [], 7: [], 11: [(642, 57)], 13: [], 17: [(79, 49)], 19: [], 23: [(82, 72)],
                 29: [(111, 97)], 31: [(389, 108)]}
CARNIELLI_NEG = {3: [(-22, 5)], 5: [(-57, 10)], 7: [], 11: [], 13: [(-147, 36)], 17: [], 19: [(-954, 117)],
                 23: [(-800, 75)], 29: [], 31: []}


def shape(p, x0):
    """(length, k1 = division steps, k2 = multiplication steps, lambda, m = number of runs of divisions) of the cycle
    through x0; m = 1 means one run of multiplications followed by one run of divisions (a 1-cycle in the sense of
    Steiner 1977 and Simons-de Weger 2005)."""
    pat = []; x = x0
    while True:
        pat.append(x % p == 0)
        x = C(p, x)
        if x == x0:
            break
    k1 = sum(pat); k2 = len(pat) - k1
    m = sum(1 for j in range(len(pat)) if pat[j] and not pat[j - 1])
    return k1 + k2, k1, k2, (1 / p) ** k1 * ((p + 1) / p) ** k2, m


def survey_neg(p, N):
    """Cycles reached from the starts -N..-1, labelled by their element of least absolute value: {label: (length, share)}."""
    term = {}; cyc = {}; share = {}
    lim = 50 * N
    for n in range(1, N + 1):
        x = -n
        if x in term:
            share[term[x]] = share.get(term[x], 0) + 1
            continue
        path = []; pos = {}
        while x not in term and x not in pos:
            pos[x] = len(path); path.append(x); x = C(p, x)
            if len(path) > 10 ** 6:
                return None
        if x in term:
            m = term[x]
        else:
            c = path[pos[x]:]; m = max(c); cyc[m] = len(c)
        for z in path:
            if -z <= lim:
                term[z] = m
        share[m] = share.get(m, 0) + 1
    return {m: (cyc[m], share[m] / N) for m in sorted(share, reverse=True)}


def run_equation(p, x0):
    """For a cycle through x0 with one run of divisions: (x, k1, k2, D, c_w), read from the start x of its multiplication run."""
    P = p + 1
    cyc = [x0]; z = C(p, x0)
    while z != x0:
        cyc.append(z); z = C(p, z)
    div = [x % p == 0 for x in cyc]
    starts = [j for j in range(len(cyc)) if not div[j] and div[j - 1]]
    assert len(starts) == 1
    x = cyc[starts[0]]
    digits = []; y = x
    while y % p:
        digits.append(y % p); y = C(p, y)
    k2 = len(digits); k1 = len(cyc) - k2
    c = sum(P ** (k2 - 1 - j) * p ** j * (p - digits[j]) for j in range(k2))
    return x, k1, k2, p ** (k1 + k2) - P ** k2, c


def steiner_check(K2=2000):
    """p = 2: primitive shapes (k1, k2), k2 <= K2, |2^(k1+k2) - 3^k2| <= 3^k2 - 2^k2, whose single candidate is divisible."""
    pos = []; neg = []
    for k2 in range(1, K2 + 1):
        c = 3 ** k2 - 2 ** k2
        k1 = 0 if k2 == 1 else 1           # (0, k2) with k2 > 1 is the fixed point -1 repeated
        while True:
            D = 2 ** (k1 + k2) - 3 ** k2
            if D > c:
                break
            if D != 0 and abs(D) <= c and c % D == 0:
                (pos if D > 0 else neg).append((k1, k2, c // D))
            k1 += 1
    return pos, neg


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
    ok_shape = True
    for p in table:
        extra = [(m, L) for m, (L, s_) in table[p].items() if m != 1]
        ok_shape &= extra == CARNIELLI_POS[p] and table[p][1][0] == p
        sh = [(m,) + shape(p, m)[1:3] + (round(shape(p, m)[3], 4), shape(p, m)[4]) for m in table[p]]
        print(f'     p = {p}: cycle shapes (minimum, k1, k2, lambda, runs of divisions) {sh}')
        ok_shape &= all(x[3] < 1 and x[4] == (2 if (p, x[0]) == (11, 642) else 1) for x in sh)
    check(ok_shape, "the cycles agree with Carnielli's Table 1 (positive side; trivial cycle 1..p of length p); every cycle has "
          'lambda = (1/p)^k1 ((p+1)/p)^k2 < 1; every cycle has a single run of divisions (a 1-cycle) except the 642-cycle of C_11, '
          'which has two (a 2-cycle) (FINITE-EXACT; cycles KNOWN)')
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
    print('(S4) the sign-conjugate side: cycles of C_p on the negative integers (= cycles of F_p on the positive integers), '
          'starts -10^5..-1')
    ok4 = True; negtab = {}
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        tab = survey_neg(p, 10 ** 5); negtab[p] = tab
        if tab is None:
            ok4 = False; print(f'     p = {p}: an orbit ran over 10^6 steps'); continue
        elem = sorted((m for m in tab if tab[m][0] == 1), reverse=True)
        extra = [(m, L) for m, (L, s_) in tab.items() if L > 1]
        ok4 &= elem == list(range(-1, -p, -1)) and extra == CARNIELLI_NEG[p]
        ok4 &= all(shape(p, m)[4] == 1 and shape(p, m)[3] > 1 for m, (L, s_) in tab.items() if L > 1)
        print(f'     p = {p}: fixed points {elem[0]}..{elem[-1]}; extra cycles (least |x|, length, share, k1, k2, lambda, runs): '
              f'{[(m, L, round(s_, 4)) + shape(p, m)[1:3] + (round(shape(p, m)[3], 4), shape(p, m)[4]) for m, (L, s_) in tab.items() if L > 1]}')
    check(ok4, "on the negative integers C_p has the fixed points -1..-(p-1) and exactly the extra cycles of Carnielli's Table 1 "
          '(p = 3: -22; 5: -57; 13: -147; 19: -954; 23: -800; none for 7, 11, 17, 29, 31), for starts down to -10^5; so F_p, '
          'conjugate to C_p, has p - 1 fixed points on the positive integers; every extra negative cycle has lambda > 1 and a '
          'single run of divisions (FINITE-EXACT; cycles KNOWN)')
    print('(S5) the 1-cycle equation x (p^(k1+k2) - (p+1)^k2) = c_w, and the single candidate for p = 2')
    ok5 = True; rows5 = []
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        ms = [m for m in table[p] if m != 1 and shape(p, m)[4] == 1]
        ms += [m for m, (L, s_) in negtab[p].items() if L > 1 and shape(p, m)[4] == 1]
        for m in ms:
            x, k1, k2, D, c = run_equation(p, m)
            ok5 &= x * D == c
            rows5.append((p, m, x, k1, k2))
    print(f'     extra 1-cycles (p, label, start of the multiplication run, k1, k2): {rows5}')
    pos, neg = steiner_check()
    print(f'     p = 2, k2 <= 2000: divisible single candidates, positive side {pos}, negative side {neg}')
    check(ok5 and len(rows5) == 10 and pos == [(1, 1, 1)] and neg == [(0, 1, -1), (1, 2, -5)],
          'every extra 1-cycle of C_p (p <= 31, both signs; 10 cycles) satisfies x (p^(k1+k2) - (p+1)^k2) = c_w for its run digits; for '
          'p = 2 the single candidate 3^k2 - 2^k2 is divisible only for the trivial cycle (positive side, as Steiner proved) and for '
          'the cycles -1, -5 (negative side), k2 <= 2000 (FINITE-EXACT)')
    naive = {}
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        P = p + 1; e_pos = e_neg = 0.0
        for k1 in range(1, 12):
            for k2 in range(1, 60 * p):
                D = p ** (k1 + k2) - P ** k2
                if D == 0 or abs(D) > (p - 1) * (P ** k2 - p ** k2):
                    continue
                if D > 0:
                    e_pos += (p - 1) ** k2 / D
                else:
                    e_neg += (p - 1) ** k2 / -D
        naive[p] = (round(e_pos, 2), round(e_neg, 2))
    print(f'     naive expected 1-cycles with a division (positive, negative; trivial cycle included): {naive}')
    t_small = [0.0, 0.0]; t_tail = [0.0, 0.0]
    for k2 in range(1, 2001):
        c = 3 ** k2 - 2 ** k2
        k1 = 0 if k2 == 1 else 1
        while True:
            D = 2 ** (k1 + k2) - 3 ** k2
            if D > c:
                break
            if abs(D) > 1:
                (t_small if k2 <= 5 else t_tail)[0 if D > 0 else 1] += 1 / abs(D)
            k1 += 1
    print(f'     p = 2, beyond the |D| = 1 shapes (the cycles 1, -1, -5): expected (positive, negative) from k2 <= 5: '
          f'({t_small[0]:.3f}, {t_small[1]:.3f}); from 5 < k2 <= 2000: ({t_tail[0]:.2e}, {t_tail[1]:.2e})')
    check(all(0.5 < naive[p][0] < 1.2 for p in naive) and t_tail[0] < 0.05 and t_tail[1] < 0.05,
          'the naive count gives 0.5-1.2 expected positive 1-cycles per base p = 3..31 (the same order as found), while for p = 2 '
          'the expectation beyond the three |D| = 1 shapes sits in shapes with k2 <= 5 (checked exactly) and the tail is below 0.05 '
          '(HEURISTIC count, exact arithmetic)')
    extra_any = sorted({p for p in table if len(table[p]) > 1} |
                       {p for p in negtab if any(L > 1 for L, s_ in negtab[p].values())})
    print(f'     primes p <= 31 (p >= 3) with an extra cycle of either sign: {extra_any}')
    check(extra_any == [3, 5, 11, 13, 17, 19, 23, 29, 31],
          'among the primes 3 <= p <= 31, only p = 7 has no extra cycle of either sign in these ranges (FINITE-EXACT; Carnielli: '
          'd = 7 is the first d with none)')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
