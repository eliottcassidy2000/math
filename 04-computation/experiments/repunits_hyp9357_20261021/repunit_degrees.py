"""The repunit degree law (session opus-2026-10-21-S51; THM-4749).  Note: 05-knowledge/results/repunit_rarity_hyp9357_20261021.md.

KNOWN (Session24, 05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md):
for fixed odd h, the number N(n,h) of n-vertex tournament classes with h Hamiltonian paths is, for n >= floor(3 log_5 h),
a polynomial in n of degree w(h) = the largest number of nontrivial strong components whose path counts multiply to h.
With the strong value set V = {odd d >= 3} - {7, 21} (FINITE-EXACT for d <= 80,265: THM-4097, 4102, 4104 and S50),
w(h) = the longest factorization of h into factors in V.

This script
  (W) computes w(2^k - 1) and w(2^k + 1) for 2 <= k <= 40 from the cyclotomic factorizations; for k <= 16 every
      prime factor is below 80,265, so the values are FINITE-EXACT, and beyond they assume V-membership of the factors;
  (P) checks the degree law on the exact counts: N(n,h) for 1 <= n <= 40 and every h = 2^k +- 1 < 2100, by sequences of
      strong components from the complete strong inventory below 2100; the w(h)-th finite difference is a positive
      constant and the next one vanishes from the onset on;
  (A) evaluates the cyclotomic corollaries: odd k (Phi_2(2) = 3 divides 2^k + 1), 3 | k (Phi_3(2) = 7 divides 2^k - 1),
      k = 2^j (Fermat products), and the asymptotic sign of N(n, 2^k - 1) - N(n, 2^k + 1).
Prints ALL CHECKS PASSED."""
import os, sys, time, math
from collections import Counter, defaultdict
from functools import lru_cache
from multiprocessing import Pool
from sympy import factorint, isprime, divisors

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'h_spectrum_multiplicities_20261021', 'audit_S50'))
import alib as L                                       # noqa: E402
from b3_ear_census_625 import grow                     # noqa: E402

FAIL = []
T0 = time.time()
BOUND = 80265                     # every odd d in [3, BOUND] except 7, 21 is a strong count (FINITE-EXACT)


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


def in_V(d):
    return d % 2 == 1 and d >= 3 and d not in (7, 21)


@lru_cache(maxsize=None)
def w(h):
    """longest factorization of h into factors in V (0 if h = 1; -inf if none)"""
    if h == 1:
        return 0
    best = -10 ** 9
    for d in divisors(h):
        if d > 1 and in_V(d):
            r = w(h // d)
            if r >= 0:
                best = max(best, r + 1)
    return best


def certified(h):
    return all(p <= BOUND for p in factorint(h))


if __name__ == '__main__':
    print('(W) w(2^k - 1) and w(2^k + 1)')
    rows = []
    for k in range(2, 41):
        a, b = 2 ** k - 1, 2 ** k + 1
        wa, wb = w(a), w(b)
        cert = certified(a) and certified(b)
        rows.append((k, wa, wb, cert))
        fa = '*'.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(factorint(a).items()))
        fb = '*'.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(factorint(b).items()))
        sign = '<' if wa < wb else '>' if wa > wb else '='
        print(f'     k = {k:2d}: w(2^k-1) = {wa:2d} [{fa}]   {sign}   w(2^k+1) = {wb:2d} [{fb}]   {"FINITE-EXACT" if cert else "conditional on V-membership of large factors"}')
    check(all(c for k, _, _, c in rows if k <= 16), 'for k <= 16 every prime factor of 2^k +- 1 is below 80,265 (degrees certified)')
    check(w(7) < 0 and w(21) < 0 and w(63) == 1 and w(511) == 1 and w(31) == 1 and w(127) == 1 and w(8191) == 1,
          '7, 21 have no factorization into V; w = 1 for 31, 63, 127, 511, 8191 (63 = 3*21 = 9*7 and 511 = 7*73 are blocked)')
    check(all((2 ** k + 1) % 3 == 0 and (2 ** k - 1) % 3 != 0 for k in range(3, 41, 2)), 'odd k: Phi_2(2) = 3 = H(C_3) divides 2^k + 1 and not 2^k - 1')
    check(all(in_V((2 ** k + 1) // 3) and (2 ** k + 1) // 3 <= BOUND and w(2 ** k + 1) >= 2 for k in range(3, 18, 2)),
          'odd 3 <= k <= 17: (2^k + 1)/3 is a strong value below 80,265, so w(2^k + 1) >= 2')
    pb = [k for k in range(3, 18, 2) if isprime(2 ** k - 1) or w(2 ** k - 1) <= 1]
    check(pb == [3, 5, 7, 9, 13, 17] and all(w(2 ** k - 1) <= 1 < 2 <= w(2 ** k + 1) for k in pb),
          f'odd k <= 17 with 2^k - 1 prime or blocked: {pb}; there deg N(., 2^k - 1) <= 1 < 2 <= deg N(., 2^k + 1)')
    odd = [(k, wa, wb) for k, wa, wb, c in rows if k % 2 == 1 and k >= 5]
    ties = [k for k, wa, wb in odd if wa == wb]
    check(not any(wa > wb for k, wa, wb in odd) and ties == [11, 23, 29],
          f'odd 5 <= k <= 40 (k > 17 conditional): w(2^k - 1) > w(2^k + 1) never occurs; ties at k = {ties} (then the leading coefficients decide)')
    pw = [(k, wa, wb) for k, wa, wb, c in rows if k in (4, 8, 16)]
    check(all(wa > wb for k, wa, wb in pw), f'k = 2^j (4, 8, 16): 2^k - 1 is a product of distinct Fermat numbers, 2^k + 1 a Fermat prime: {pw}')

    print('(P) the degree law on exact counts (h = 2^k +- 1 < 2100, n <= 40)')
    inv = {3: {L.canon(L.chain((2,))): 3}}
    LIM = 2100
    with Pool(12) as pool:
        n = 4
        while True:
            par = [c for c, h in inv[n - 1].items() if h <= LIM]
            got, _ = grow(pool, n, par, cap=LIM) if par else ({}, None)
            got = {c: h for c, h in got.items() if h <= LIM}
            if not got:
                break
            inv[n] = got
            n += 1
    print(f'     strong inventory below {LIM}: {sum(len(v) for v in inv.values())} classes [{time.time() - T0:.0f}s]')
    targets = sorted({2 ** k + e for k in range(2, 12) for e in (-1, 1)} - {1})
    targets = [h for h in targets if h < LIM]
    # strong classes whose H divides some target suffice
    need = set()
    for h in targets:
        need |= set(divisors(h))
    s = defaultdict(Counter); s[1][1] = 1
    for m, dd in inv.items():
        for c, h in dd.items():
            if h in need:
                s[m][h] += 1
    NM = 40
    F = [Counter() for _ in range(NM + 1)]
    F[0][1] = 1
    for n in range(1, NM + 1):
        for m in range(1, n + 1):
            for d, cnt in s[m].items():
                for h, f in F[n - m].items():
                    if h * d in need:
                        F[n][h * d] += cnt * f
    allok = True
    for h in targets:
        seq = [F[n][h] for n in range(NM + 1)]
        deg = w(h)
        onset = math.floor(3 * math.log(h, 5)) if h > 1 else 0
        diffs = seq[onset:]
        for _ in range(deg):
            diffs = [diffs[i + 1] - diffs[i] for i in range(len(diffs) - 1)]
        ok = len(set(diffs)) == 1 and diffs[0] > 0 if deg >= 0 else all(x == 0 for x in seq)
        allok &= ok
        print(f'     h = {h}: w = {deg}; N(n,h) for n = 8..14: {seq[8:15]}; {deg}-th difference from n = {onset}: {sorted(set(diffs))} {"ok" if ok else "FAIL"}')
    check(allok, 'every N(n, 2^k +- 1) (h < 2100) is a polynomial of degree w(h) from n = floor(3 log_5 h) on (Session24 law, exact data)')

    print('(A) asymptotic sign of N(n, 2^k - 1) - N(n, 2^k + 1)')
    for k, wa, wb, c in rows:
        if 4 <= k <= 16:
            print(f'     k = {k:2d}: {"2^k - 1 rarer" if wa < wb else "2^k - 1 commoner" if wa > wb else "same degree (leading coefficients decide)"} as n -> infinity')
    check(True, 'asymptotic signs listed')

    print()
    print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
