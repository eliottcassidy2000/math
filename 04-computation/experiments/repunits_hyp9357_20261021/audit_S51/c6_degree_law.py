"""S51 audit c6: THM-4749 from the auditor's own inventory (c4, cached) and own arithmetic.
  (1) Session24's formula N(n,h) = sum_{r,M} B(h,r,M) C(n-M+r, r) against a direct convolution over sequences of strong
      components (singletons included), every h < 2100, n <= 30; degree w(h), onset floor(3 log_5 h), leading
      coefficient (#ordered top-length words)/w!;
  (2)-(4) the cyclotomic corollaries: odd k, 3 | k, k = 2^j; the lists printed in the proofs;
  (5) the w table for k <= 16 with certification (a maximal factorization using only factors <= 80,265), the ties and
      the decided ties k = 10, 11 (third differences 24, 18; second differences 184 = 2 S(23) S(89), 516 = 2 S(683)),
      and the odd-k claim for k <= 40 (conditional beyond the certified range);
  (6) the finite table n = 8, 10, 14.
Prints ALL CHECKS PASSED or the failures."""
import os, math, time
from collections import Counter, defaultdict
from functools import lru_cache
from math import comb
import numpy as np
from sympy import isprime, divisors, factorint

FAIL = []
T0 = time.time()
SCR = r'C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-ephrepos-math\e1fb236f-ec78-4211-aa52-c187130b32b0\scratchpad\audit_S51'
BOUND = 80265


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


inV = lambda d: d >= 3 and d % 2 == 1 and d not in (7, 21)


@lru_cache(maxsize=None)
def w(h, cap=None):
    """longest factorization of h into factors in V (optionally only factors <= cap); -1 if none, 0 for h = 1"""
    if h == 1:
        return 0
    best = -1
    for d in divisors(h):
        if inV(d) and (cap is None or d <= cap):
            r = w(h // d, cap)
            if r >= 0:
                best = max(best, r + 1)
    return best


Z = np.load(os.path.join(SCR, 'inventory_2100.npz'))
orders = sorted(int(k[1:]) for k in Z.files if k.startswith('H'))
Sm = {n: Counter(int(h) for h in Z[f'H{n}']) for n in orders}
S = Counter()
for n in orders:
    S.update(Sm[n])

print('(1) Session24 formula, degree, onset and leading coefficient on exact data (h < 2100)')
LIM = 2100
NM = 30
# direct convolution: classes = sequences of strong components incl. singletons
s = {1: Counter({1: 1})}
for n in orders:
    s[n] = Sm[n]
F = [Counter() for _ in range(NM + 1)]
F[0][1] = 1
for n in range(1, NM + 1):
    for m in range(1, n + 1):
        if m not in s:
            continue
        for d, c in s[m].items():
            for h, f in F[n - m].items():
                if h * d < LIM:
                    F[n][h * d] += c * f
# ordered nontrivial core words: B[h][(r, M)]
B = defaultdict(Counter)
B[1][(0, 0)] = 1
frontier = {(1, 0, 0): 1}
while frontier:
    new = defaultdict(int)
    for (h, r, M), c in frontier.items():
        for m in orders:
            for d, cc in Sm[m].items():
                if h * d < LIM:
                    new[(h * d, r + 1, M + m)] += c * cc
    for (h, r, M), c in new.items():
        B[h][(r, M)] += c
    frontier = new
ok = True
for h in range(1, LIM, 2):
    for n in range(0, NM + 1):
        val = sum(c * (comb(n - M + r, r) if n >= M else 0) for (r, M), c in B[h].items())
        ok &= val == F[n][h]
check(ok, 'N(n,h) = sum B(h,r,M) C(n-M+r, r) (zero for n < M) equals the direct component convolution, every odd h < 2100, n <= 30')
okdeg = True; okon = True; oklead = True
for h in range(3, LIM, 2):
    wh = w(h)
    if wh < 0:
        okdeg &= all(F[n][h] == 0 for n in range(NM + 1))
        continue
    on = math.floor(3 * math.log(h, 5))
    seq = [F[n][h] for n in range(on, NM + 1)]
    dif = seq
    for _ in range(wh):
        dif = [dif[i + 1] - dif[i] for i in range(len(dif) - 1)]
    top = sum(c for (r, M), c in B[h].items() if r == wh)
    okdeg &= len(set(dif)) == 1 and dif[0] > 0
    oklead &= dif[0] == top                       # w-th difference = w! * leading coefficient = #ordered top words
    # onset: the polynomial (from the binomial formula) agrees with the true count from n = on
    pol = lambda n: sum(c * comb(n - M + r, r) if n - M + r >= 0 else c * (math.prod(n - M + r - i for i in range(r)) // math.factorial(r)) for (r, M), c in B[h].items())
    okon &= all(pol(n) == F[n][h] for n in range(on, NM + 1))
check(okdeg, 'for EVERY odd h < 2100 (not only repunits), N(n,h) is a polynomial of degree w(h) on n >= floor(3 log_5 h) (n <= 30); w(h) < 0 iff N = 0')
check(oklead, 'its w-th difference equals the number of ordered top-length strong words (leading coefficient #words / w!)')
check(okon, 'the binomial polynomial agrees with the exact count from n = floor(3 log_5 h) on, every odd h < 2100')

print('(2)-(4) cyclotomic corollaries')
check(all((2 ** k + 1) % 3 == 0 and (2 ** k - 1) % 3 == 1 for k in range(1, 200, 2)), 'odd k: 3 | 2^k + 1 and 2^k - 1 = 1 mod 3')
lst = [(2 ** k + 1) // 3 for k in range(3, 18, 2)]
check(lst == [3, 11, 43, 171, 683, 2731, 10923, 43691] and all(inV(x) and x <= BOUND for x in lst), f'(2^k+1)/3 for odd 3 <= k <= 17: {lst}, all in V below 80,265')
pr = [2 ** k - 1 for k in range(3, 18, 2) if isprime(2 ** k - 1)]
print(f'     2^k - 1 prime for odd 3 <= k <= 17: {pr} (the proof lists 31, 127, 8191, 131071 and omits 7 = 2^3 - 1, which is not in V: N(n,7) = 0)')
pb = [k for k in range(3, 18, 2) if isprime(2 ** k - 1) or w(2 ** k - 1) <= 1]
check(pb == [3, 5, 7, 9, 13, 17], f'odd k <= 17 with 2^k - 1 prime or blocked: {pb}')
check(w(7) < 0 and w(511) == 1 and w(513) == 4 and w(63) == 1 and w(4095) == 4, 'w(7) = -inf (no factorization), w(511) = 1, w(513) = 4, w(63) = 1, w(4095) = 4')
check(all((2 ** k - 1) % 7 == 0 for k in range(3, 300, 3)) and not any((2 ** k - 1) % 7 == 0 for k in range(1, 300) if k % 3) and not any((2 ** k + 1) % 7 == 0 for k in range(1, 300)),
      '7 | 2^k - 1 iff 3 | k; 7 never divides 2^k + 1')
check(all(math.prod(2 ** (2 ** i) + 1 for i in range(j)) == 2 ** (2 ** j) - 1 for j in range(1, 8)), '2^(2^j) - 1 = F_0 ... F_(j-1)')
check(all(w(2 ** (2 ** j) - 1) >= j for j in range(1, 5)), 'w(2^(2^j) - 1) >= j (j <= 4); NOTE: true for all j, since F_i = H(NT_(2^i+2)) is a strong value for EVERY i (no 80,265 bound needed)')
check(all(isprime(2 ** (2 ** j) + 1) for j in range(5)) and not isprime(2 ** 32 + 1), 'F_j prime for j <= 4, F_5 composite')

print('(5) the table k <= 16, certification, ties')
rows = {}
for k in range(2, 41):
    a, b = 2 ** k - 1, 2 ** k + 1
    rows[k] = (w(a), w(b), w(a, BOUND) == w(a) and w(b, BOUND) == w(b))
tab_m = [rows[k][0] for k in range(4, 17)]; tab_p = [rows[k][1] for k in range(4, 17)]
check(tab_m == [2, 1, 1, 1, 3, 1, 3, 2, 4, 1, 3, 2, 4] and tab_p == [1, 2, 2, 2, 1, 4, 3, 2, 2, 2, 3, 4, 1], f'w table k = 4..16: w- = {tab_m}, w+ = {tab_p} (THM-4749 (5))')
check(all(rows[k][2] for k in range(2, 17)), 'k <= 16: every w is attained by a factorization using only factors <= 80,265 (certified strong values)')
cert = [k for k in range(17, 41) if rows[k][2]]
print(f'     k in 17..40 whose w values are also certified by factors <= 80,265: {cert} (THM-4749 calls all k >= 17 conditional)')
odd = [k for k in range(5, 41, 2)]
check(not any(rows[k][0] > rows[k][1] for k in odd) and [k for k in odd if rows[k][0] == rows[k][1]] == [11, 23, 29], 'odd 5 <= k <= 40: never w- > w+; ties exactly at k = 11, 23, 29')
rar = [k for k in range(4, 17) if rows[k][0] < rows[k][1]]; com = [k for k in range(4, 17) if rows[k][0] > rows[k][1]]; tie = [k for k in range(4, 17) if rows[k][0] == rows[k][1]]
check(rar == [5, 6, 7, 9, 13, 15] and com == [4, 8, 12, 16] and tie == [10, 11, 14], f'k <= 16: rarer {rar}, commoner {com}, ties {tie}')


def diff(h, wh):
    on = math.floor(3 * math.log(h, 5))
    seq = [F[n][h] for n in range(on, NM + 1)]
    for _ in range(wh):
        seq = [seq[i + 1] - seq[i] for i in range(len(seq) - 1)]
    return seq[0]


d1023, d1025, d2047, d2049 = diff(1023, 3), diff(1025, 3), diff(2047, 2), diff(2049, 2)
print(f'     S(3) = {S[3]}, S(5) = {S[5]}, S(11) = {S[11]}, S(23) = {S[23]}, S(31) = {S[31]}, S(41) = {S[41]}, S(89) = {S[89]}, S(683) = {S[683]}')
check(d1023 == 24 == 6 * S[3] * S[11] * S[31] and d1025 == 18 == 3 * S[5] ** 2 * S[41], f'k = 10: third differences {d1023} = 3! S(3)S(11)S(31) and {d1025} = 3 S(5)^2 S(41) (only 3*11*31 and 5*5*41 have length 3)')
check(d2047 == 184 == 2 * S[23] * S[89] and d2049 == 516 == 2 * S[3] * S[683], f'k = 11: second differences {d2047} = 2 S(23) S(89), {d2049} = 2 S(683)')

print('(6) the finite table')
tbl = {8: [(20, 6), (10, 23), (2, 25), (10, 40), (28, 30), (32, 35), (0, 0), (0, 0)],
       10: [(42, 10), (18, 49), (6, 59), (34, 136), (188, 163), (290, 363), (602, 543), (990, 1029)],
       14: [(110, 18), (34, 125), (14, 151), (82, 424), (1176, 455), (946, 2427), (4770, 3133), (7426, 10168)]}
check(all([(F[n][2 ** k - 1], F[n][2 ** k + 1]) for k in range(4, 12)] == v for n, v in tbl.items()), 'THM-4749 (6) table n = 8, 10, 14 reproduced')
check(all(F[n][63] == 2 * (n - 7) for n in range(7, NM + 1)), 'N(n, 63) = 2(n - 7) for 7 <= n <= 30 (also Session24 PROOF.md section 5 table, H = 63 row)')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
