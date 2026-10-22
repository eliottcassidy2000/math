"""S51 audit c5b: wording variants of the second pre-registration and the test-range composition.
The G-section baseline sentence defines the generic remainder as "no exact cover and trivial Aut", while G2 and the
session script use "no exact cover and 3 not dividing |Aut|".  This script evaluates G2-G3 under both readings, prints
the training baseline under both (the prereg quotes +0.049 +- 0.016), and the share of the test-range classes that
have order <= 9 (the order-9 slice, whose class counts N(9, d) = S_9(d) for d > 661 were computed before the
pre-registration by the session's scratch rep1.py and are in the public S50 census).
Uses the auditor's cached inventory (c4) and the auditor's own cover routine (imported from c5's definitions by
re-implementation here).  Prints ALL CHECKS PASSED or the failures."""
import os, math, time
from collections import Counter
import numpy as np
import numba as nb
from sympy import isprime
import alib51 as L

FAIL = []
T0 = time.time()
SCR = r'C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-ephrepos-math\e1fb236f-ec78-4211-aa52-c187130b32b0\scratchpad\audit_S51'


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


@nb.njit(cache=True)
def popc(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@nb.njit(cache=True)
def has_cover(n, adj):
    Hs = L.ham_subsets(n, adj)
    full = (1 << n) - 1
    H = Hs[full]
    for v in range(n):
        rest = full ^ (1 << v)
        low = rest & -rest
        S_ = rest
        while True:
            S_ = (S_ - 1) & rest
            if S_ == 0:
                break
            if not (S_ & low) or popc(S_) < 2 or popc(rest ^ S_) < 2:
                continue
            a = S_ | (1 << v); b = (rest ^ S_) | (1 << v)
            if Hs[a] * Hs[b] == H and L.is_strong_mask(n, adj, a) and L.is_strong_mask(n, adj, b):
                return 1
    return 0


@nb.njit(parallel=True, cache=True)
def cov_all(n, A):
    out = np.zeros(A.shape[0], dtype=np.int64)
    for t in nb.prange(A.shape[0]):
        out[t] = has_cover(n, A[t])
    return out


def F5(d):
    Ld = math.log(d)
    return [1.0, Ld, Ld ** 2, Ld ** 3, d % 3 == 0, d % 9 == 0, d % 5 == 0, d % 8 == 1, d % 8 == 7]


def ols(cnt, ds):
    ds = [d for d in ds if cnt[d] > 0]
    X = np.array([F5(d) for d in ds], dtype=float); y = np.log(np.array([cnt[d] for d in ds], dtype=float))
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    s2 = float(res @ res) / (len(ds) - X.shape[1])
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    return dict(zip(['3|d', '9|d', '5|d', '1 mod 8', '7 mod 8'], beta[4:])), dict(zip(['3|d', '9|d', '5|d', '1 mod 8', '7 mod 8'], se[4:]))


Z = np.load(os.path.join(SCR, 'inventory_2100.npz'))
orders = sorted(int(k[1:]) for k in Z.files if k.startswith('H'))
S = Counter(); G3n = Counter(); G1a = Counter(); byord = Counter()
for n in orders:
    H, aut, A = Z[f'H{n}'], Z[f'aut{n}'], Z[f'A{n}']
    cv = cov_all(n, A)
    for h, a, c in zip(H, aut, cv):
        h = int(h); S[h] += 1
        if 1101 <= h <= 2099:
            byord[n] += 1
        if not c:
            if a % 3:
                G3n[h] += 1
            if a == 1:
                G1a[h] += 1
TR = list(range(301, 1099, 2)); TE = list(range(1101, 2100, 2))
for name, cnt in (('no cover, 3 not | Aut (script)', G3n), ('no cover, |Aut| = 1 (prereg baseline wording)', G1a)):
    btr, setr = ols(cnt, TR); bte, sete = ols(cnt, TE)
    print(f'     {name}: training 3|d {btr["3|d"]:+.3f}, 1 mod 8 {btr["1 mod 8"]:+.3f} +- {setr["1 mod 8"]:.3f}; '
          f'test 3|d {bte["3|d"]:+.3f}, 1 mod 8 {bte["1 mod 8"]:+.3f}, 7 mod 8 {bte["7 mod 8"]:+.3f}')
    g2 = abs(bte['3|d']) < 0.04; g3 = bte['1 mod 8'] > 0 and abs(bte['7 mod 8']) < 0.03
    check(g2 and g3, f'G2 and G3 pass under the reading "{name}"')
tot = sum(byord.values())
print(f'     test-range classes by order: {dict(sorted(byord.items()))}; order <= 9 share {sum(v for k, v in byord.items() if k <= 9) / tot:.3f}')
check(byord[9] / tot > 0.1, 'the order-9 slice is a non-negligible part (>10%) of the test-range classes')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
