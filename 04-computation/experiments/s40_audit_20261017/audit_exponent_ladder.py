"""S40 audit: independent exponent ladder of imaginary quadratic class groups for -X <= D < 0.

Independent of the session's code: class numbers by Cohen-style reduced-form counting (audit_common.h_forms), forms
listed by the same loop, composition by Cohen's Algorithm 5.4.7 (Shanks), reduction by Cohen's Algorithm 5.4.2, exact
test f^e = 1 for every reduced form (binary powering). For e = 2 also the direct criterion 'every reduced form is
ambiguous (b = 0, b = a or a = c)'. Group structure at the largest members by counting elements of each order.
Usage: python audit_exponent_ladder.py X
"""
import sys, os, json, math, time
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_common import h_forms, fundamental_neg

EXPS = (1, 2, 3, 4, 5, 7, 9)


@njit(cache=True)
def red(a, b, c):
    while True:
        if not (-a < b and b <= a):
            q = b // (2 * a)
            r = b - 2 * a * q
            if r > a:
                r -= 2 * a
                q += 1
            c = c - ((b + r) // 2) * q
            b = r
        if a > c:
            b = -b
            a, c = c, a
            continue
        if a == c and b < 0:
            b = -b
        return a, b, c


@njit(cache=True)
def xgcd(a, b):
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b != 0:
        q = a // b
        a, b = b, a - q * b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return x0, y0, a


@njit(cache=True)
def comp(a1, b1, c1, a2, b2, c2, D):
    if a1 > a2:
        a1, b1, c1, a2, b2, c2 = a2, b2, c2, a1, b1, c1
    s = (b1 + b2) // 2
    n = b2 - s
    if a2 % a1 == 0:
        y1 = 0
        d = a1
    else:
        u, v, d = xgcd(a2, a1)       # u a2 + v a1 = d
        y1 = u
    if s % d == 0:
        y2 = -1
        x2 = 0
        d1 = d
    else:
        x2, y2, d1 = xgcd(s, d)      # x2 s + y2 d = d1
        y2 = -y2
    v1 = a1 // d1
    v2 = a2 // d1
    r = (y1 * y2 * n - x2 * c2) % v1
    b3 = b2 + 2 * v2 * r
    a3 = v1 * v2
    num = b3 * b3 - D
    if num % (4 * a3) != 0:
        return -1, 0, 0              # signals a bug
    c3 = num // (4 * a3)
    return red(a3, b3, c3)


@njit(cache=True)
def pw(a, b, c, e, D):
    b0 = D & 1
    ra, rb, rc = 1, b0, (b0 * b0 - D) // 4
    while e > 0:
        if e & 1:
            ra, rb, rc = comp(ra, rb, rc, a, b, c, D)
        e >>= 1
        if e > 0:
            a, b, c = comp(a, b, c, a, b, c, D)
    return ra, rb, rc


@njit(cache=True)
def forms(D, h):
    n = -D
    F = np.zeros((h, 3), dtype=np.int64)
    k = 0
    b = n % 2
    while 3 * b * b <= n:
        q = (b * b + n) // 4
        a = b if b > 1 else 1
        while a * a <= q:
            if q % a == 0:
                c = q // a
                g = a
                x = b
                while x:
                    g, x = x, g % x
                x = c
                while x:
                    g, x = x, g % x
                if g == 1:
                    F[k, 0] = a; F[k, 1] = b; F[k, 2] = c
                    k += 1
                    if not (a == b or a == c or b == 0):
                        F[k, 0] = a; F[k, 1] = -b; F[k, 2] = c
                        k += 1
            a += 1
        b += 2
    return F, k


@njit(cache=True)
def exp_divides(F, h, e, D):
    b0 = D & 1
    i0 = (1, b0, (b0 * b0 - D) // 4)
    for i in range(h - 1, -1, -1):          # large forms first
        a, b, c = pw(F[i, 0], F[i, 1], F[i, 2], e, D)
        if a < 0:
            return -1
        if a != i0[0] or b != i0[1] or c != i0[2]:
            return 0
    return 1


@njit(cache=True)
def all_ambiguous(F, h):
    for i in range(h):
        a, b, c = F[i, 0], F[i, 1], F[i, 2]
        if not (b == 0 or b == a or a == c):
            return False
    return True


@njit(cache=True)
def order_counts(F, h, D, maxk):
    """cnt[k] = number of classes g with g^k = 1, k = 1..maxk."""
    b0 = D & 1
    cnt = np.zeros(maxk + 1, dtype=np.int64)
    for k in range(1, maxk + 1):
        for i in range(h):
            a, b, c = pw(F[i, 0], F[i, 1], F[i, 2], k, D)
            if a == 1 and b == b0:
                cnt[k] += 1
    return cnt


def prime_factors(m):
    f = set()
    d = 2
    while d * d <= m:
        while m % d == 0:
            f.add(d); m //= d
        d += 1
    if m > 1:
        f.add(m)
    return f


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 300000
    t0 = time.time()
    lists = {e: [] for e in EXPS}
    bugs = 0
    amb_mismatch = 0
    hist = {}
    idoneal = []
    for n in range(3, X + 1):
        if n % 4 not in (0, 3):
            continue
        D = -n
        h = h_forms(D)
        pf = prime_factors(h) if h > 1 else set()
        cand = [e for e in EXPS if pf <= prime_factors(e) and (e % 2 == 0 or h % 2 == 1)]
        if not cand:
            continue
        F, k = forms(D, h)
        assert k == h
        okset = []
        for e in cand:
            r = exp_divides(F, h, e, D)
            if r < 0:
                bugs += 1
            if r == 1:
                okset.append(e)
        if 2 in cand:
            amb = all_ambiguous(F, h)
            if amb != (2 in okset):
                amb_mismatch += 1
        if h == 1:
            okset = list(EXPS)
        for e in okset:
            lists[e].append((n, h, fundamental_neg(D)))
        if 2 in okset and n % 4 == 0:
            idoneal.append(n // 4)
            hist[h] = hist.get(h, 0) + 1
    print('X = %d; composition bugs %d; e=2 ambiguous-criterion mismatches %d; %.1f s' % (X, bugs, amb_mismatch, time.time() - t0))
    known65 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 18, 21, 22, 24, 25, 28, 30, 33, 37, 40, 42, 45, 48, 57, 58,
               60, 70, 72, 78, 85, 88, 93, 102, 105, 112, 120, 130, 133, 165, 168, 177, 190, 210, 232, 240, 253, 273, 280,
               312, 330, 345, 357, 385, 408, 462, 520, 760, 840, 1320, 1365, 1848]
    print('idoneal (D = -4n, exponent | 2): %d, equals Euler 65: %s; class-number histogram %s'
          % (len(idoneal), sorted(idoneal) == known65, dict(sorted(hist.items()))))
    summ = {}
    for e in EXPS:
        L = lists[e]
        fund = [t for t in L if t[2]]
        big = max(L) if L else None
        bigf = max(fund) if fund else None
        summ[e] = dict(count=len(L), fundamental=len(fund), largest=big, largest_fundamental=bigf,
                       fundamental_list=[t[0] for t in fund], all_list=[t[0] for t in L])
        print('e=%d: %d discriminants, %d fundamental, largest %s, largest fundamental %s' % (e, len(L), len(fund), big, bigf),
              flush=True)
    # group structures at the largest members
    for e in EXPS:
        for t in {summ[e]['largest'], summ[e]['largest_fundamental']}:
            if t is None or t[1] == 1:
                continue
            n, h = t[0], t[1]
            F, k = forms(-n, h)
            cnt = order_counts(F, h, -n, e)
            print('   e=%d |D|=%d h=%d: #{g : g^k = 1} for k=1..%d: %s' % (e, n, h, e, cnt[1:].tolist()))
    json.dump(dict(X=X, idoneal=sorted(idoneal), hist=hist, ladder={str(e): summ[e] for e in EXPS}),
              open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'audit_exponent_ladder_%d.json' % X), 'w'))


if __name__ == '__main__':
    main()
