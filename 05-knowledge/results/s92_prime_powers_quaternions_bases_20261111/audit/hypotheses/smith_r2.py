#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9666 (C). For B(u, v) = u^T G v with u, v independent and uniform in (Z/p^k)^2:
# does the class law depend only on the Smith normal form of G over Z_p, and is it the additive convolution of
# the laws of the scaled tables p^(e_i) a_i b_i ?
# Full brute force (no theory used in the counting):
#   modulo 9 : ALL 2 x 2 matrices with entries in 0..8 (6561 matrices), all 9^4 pairs (u, v);
#   modulo 27: all matrices with entries in -4..4 (6561 matrices), all 27^4 pairs (u, v);
#   modulo 25: all matrices with entries in -3..3 (2401 matrices), all 25^4 pairs;
#   modulo 8 : all matrices with entries in 0..7 (4096 matrices), all 8^4 pairs  (the prime 2: same statement).
# Smith form over Z_p of an integer matrix: e1 = min valuation of the entries, e1 + e2 = v_p(det) (infinite if 0),
# truncated at k.  Reference law: counts of p^e1 a1 b1 + p^e2 a2 b2 = c over (a1, b1, a2, b2) modulo p^k.
import sys, time
import numpy as np
from numba import njit, prange


@njit(parallel=True, cache=True)
def laws(mats, q):
    nm = mats.shape[0]
    out = np.zeros((nm, q), dtype=np.int64)
    for i in prange(nm):
        g11 = mats[i, 0] % q
        g12 = mats[i, 1] % q
        g21 = mats[i, 2] % q
        g22 = mats[i, 3] % q
        for u1 in range(q):
            for u2 in range(q):
                w1 = (u1 * g11 + u2 * g21) % q
                w2 = (u1 * g12 + u2 * g22) % q
                for v1 in range(q):
                    t = w1 * v1
                    for v2 in range(q):
                        out[i, (t + w2 * v2) % q] += 1
    return out


def vp(n, p, cap):
    if n == 0:
        return cap
    v = 0
    while n % p == 0 and v < cap:
        n //= p
        v += 1
    return min(v, cap)


def ref_law(p, k, e1, e2):
    q = p ** k
    s1 = (p ** e1) % q if e1 < k else 0
    s2 = (p ** e2) % q if e2 < k else 0
    a = np.arange(q)
    prod = (a[:, None] * a[None, :]) % q
    t1 = np.bincount(((s1 * prod) % q).ravel(), minlength=q)
    t2 = np.bincount(((s2 * prod) % q).ravel(), minlength=q)
    out = np.zeros(q, dtype=np.int64)
    for c1 in range(q):
        if t1[c1]:
            out[(c1 + np.arange(q)) % q] += t1[c1] * t2
    return out


def run(p, k, lo, hi):
    q = p ** k
    rng = range(lo, hi + 1)
    mats = np.array([(a, b, c, d) for a in rng for b in rng for c in rng for d in rng], dtype=np.int64)
    t0 = time.time()
    L = laws(mats, q)
    groups = {}
    for i in range(mats.shape[0]):
        a, b, c, d = (int(x) for x in mats[i])
        e1 = min(vp(a, p, k), vp(b, p, k), vp(c, p, k), vp(d, p, k))
        det = a * d - b * c
        # v_p(det) can exceed 2k - truncate the second divisor at k
        vdet = vp(det, p, 10 ** 6) if det != 0 else 10 ** 6
        e2 = min(vdet - e1, k) if e1 < k else k
        groups.setdefault((e1, e2), []).append(i)
    print("modulo %d = %d^%d, entries in [%d, %d]: %d matrices, %d Smith classes (e1, e2 truncated at %d)   %.1f s" % (
        q, p, k, lo, hi, mats.shape[0], len(groups), k, time.time() - t0))
    allok = True
    for key in sorted(groups):
        idx = groups[key]
        first = L[idx[0]]
        same = all(np.array_equal(L[i], first) for i in idx)
        ref = ref_law(p, k, key[0], key[1])
        eq = np.array_equal(first, ref)
        allok &= same and eq
        # the law as class factors by valuation (q * P(c)), for the record
        fac = {}
        for c in range(q):
            v = vp(c, p, k)
            fac.setdefault(v, set()).add(int(first[c]))
        desc = "; ".join("v=%d: %s" % (v, "/".join("%.4f" % (x * q / q ** 4) for x in sorted(s))) for v, s in sorted(fac.items()))
        print("   Smith (p^%d, p^%d)%s: %5d matrices; one law for all: %s; equals the convolution of the scaled tables: %s;  q*P by valuation: %s" % (
            key[0], key[1], " [p^%d = 0]" % k if k in key else "", len(idx), same, eq, desc))
    return allok


def main():
    print("# smith_r2.py (auditor hypotheses; own code)")
    ok = True
    ok &= run(3, 2, 0, 8)
    ok &= run(3, 3, -4, 4)
    ok &= run(5, 2, -3, 3)
    ok &= run(2, 3, 0, 7)
    # the two instances named in HYP-9666 (C)
    for name, G, p, k in (("[[11,6],[6,3]] (det -3)", (11, 6, 6, 3), 3, 2), ("[[11,6],[6,3]] (det -3)", (11, 6, 6, 3), 3, 3),
                          ("Session261 [[149,72],[72,37]] (det 329)", (149, 72, 72, 37), 7, 1),
                          ("Session261 [[149,72],[72,37]] (det 329)", (149, 72, 72, 37), 7, 2),
                          ("Session261 [[149,72],[72,37]] (det 329)", (149, 72, 72, 37), 47, 1)):
        q = p ** k
        L = laws(np.array([G], dtype=np.int64), q)[0]
        ref = ref_law(p, k, 0, 1)
        print("%s modulo %d: law = law of a1 b1 + %d a2 b2: %s" % (name, q, p, np.array_equal(L, ref)))
        ok &= np.array_equal(L, ref)
    print("ALL AGREE" if ok else "DISAGREEMENT FOUND")


if __name__ == "__main__":
    main()
