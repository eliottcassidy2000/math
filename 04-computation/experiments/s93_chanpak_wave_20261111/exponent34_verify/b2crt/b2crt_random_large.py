#!/usr/bin/env python3
"""
b2crt_random_large.py -- the CRT form of (B2) at sizes where exhaustive enumeration is impossible
(verifier lane b2crt).  Exact integer / rational arithmetic.

For random (s1, s2) with s1 - s2 in {0, 1}, s2 <= SMAX, random P_u = (a, b) (b odd, |a|, |b| < 2^s1, a != 0)
and random box points (a', b') (b' odd, |a'|, |b'| < 2^s2), biased towards points of SMALL height
(the points are generated on a prescribed primitive direction when possible):
  (i)  from the point: direction (k, k'), e = 2^ds k' a' - k a; check 2^s2 e = k b' - k' b,
       |e| < |k| + 2^ds k', e = e0 (mod m) where e0 is computed from the two congruences by CRT,
       and e0/m = -sgn(k) a/(2^ds k') - b inv(2^s1;|k|)/|k| (mod 1).
  (ii) from the class: every e = e0 (mod m), |e| < |k| + 2^ds k' gives integers (a', b'); if they lie in the
       box then their direction is (k, k'); the original point is among them; the number of such e is
       <= 2(|k| + 2^ds k')/m + 1 <= 6M/m + 1.
The same for Deg' (P_v fixed, (a, b) varying) with e = k a - 2^ds k' a', m' = |k| k',
e0/m' = -sgn(k) (2^ds a'/k + b' inv(2^s2; k')/k') (mod 1).

usage: python b2crt_random_large.py NTRIALS SMAX SEED
"""
import sys
import random
from fractions import Fraction
from math import gcd


def sgn(x):
    return 1 if x > 0 else -1


def inv_mod(x, n):
    return 0 if n == 1 else pow(x, -1, n)


def direction(X, Xi):
    g = gcd(abs(X), abs(Xi))
    return (X // g) * sgn(Xi), abs(Xi) // g


def crt(r1, m1, r2, m2):
    """the e in [0, m1 m2) with e = r1 (mod m1), e = r2 (mod m2); gcd(m1, m2) = 1"""
    t = ((r2 - r1) * inv_mod(m1, m2)) % m2
    return (r1 + t * m1) % (m1 * m2)


def is_int(fr):
    return fr.denominator == 1


def one_deg(rng, s1, s2, stats):
    ds = s1 - s2
    B1, B2 = 1 << s1, 1 << s2
    a = rng.randrange(1, B1) * rng.choice((1, -1))
    b = (2 * rng.randrange(0, B1 // 2) + 1) * rng.choice((1, -1))
    # choose a point: half of the time uniformly, half of the time on a short direction
    if rng.random() < 0.5:
        ap = rng.randrange(-B2 + 1, B2)
        bp = (2 * rng.randrange(0, B2 // 2) + 1) * rng.choice((1, -1))
    else:
        # short direction: pick odd coprime (k, k') small, then try to hit the class
        Kb = 1 << rng.randrange(1, s2 + 2)
        while True:
            k = (2 * rng.randrange(0, Kb) + 1) * rng.choice((1, -1))
            kp = 2 * rng.randrange(0, Kb) + 1
            if gcd(abs(k), kp) == 1:
                break
        m1, m2 = (1 << ds) * kp, abs(k)
        e0 = crt((-k * a) % m1, m1, (-kp * b * inv_mod(B2, m2)) % m2, m2)
        m = m1 * m2
        T = abs(k) + m1
        cand = []
        e = -T + 1 + ((e0 - (-T + 1)) % m)
        while e < T:
            ap, bp = (e + k * a) // m1, (B2 * e + kp * b) // k
            assert (e + k * a) % m1 == 0 and (B2 * e + kp * b) % k == 0
            if abs(ap) < B2 and abs(bp) < B2:
                cand.append((ap, bp))
            e += m
        stats["short_dirs"] += 1
        if not cand:
            stats["short_dirs_empty"] += 1
            return
        ap, bp = rng.choice(cand)
    X, Xi = ap * B1 + b, a * B2 + bp
    k, kp = direction(X, Xi)
    assert k % 2 != 0 and kp % 2 == 1
    m1, m2 = (1 << ds) * kp, abs(k)
    m = m1 * m2
    assert gcd(m1, m2) == 1
    e = (1 << ds) * kp * ap - k * a
    T = abs(k) + m1
    assert k * bp - kp * b == B2 * e
    assert abs(e) < T, ("|e| bound", s1, s2, a, b, ap, bp)
    e0 = crt((-k * a) % m1, m1, (-kp * b * inv_mod(B2, m2)) % m2, m2)
    assert (e - e0) % m == 0, ("class", s1, s2, a, b, ap, bp)
    rhs = -sgn(k) * Fraction(a, m1) - Fraction(b * inv_mod(B1, m2), m2)
    assert is_int(Fraction(e0, m) - rhs), ("formula", s1, s2, a, b, k, kp)
    # (ii)
    ee = -T + 1 + ((e0 - (-T + 1)) % m)
    pts = []
    ne = 0
    while ee < T:
        ne += 1
        assert (ee + k * a) % m1 == 0 and (B2 * ee + kp * b) % k == 0
        a2, b2 = (ee + k * a) // m1, (B2 * ee + kp * b) // k
        if abs(a2) < B2 and abs(b2) < B2:
            assert direction(a2 * B1 + b, a * B2 + b2) == (k, kp), ("direction of the image", s1, s2, a, b, k, kp, ee)
            pts.append((a2, b2))
        ee += m
    assert (ap, bp) in pts
    assert len(set(pts)) == len(pts)
    M = max(abs(k), kp)
    assert ne * m <= 2 * T + m and len(pts) * m <= 6 * M + m
    stats["deg"] += 1
    stats["deg_maxpts"] = max(stats["deg_maxpts"], len(pts))
    if M <= B2:
        stats["deg_short"] += 1


def one_degp(rng, s1, s2, stats):
    ds = s1 - s2
    B1, B2 = 1 << s1, 1 << s2
    ap = rng.randrange(1, B2) * rng.choice((1, -1))
    bp = (2 * rng.randrange(0, B2 // 2) + 1) * rng.choice((1, -1))
    if rng.random() < 0.5:
        a = rng.randrange(-B1 + 1, B1)
        b = (2 * rng.randrange(0, B1 // 2) + 1) * rng.choice((1, -1))
    else:
        Kb = 1 << rng.randrange(1, s2 + 2)
        while True:
            k = (2 * rng.randrange(0, Kb) + 1) * rng.choice((1, -1))
            kp = 2 * rng.randrange(0, Kb) + 1
            if gcd(abs(k), kp) == 1:
                break
        m1, m2 = abs(k), kp
        e0 = crt((-(1 << ds) * kp * ap) % m1, m1, (-k * bp * inv_mod(B2, m2)) % m2, m2)
        m = m1 * m2
        T = abs(k) + (1 << ds) * kp
        cand = []
        e = -T + 1 + ((e0 - (-T + 1)) % m)
        while e < T:
            assert (e + (1 << ds) * kp * ap) % k == 0 and (k * bp + B2 * e) % kp == 0
            a, b = (e + (1 << ds) * kp * ap) // k, (k * bp + B2 * e) // kp
            if abs(a) < B1 and abs(b) < B1:
                cand.append((a, b))
            e += m
        stats["short_dirs"] += 1
        if not cand:
            stats["short_dirs_empty"] += 1
            return
        a, b = rng.choice(cand)
    X, Xi = ap * B1 + b, a * B2 + bp
    k, kp = direction(X, Xi)
    m1, m2 = abs(k), kp
    m = m1 * m2
    assert gcd(m1, m2) == 1
    e = k * a - (1 << ds) * kp * ap
    T = abs(k) + (1 << ds) * kp
    assert kp * b == k * bp + B2 * e
    assert abs(e) < T, ("|e| bound'", s1, s2, ap, bp, a, b)
    e0 = crt((-(1 << ds) * kp * ap) % m1, m1, (-k * bp * inv_mod(B2, m2)) % m2, m2)
    assert (e - e0) % m == 0, ("class'", s1, s2, ap, bp, a, b)
    rhs = -sgn(k) * (Fraction((1 << ds) * ap, k) + Fraction(bp * inv_mod(B2, kp), kp))
    assert is_int(Fraction(e0, m) - rhs), ("formula'", s1, s2, ap, bp, k, kp)
    ee = -T + 1 + ((e0 - (-T + 1)) % m)
    pts = []
    ne = 0
    while ee < T:
        ne += 1
        assert (ee + (1 << ds) * kp * ap) % k == 0 and (k * bp + B2 * ee) % kp == 0
        a2, b2 = (ee + (1 << ds) * kp * ap) // k, (k * bp + B2 * ee) // kp
        if abs(a2) < B1 and abs(b2) < B1:
            assert direction(ap * B1 + b2, a2 * B2 + bp) == (k, kp), ("direction of the image'",)
            pts.append((a2, b2))
        ee += m
    assert (a, b) in pts
    assert len(set(pts)) == len(pts)
    M = max(abs(k), kp)
    assert ne * m <= 2 * T + m and len(pts) * m <= 6 * M + m
    stats["degp"] += 1
    stats["degp_maxpts"] = max(stats["degp_maxpts"], len(pts))
    if M <= B2:
        stats["degp_short"] += 1


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    smax = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20261111
    rng = random.Random(seed)
    stats = dict(deg=0, degp=0, deg_short=0, degp_short=0, deg_maxpts=0, degp_maxpts=0, short_dirs=0, short_dirs_empty=0)
    for _ in range(n):
        s2 = rng.randrange(2, smax + 1)
        s1 = s2 + rng.randrange(0, 2)
        one_deg(rng, s1, s2, stats)
        one_degp(rng, s1, s2, stats)
    print(f"trials {n}, s2 in [2, {smax}], seed {seed}")
    print(f"Deg : {stats['deg']} points verified ({stats['deg_short']} of height <= 2^s2), max points on a direction {stats['deg_maxpts']}")
    print(f"Deg': {stats['degp']} points verified ({stats['degp_short']} of height <= 2^s2), max points on a direction {stats['degp_maxpts']}")
    print(f"prescribed short directions tried {stats['short_dirs']}, empty {stats['short_dirs_empty']}")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
