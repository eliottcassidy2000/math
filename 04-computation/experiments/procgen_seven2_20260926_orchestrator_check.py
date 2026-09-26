#!/usr/bin/env python3
"""Orchestrator audit of lane `seven2`, written from the note's statements and
from the certificate lemmas G1, G2, L of the floor note (THM-4486); the lane's
scripts and checkers were not read.

Setting: level k, odd q, N = 2^k, H = 2^(k-1). Pairs P in Z/H with lifts P, P+H.
Options: even node t -> pair t/2; odd node t -> pairs (q t + 1)/2 mod H (sign +)
and (q t - 1)/2 mod H (sign -). Threshold F = fn/fd, weights e(t) = fd - fn
(odd t), -fn (even t).

  1. Independent game solver (numpy, least fixed points on pairs):
     - upper (Min, Lemma L): u(P) = max over lifts t of max(0, e(t) + min over options u);
     - lower (Max, Lemma G2): g(P) = max(0, min over lifts t of (max over options g) - e(t)).
     Both stabilise at F = rho*(q,k) => rho*(q,k) = F exactly. Checked for
     q = 7, k = 8..21 against the published table, and q = 5, k = 13..16
     (5/12 at k = 14, 2/5 at k = 15: "why 15").
  2. The explicit q = 5 rule (MH flipped on 46 classes) has rho_max = 2/5 at
     k = 15, 16, 17, and the explicit q = 7 rule (18 classes) has rho_max = 2/5
     at k = 10, 12, 14, 16: fixed-strategy Kleene iteration at F = 2/5 plus a
     rational periodic point of density 2/5 (Lemma C direction (a)).
  3. Lemma MH: rho_max(MH) = 1/2 for q = 5, 7, 9, 11 at k = 8, 12, 16; the
     S_inf periodic points: exactly 2^n of period dividing n (n <= 8), all
     rational with |x| <= 1/(q-4), all valuations 2.
  4. Lemma C(b): for random strategies at k = 6, 7 (q = 7), every closed walk
     of length <= 9 of G_sigma comes from a rational periodic point with the
     walk's residues.
  5. Lemma F (q = 7) and its q = 5 analogue on random 2-adic integers
     (mod 2^200), and the gaining pattern = classes 11, 21 mod 32 (q = 7).
"""
import math, random, itertools
from fractions import Fraction
import numpy as np


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg, flush=True)


# ------------------------------------------------------------------ arena
def arena(q, k):
    N = 1 << k
    H = N >> 1
    t = np.arange(N, dtype=np.int64)
    odd = (t & 1).astype(bool)
    op = np.where(odd, ((q * t + 1) // 2) % H, t // 2).astype(np.int64)
    om = np.where(odd, ((q * t - 1) // 2) % H, t // 2).astype(np.int64)
    return N, H, odd, op, om


def solve_upper(q, k, fn, fd, maxsweeps=200000):
    N, H, odd, op, om = arena(q, k)
    e = np.where(odd, fd - fn, -fn).astype(np.int64)
    u = np.zeros(H, dtype=np.int64)
    for sweep in range(maxsweeps):
        psi = np.maximum(0, e + np.minimum(u[op], u[om]))
        un = np.maximum(psi[:H], psi[H:])
        if np.array_equal(un, u):
            return True, sweep, int(u.max())
        u = un
        if u.max() > 10 ** 9:
            return False, sweep, int(u.max())
    return False, maxsweeps, int(u.max())


def solve_lower(q, k, fn, fd, maxsweeps=200000):
    N, H, odd, op, om = arena(q, k)
    e = np.where(odd, fd - fn, -fn).astype(np.int64)
    g = np.zeros(H, dtype=np.int64)
    for sweep in range(maxsweeps):
        m = np.maximum(g[op], g[om]) - e
        gn = np.maximum(0, np.minimum(m[:H], m[H:]))
        if np.array_equal(gn, g):
            return True, sweep, int(g.max())
        g = gn
        if g.max() > 10 ** 9:
            return False, sweep, int(g.max())
    return False, maxsweeps, int(g.max())


def verify_upper_cert(q, k, fn, fd, u):
    """exact re-check of Lemma G1 from a pair potential u: psi(t) = max(0, e + min option u); sigma = argmin;
    need psi(t') + e(t) <= psi(t) for both lifts t' of sigma(t)."""
    N, H, odd, op, om = arena(q, k)
    e = np.where(odd, fd - fn, -fn).astype(np.int64)
    psi = np.maximum(0, e + np.minimum(u[op], u[om]))
    ch = np.where(u[op] <= u[om], op, om)
    return bool(np.all(np.maximum(psi[ch], psi[ch + H]) + e <= psi))


table7 = {}
for kk, v in [(7, (3, 7)), (10, (2, 5)), (14, (15, 38)), (16, (7, 18)), (18, (19, 49)), (19, (13, 34)),
              (20, (34, 89)), (21, (8, 21))]:
    table7[kk] = v
vals = {}
cur = None
for kk in range(7, 22):
    if kk in table7:
        cur = table7[kk]
    vals[kk] = cur
for kk in range(8, 22):
    fn, fd = vals[kk]
    okU, sU, mU = solve_upper(7, kk, fn, fd)
    okL, sL, mL = solve_lower(7, kk, fn, fd)
    assert okU and okL, (kk, fn, fd, okU, okL)
    print(f"    q=7 k={kk}: rho* = {fn}/{fd}: upper lfp in {sU} sweeps (max {mU}), lower lfp in {sL} sweeps (max {mL})", flush=True)
check(True, "independent game solver: rho*(7,k) equals the published values at every k = 8..21 (both least fixed points finite at F = rho*)")
for kk, (fn, fd) in [(13, (5, 12)), (14, (5, 12)), (15, (2, 5)), (16, (2, 5))]:
    okU, _, _ = solve_upper(5, kk, fn, fd)
    okL, _, _ = solve_lower(5, kk, fn, fd)
    assert okU and okL, (kk,)
check(True, "q = 5: rho*(5,13) = rho*(5,14) = 5/12 > 2/5 and rho*(5,15) = rho*(5,16) = 2/5 exactly (both least fixed points finite): level 15 is needed")

# explicit certificate re-check at one level (exact integers)
N, H, odd, op, om = arena(7, 16)
fn, fd = 7, 18
e = np.where(odd, fd - fn, -fn).astype(np.int64)
u = np.zeros(H, dtype=np.int64)
while True:
    psi = np.maximum(0, e + np.minimum(u[op], u[om]))
    un = np.maximum(psi[:H], psi[H:])
    if np.array_equal(un, u):
        break
    u = un
check(verify_upper_cert(7, 16, 7, 18, u), "the q = 7, k = 16 least Min potential is an exact Lemma G1 certificate for 7/18")


# ------------------------------------------------------------------ fixed strategies
def mh_sign(q, x):
    """sign s with q x + s = 0 mod 4 (x odd)."""
    return 1 if (q * x + 1) % 4 == 0 else -1


def rule_sign(q, flips, x):
    s = mh_sign(q, x)
    for (mod, cls) in flips:
        if x % mod == cls:
            return -s
    return s


def strategy_array(q, k, flips):
    N = 1 << k
    sg = np.zeros(N, dtype=np.int64)
    for x in range(1, N, 2):
        sg[x] = rule_sign(q, flips, x)
    return sg


def rho_le(q, k, sg, fn, fd, maxsweeps=100000):
    N = 1 << k
    H = N >> 1
    t = np.arange(N, dtype=np.int64)
    odd = (t & 1).astype(bool)
    P = np.where(odd, ((q * t + sg) // 2) % H, t // 2)
    e = np.where(odd, fd - fn, -fn).astype(np.int64)
    psi = np.zeros(N, dtype=np.int64)
    for _ in range(maxsweeps):
        pn = np.maximum(0, e + np.maximum(psi[P], psi[P + H]))
        if np.array_equal(pn, psi):
            return True
        psi = pn
        if psi.max() > 10 ** 8:
            return False
    return False


def periodic_density(q, flips, x0, maxsteps=200):
    """iterate T_sigma on the rational x0 (sigma read from the 2-adic residue); return (odd steps, period) if periodic."""
    x = Fraction(x0)
    odds = 0
    for step in range(1, maxsteps + 1):
        num, den = x.numerator, x.denominator
        assert den % 2 == 1
        if num % 2 == 0:
            x = x / 2
        else:
            # residue of x mod 2^40 decides the sign (rules read <= 15 bits)
            r = (num * pow(den, -1, 1 << 40)) % (1 << 40)
            s = rule_sign(q, flips, r)
            x = (q * x + s) / 2
            odds += 1
        if x == Fraction(x0):
            return odds, step
    return None


flips5 = []
for mod_e, cl in [(6, [31, 33]), (8, [57, 199]), (9, [1, 65, 111, 401, 447, 511]), (10, [7, 71, 953, 1017]),
                  (11, [63, 191, 327, 337, 761, 1287, 1711, 1721, 1857, 1985]), (12, [575, 879, 1345, 2751, 3217, 3521]),
                  (13, [2735, 3921, 4271, 5457]), (14, [1873, 2927, 3409, 5999, 10385, 12975, 13457, 14511]),
                  (15, [7023, 10095, 22673, 25745])]:
    for c in cl:
        flips5.append((1 << mod_e, c))
assert len(flips5) == 46
# disjoint and negation-closed
M15 = 1 << 15
cover = {}
for (mod, c) in flips5:
    for r in range(c, M15, mod):
        assert r not in cover, ("overlap", mod, c)
        cover[r] = (mod, c)
    assert any(m2 == mod and (c + c2) % mod == 0 for (m2, c2) in flips5), ("negation", mod, c)
for kk in (15, 16, 17):
    sg = strategy_array(5, kk, flips5)
    assert rho_le(5, kk, sg, 2, 5), kk
pd = periodic_density(5, flips5, 1)
check(pd == (2, 5), f"q = 5 explicit rule (46 flipped classes, disjoint, negation-closed): rho_max <= 2/5 at k = 15, 16, 17 and the periodic point 1 has density {pd[0]}/{pd[1]} (the sporadic cycle 1,3,8,4,2)")

flips7 = [(32, 11), (32, 21)] + [(256, c) for c in (35, 67, 91, 165, 189, 221)] + \
         [(1024, c) for c in (93, 157, 349, 381, 413, 611, 643, 675, 867, 931)]
assert len(flips7) == 18
for kk in (10, 12, 14, 16):
    sg = strategy_array(7, kk, flips7)
    assert rho_le(7, kk, sg, 2, 5), kk
pd7 = periodic_density(7, flips7, Fraction(-9, 17))
check(pd7 is not None and Fraction(pd7[0], pd7[1]) == Fraction(2, 5),
      f"q = 7 explicit rule (18 flipped classes): rho_max <= 2/5 at k = 10, 12, 14, 16; -9/17 is a periodic point of density {pd7}")
check(all(x % 32 in (11, 21) for x in (pow(3, -1, 1 << 20), (-pow(3, -1, 1 << 20)) % (1 << 20))), "the classes of 1/3 and -1/3 mod 32 are 11 and 21")

# ------------------------------------------------------------------ Lemma MH
for q in (5, 7, 9, 11):
    for kk in (8, 12, 16):
        N = 1 << kk
        sg = np.zeros(N, dtype=np.int64)
        for x in range(1, N, 2):
            sg[x] = mh_sign(q, x)
        assert rho_le(q, kk, sg, 1, 2), (q, kk)
    x0 = Fraction(-1, q - 4)
    assert periodic_density(q, [], x0) == (1, 2) or periodic_density(q, [], Fraction(1, q - 4)) == (1, 2), q
check(True, "Lemma MH: rho_max(MH) <= 1/2 at k = 8, 12, 16 and a periodic point +-1/(q-4) of density 1/2, for q = 5, 7, 9, 11")


def accel_mh(q, x):
    s = mh_sign(q, x.numerator * pow(x.denominator, -1, 1 << 60) % (1 << 60))
    y = q * x + s
    v = 0
    while y.numerator % 2 == 0:
        y /= 2
        v += 1
    return s, v, y


for q in (5, 7):
    for n in range(1, 9):
        pts = set()
        for word in itertools.product((1, -1), repeat=n):
            # fixed point of g_{s_1} o ... o g_{s_n}, g_s(y) = (4y - s)/q: solve x = G(x) exactly
            a, b = Fraction(1), Fraction(0)          # G(x) = a x + b
            for s in reversed(word):
                a, b = a * Fraction(4, q), (b * 4 - s) / q
            x = b / (1 - a)
            # x must be the MH-periodic point with signs word and valuations 2
            y = x
            for s in word:
                s2, v, y = accel_mh(q, y)
                assert v == 2 and s2 == s, (q, n, word)
            assert y == x, (q, n, word)
            assert abs(x) <= Fraction(1, q - 4)
            pts.add(x)
        assert len(pts) == 2 ** n, (q, n, len(pts))
check(True, "S_inf: exactly 2^n periodic points of period dividing n (n = 1..8, q = 5, 7), all with valuations 2 and |x| <= 1/(q-4)")

# ------------------------------------------------------------------ Lemma C(b)
random.seed(11)
walks = 0
for kk in (6, 7):
    N = 1 << kk
    H = N >> 1
    for trial in range(3):
        sig = {x: random.choice((1, -1)) for x in range(1, N, 2)}

        def succ(w):
            P = (w // 2) if w % 2 == 0 else ((7 * w + sig[w]) // 2) % H
            return (P, P + H)
        for start in range(N):
            stack = [(start, [start])]
            while stack:
                w, path = stack.pop()
                if len(path) > 9:
                    continue
                for nx in succ(w):
                    if nx == start:
                        # closed walk path -> start; build the periodic point
                        a, b = Fraction(1), Fraction(0)     # x_0 = a x_p + b composition of inverse branches
                        for wi in reversed(path):
                            if wi % 2 == 0:
                                a, b = 2 * a, 2 * b
                            else:
                                a, b = 2 * a / 7, (2 * b - sig[wi]) / 7
                        x0 = b / (1 - a)
                        x = x0
                        for wi in path:
                            r = x.numerator * pow(x.denominator, -1, N) % N
                            assert r == wi, (kk, path, x0)
                            x = x / 2 if wi % 2 == 0 else (7 * x + sig[wi]) / 2
                        assert x == x0
                        walks += 1
                    else:
                        stack.append((nx, path + [nx]))
check(walks > 1000, f"Lemma C(b): {walks} closed walks (length <= 9, k = 6, 7, random strategies, q = 7) each come from a rational periodic point with exactly the walk's residues")

# ------------------------------------------------------------------ Lemma F
random.seed(5)
MOD = 1 << 200


def v2(z):
    return (z & -z).bit_length() - 1


def itin(q, x, n):
    out = []
    for _ in range(n):
        s = 1 if (q * x + 1) % 4 == 0 else -1
        z = q * x + s
        v = v2(z)
        out.append((s, v))
        x = z >> v
    return out


gain_ok = 0
for _ in range(100000):
    x = random.randrange(1, MOD, 2)
    (s1, v1), (s2, v2_), (s3, v3) = itin(7, x, 3)
    y = (7 * x - s1) // 2
    sy = 1 if (7 * y + 1) % 4 == 0 else -1
    vy = v2(7 * y + sy)
    if v1 == 2 and s2 != s1:
        pred = 2
    elif v1 == 2 and s2 == s1 and v2_ >= 3:
        pred = 3
    elif v1 == 2 and s2 == s1 and v2_ == 2 and s3 != s1:
        pred = 4
    elif v1 == 2 and s2 == s1 and v2_ == 2 and s3 == s1:
        pred = 5 if v3 == 2 else (6 if v3 >= 4 else ">=7")
    elif v1 == 3:
        pred = 2
    elif v1 == 4 and s2 != s1:
        pred = 3 + v2_
    elif v1 == 4 and s2 == s1:
        pred = 4
    else:
        pred = 3
    assert (vy >= 7) if pred == ">=7" else (vy == pred), (x % 4096, pred, vy)
    gain = (1 + vy) - (v1 + v2_)
    gaining = (v1 == 2 and v2_ == 2 and s1 == s2)
    neutral = (v1 == 4 and s2 == -s1)
    assert (gain >= 1) == gaining and (gain == 0) == neutral, (gain, v1, v2_, s1, s2)
    if gaining:
        assert x % 32 in (11, 21)
    gain_ok += 1
check(True, f"Lemma F (q = 7): the valuation table and the gain classification hold on {gain_ok} random 2-adic integers; the gaining pattern lies in classes 11, 21 mod 32")
for _ in range(100000):
    x = random.randrange(1, MOD, 2)
    (s1, v1), (s2, v2_), (s3, v3) = itin(5, x, 3)
    y = (5 * x - s1) // 2
    sy = 1 if (5 * y + 1) % 4 == 0 else -1
    vy = v2(5 * y + sy)
    if v1 == 2 and s2 != s1:
        ok = vy == 2
    elif v1 == 2 and s2 == s1 and v2_ >= 3:
        ok = vy == 3
    elif v1 == 2 and s2 == s1 and v2_ == 2 and s3 == s1:
        ok = vy == 4
    elif v1 == 2 and s2 == s1 and v2_ == 2 and s3 != s1:
        ok = vy >= 5
    elif v1 == 3 and s2 != s1:
        ok = vy == 2 + v2_
    elif v1 == 3 and s2 == s1:
        ok = vy == 3
    else:
        ok = vy == 2
    assert ok, (x % 4096, v1, v2_, s1, s2, s3, vy)
check(True, "Lemma F, q = 5 analogue (remark c): the valuation table holds on 100000 random 2-adic integers")
