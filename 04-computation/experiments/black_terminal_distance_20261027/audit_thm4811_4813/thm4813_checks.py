"""Independent auditor checks for THM-4813 (A), (B), (C) and THM-4811 (A) on random tournaments.

Written from the theorem statements only (no import of the repository scripts).
"""
from fractions import Fraction as Fr
import math
import random
import sys

import mpmath
import numpy as np

mpmath.mp.dps = 60
out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


# ----------------------------------------------------------------------------------------------------------------
# 2-adic helpers on rationals with odd denominator
def v2(x):
    """2-adic valuation of a nonzero rational."""
    x = Fr(x)
    p, q = x.numerator, x.denominator
    assert q % 2 == 1 and p != 0
    return (p & -p).bit_length() - 1


def mod8(x):
    x = Fr(x)
    return (x.numerator * pow(x.denominator, -1, 8)) % 8


def is_odd(x):
    return v2(x) == 0


def U(x):
    """Syracuse map on odd 2-adic rationals."""
    z = 3 * Fr(x) + 1
    return z / (2 ** v2(z))


def B(x):
    return 4 * Fr(x) + 1


def rho_strip(x):
    """primitive root by stripping B-layers: while x = 5 mod 8, x -> (x-1)/4"""
    x = Fr(x)
    j = 0
    while mod8(x) == 5:
        x = (x - 1) / 4
        j += 1
    return x, j


def rho_val(x):
    """primitive root by valuations: 3x+1 = 4^j (3y+1) with v2(3y+1) in {1,2}"""
    x = Fr(x)
    w = v2(3 * x + 1)
    j = (w - 1) // 2
    y = ((3 * x + 1) / Fr(4) ** j - 1) / 3
    return y, j


def a_of(s):
    return 2 if mod8(s) == 1 else 1


def native(s):
    """native compressed seed map on a primitive odd s: returns (y, a, j, t)"""
    s = Fr(s)
    assert is_odd(s) and mod8(s) != 5
    a = a_of(s)
    z = (3 * s + 1) / 2 ** a
    assert is_odd(z), (s, z)
    y, j = rho_strip(z)
    y2, j2 = rho_val(z)
    assert (y, j) == (y2, j2), (s, y, j, y2, j2)
    return y, a, j, a + 2 * j


# ---------------- A1 ------------------------------------------------------------------------------------------
def check_A1(N=2 * 10 ** 6):
    bad = 0
    cnt = 0
    for s in range(1, N, 2):
        if s % 8 == 5:
            continue
        a = 2 if s % 8 == 1 else 1
        z = (3 * s + 1) >> a
        assert z % 2 == 1
        # valuation route for the root
        w = ((3 * z + 1) & -(3 * z + 1)).bit_length() - 1
        j = (w - 1) // 2
        Y = (3 * z + 1) >> (2 * j)
        assert (Y - 1) % 3 == 0
        y = (Y - 1) // 3
        t = a + 2 * j
        X = 3 * s + 1
        cnt += 1
        if (1 << t) * Y != 3 * X + (1 << a):
            bad += 1
        if y % 8 == 5 or y % 2 == 0:
            bad += 1
        # Y = (3X/2^a + 1)/4^j
        if Fr(Fr(3 * X, 1 << a) + 1, 4 ** j) != Y:
            bad += 1
    say(f"A1: 2^t Y = 3X + 2^a, Y = (3X/2^a+1)/4^j, P(s) primitive, for {cnt} primitive s < {N} (root via valuations): failures {bad}")
    bad = 0
    for x in range(-10001, 10002, 2):
        X = 3 * x + 1
        if 3 * (4 * x + 1) + 1 != 4 * X:
            bad += 1
        if X != 0:
            Ux = U(x)
            if 3 * Ux + 1 != Fr(3 * X, 2 ** v2(X)) + 1:
                bad += 1
            if x != -1 and 4 * x + 1 != 0 and U(4 * x + 1) != Ux:
                bad += 1
    say(f"A1: 3B(x)+1 = 4(3x+1), 3U(x)+1 = 3X/2^v+1, U(B(x)) = U(x) for odd |x| <= 10001: failures {bad}")


# ---------------- A2 ------------------------------------------------------------------------------------------
def check_A2():
    bad = 0
    rng = random.Random(4813)
    for s in list(range(1, 20001, 2)) + [rng.randrange(1, 10 ** 30, 2) for _ in range(2000)]:
        if s % 8 == 5:
            continue
        y, a, j, t = native(s)
        X, Y = 3 * Fr(s) + 1, 3 * y + 1
        b = 3 * X / (2 ** t * Y)
        lam = (1 - b) / (1 + b)
        if lam != Fr(2 ** a) / (6 * X + 2 ** a):
            bad += 1
        if 1 / b != 1 + Fr(2 ** a) / (3 * X):
            bad += 1
        if (1 + lam) / (1 - lam) != 1 / b:          # 2 artanh(lam) = -log b
            bad += 1
        if 2 ** t * Y / (3 * X) != 1 / b:            # = exp(t log 2 - log 3 + log Y - log X)
            bad += 1
    say(f"A2: lambda = 2^a/(6X+2^a), 1/b = 1 + 2^a/(3X) = 2^t Y/(3X), (1+lam)/(1-lam) = 1/b (exact rationals): failures {bad}")
    # carry formula along random words (native orbits of random seeds, arbitrary length)
    bad = 0
    for trial in range(3000):
        s = rng.randrange(1, 10 ** 12, 2)
        while s % 8 == 5:
            s = rng.randrange(1, 10 ** 12, 2)
        d = rng.randrange(1, 60)
        X0 = 3 * s + 1
        T = 0
        terms = []
        cur = Fr(s)
        prodb = Fr(1)
        for i in range(1, d + 1):
            y, a, j, t = native(cur)
            X, Y = 3 * cur + 1, 3 * y + 1
            terms.append((T, a))
            prodb *= 3 * X / (2 ** t * Y)
            T += t
            cur = y
            if cur == 1:
                d = i
                break
        Cd = sum(3 ** (d - i) * 2 ** (Ti + ai) for i, (Ti, ai) in enumerate(terms, start=1))
        Xd = 3 * cur + 1
        if 2 ** T * Xd != 3 ** d * X0 + Cd or Cd <= 0:
            bad += 1
        if prodb != Fr(3 ** d * X0) / (2 ** T * Xd):
            bad += 1
    say(f"A2: carry 2^T X_d = 3^d X_0 + C_d with C_d = sum 3^(d-i) 2^(T_(i-1)+a_i), prod b = 3^d X_0/(2^T X_d), 3000 random words: failures {bad}")


# ---------------- A3: cycles ----------------------------------------------------------------------------------
def syracuse_cycles_k(k, bound):
    """integer cycles of n -> (3n+k)/2^v (odd n, |n| <= bound); returned as rational Syracuse cycles {n/k}"""
    cycles = set()
    for n0 in range(-bound | 1, bound + 1, 2):
        if n0 % 2 == 0:
            continue
        if math.gcd(n0, k) != 1:
            continue
        x = Fr(n0, k)
        if 3 * x + 1 == 0:
            continue
        seen = {}
        cur = x
        steps = 0
        while cur not in seen and abs(cur.numerator) < 10 ** 8 and steps < 2000:
            seen[cur] = steps
            cur = U(cur)
            steps += 1
            if 3 * cur + 1 == 0:
                break
        if cur in seen:
            # extract cycle
            cyc = []
            c = cur
            while True:
                cyc.append(c)
                c = U(c)
                if c == cur:
                    break
            m = min(range(len(cyc)), key=lambda i: cyc[i])
            cycles.add(tuple(cyc[m:] + cyc[:m]))
    return cycles


def native_cycles_k(k, bound):
    """cycles of the native map P on primitive odd rationals n/k (|n| <= bound)"""
    cycles = set()
    for n0 in range(-bound | 1, bound + 1, 2):
        if math.gcd(n0, k) != 1:
            continue
        x = Fr(n0, k)
        if mod8(x) == 5 or 3 * x + 1 == 0:
            continue
        seen = {}
        cur = x
        steps = 0
        while cur not in seen and abs(cur.numerator) < 10 ** 8 and steps < 2000:
            seen[cur] = steps
            cur = native(cur)[0]
            steps += 1
            if 3 * cur + 1 == 0:
                break
        if cur in seen:
            cyc = []
            c = cur
            while True:
                cyc.append(c)
                c = native(c)[0]
                if c == cur:
                    break
            m = min(range(len(cyc)), key=lambda i: cyc[i])
            cycles.add(tuple(cyc[m:] + cyc[:m]))
    return cycles


def canon_cycle(c):
    m = min(range(len(c)), key=lambda i: c[i])
    return tuple(c[m:] + c[:m])


def check_A3():
    ks = [1, 5, 7, 11, 13, 17, 19, 23, 25, 29, 31, 35, 37, 41, 43, 47]
    total_s = total_n = 0
    bad = 0
    normal = []
    neg_ineq = []
    for k in ks:
        bound = 3000 if k > 1 else 20000
        S = syracuse_cycles_k(k, bound)
        Nc = native_cycles_k(k, bound)
        # map native -> Syracuse by z_i = U(y_(i-1)); Syracuse -> native by rho
        img = set()
        for cyc in Nc:
            d = len(cyc)
            T = 0
            for y in cyc:
                T += native(y)[3]
            z = canon_cycle([U(y) for y in cyc])
            # valuations
            vs = sum(v2(3 * zz + 1) for zz in z)
            if len(set(z)) != d or vs != T:
                bad += 1
            # it is a Syracuse cycle
            if any(U(z[i]) != z[(i + 1) % d] for i in range(d)):
                bad += 1
            img.add(z)
            # forcing identities on the cycle
            lam = [mpmath.mpf(2 ** native(y)[1]) / (6 * (3 * mpmath.mpf(y.numerator) / y.denominator + 1) + 2 ** native(y)[1]) for y in cyc]
            S2 = sum(2 * mpmath.atanh(l) for l in lam)
            S1 = sum(mpmath.atanh(l) for l in lam)
            S0 = sum(lam)
            D = T * mpmath.log(2) - d * mpmath.log(3)
            if abs(S2 - D) > mpmath.mpf(10) ** -50:
                bad += 1
            if not (S0 < S1):
                if all(l > 0 for l in lam):
                    bad += 1
                else:
                    neg_ineq.append((k, d, T))
            # factor-by-factor Steiner-Eliahou: 1 + 2^a/(3X) at y equals 1 + 1/(3 U(y))
            for y in cyc:
                a = native(y)[1]
                if 1 + Fr(2 ** a) / (3 * (3 * y + 1)) != 1 + 1 / (3 * U(y)):
                    bad += 1
            if abs(float(D)) > 0:
                normal.append((k, d, T, float(S0), float(S1), float(S2), float(D)))
        rhoimg = set()
        for cyc in S:
            r = [rho_strip(z)[0] for z in cyc]
            rhoimg.add(canon_cycle(r))
        if img != set(S) or rhoimg != set(Nc):
            bad += 1
            say(f"  MISMATCH k={k}: native {len(Nc)} syracuse {len(S)}")
        total_s += len(S)
        total_n += len(Nc)
    say(f"A3: rational cycles x = n/k, k in {ks}: {total_s} Syracuse cycles, {total_n} native cycles; bijection z_i = U(y_(i-1)) / rho,"
        f" same (d, T), sum 2artanh(lam) = T log2 - d log3 at 50 digits, factorwise Steiner-Eliahou: failures {bad}")
    say(f"  (cycles with some negative lambda, where sum lam < sum artanh lam fails as expected since artanh(x) < x for x < 0: {len(neg_ineq)};"
        f" positive-lambda cycles all satisfy sum lam < sum artanh lam = (T log2 - d log3)/2)")
    say("  sample (k, d, T, sum lam, sum artanh lam, sum 2artanh lam, T log2 - d log3):")
    for row in sorted(set(normal))[:8]:
        say("   ", tuple(round(v, 9) if isinstance(v, float) else v for v in row))
    # trivial cycle
    y, a, j, t = native(1)
    lam = Fr(2 ** a, 6 * 4 + 2 ** a)
    say(f"A3: trivial cycle P(1) = {y}, a = {a}, t = {t}, lambda = {lam}, b = {Fr(12, 2 ** t * 4)}; (1+lam)/(1-lam) = {(1 + lam) / (1 - lam)};"
        f" so 2 artanh(1/7) = log(4/3) = T log2 - d log3 with (d,T) = (1,2); artanh(1/7) = {mpmath.nstr(mpmath.atanh(mpmath.mpf(1)/7), 12)}")


def check_window():
    # R2 per-clock bound on exterior seeds s >= 1009 and the artanh ratio
    worst = Fr(0)
    for s in range(1009, 400001, 2):
        if s % 8 == 5:
            continue
        y, a, j, t = native(s)
        lam = Fr(2 ** a, 6 * (3 * s + 1) + 2 ** a)
        r = lam / t
        if r > worst:
            worst = r
            arg = s
    lmax = Fr(2, 9 * 1009 + 5)
    ratio = 2 * mpmath.atanh(mpmath.mpf(lmax.numerator) / lmax.denominator) / (mpmath.mpf(lmax.numerator) / lmax.denominator)
    say(f"A3 window: max lambda/t over primitive 1009 <= s <= 400001 is {worst} = 1/{1/worst} at s = {arg} (packet R2: 1/9139 at 1015);"
        f" max exterior lambda <= 2/(9*1009+5), 2artanh(l)/l <= {mpmath.nstr(ratio, 15)} < 2 + 1e-6")


def check_A4():
    rows = []
    bad = 0
    for a in (1, 2):
        for j in range(0, 6):
            t = a + 2 * j
            Xs = Fr(2 ** a, 2 ** t - 3)
            if (3 * Xs + 2 ** a) / 2 ** t != Xs:
                bad += 1
            ss = (Xs - 1) / 3
            if ss != Fr(2 ** a - 2 ** t + 3, 3 * (2 ** t - 3)):
                bad += 1
            rows.append((a, j, t, ss))
    anchors = {a: Fr(1, 2 ** a - 3) for a in range(1, 11)}
    hits = {(a, j): [b for b, c in anchors.items() if c == ss] for (a, j, t, ss) in rows}
    say(f"A4: native-edge fixed points (a, j, t, s*): {[(a, j, t, str(ss)) for (a, j, t, ss) in rows]}; failures {bad}")
    say(f"A4: which THM-4600 anchors c_b (b = 1..10) are native-edge fixed points: {[(k, v) for k, v in hits.items() if v]}"
        f" -> only c_1 = -1 and c_2 = 1 (j = 0); c_3 = 1/5, c_4 = 1/13, ... are fixed points of Syracuse letters X -> 3X/2^b + 1, not of native edges")
    # c_b are fixed points of the Syracuse letters in numerator coordinate
    ok = all((3 * (1 + 3 * c) / 2 ** b + 1) == 1 + 3 * c for b, c in anchors.items())
    say(f"A4: X* = 1 + 3c_b is the fixed point of X -> 3X/2^b + 1 for b = 1..10: {ok}")


# ---------------- B ---------------------------------------------------------------------------------------------
def polymul_mod_phi(p1, p2, q):
    """multiply polynomials in zeta (coeff lists length q) modulo zeta^q - 1, then reduce mod Phi_q (q prime)"""
    r = [0] * q
    for i, c in enumerate(p1):
        if c:
            for j, d in enumerate(p2):
                if d:
                    r[(i + j) % q] += c * d
    # reduce modulo 1 + z + ... + z^(q-1): subtract r[q-1] * Phi
    c = r[q - 1]
    return [x - c for x in r[:q - 1]] + [0]


def is_qr(x, q):
    return pow(x % q, (q - 1) // 2, q) == 1


def check_B():
    bad = 0
    for q in [3, 7, 11, 19, 23, 31, 43, 47]:
        S = [s for s in range(1, q) if is_qr(s, q)]
        assert all((-s) % q not in S for s in S)
        g = [0] * q
        for s in S:
            g[s] += 1
            g[(-s) % q] -= 1
        g2 = polymul_mod_phi(g, g, q)
        if g2 != [-q] + [0] * (q - 1):
            bad += 1
        # matrix check
        K = np.zeros((q, q), dtype=np.int64)
        for u in range(q):
            for s in S:
                K[u, (u + s) % q] = 1
                K[(u + s) % q, u] = -1
        if not np.array_equal(K @ K.T, q * np.eye(q, dtype=np.int64) - np.ones((q, q), dtype=np.int64)):
            bad += 1
        # K^2 = -q on Lambda_q: K^2 x = -q x for x with sum 0
        Lb = np.zeros((q, q - 1), dtype=np.int64)
        for i in range(q - 1):
            Lb[i, i] = 1
            Lb[q - 1, i] = -1
        if not np.array_equal(K @ K @ Lb, -q * Lb):
            bad += 1
        # K equals sum_{s in S} (R^s - R^{-s}) with (Rx)_i = x_(i+1)
        R = np.zeros((q, q), dtype=np.int64)
        for i in range(q):
            R[i, (i + 1) % q] = 1
        K2 = sum(np.linalg.matrix_power(R, s) - np.linalg.matrix_power(R, q - s) for s in S)
        if not np.array_equal(K, K2):
            bad += 1
    say(f"B1: Paley q in [3,7,11,19,23,31,43,47]: G^2 = -q in Z[zeta], K K^T = qI - J, K^2 = -q on Lambda_q, K = sum (R^s - R^-s): failures {bad}")
    # nuclear norm formula for general circulant tournaments
    bad = 0
    rng = random.Random(1)
    for q in [5, 7, 9, 11, 13, 15, 17, 21]:
        for trial in range(20):
            S = []
            for s in range(1, (q - 1) // 2 + 1):
                S.append(s if rng.random() < 0.5 else q - s)
            K = np.zeros((q, q))
            for u in range(q):
                for s in S:
                    K[u, (u + s) % q] = 1
                    K[(u + s) % q, u] = -1
            nuc = np.linalg.svd(K, compute_uv=False).sum()
            form = sum(abs(sum(2 * math.sin(2 * math.pi * j * s / q) for s in S)) for j in range(1, q))
            if abs(nuc - form) > 1e-8:
                bad += 1
    say(f"B1: ||K_C(S)||_* = sum_j |sum_s 2 sin(2 pi j s/q)| on 160 random circulant tournaments (q = 5..21, incl. composite q): failures {bad}")
    # B2
    def iota(a, b):
        return np.array([a, -b, b - a], dtype=np.int64)
    R = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=np.int64)    # (Rx)_i = x_(i+1)
    K = R - R @ R                                                     # C3 with S = {1}: u -> u+1
    bad = 0
    for (a, b) in [(1, 0), (0, 1)] + [(rng.randrange(-50, 50), rng.randrange(-50, 50)) for _ in range(200)]:
        # (2 + w)(a + b w) = (2a - b) + (a + b) w
        root_rows = np.array([2 * a - b, -a - b, -a + 2 * b], dtype=np.int64)     # 2 Re(w^r A), r = 0, 1, 2
        w = complex(-0.5, math.sqrt(3) / 2)
        A = a + b * w
        numer = np.array([round(2 * (w ** r * A).real) for r in range(3)])
        if not np.array_equal(root_rows, numer):
            bad += 1
        if not np.array_equal(iota(2 * a - b, a + b), root_rows):
            bad += 1
        if not np.array_equal(-R @ K @ iota(a, b), root_rows):
            bad += 1
        # iota intertwines multiplication by w with R and by 1 + 2w with K
        if not np.array_equal(R @ iota(a, b), iota(-b, a - b)):            # w(a + b w) = -b + (a - b) w
            bad += 1
        if not np.array_equal(K @ iota(a, b), iota(a - 2 * b, 2 * a - b)):  # (1+2w)(a+bw) = (a - 2b) + (2a - b) w
            bad += 1
    say(f"B2: ROOT rows (2a-b, -a-b, -a+2b) = 2Re(w^r A) = iota((2+w)(a+bw)) = -R K iota(a+bw); iota intertwines w <-> R, 1+2w <-> K: failures {bad};"
        f" -RK iota(1) = {list(-R @ K @ iota(1, 0))}")
    # eigenvalues: shared modulus, different eigenvalues
    ev_root = np.linalg.eigvals((-R @ K)[:, :])
    ev_K = np.linalg.eigvals(K.astype(float))
    say(f"B2: eigenvalues of -RK on C^3: {np.round(ev_root, 6)}; of K: {np.round(ev_K, 6)} (same moduli sqrt3 on Lambda_3, different arguments)")


# ---------------- C ---------------------------------------------------------------------------------------------
def check_C():
    rows = []
    cross = None
    for q in [p for p in range(7, 400) if all(p % d for d in range(2, int(p ** 0.5) + 1)) and p % 4 == 3]:
        nuc_T = sum(abs(1 / math.tan((2 * k - 1) * math.pi / (2 * q))) for k in range(1, q + 1))
        spec = (q * (q - 1) - math.sqrt(q) * nuc_T) / 4
        rear = (q * q - 1) / 8
        rows.append((q, round(spec, 1), rear))
        if cross is None and spec > rear:
            cross = q
    later_fail = [r for r in rows if r[0] > (cross or 0) and r[1] <= r[2]]
    say(f"C: spectral FAS lower bound vs (q^2-1)/8: {rows[:8]}; first q with spectral > rearrangement: {cross}; later failures up to 400: {later_fail}")
    # transitive skew spectrum and Paley op norm
    for q in [7, 11, 19]:
        T = np.zeros((q, q))
        for i in range(q):
            for j in range(i + 1, q):
                T[i, j] = 1
                T[j, i] = -1
        op = np.linalg.norm(T, 2)
        nucT = np.linalg.svd(T, compute_uv=False).sum()
        say(f"C: q={q}: ||K_T||_op = {op:.6f} vs cot(pi/2q) = {1 / math.tan(math.pi / (2 * q)):.6f}; ||K_T||_* = {nucT:.6f} vs sum|cot| = {sum(abs(1 / math.tan((2 * k - 1) * math.pi / (2 * q))) for k in range(1, q + 1)):.6f};"
            f" move bound (cot(pi/2q) - sqrt q)/(4cos(pi/(q+1))) = {(1 / math.tan(math.pi / (2 * q)) - math.sqrt(q)) / (4 * math.cos(math.pi / (q + 1))):.3f}, q/(2pi) = {q / (2 * math.pi):.3f}")
    # Lipschitz constants on random tournaments with random Hamiltonian paths (also the dipole law for THM-4811 (A))
    rng = np.random.default_rng(4811)
    bad = {"K": 0, "d": 0, "cut": 0, "c3": 0, "op": 0, "nuc": 0, "pathspec": 0}
    trials = 0
    worst_op = 0.0
    worst_nuc = 0.0
    for n in list(range(4, 13)) + [16, 20, 25, 31, 40]:
        for _ in range(150):
            A = np.triu((rng.random((n, n)) < rng.random()).astype(int), 1)
            K = A - A.T
            K = K + 0
            # tournament: K_ij = +1 if i->j; fill: pairs with A=0 above diagonal mean j->i
            Kt = np.zeros((n, n), dtype=int)
            for i in range(n):
                for j in range(i + 1, n):
                    if A[i, j]:
                        Kt[i, j], Kt[j, i] = 1, -1
                    else:
                        Kt[i, j], Kt[j, i] = -1, 1
            K = Kt
            # Hamiltonian path by insertion (random order)
            order = list(rng.permutation(n))
            path = [order[0]]
            for v in order[1:]:
                # insert v: before first u with v -> u
                pos = len(path)
                for idx, u in enumerate(path):
                    if K[v, u] == 1:
                        pos = idx
                        break
                path.insert(pos, v)
            assert all(K[path[i], path[i + 1]] == 1 for i in range(n - 1))
            # random rotation is not needed; compute move
            Kp = -K.copy()
            for i in range(n - 1):
                u, v = path[i], path[i + 1]
                Kp[u, v], Kp[v, u] = 1, -1
            AQ = np.zeros((n, n), dtype=int)
            for i in range(n - 1):
                u, v = path[i], path[i + 1]
                AQ[u, v], AQ[v, u] = 1, -1
            trials += 1
            if not np.array_equal(Kp, -K + 2 * AQ):
                bad["K"] += 1
            d = K.sum(axis=1) / 2
            dp = Kp.sum(axis=1) / 2
            e = np.zeros(n)
            e[path[0]] += 1
            e[path[-1]] -= 1
            if not np.allclose(dp, -d + e):
                bad["d"] += 1
            for _r in range(30):
                X = rng.random(n) < 0.5
                if X.all() or (~X).all():
                    continue
                f = lambda M, Y: int(((M == 1) & np.outer(Y, ~Y)).sum())
                if f(Kp, X) != f(K, ~X) + int(X[path[0]]) - int(X[path[-1]]):
                    bad["cut"] += 1
            s = (K == 1).sum(axis=1)
            sp = (Kp == 1).sum(axis=1)
            c3 = lambda sc: math.comb(n, 3) - sum(math.comb(int(x), 2) for x in sc)
            if c3(sp) - c3(s) != s[path[0]] - s[path[-1]] - 1:
                bad["c3"] += 1
            sv = np.linalg.svd(K.astype(float), compute_uv=False)
            svp = np.linalg.svd(Kp.astype(float), compute_uv=False)
            Lop = 4 * math.cos(math.pi / (n + 1))
            Lnuc = 2 * sum(abs(2 * math.cos(k * math.pi / (n + 1))) for k in range(1, n + 1))
            dop = abs(sv[0] - svp[0])
            dnuc = abs(sv.sum() - svp.sum())
            worst_op = max(worst_op, dop / Lop)
            worst_nuc = max(worst_nuc, dnuc / Lnuc)
            if dop > Lop + 1e-9:
                bad["op"] += 1
            if dnuc > Lnuc + 1e-9:
                bad["nuc"] += 1
            ev = np.sort(np.abs(np.linalg.eigvals(AQ.astype(float))))
            ex = np.sort(np.abs([2 * math.cos(k * math.pi / (n + 1)) for k in range(1, n + 1)]))
            if not np.allclose(ev, ex, atol=1e-7):
                bad["pathspec"] += 1
    say(f"THM-4811 (A) + THM-4813 (C) on {trials} random (tournament, Hamiltonian path) pairs, n = 4..12, 16, 20, 25, 31, 40: failures {bad};"
        f" max observed |d op|/4cos(pi/(n+1)) = {worst_op:.3f}, max |d nuc|/(2 sum|2cos|) = {worst_nuc:.3f}")


if __name__ == "__main__":
    check_A1()
    check_A2()
    check_A3()
    check_window()
    check_A4()
    check_B()
    check_C()
