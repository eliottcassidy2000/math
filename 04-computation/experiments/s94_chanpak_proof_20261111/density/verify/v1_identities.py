# Verifier (lane density, S94). Independent exact checks of P1-P4 of the lane's PROOFS.txt.
# Independence: eta is computed from the g-recursion of THM-4890 (g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v mod 2)),
# eta(h) = g(2h), q(h) = (2h+1-eta(h))/2; the lane used the q-recursion.  General (A,B), not only F1 rows.
import sys, random
from math import gcd
from functools import lru_cache
sys.setrecursionlimit(10000)
@lru_cache(maxsize=None)
def g(v):
    if v == 0: return 1
    if v & 1: return g(v >> 1)
    u = v >> 1
    return g(u) + u + (u & 1)
def eta(h): return g(2 * h)
def q(h):
    e = 2 * h + 1 - eta(h)
    assert e % 2 == 0
    return e // 2
def qrec(h):
    if h == 0: return 0
    if h & 1: return qrec(h >> 1) + ((h >> 1) + 1) // 2
    return qrec(h >> 1)
def rev(x, w):
    r = 0
    for i in range(w):
        if (x >> i) & 1: r |= 1 << (w - 1 - i)
    return r
def wt(d): return 1 if d == 1 else 1 << (d - 2)
ok = True
# 0. q from g agrees with the q-recursion; bounds 0 <= q <= (h-1)/2, h+2 <= eta <= 2h+1
bad = 0
for h in range(1, 1 << 14):
    if q(h) != qrec(h): bad += 1
    if not (0 <= 2 * q(h) <= h - 1): bad += 1
    if not (h + 2 <= eta(h) <= 2 * h + 1): bad += 1
print("0. q via g-recursion = q-recursion, size bounds, h < 2^14: violations", bad); ok &= bad == 0
# P1. readout form on random large arguments
rnd = random.Random(20261111)
bad = 0; cnt = 0
for _ in range(20000):
    h1 = rnd.randrange(0, 1 << rnd.randrange(1, 40)); h2 = rnd.randrange(0, 1 << rnd.randrange(1, 40)); s = rnd.randrange(1, 30)
    x1, y1, x2, y2 = eta(h1), 2 * h1 + 1, eta(h2), 2 * h2 + 1
    A, B = x1 + (1 << s) * y1, x1
    if x1 * x2 + (1 << s) * y1 * y2 != A * (2 * h2 + 1) - 2 * B * q(h2): bad += 1
    if gcd(A, B) != gcd(x1, y1): bad += 1
    cnt += 1
print("P1. readout identity and gcd(A,B)=gcd(x1,y1), random h1,h2 < 2^40, s < 30: instances", cnt, "violations", bad); ok &= bad == 0
# P2. general (A,B): divisibility, localisation, injectivity threshold; sharpness of the threshold
def collisions(A, B, n):
    d = {}
    for h in range(1 << (n - 1), 1 << n):
        d.setdefault(A * (2 * h + 1) - 2 * B * q(h), []).append(h)
    return [L for L in d.values() if len(L) > 1]
bad = 0; npairs = 0; tested = 0; inj_fail = 0; thr_cases = 0
for n in range(4, 13):
    for _ in range(250):
        A = rnd.randrange(1, 1 << (n + 1)); B = rnd.randrange(1, 1 << (n + 1))
        d = gcd(A, B); tested += 1
        cl = collisions(A, B, n)
        if A >= (1 << (n - 1)) * d:
            thr_cases += 1
            if cl: inj_fail += 1
        for L in cl:
            for h in L:
                for hp in L:
                    if h == hp: continue
                    npairs += 1
                    if (h - hp) % (B // d): bad += 1
                    if (q(h) - q(hp)) % (A // d): bad += 1
                    if q(h) == q(hp): bad += 1
                    if A // d > (1 << (n - 1)) - 1: bad += 1
                    if abs(h - hp) * A > B * ((1 << (n - 1)) - 1): bad += 1
print("P2. random (A,B) < 2^(n+1), n = 4..12: readouts %d, colliding ordered pairs %d, violations %d; readouts with A >= 2^(n-1) gcd: %d, of which non-injective: %d"
      % (tested, npairs, bad, thr_cases, inj_fail)); ok &= (bad == 0 and inj_fail == 0)
# sharpness: largest A/d (d = gcd) for which a collision exists on n-bit words, B ranging freely
for n in range(4, 11):
    best = 0; arg = None
    for A in range(1, 1 << (n - 1)):
        for B in range(1, 4 * A + 2):
            if gcd(A, B) != 1: continue
            if A <= best: continue
            if collisions(A, B, n):
                best = A; arg = (A, B)
    print("   sharpness n=%d: largest coprime A with a collision (B <= 4A+1): %s  (lemma: A <= 2^(n-1)-1 = %d)" % (n, arg, (1 << (n - 1)) - 1))
# P2 corollary on F1 rows: n2 <= sigma + n1 + 1 injective if gcd(x1,y1)=1; what happens at n2 = sigma + n1 + 2
bad = 0; rows = 0; beyond = 0; beyond_coll = 0
for n1 in range(1, 7):
    for h1 in range(1 << (n1 - 1), 1 << n1):
        x1, y1 = eta(h1), 2 * h1 + 1
        if gcd(x1, y1) != 1: continue
        for s in (1, 2, 3):
            A, B = x1 + (1 << s) * y1, x1
            for n2 in range(1, s + n1 + 2):
                rows += 1
                if collisions(A, B, n2): bad += 1
            if s + n1 + 2 <= 13:
                beyond += 1
                if collisions(A, B, s + n1 + 2): beyond_coll += 1
print("P2 corollary: coprime rows h1 (n1 <= 6), sigma = 1..3, all n2 <= sigma+n1+1: rows %d, non-injective %d; at n2 = sigma+n1+2: %d rows, %d non-injective"
      % (rows, bad, beyond, beyond_coll)); ok &= bad == 0
# h1 = 0 row: injective for n2 <= sigma + 1
bad = 0
for s in range(1, 10):
    for n2 in range(1, s + 2):
        if collisions(1 + (1 << s), 1, n2): bad += 1
print("P2 corollary, h1 = 0 (A = 2^sigma + 1, B = 1), n2 <= sigma + 1, sigma <= 9: non-injective", bad); ok &= bad == 0
# P3. affine fibre lemma on random pairs of wide words, and (iii) exhaustively for w <= 8
bad = 0; cnt = 0
for _ in range(4000):
    w = rnd.randrange(2, 48)
    x = rnd.randrange(0, 1 << w); xp = rnd.randrange(0, 1 << w)
    eps = [((x >> i) & 1) - ((xp >> i) & 1) for i in range(w)]
    xplus = sum(1 << i for i in range(w) if eps[i] == 1); xminus = sum(1 << i for i in range(w) if eps[i] == -1)
    lin = 0
    for i in range(w):
        if eps[i] == 0 and (x >> i) & 1:
            lin += sum(eps[j] * wt(abs(i - j)) for j in range(w) if j != i and eps[j])
    if q(x) - q(xp) != q(xplus) - q(xminus) + lin: bad += 1
    cnt += 1
print("P3. affine fibre lemma, random pairs of width < 48: instances", cnt, "violations", bad); ok &= bad == 0
bad = 0; cnt = 0
import itertools
for w in range(1, 9):
    for eps in itertools.product((-1, 0, 1), repeat=w):
        X = sum(eps[i] << i for i in range(w)); Xi = sum(eps[i] << (w - 1 - i) for i in range(w))
        for i in range(w):
            if eps[i]: continue
            L = sum(eps[j] * wt(abs(i - j)) for j in range(w) if j != i)
            Xlt = sum(eps[j] << j for j in range(i)); Xigt = sum(eps[j] << (w - 1 - j) for j in range(i + 1, w))
            e1 = eps[i + 1] if i + 1 < w else 0; e0 = eps[i - 1] if i >= 1 else 0
            # 4 L_i = (X - X_<i)/2^i + (Xi - Xi_>i)/2^(w-1-i) + 2(e1 + e0); the divisions are exact
            if (X - Xlt) % (1 << i) or (Xi - Xigt) % (1 << (w - 1 - i)): bad += 1
            if 4 * L != (X - Xlt) // (1 << i) + (Xi - Xigt) // (1 << (w - 1 - i)) + 2 * (e1 + e0): bad += 1
            cnt += 1
print("P3 (iii). formula for 4 L_i, all signed strings of width <= 8: instances", cnt, "violations", bad); ok &= bad == 0
# P4. general A >= B >= 1
def maxrun(h, n):
    b = [(h >> i) & 1 for i in range(n)]
    m = 1; r = 1
    for i in range(1, n):
        r = r + 1 if b[i] == b[i - 1] else 1
        m = max(m, r)
    return m
from math import log2
bad_eq = bad_a = bad_b = bad_c = bad_d = 0; npairs = 0; rd = 0; maxdev = 0.0
for n in range(6, 13):
    for _ in range(60):
        B = rnd.randrange(1, 1 << rnd.randrange(1, n + 1)); A = B + rnd.randrange(0, 1 << rnd.randrange(1, n + 2))
        dg = gcd(A, B); rd += 1
        d = {}
        for h in range(1 << (n - 1), 1 << n):
            d.setdefault(A * (2 * h + 1) - 2 * B * q(h), []).append(h)
        sols = {}
        for N, L in d.items():
            # (c) prefix classes
            if len(L) > 1:
                for t in range(1, n + 1):
                    k = len(set(h >> (n - t) for h in L))
                    if k * A > 2 * A + (1 << (t - 1)) * B: bad_c += 1
            for h in L:
                for hp in L:
                    if h == hp: continue
                    npairs += 1
                    dd = h ^ hp; top = dd.bit_length() - 1; s = (dd & -dd).bit_length() - 1
                    w = top - s + 1; tp = n - 1 - top
                    P = h >> (top + 1); Z = h & ((1 << s) - 1)
                    x = (h >> s) & ((1 << w) - 1); xp = (hp >> s) & ((1 << w) - 1)
                    c = (1 << (s + 1)) - rev(Z, s)
                    X = x - xp; Xi = rev(x, w) - rev(xp, w)
                    gam = 2 * (q(x) - q(xp)) + (P & 1) * (((x >> (w - 1)) & 1) - ((xp >> (w - 1)) & 1)) + (c & 1) * ((x & 1) - (xp & 1))
                    if (B * c + (A - B) * (1 << (s + 1))) * X - B * P * Xi != B * gam: bad_eq += 1
                    if abs(gam) > (1 << w): bad_eq += 1
                    if not (X & 1 and Xi & 1): bad_eq += 1
                    if B & 1 and X % (B // dg): bad_b += 1
                    key = (tp, s, x, xp, P & 1, c & 1)
                    sols.setdefault(key, set()).add((c, P))
                    R = max(maxrun(h, n), maxrun(hp, n))
                    dev = abs(tp - s - log2(A / B))
                    if dev > R + 4 + 1e-9: bad_d += 1
                    maxdev = max(maxdev, dev - R)
        for (tp, s, x, xp, p1, c1), S in sols.items():
            w = n - tp - s
            X = x - xp; Xi = rev(x, w) - rev(xp, w); gg = gcd(X, Xi)
            # bound: 1 + g * min(2^(t'-1)/|X|, 2^s/|Xi|)
            lim = 1 + min(gg * (1 << (tp - 1)) / abs(X), gg * (1 << s) / abs(Xi))
            if len(S) > lim + 1e-9: bad_a += 1
print("P4. random A >= B >= 1, n = 6..12: readouts %d, colliding ordered pairs %d; violations: equation/|gamma|/oddness %d, (a) %d, (b) %d, (c) %d, (d) %d; max of |t'-s-log2(A/B)| - R = %.3f"
      % (rd, npairs, bad_eq, bad_a, bad_b, bad_c, bad_d, maxdev)); ok &= (bad_eq + bad_a + bad_b + bad_c + bad_d == 0)
print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")
