#!/usr/bin/env python3
"""Verification of the corrected S92 note (README.md): small exact statements, recomputed with independent code.

Every block prints what it checks and the README passage it belongs to.  Exit status 1 on any failure.
No file of the session or of the auditors is imported.
"""
import itertools
import math
import sys
from fractions import Fraction

import numpy as np

OK = True


def check(cond, msg):
    global OK
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        OK = False


# ----------------------------------------------------------------------------------------------------------------
# states
# ----------------------------------------------------------------------------------------------------------------
def q_of(h):
    """q(2h) = q(h), q(2h + 1) = q(h) + ceil(h/2)."""
    if h == 0:
        return 0
    k = h >> 1
    return q_of(k) + ((k + 1) >> 1 if h & 1 else 0)


def psi(h):
    return h - q_of(h)


def odd_part(n):
    n = abs(n)
    while n % 2 == 0:
        n //= 2
    return n


def chi4(m):
    return 1 if m % 4 == 1 else -1


# ----------------------------------------------------------------------------------------------------------------
print("A. one-apex tournaments and joins (section 5): H(R(h)) = 2 psi(h) + 1, H(J(h1, h2, s)) = x1 x2 + 2^s y1 y2 - [s = 0]")


def one_apex(h):
    """transitive c_1 -> ... -> c_k, k = bit length of 2h; apex u -> c_i iff the i-th digit of 2h (msb first) is 1."""
    u = 2 * h
    k = u.bit_length()
    n = k + 1
    adj = [[0] * n for _ in range(n)]
    for i in range(k):
        for j in range(i + 1, k):
            adj[i][j] = 1
        if (u >> (k - 1 - i)) & 1:
            adj[k][i] = 1
        else:
            adj[i][k] = 1
    return adj, n


def ham(adj, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if c:
                av = adj[v]
                for w in range(n):
                    if av[w] and not (S >> w) & 1:
                        dp[S | (1 << w)][w] += c
    return sum(dp[(1 << n) - 1])


def join(h1, h2, s):
    a1, n1 = one_apex(h1)
    a2, n2 = one_apex(h2)
    n = n1 + s + n2
    adj = [[0] * n for _ in range(n)]
    for i in range(n1):
        for j in range(n1):
            adj[i][j] = a1[i][j]
    for i in range(n2):
        for j in range(n2):
            adj[n1 + s + i][n1 + s + j] = a2[j][i]
    blocks = [range(0, n1)] + [range(n1 + t, n1 + t + 1) for t in range(s)] + [range(n1 + s, n)]
    for bi in range(len(blocks)):
        for bj in range(bi + 1, len(blocks)):
            for i in blocks[bi]:
                for j in blocks[bj]:
                    adj[i][j] = 1
    ap1, ap2 = n1 - 1, n - 1
    adj[ap1][ap2] = 0
    adj[ap2][ap1] = 1
    return adj, n


good = all(ham(*one_apex(h)) == 2 * psi(h) + 1 for h in range(64))
check(good, "H(R(h)) = 2 psi(h) + 1 for h < 64")
nj = 0
good = True
for h1 in range(8):
    for h2 in range(8):
        for s in range(4):
            adj, n = join(h1, h2, s)
            if n > 11:
                continue
            x1, y1, x2, y2 = 2 * psi(h1) + 1, 2 * h1 + 1, 2 * psi(h2) + 1, 2 * h2 + 1
            good &= ham(adj, n) == x1 * x2 + (y1 * y2 << s) - (1 if s == 0 else 0)
            nj += 1
check(good, f"H(J(h1, h2, s)) = x1 x2 + 2^s y1 y2 - [s = 0] on {nj} joins with at most 11 vertices (Hamiltonian paths by DP)")

# digit bound, the family of 10,400 joins
fam = [(h1, h2, s) for h1 in range(64) for h2 in range(h1, 64) for s in range(5)]
check(len(fam) == 10400, f"the family h1 <= h2 < 64, s <= 4 has {len(fam)} joins (README: 10,400)")


def alphas(h1, h2, s):
    c3 = ((2 * h1 + 1) * (2 * h2 + 1) << s) - (1 if s == 0 else 0)
    return psi(h1) + psi(h2) + c3 // 2, psi(h1) * psi(h2)


shares = {}
for b in (2, 3, 10, 100):
    cntb = sum(1 for t in fam if max(alphas(*t)) <= b - 1)
    shares[b] = (cntb, cntb / len(fam))
print("    joins with digits at most b - 1:", {b: (c, round(sh, 4)) for b, (c, sh) in shares.items()})
check([round(shares[b][1], 4) for b in (2, 3, 10)] == [0.0002, 0.0004, 0.0014] and round(shares[100][1], 3) == 0.022,
      "shares 0.0002, 0.0004, 0.0014, 0.022 for b = 2, 3, 10, 100")
vecs17 = sorted({alphas(*t) for t in fam if 1 + 2 * alphas(*t)[0] + 4 * alphas(*t)[1] == 17})
check(vecs17 == [(6, 1), (8, 0)], f"17 is the count of joins with (alpha_1, alpha_2) in {vecs17}")
good = True
for t in fam[:3000]:
    a1, a2 = alphas(*t)
    for b in (max(a1, a2) + 1, max(a1, a2) + 7):
        if b < 2:
            continue
        Nb = (1 + b * psi(t[0])) * (1 + b * psi(t[1])) + b * ((((2 * t[0] + 1) * (2 * t[1] + 1)) << t[2]) - (t[2] == 0)) // 2
        good &= (Nb % b, (Nb // b) % b, Nb // (b * b)) == (1, a1, a2) and Nb == 1 + a1 * b + a2 * b * b
check(good, "for b > max(alpha_1, alpha_2) the base-b numeral of N_b = x_b(h1) x_b(h2) + b (2^s y1 y2 - [s = 0])/2 is (alpha_2)(alpha_1)(1)")
# balanced states
size = 1 << 20
qs = np.zeros(size, np.int64)
for k in range(20):
    lo, hi = 1 << k, 1 << (k + 1)
    idx = np.arange(lo, hi, dtype=np.int64)
    half = idx >> 1
    qs[lo:hi] = qs[half] + np.where(idx & 1, (half + 1) >> 1, 0)
hs = np.arange(size, dtype=np.int64)
ps = hs - qs
check(all(int(ps[h]) == psi(h) for h in (0, 1, 2, 3, 5, 6, 45, 117, 309, 9045, 62805, 777777)), "vectorised psi = recursive psi")
bal2 = hs[2 * ps == 2 * hs].tolist()
check(all(v & (v - 1) == 0 for v in bal2), "b = 2: balanced states x = y are h = 0 and the powers of 2 (below 2^20)")
bal3 = hs[3 * ps == 2 * hs].tolist()
check(bal3[:8] == [0, 3, 45, 117, 309, 9045, 62805, 111957], f"b = 3: balanced states below 2^20: {bal3}")
check(all(hs[b * ps == 2 * hs].tolist() == [0] for b in (4, 5, 6, 10)) and bool(np.all(2 * ps[1:] > hs[1:])),
      "b >= 4: only h = 0 (psi(h) > h/2 for 1 <= h < 2^20)")
check((psi(3), 2 * 3 + 1, 1 + 3 * psi(3), 1 + 10 * psi(3)) == (2, 7, 7, 21), "h = 3: psi = 2, y = 7, x_3 = 7, x_10 = 21")
check((2 * 6 + 1, 1 + 5 * psi(6)) == (13, 26), "h = 6: y = 13, x_5 = 26 = 2y")

# ----------------------------------------------------------------------------------------------------------------
print("B. the signed coordinate (section 6)")
R = 5000
lab = {n: abs(4 * n + 1) for n in range(-R, R + 1)}
check(sorted(lab.values()) == list(range(1, 4 * R + 2, 2)),
      "n -> |4n + 1| maps [-5000, 5000] bijectively onto the odd numbers 1..20001")
check(all((chi4(m) * m - 1) % 4 == 0 and abs(4 * ((chi4(m) * m - 1) // 4) + 1) == m for m in range(1, 20001, 2)),
      "every odd m >= 1 is |4n + 1| for n = (chi_4(m) m - 1)/4 (m < 20000)")
good = True
for n in range(-R, R + 1):
    ms = 4 * n + 1
    if n % 2:
        good &= 4 * (3 * n + 1) + 1 == 3 * ms + 2
    else:
        good &= 4 * (n // 2) + 1 == (ms + 1) // 2 and (ms + 1) % 2 == 0
check(good, "odd step m* -> 3 m* + 2, even step m* -> (m* + 1)/2 for m* = 4n + 1")
check(all((4 * n + 1) * (4 * k + 1) == 4 * (4 * n * k + n + k) + 1 for n in range(-40, 41) for k in range(-40, 41)),
      "product of signed labels: n * k = 4nk + n + k")
check(all(3 * abs(4 * n + 1) == abs(4 * (-3 * n - 1) + 1) for n in range(-R, R + 1)), "tripling a count is n -> -3n - 1")
check(all((n >= 0) == (chi4(abs(4 * n + 1)) == 1) for n in range(-R, R + 1)), "the Collatz side (sign of n, 0 positive) is chi_4 of the count")
# Session266 coordinates: h = 4x + s, n = s x
check(all(abs(4 * (s * x) + 1) == 4 * x + s for x in range(1, 500) for s in (1, -1)), "with h = 4x + eps (Session266): n = eps x gives |4n + 1| = h")
check(all((4 * x + e) * 3 == 4 * (3 * x + e) - e for x in range(1, 500) for e in (1, -1)), "tripling: J(x, eps) = (3x + eps, -eps)")
check(all(4 * (3 * x + e) + e == 3 * (4 * x + e) + 2 * e for x in range(1, 500, 2) for e in (1, -1))
      and all(4 * (x // 2) + e == ((4 * x + e) + e) // 2 for x in range(2, 500, 2) for e in (1, -1)),
      "central-modular packet: m -> 3m + 2 eps (x odd), m -> (m + eps)/2 (x even)")


def collatz(n):
    return 3 * n + 1 if n % 2 else n // 2


def cycle_from(n0):
    c = [n0]
    n = collatz(n0)
    while n != n0:
        c.append(n)
        n = collatz(n)
    return c


cycles = [cycle_from(0), cycle_from(1), cycle_from(-1), cycle_from(-5), cycle_from(-17)]
labels = [[abs(4 * n + 1) for n in c] for c in cycles]
print("    cycles:", [(c[0], len(c)) for c in cycles])
print("    labels:", labels[:4], "and", len(labels[4]), "labels of the cycle of -17")
flat = [m for L in labels for m in L]
check([len(c) for c in cycles] == [1, 3, 2, 5, 18] and len(flat) == 29 and len(set(flat)) == 29, "five known cycles: 29 arcs, 29 distinct labels")
check(sorted(labels[1]) == [5, 9, 17] and sorted(labels[2]) == [3, 7] and sorted(labels[3]) == [19, 27, 39, 55, 79],
      "labels {5, 17, 9}, {3, 7}, {19, 55, 27, 79, 39}")
check([m for m in flat if m in (7, 21)] == [7] and max(flat) == 1087,
      "among the 29 labels the only one of the two known non-values 7, 21 is 7 (largest label 1087 < 2^32: THM-4875 covers the others)")
check(abs(4 * (-2) + 1) == 7 and collatz(-2) == -1 and abs(4 * 5 + 1) == 21 and collatz(5) == 16, "7 = arc -2 -> -1, 21 = arc 5 -> 16")
check(all(collatz((4 ** k - 1) // 3) == 4 ** k and 4 * ((4 ** k - 1) // 3) + 1 == (4 ** (k + 1) - 1) // 3 for k in range(1, 12)),
      "trunk arcs (4^k - 1)/3 -> 4^k have labels (4^(k+1) - 1)/3 = 5, 21, 85, ...")
check([abs(4 * n + 1) for n in (0, -1, 1)] == [1, 3, 5] and [collatz(n) for n in (0, -1, 1)] == [0, -2, 4],
      "1, 3, 5 are the labels of the arcs 0 -> 0, -1 -> -2, 1 -> 4")
# favoured classes
good = True
for p in (3, 5, 7, 11, 13, 17, 19, 23):
    i4, i2 = pow(4, -1, p), pow(2, -1, p)
    for n in range(-600, 601):
        N = abs(4 * n + 1)
        good &= (N % p == 0) == ((n + i4) % p == 0)
        if n >= 0:
            good &= ((N + 1) % p == 0) == ((n + i2) % p == 0)
        else:
            good &= ((N + 1) % p == 0) == (n % p == 0)
check(good, "p | N iff n = -1/4; N = -1 (mod p) iff n = -1/2 (n >= 0) or n = 0 (n < 0)   (p <= 23, |n| <= 600)")
check(all(odd_part(abs(4 * n + 1) + 1) == (2 * n + 1 if n >= 0 else odd_part(n)) for n in range(-3000, 3001)),
      "argument of the series of sheet 0 (odd part of N + 1): 2n + 1 for n >= 0, the odd part of n for n < 0")
check(Fraction(-1, 2) == 3 * Fraction(-1, 2) + 1 and Fraction(1, 2) == 3 * Fraction(1, 2) - 1, "fixed points: -1/2 of 3n + 1, 1/2 of 3n - 1")

# ----------------------------------------------------------------------------------------------------------------
print("C. lines (section 6, last bullet)")


def Tm(n):
    return odd_part(3 * n - 1)


def Tp(n):
    return odd_part(3 * n + 1)


def cm(N):
    return odd_part(N - 1)


def cp(N):
    return odd_part(N + 1)


rng = range(3, 100001, 2)
check(all(cm(3 * N) == Tm(N) and cp(3 * N) == Tp(N) for N in rng), "c_-(3N) = T_-(N), c_+(3N) = T_+(N) (definitions), odd 3 <= N < 10^5")
check(all((cm(3 * N) == Tp(cm(N))) == (N % 4 == 3) for N in rng), "c_-(3N) = T_+(c_-(N)) iff N = 3 (mod 4)")
check(all((cp(3 * N) == Tm(cp(N))) == (N % 4 == 1) for N in rng), "c_+(3N) = T_-(c_+(N)) iff N = 1 (mod 4)")
check(all(3 * (2 * a + 1) - 1 == 2 * (3 * a + 1) and (2 * a + 1) * (2 * b + 1) == 2 * (a + b + 2 * a * b) + 1 for a in range(60) for b in range(60)),
      "star law a * b = a + b + 2ab on (H - 1)/2, 1 * a = 3a + 1 and 3N - 1 = 2 (3a + 1) for N = 2a + 1")

# ----------------------------------------------------------------------------------------------------------------
print("D. fan primes (section 6)")


def primes_below(n):
    s = bytearray([1]) * n
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n) if s[i]]


xs_all = 2 * ps + 1
ys_all = 2 * hs + 1
fan_first = []
table = {}
for p in primes_below(400)[1:]:
    o = 1
    while pow(2, o, p) != 1:
        o += 1
    div = np.nonzero((xs_all[1:] % p == 0) & (ys_all[1:] % p == 0))[0]
    first = int(div[0]) + 1 if len(div) else None
    fan = (1 << (o // 2 - 1)) if o % 2 == 0 else None
    table[p] = (o, fan, first)
    if fan is not None:
        if not (2 ** (o // 2) + 1) % p == 0:
            check(False, f"p = {p}: 2^(ord/2) + 1 not divisible")
        if fan < size and not (xs_all[fan] == ys_all[fan] and xs_all[fan] % p == 0):
            check(False, f"p = {p}: the fan is not a balanced p-divisible state")
    if first is None:
        check(False, f"p = {p}: no p-divisible state below 2^20")
    if fan is not None and first == fan:
        fan_first.append(p)
print("    (p, ord_p(2), fan, first p-divisible state) for p < 60:", [(p,) + table[p] for p in primes_below(60)[1:]])
check([p for p in fan_first if p < 200] == [3, 5, 11, 13, 17, 43, 113], f"p < 200 with the fan as first p-divisible state: {[p for p in fan_first if p < 200]}")
check([p for p in fan_first if p >= 200] == [241, 257, 331], f"200 < p < 400: {[p for p in fan_first if p >= 200]} (README: the next are 241, 257, 331)")
check(all((table[p][1] is not None) == (table[p][0] % 2 == 0) for p in table) and all(
    any(pow(2, t, p) == p - 1 for t in range(1, p)) == (table[p][0] % 2 == 0) for p in table),
      "p | 2^t + 1 for some t iff ord_p(2) is even (p < 400)")

# ----------------------------------------------------------------------------------------------------------------
print("E. tournaments of order at most 6 (sections 7, 8)")


def tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    return pairs


def analyse(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pidx = {pr: k for k, pr in enumerate(pairs)}
    perms = list(itertools.permutations(range(n)))
    npairs = len(pairs)
    full = (1 << npairs) - 1
    seen = {}
    classes = []
    for code in range(1 << npairs):
        if code in seen:
            continue
        adj = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (code >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        orbit = set()
        aut = 0
        for pm in perms:
            c2 = 0
            for k, (i, j) in enumerate(pairs):
                a, b = pm[i], pm[j]
                if adj[i][j]:
                    src, dst = a, b
                else:
                    src, dst = b, a
                if src < dst:
                    c2 |= 1 << pidx[(src, dst)]
            orbit.add(c2)
            if c2 == code:
                aut += 1
        idx = len(classes)
        for c2 in orbit:
            seen[c2] = idx
        H = ham(adj, n)
        c3 = sum(1 for a, b, c in itertools.combinations(range(n), 3)
                 if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[b][a] and adj[c][b] and adj[a][c]))
        sc = sorted(sum(r) for r in adj)
        strong = all(sum(sc[:k]) > k * (k - 1) // 2 for k in range(1, n))
        classes.append(dict(code=code, H=H, c3=c3, aut=aut, strong=strong, size=len(orbit), conv=(full ^ code), adj=adj))
    for c in classes:
        c["self_converse"] = seen[c["conv"]] == seen[c["code"]]
    return classes, seen


def arc_usage(adj, n):
    use = {}
    for pm in itertools.permutations(range(n)):
        if all(adj[pm[i]][pm[i + 1]] for i in range(n - 1)):
            for i in range(n - 1):
                use[(pm[i], pm[i + 1])] = use.get((pm[i], pm[i + 1]), 0) + 1
    return sorted(use.values(), reverse=True)


nclasses, nstrong, nsc, Hsets = {}, {}, {}, {}
cls = {}
for n in range(1, 7):
    cl, seen = analyse(n)
    cls[n] = cl
    nclasses[n] = len(cl)
    nstrong[n] = sum(c["strong"] for c in cl)
    nsc[n] = sum(c["self_converse"] for c in cl)
    Hsets[n] = sorted({c["H"] for c in cl})
print("    classes:", nclasses, " strong:", nstrong, " self-converse:", nsc)
check([nclasses[n] for n in range(1, 7)] == [1, 1, 2, 4, 12, 56], "classes 1, 1, 2, 4, 12, 56 at orders 1 to 6")
check([nstrong[n] for n in range(1, 7)] == [1, 0, 1, 1, 6, 35], "strong classes 1, 0, 1, 1, 6, 35")
check([nsc[n] for n in range(3, 7)] == [2, 2, 8, 12], "self-converse classes 2, 2, 8, 12 at orders 3 to 6 (HYP-8205); 8 of 12 at order 5")
# B = 1/(1 - S)
S = [0, 1, 0, 1, 1, 6, 35]


def inv_series(S, N):
    B = [1] + [0] * N
    for n in range(1, N + 1):
        B[n] = sum(S[k] * B[n - k] for k in range(1, min(n, len(S) - 1) + 1))
    return B


check(inv_series(S, 6) == [1, 1, 1, 2, 4, 12, 56], "B(x) = 1/(1 - S(x)) with S = x + x^3 + x^4 + 6x^5 + 35x^6 gives 1, 1, 2, 4, 12, 56")
check(inv_series(S[:6], 6)[6] == 21, "S truncated at 6x^5 gives 21 at order 6")
c4 = cls[4]
info4 = sorted((c["H"], c["c3"], c["self_converse"], c["strong"], c["size"]) for c in c4)
print("    order 4 (H, 3-cycles, self-converse, strong, labelled copies):", info4)
check(info4 == [(1, 0, True, False, 24), (3, 1, False, False, 8), (3, 1, False, False, 8), (5, 2, True, True, 24)],
      "order 4: H = 1 + 2 c_3; the two classes with H = 3 are a converse pair, the others self-converse; 48 labelled self-converse (H = 1, 5), 16 others (H = 3)")
check(all((c["self_converse"]) == (c["H"] % 4 == 1) == (c["c3"] % 2 == 0) == ((chi4(c["H"]) * c["H"] - 1) // 4 >= 0) for c in c4),
      "order 4: self-converse iff H = 1 (mod 4) iff n >= 0 iff the number of 3-cycles is even")
prof4 = sorted(arc_usage(c["adj"], 4) for c in c4)
check(prof4 == sorted([[1, 1, 1], [2, 2, 2, 1, 1, 1], [2, 2, 2, 1, 1, 1], [4, 3, 3, 3, 1, 1]]), f"order 4 arc-usage profiles: {prof4}")
s5 = sorted((c["H"], c["aut"], c["self_converse"]) for c in cls[5] if c["strong"])
print("    strong classes of order 5 (H, |Aut|, self-converse):", s5)
check([t[0] for t in s5] == [9, 9, 11, 13, 15, 15] and sorted(t[1] for t in s5) == [1, 1, 1, 1, 3, 5] and all(t[2] for t in s5),
      "the six strong classes of order 5: H = 9, 9, 11, 13, 15, 15, |Aut| = 1, 1, 1, 1, 3, 5, all self-converse")
check([(t[0], t[1]) for t in s5][-2:] == [(15, 3), (15, 5)], "the automorphism orders 3 and 5 belong to the two classes with H = 15")
ns5 = sorted((c["H"], c["self_converse"]) for c in cls[5] if not c["strong"])
check(ns5 == [(1, True), (3, False), (3, False), (3, True), (5, False), (5, False)], f"the six non-strong classes of order 5: {ns5} (two converse pairs)")
check(all(c["H"] % 4 == 1 for n in (4, 6) for c in cls[n] if c["self_converse"]), "even order 4, 6: self-converse implies H = 1 (mod 4)")
check(any(c["H"] % 4 == 1 and not c["self_converse"] for c in cls[6]), "order 6: the converse fails (a class with H = 1 (mod 4) that is not self-converse)")
ints = {n: sorted((chi4(H) * H - 1) // 4 for H in Hsets[n]) for n in (4, 5, 6)}
print("    H at order 5:", Hsets[5])
print("    H at order 6:", Hsets[6])
check(ints[4] == [-1, 0, 1], "order 4: integers n = -1, 0, 1")
check(ints[5] == [n for n in range(-4, 4) if n != -2], "order 5: [-4, 3] without -2")
check(ints[6] == [n for n in range(-11, 12) if n not in (-10, -9, -2, 5)], "order 6: [-11, 11] without -10, -9, -2, 5")
c3cls = [c for c in cls[3] if c["H"] == 3][0]
# Pi_2: ordered pairs of disjoint possibly empty paths covering the tournament
adj3 = c3cls["adj"]


def path_count(adj, verts):
    if not verts:
        return 1
    return sum(1 for pm in itertools.permutations(verts) if all(adj[pm[i]][pm[i + 1]] for i in range(len(pm) - 1)))


pi2 = sum(path_count(adj3, [v for v in range(3) if (mask >> v) & 1]) * path_count(adj3, [v for v in range(3) if not (mask >> v) & 1]) for mask in range(8))
check(pi2 == 12 and Fraction(pi2, 3) == 4 and Fraction(1 + 3, 3) == Fraction(4, 3), "Pi_2(C_3) = 12, Pi_2/H = 4; sigma_1(3)/3 = 4/3")
check(math.factorial(7) // 21 == 240 == 16 * 15, "7!/21 = 240 = 16 x 15")

# ----------------------------------------------------------------------------------------------------------------
print("F. the eight extensions of F_2^2 by F_2 (section 3 table)")
V = [(0, 0), (0, 1), (1, 0), (1, 1)]
nz = V[1:]
cocycles = []
for bits in itertools.product((0, 1), repeat=9):
    f = {(a, b): 0 for a in V for b in V}
    for k, (a, b) in enumerate([(a, b) for a in nz for b in nz]):
        f[(a, b)] = bits[k]
    add = lambda a, b: ((a[0] + b[0]) % 2, (a[1] + b[1]) % 2)
    if all((f[(b, c)] + f[(a, add(b, c))] + f[(add(a, b), c)] + f[(a, b)]) % 2 == 0 for a in V for b in V for c in V):
        cocycles.append(f)
check(len(cocycles) == 16, f"normalized 2-cocycles: {len(cocycles)}")
bygroup = {}
for f in cocycles:
    Sset = tuple(a for a in nz if f[(a, a)] == 1)
    comm = all(f[(a, b)] == f[(b, a)] for a in nz for b in nz)
    elems = [(a, z) for a in V for z in (0, 1)]
    mul = lambda u, v: (((u[0][0] + v[0][0]) % 2, (u[0][1] + v[0][1]) % 2), (u[1] + v[1] + f[(u[0], v[0])]) % 2)
    order4 = sum(1 for e in elems if mul(e, e) != ((0, 0), 0))
    grp = {(0, True): "C2^3", (4, True): "C4 x C2", (2, False): "D4", (6, False): "Q8"}[(order4, comm)]
    bygroup.setdefault((len(Sset), grp, comm), set()).add(Sset)
summary = sorted((k, len(v)) for k, v in bygroup.items())
print("    (#units of square -1, group, commutative) -> number of classes (subsets):", summary)
check(summary == [((0, "C2^3", True), 1), ((1, "D4", False), 3), ((2, "C4 x C2", True), 3), ((3, "Q8", False), 1)],
      "0: C_2^3 (1 class), 1: D_4 (3), 2: C_4 x C_2 (3), 3: Q_8 (1); commutative iff the number of units of square -1 is even")
check(sum(1 for f in cocycles if not all(f[(a, b)] == f[(b, a)] for a in nz for b in nz)) == 8, "8 cocycles with non-trivial commutator, 8 with trivial (four commutative classes)")
from collections import Counter
perS = Counter(tuple(a for a in nz if f[(a, a)] == 1) for f in cocycles)
check(len(perS) == 8 and set(perS.values()) == {2}, "every subset of the three units is the square set of exactly two cocycles (8 classes)")

# ----------------------------------------------------------------------------------------------------------------
print("G. norm counts on the orders (section 3, 'Norms, and on which orders')")
NMAX = 400


def r2(n):
    c = 0
    a = 0
    while a * a <= n:
        b2 = n - a * a
        b = math.isqrt(b2)
        if b * b == b2:
            c += (1 if a == 0 else 2) * (1 if b == 0 else 2)
        a += 1
    return c


R2 = [r2(n) for n in range(NMAX + 1)]


def sigma1(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)


add_conv = [sum(R2[a] * R2[n - a] for a in range(n + 1)) for n in range(NMAX + 1)]
mul_conv = [0] + [sum(R2[a] * R2[n // a] for a in range(1, n + 1) if n % a == 0) for n in range(1, NMAX + 1)]
check(all(add_conv[n] == 8 * sigma1(n) for n in range(1, NMAX + 1, 2)), "Lipschitz order: r_2 additively convolved with itself is 8 sigma_1(n) for odd n <= 400")
check(mul_conv[1:9] == [16, 32, 0, 48, 64, 0, 0, 64], f"Z[i] x Z[i]: multiplicative convolution at n = 1..8: {mul_conv[1:9]}")
# direct counts on Z[i] x Z[i] and on the lattice z1 = z2 (mod 2)
gauss = {}
B0 = 21
for a in range(-B0, B0 + 1):
    for b in range(-B0, B0 + 1):
        gauss.setdefault(a * a + b * b, []).append((a, b))
NDIR = 120
full_direct, sub_direct = [0] * (NDIR + 1), [0] * (NDIR + 1)
for n in range(1, NDIR + 1):
    for d in range(1, n + 1):
        if n % d:
            continue
        for (a, b) in gauss.get(d, []):
            for (c, e) in gauss.get(n // d, []):
                full_direct[n] += 1
                if (a - c) % 2 == 0 and (b - e) % 2 == 0:
                    sub_direct[n] += 1
check(full_direct[1:] == mul_conv[1:NDIR + 1], "direct count of |z1|^2 |z2|^2 = n on Z[i] x Z[i] equals the multiplicative convolution (n <= 120)")
check(sub_direct[1:9] == [8, 0, 0, 16, 32, 0, 0, 0], f"lattice z1 = z2 (mod 2), n = 1..8: {sub_direct[1:9]}")
check(all(2 * sub_direct[n] == full_direct[n] for n in range(1, NDIR + 1, 2)), "sign-table lattice: half of the count at every odd n <= 120")
check(all(sub_direct[n] == 0 for n in range(2, NDIR + 1, 4)), "sign-table lattice: 0 at n = 2 (mod 4)")
check(all(full_direct[n] == 0 and sub_direct[n] == 0 for n in range(3, NDIR + 1, 4)) and all(mul_conv[n] == 0 for n in range(3, NMAX + 1, 4)),
      "both counts vanish at n = 3 (mod 4)")
# the lattice z1 = z2 (mod 2) is the span of 1, i, j, ij: index 4, units 8 against 16
units_full = full_direct[1]
units_sub = sub_direct[1]
check((units_full, units_sub) == (16, 8), "16 units in Z[i] x Z[i], 8 in the span of the sign table")
# Hurwitz order: 24 sigma_1 at odd n
hur = [0] * 61
for a2 in range(-16, 17):
    for b2 in range(-16, 17):
        if (a2 - b2) % 2:
            continue
        for c2 in range(-16, 17):
            if (a2 - c2) % 2:
                continue
            for d2 in range(-16, 17):
                if (a2 - d2) % 2:
                    continue
                s4 = a2 * a2 + b2 * b2 + c2 * c2 + d2 * d2
                if s4 % 4 == 0 and 0 < s4 // 4 <= 60:
                    hur[s4 // 4] += 1
check(all(hur[n] == 24 * sigma1(n) for n in range(1, 60, 2)), "Hurwitz order: 24 sigma_1(n) elements of norm n for odd n < 60")


# norm form on the span of the units of M_2: x0^2 + x1^2 - x2^2 - x3^2 against the seam a1 a2 + b1 b2
def class_counts(fn, m):
    cnt = [0] * m
    for v in itertools.product(range(m), repeat=4):
        cnt[fn(*v) % m] += 1
    return cnt


seam = lambda a, b, c, d: a * b + c * d
split = lambda a, b, c, d: a * a + b * b - c * c - d * d
four = lambda a, b, c, d: a * a + b * b + c * c + d * d
eq = {m: class_counts(seam, m) == class_counts(split, m) for m in (2, 3, 4, 5, 7, 9)}
print("    x0^2 + x1^2 - x2^2 - x3^2 has the class counts of the seam:", eq)
check(eq == {2: False, 3: True, 4: False, 5: True, 7: True, 9: True}, "the span form has the seam law at the odd primes tested (3, 9, 5, 7) and not at 2, 4")
check(all(class_counts(seam, m) == class_counts(four, m) for m in (3, 5, 7, 9)) and class_counts(seam, 2) != class_counts(four, 2),
      "four squares = seam modulo 3, 5, 7, 9; different modulo 2 (8, 8 against 10, 6)")

# ----------------------------------------------------------------------------------------------------------------
print("H. local algebra (sections 1, 2)")


def Lm(m, e):
    tot = Fraction(0)
    g = math.gcd(m, e)
    for d in range(1, g + 1):
        if g % d == 0:
            t = Fraction(1, d)
            md = m // d
            pp = 2
            x = md
            while x > 1:
                if x % pp == 0:
                    t *= 1 - Fraction(1, pp * pp)
                    while x % pp == 0:
                        x //= pp
                pp += 1
            tot += t
    return tot


def table_counts(m):
    Nm = [0] * m
    for a in range(m):
        for b in range(m):
            Nm[a * b % m] += 1
    return Nm


good = True
goodN = True
for m in (9, 15, 21, 25, 27, 45, 49, 63, 75, 81, 105, 135):
    Nm = table_counts(m)
    goodN &= all(Nm[t] == sum(d * sum(1 for k in range(1, m // d + 1) if math.gcd(k, m // d) == 1) for d in range(1, m + 1) if math.gcd(t, m) % d == 0)
                 for t in range(m))
    sq = [sum(Nm[a] * Nm[(c - a) % m] for a in range(m)) for c in range(m)]
    good &= all(sq[c] == m ** 3 * Lm(m, c) for c in range(m))
    sq2 = [sum(Nm[a] * Nm[((c - a) * pow(2, -1, m)) % m] for a in range(m)) for c in range(m)]
    good &= sq2 == sq
check(goodN, "N_m(t) = sum_{d | gcd(t, m)} d phi(m/d) for twelve odd m <= 135")
check(good, "(N_m * N_m)(c) = m^3 L_m(c), also with one table multiplied by the unit 2 (twelve odd m <= 135)")
good = True
for p, k in ((3, 2), (3, 3), (3, 4), (5, 2), (5, 3), (7, 2), (11, 2)):
    m = p ** k
    for e in range(m):
        v = 0
        while v < k and e % p ** (v + 1) == 0:
            v += 1
        if v < k:
            good &= Lm(m, e) == (1 - Fraction(1, p * p)) * Fraction(sum(p ** i for i in range(v + 1)), p ** v)
        else:
            good &= Lm(m, e) == 1 + Fraction(1, p) - Fraction(1, p ** (k + 1))
check(good, "L_(p^k) = (1 - p^-2) sigma_1(p^v)/p^v at v < k, 1 + 1/p - p^(-k-1) at the class 0")
check(Lm(9, 0) == Fraction(35, 27) and abs(float(Lm(9, 0)) - 1.296) < 5e-4, "L_9(0) = 35/27 = 1.296")
check(all(Lm(a * b, e) == Lm(a, e) * Lm(b, e) for a, b in ((9, 5), (9, 25), (3, 49), (27, 7), (25, 49)) for e in range(a * b)), "L_m is multiplicative")
good = True
for m in (9, 27, 45, 75, 225, 315):
    cd = {d: Fraction(1, d * d) * math.prod(1 - Fraction(1, p * p) for p in primes_below(m + 1) if (m // d) % p == 0) for d in range(1, m + 1) if m % d == 0}
    good &= sum(cd.values()) == 1
    # share of residue vectors with gcd(a, b, m) = d
    for d in cd:
        cntv = sum(1 for a in range(m) for b in range(m) if math.gcd(math.gcd(a, b), m) == d) if m <= 75 else None
        if cntv is not None:
            good &= Fraction(cntv, m * m) == cd[d]
    good &= all(Lm(m, e) == sum(d * cd[d] for d in cd if math.gcd(m, e) % d == 0) for e in range(m))
check(good, "c_d = d^-2 prod_{p | m/d} (1 - p^-2) is the share of the residue vectors with gcd(a, b, m) = d, sum c_d = 1, L_m(e) = sum_{d | gcd(m, e)} d c_d")
# the series
good = True
for n in (1, 3, 5, 7, 9, 15, 21, 45, 105, 315, 1001):
    prod = 1.0
    for p in primes_below(200000)[1:]:
        v = 0
        x = n
        while x % p == 0:
            x //= p
            v += 1
        prod *= (1 - 1 / p ** 2) * sum(p ** i for i in range(v + 1)) / p ** v
    good &= abs(prod - (8 / math.pi ** 2) * sigma1(n) / n) < 2e-5
    good &= Fraction(sigma1(n), n) == sum(Fraction(1, d) for d in range(1, n + 1) if n % d == 0)
check(good, "prod over odd p of (1 - p^-2) sigma_1(p^v)/p^v = (8/pi^2) sigma_1(n)/n (primes < 200000), and sigma_1(n)/n = sum_{d | n} 1/d")
check(abs((6 / math.pi ** 2) / (3 / 4) - 8 / math.pi ** 2) < 1e-15 and abs((8 / math.pi ** 2) / (1 - 1 / 9) - 9 / math.pi ** 2) < 1e-15,
      "(6/pi^2)/(3/4) = 8/pi^2 (Ingham's coefficient without the prime 2); (8/pi^2)/(1 - 1/9) = 9/pi^2 (bases 3, 6)")
check(all(8 * sigma1(n) == round(math.pi ** 2 * n * (8 / math.pi ** 2) * sigma1(n) / n) for n in range(1, 400, 2)), "r_4(n) = 8 sigma_1(n) = pi^2 n * S for odd n")

# ----------------------------------------------------------------------------------------------------------------
print("I. E8 and the rank-4 law (section 4)")
C = [[0] * 8 for _ in range(8)]
for i in range(8):
    C[i][i] = 2
for i, j in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)]:
    C[i][j] = C[j][i] = -1
detC = round(np.linalg.det(np.array(C, dtype=float)))
check(detC == 1, "the E8 Cartan matrix has determinant 1")


def e8_counts(m):
    cnt = [0] * m
    for v in itertools.product(range(m), repeat=8):
        val = sum(v[i] * v[i] for i in range(8))
        for i, j in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)]:
            val -= v[i] * v[j]
        cnt[val % m] += 1
    return cnt


def seam4_counts(m):
    c2 = [0] * m
    for a in range(m):
        for b in range(m):
            c2[a * b % m] += 1
    cur = c2[:]
    for _ in range(3):
        nxt = [0] * m
        for x in range(m):
            if cur[x]:
                for y in range(m):
                    nxt[(x + y) % m] += cur[x] * c2[y]
        cur = nxt
    return cur


e2 = e8_counts(2)
check(e2 == [136, 120] and seam4_counts(2) == [136, 120], f"modulo 2: E8 form {e2}, rank-4 seam {seam4_counts(2)} (136 zeros among 256 vectors)")
e3 = e8_counts(3)
check(e3 == seam4_counts(3), "modulo 3: E8 form = rank-4 seam")


def lifted(N1, p, n, K):
    """N_k(c) = p^((n-1)(k-1)) (N_1(c mod p) - [p | c]) + [p^2 | c] p^n N_(k-2)(c/p^2)."""
    out = {0: [1], 1: list(N1)}
    for k in range(2, K + 1):
        m = p ** k
        row = []
        for c in range(m):
            val = p ** ((n - 1) * (k - 1)) * (N1[c % p] - (1 if c % p == 0 else 0))
            if c % (p * p) == 0:
                val += p ** n * out[k - 2][(c // (p * p)) % (p ** (k - 2))]
            row.append(val)
        out[k] = row
    return out


lift = lifted([136, 120], 2, 8, 8)
good = True
for k in range(1, 9):
    m = 2 ** k
    for c in range(m):
        v = 0
        while v < k and c % 2 ** (v + 1) == 0:
            v += 1
        r = 4
        if v < k:
            f = (1 - Fraction(1, 2 ** r)) * Fraction(sum(2 ** (3 * i) for i in range(v + 1)), 2 ** (3 * v))
        else:
            f = (1 - Fraction(1, 2 ** r)) * Fraction(sum(2 ** (3 * i) for i in range(k)), 2 ** (3 * (k - 1))) + Fraction(1, 2 ** (3 * k))
        good &= Fraction(lift[k][c], m ** 7) == f
check(good, "the lifting recursion started from 136, 120 gives the rank-4 law (1 - 2^-4) sigma_3(2^v)/2^(3v) (and the class 0) modulo 2^k, k <= 8")
check(lift[2] == seam4_counts(4) and lift[3] == seam4_counts(8), "the recursion reproduces the enumerated rank-4 seam counts modulo 4 and 8")
check(lift[2] == [17536, 15360, 17280, 15360], "modulo 4: 17536, 15360, 17280, 15360 (the E8 counts printed in bridges_exact.out S3)")
lift3 = lifted(e3, 3, 8, 3)
check(lift3[2] == seam4_counts(9), "the same recursion at p = 3 from the counts modulo 3 gives the seam counts modulo 9")
zeta = lambda r: sum(1 / k ** r for k in range(1, 200000)) + 1 / ((r - 1) * 200000 ** (r - 1))
bern = {2: Fraction(1, 6), 4: Fraction(1, 30)}
for r, want in ((2, 24), (4, 240)):
    val = (2 * math.pi) ** r / (math.factorial(r - 1) * zeta(r))
    check(abs(val - want) < 1e-3 and Fraction(2 * r) / bern[r] == want, f"r = {r}: (2 pi)^r/((r-1)! zeta(r)) = {val:.4f} = 2r/|B_r| = {want}")
check(abs((8 * math.pi ** 4 / 3) / (math.pi ** 4 / 90) - 240) < 1e-9 and abs((2 * math.pi) ** 4 / 6 - 8 * math.pi ** 4 / 3) < 1e-9, "(8 pi^4/3)/zeta(4) = 240")
check(1 + Fraction(7 * 2, 2 * 27) == Fraction(34, 27) and 1 + Fraction(1, 27) - Fraction(1, 81) == Fraction(83, 81)
      and abs(float(1 + Fraction(1, 125) - Fraction(1, 625)) - 1.0064) < 1e-12,
      "THM-4916 factor at 3 is 34/27; uniform rank-4 value 83/81 at 3 and 1.0064 at 5")
check(all(Fraction(seam4_counts(p)[0], p ** 7) == 1 + Fraction(1, p ** 3) - Fraction(1, p ** 4) for p in (3, 5, 7)), "rank-4 seam, class 0 modulo p: factor 1 + p^-3 - p^-4 (p = 3, 5, 7)")

# ----------------------------------------------------------------------------------------------------------------
print("J. Ingham's sum (section 2, item 2)")
NI = 20400
dcount = np.zeros(NI + 1, np.int64)
for d in range(1, NI + 1):
    dcount[d::d] += 1
sig = np.zeros(NI + 1, np.int64)
for d in range(1, NI + 1):
    sig[d::d] += d
Dn, sn, ratio = [], [], []
for n in range(20000, 20400):
    D = int(np.dot(dcount[1:n], dcount[n - 1:0:-1]))
    Dn.append(D / n)
    sn.append(sig[n] / n)
    ratio.append(D / ((6 / math.pi ** 2) * sig[n] * math.log(n) ** 2))
corr = float(np.corrcoef(Dn, sn)[0, 1])
print(f"    correlation of D(n)/n with sigma_1(n)/n on [20000, 20400): {corr:.4f}; D(n)/main term: mean {np.mean(ratio):.3f}")
check(round(corr, 3) == 0.979, "correlation 0.979")
# D(n) is the number of positive solutions of a1 a2 + b1 b2 = n
good = True
for n in range(2, 41):
    direct = sum(1 for a1 in range(1, n) for a2 in range(1, n) if a1 * a2 < n for b1 in range(1, n) if (n - a1 * a2) % b1 == 0)
    good &= direct == int(np.dot(dcount[1:n], dcount[n - 1:0:-1]))
check(good, "D(n) = sum_{k < n} d(k) d(n - k) is the number of solutions of a1 a2 + b1 b2 = n in positive integers (n <= 40)")

# ----------------------------------------------------------------------------------------------------------------
print("K. unit equations (section 7, last paragraph)")
sols = [d for d in range(1, 60) if 2 ** d - 3 ** (d - 1) == 1]
check(sols == [1, 2], f"2^delta - 3^(delta-1) = 1 exactly for delta in {sols} (delta < 60): 2 - 1, 4 - 3")
check((2 ** 2 - 3 ** 1, 2 ** 1 - 3 ** 1, 2 ** 3 - 3 ** 2, 2 ** 11 - 3 ** 7) == (1, -1, -1, -139), "D = 2^L - 3^r = 1, -1, -1, -139 for (L, r) = (2, 1), (1, 1), (3, 2), (11, 7)")


def shortcut_cycle(n0):
    L = r = 0
    n = n0
    while True:
        if n % 2:
            n = (3 * n + 1) // 2
            r += 1
        else:
            n //= 2
        L += 1
        if n == n0:
            return L, r


check([shortcut_cycle(n) for n in (1, -1, -5, -17)] == [(2, 1), (1, 1), (3, 2), (11, 7)], "(L, r) of the shortcut cycles of 1, -1, -5, -17")
check({(4, 3), (2, 3), (8, 9)} & {(2, 1), (4, 3)} == {(4, 3)}, "only 4 - 3 is common to the cycle units {4 - 3, 2 - 3, 8 - 9} and the criterion {2 - 1, 4 - 3}")

print("ALL CHECKS PASSED" if OK else "SOME CHECKS FAILED")
sys.exit(0 if OK else 1)
