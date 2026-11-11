#!/usr/bin/env python3
"""Verifier of THM-5134 (after the audit corrections): exact checks of the rewritten statements and proofs.

Own code, written from the text of the theorem file (not from the session or audit scripts).
Sections:
  E1  recursions of psi, q; q(h) < h/2; x, y
  E2  the tournaments R(h), J(h1, h2, s) as defined in (G): Hamiltonian paths = N, packing polynomial = formula
  E3  THM-5070 (2) dictionary: a(W), trailing ones, the bracket for every split s1 + s2 = s
  E4  THM-5085 (2) dictionary and the Session206 dictionary of (A)
  E5  c_d, the fibre count of (A), L_m: definition, both prime-power cases, multiplicativity, mean 1,
      m^3 L_m(c) = #{a1 a2 + u b1 b2 = c} for every m (even included), four squares for odd m
  E6  base b: N_b = a x2 + b' y2 + k; the value V; content; affine bijection; the two congruences of alpha_1
  E7  balanced states; early boosts
Exit status 1 on failure.
"""
import itertools
import math
import sys
from fractions import Fraction

FAIL = []


def check(name, ok, extra=""):
    print(("  ok   " if ok else "  FAIL ") + name + ((" : " + extra) if extra else ""))
    if not ok:
        FAIL.append(name)


def psi(h):
    if h == 0:
        return 0
    k = h >> 1
    if h % 2 == 0:
        return psi(k) + k
    return psi(k) + (k >> 1) + 1


def q(h):
    return h - psi(h)


def x(h):
    return 2 * psi(h) + 1


def y(h):
    return 2 * h + 1


def xb(b, h):
    return 1 + b * psi(h)


# ---------------------------------------------------------------- E1
print("E1 recursions")
ok = all(q(2 * h) == q(h) and q(2 * h + 1) == q(h) + (h + 1) // 2 for h in range(1 << 14))
check("q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2), h < 2^14", ok)
check("q(h) < h/2 (strict) for 1 <= h < 2^16", all(2 * q(h) < h for h in range(1, 1 << 16)))
check("y/2 < x <= y for h < 2^16", all(y(h) < 2 * x(h) <= 2 * y(h) for h in range(1 << 16)))
check("psi(h) >= 1 for h >= 1 (so x_b >= 1 + b)", all(psi(h) >= 1 for h in range(1, 1 << 14)))


# ---------------------------------------------------------------- E2
def one_apex(h):
    """R(h) as in (G): chain c_1 -> ... -> c_k (k = bit length of 2h), apex u (last index);
    u -> c_i iff the i-th binary digit of 2h, most significant first, is 1."""
    w = 2 * h
    k = w.bit_length()
    n = k + 1
    adj = [[0] * n for _ in range(n)]
    for i in range(k):
        for j in range(i + 1, k):
            adj[i][j] = 1
        if (w >> (k - 1 - i)) & 1:
            adj[k][i] = 1
        else:
            adj[i][k] = 1
    return adj, k  # apex index k


def join(h1, h2, s):
    """J(h1, h2, s): R(h1), then a transitive tournament on s vertices, then the converse of R(h2);
    every vertex beats every vertex of a later block, except that the apex of the last block beats the apex of
    the first."""
    a1, ap1 = one_apex(h1)
    a2, ap2 = one_apex(h2)
    n1, n2 = len(a1), len(a2)
    n = n1 + s + n2
    adj = [[0] * n for _ in range(n)]
    for i in range(n1):
        for j in range(n1):
            adj[i][j] = a1[i][j]
    for i in range(s):
        for j in range(i + 1, s):
            adj[n1 + i][n1 + j] = 1
    for i in range(n2):
        for j in range(n2):
            adj[n1 + s + i][n1 + s + j] = a2[j][i]
    blocks = [range(0, n1), range(n1, n1 + s), range(n1 + s, n)]
    for bi in range(3):
        for bj in range(bi + 1, 3):
            for i in blocks[bi]:
                for j in blocks[bj]:
                    adj[i][j] = 1
    u1, u2 = ap1, n1 + s + ap2
    adj[u1][u2] = 0
    adj[u2][u1] = 1
    return adj, n


def ham_paths(adj, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if c:
                for w in range(n):
                    if adj[v][w] and not (S >> w) & 1:
                        dp[S | (1 << w)][w] += c
    return sum(dp[(1 << n) - 1])


def odd_cycle_masks(adj, n):
    out = []
    for start in range(n):
        stack = [(start, 1 << start, 1)]
        while stack:
            v, mask, length = stack.pop()
            for w in range(start + 1, n):
                if adj[v][w] and not (mask >> w) & 1:
                    stack.append((w, mask | (1 << w), length + 1))
            if length >= 3 and length % 2 == 1 and adj[v][start]:
                out.append(mask)
    return out


def indep_poly(masks):
    alpha = [1]
    level = [(m, i) for i, m in enumerate(masks)]
    while level:
        alpha.append(len(level))
        nxt = []
        for m, i in level:
            for j in range(i + 1, len(masks)):
                if not (m & masks[j]):
                    nxt.append((m | masks[j], j))
        level = nxt
    return alpha


print("E2 tournaments of (G)")
ok = True
for h in range(0, 48):
    adj, ap = one_apex(h)
    n = len(adj)
    if n > 9:
        continue
    al = indep_poly(odd_cycle_masks(adj, n))
    al = al + [0] * (3 - len(al))
    if ham_paths(adj, n) != x(h) or al[1] != psi(h) or any(al[2:]):
        ok = False
check("R(h): #Hamiltonian paths = x(h) and I(Omega, z) = 1 + psi(h) z  (orders <= 9)", ok)
cnt = 0
ok = True
okN = True
for h1 in range(0, 12):
    for h2 in range(0, 12):
        for s in range(0, 4):
            adj, n = join(h1, h2, s)
            if n > 10:
                continue
            cnt += 1
            N = x(h1) * x(h2) + (1 << s) * y(h1) * y(h2) - (1 if s == 0 else 0)
            if ham_paths(adj, n) != N:
                okN = False
            al = indep_poly(odd_cycle_masks(adj, n))
            al = al + [0] * (3 - len(al))
            c3 = ((1 << s) * y(h1) * y(h2) - (1 if s == 0 else 0)) // 2
            if al[0] != 1 or al[1] != psi(h1) + psi(h2) + c3 or al[2] != psi(h1) * psi(h2) or any(al[3:]):
                ok = False
            for b in (2, 3, 4, 5, 6, 10):
                Nb = sum(a * b ** i for i, a in enumerate(al))
                if Nb != xb(b, h1) * xb(b, h2) + b * c3:
                    ok = False
                if Nb % b != 1 or ((Nb - 1) // b - al[1]) % b != 0:
                    ok = False
check("J(h1, h2, s) (%d joins, order <= 10): #Hamiltonian paths = x1 x2 + 2^s y1 y2 - [s = 0]" % cnt, okN)
check("J(h1, h2, s): I(Omega, b) = x_b(h1) x_b(h2) + b (2^s y1 y2 - [s = 0])/2, alpha_1, alpha_2 as stated, "
      "N_b = 1 (mod b), second digit = alpha_1 (mod b)", ok)


# ---------------------------------------------------------------- E3
def a_word(W):
    if W == 0:
        return 0
    if W % 2 == 1:
        return a_word(W // 2)
    return a_word(W // 2) + (W // 2 + 1) // 2


print("E3 THM-5070 (2) dictionary")
check("a(2h) = psi(h), h < 2^12", all(a_word(2 * h) == psi(h) for h in range(1 << 12)))
check("a(4h) = a(2h) + h and a(4h+2) = a(2h+1) + h + 1 = a(2h) + floor(h/2) + 1, h < 2^12",
      all(a_word(4 * h) == a_word(2 * h) + h and a_word(4 * h + 2) == a_word(2 * h + 1) + h + 1
          == a_word(2 * h) + h // 2 + 1 for h in range(1 << 12)))
ok = True
okb = True
for h1 in range(64):
    for s1 in range(5):
        W = (1 << s1) * y(h1) - 1
        if a_word(W) != psi(h1):
            ok = False
        cW, fW = (W + 1) // 2, 1 + W // 2
        if s1 == 0:
            if (cW, fW) != (h1, h1 + 1):
                okb = False
        else:
            if cW != (1 << (s1 - 1)) * y(h1) or fW != cW:
                okb = False
check("a(2^s1 y1 - 1) = psi(h1) (trailing ones do not change a), h1 < 64, s1 < 5", ok)
check("ceil(W/2), 1 + floor(W/2) = (h1, h1 + 1) for s1 = 0 and both 2^(s1-1) y1 for s1 >= 1", okb)
ok = True
for h1 in range(40):
    for h2 in range(40):
        for s in range(5):
            for s1 in range(s + 1):
                s2 = s - s1
                W, Wp = (1 << s1) * y(h1) - 1, (1 << s2) * y(h2) - 1
                br = ((W + 1) // 2) * (1 + Wp // 2) + (1 + W // 2) * ((Wp + 1) // 2)
                tgt = ((1 << s) * y(h1) * y(h2) - (1 if s == 0 else 0))
                if 2 * br != tgt:
                    ok = False
check("bracket of THM-5070 (2) = (2^s y1 y2 - [s = 0])/2 for every split s1 + s2 = s (h1, h2 < 40, s < 5)", ok)


# apex bridge T (.) B built from the words, compared with J by the packing polynomial and the path count
def bridge(W, kW, Wp, kWp):
    """T = TT_k + u, u beats the chain vertices marked by W (most significant bit first; kW chain vertices);
    B = TT_l + w, the chain vertices beating w marked by Wp (least significant bit first; kWp chain vertices);
    T => B with the arc u -> w reversed."""
    nT, nB = kW + 1, kWp + 1
    n = nT + nB
    adj = [[0] * n for _ in range(n)]
    for i in range(kW):
        for j in range(i + 1, kW):
            adj[i][j] = 1
        if (W >> (kW - 1 - i)) & 1:
            adj[kW][i] = 1
        else:
            adj[i][kW] = 1
    o = nT
    for i in range(kWp):
        for j in range(i + 1, kWp):
            adj[o + i][o + j] = 1
        if (Wp >> i) & 1:
            adj[o + i][o + kWp] = 1
        else:
            adj[o + kWp][o + i] = 1
    for i in range(nT):
        for j in range(nT, n):
            adj[i][j] = 1
    u, w = kW, o + kWp
    adj[u][w] = 0
    adj[w][u] = 1
    return adj, n


def canon_score(adj, n):
    """isomorphism invariant strong enough for these small cases: sorted (out-degree, sorted out-degrees of
    out-neighbours, number of 3-cycles through the vertex)."""
    out = [sum(r) for r in adj]
    sig = []
    for v in range(n):
        c3 = 0
        for a in range(n):
            if adj[v][a]:
                for b in range(n):
                    if adj[a][b] and adj[b][v]:
                        c3 += 1
        sig.append((out[v], tuple(sorted(out[w] for w in range(n) if adj[v][w])), c3))
    return tuple(sorted(sig))


ok = True
cnt = 0
for h1 in range(8):
    for h2 in range(8):
        for s in range(3):
            aj, n = join(h1, h2, s)
            if n > 9:
                continue
            for s1 in range(s + 1):
                s2 = s - s1
                k1 = (2 * h1).bit_length() + s1
                k2 = (2 * h2).bit_length() + s2
                ab, nb = bridge((1 << s1) * y(h1) - 1, k1, (1 << s2) * y(h2) - 1, k2)
                cnt += 1
                if nb != n or ham_paths(ab, nb) != ham_paths(aj, n) or canon_score(ab, nb) != canon_score(aj, n) \
                        or indep_poly(odd_cycle_masks(ab, nb)) != indep_poly(odd_cycle_masks(aj, n)):
                    ok = False
check("J(h1, h2, s) and T(W) (.) B(W') agree (order, score invariant, path count, packing polynomial) for every "
      "split, %d comparisons" % cnt, ok)

# ---------------------------------------------------------------- E4
print("E4 THM-5085 (2) and Session206 dictionaries")
ok = True
for b in (2, 3, 4, 5, 6, 7, 10):
    for h1 in range(40):
        for h2 in range(40):
            for s in range(5):
                e0 = 1 if s == 0 else 0
                Nb = xb(b, h1) * xb(b, h2) + b * ((1 << s) * y(h1) * y(h2) - e0) // 2
                h0 = Fraction(xb(b, h1)) + Fraction(b, 2) * ((1 << s) * y(h1) - e0)
                R0 = 2 * (xb(b, h1) + (1 << s) * y(h1))
                hu = xb(b, h1)
                val = h0 + Fraction(b, 2) * (R0 * h2 - hu * 2 * q(h2))
                if val != Nb:
                    ok = False
check("N_b = h0 + (b/2)[R0 h2 - h_u Q(h2)], (h0, R0, h_u) as in the status and in (G), Q = 2q", ok)
ok = True
okg = True
for h1 in range(60):
    for h2 in range(60):
        for s in range(5):
            e0 = 1 if s == 0 else 0
            N = x(h1) * x(h2) + (1 << s) * y(h1) * y(h2) - e0
            hS = x(h1) + (1 << s) * y(h1) - e0
            RS = 2 * (x(h1) + (1 << s) * y(h1))
            dS = x(h1)
            if N != hS + RS * h2 - dS * 2 * q(h2):
                ok = False
            for m in (3, 9, 15, 21, 25, 27, 45, 49, 63, 75, 81, 105):
                if math.gcd(math.gcd(RS, dS), m) != math.gcd(math.gcd(x(h1), y(h1)), m):
                    okg = False
check("(A): N = h_S + R_S h2 - d_S Q(h2) with (h_S, R_S, d_S) = (x1 + 2^s y1 - [s=0], 2(x1 + 2^s y1), x1)", ok)
check("(A): gcd(R_S, d_S, m) = gcd(x1, y1, m) for 12 odd moduli", okg)

# ---------------------------------------------------------------- E5
print("E5 c_d, fibres, L_m")


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def prime_factors(n):
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out


def c_d(m, d):
    r = Fraction(1, d * d)
    for p in prime_factors(m // d):
        r *= 1 - Fraction(1, p * p)
    return r


def L(m, e):
    g = math.gcd(m, e % m) if e % m else m
    return sum(d * c_d(m, d) for d in divisors(g))


ok = True
for m in (9, 15, 21, 25, 27, 45, 49, 63, 75, 81, 105):
    cnt = {d: 0 for d in divisors(m)}
    for al in range(m):
        for be in range(m):
            cnt[math.gcd(math.gcd(al, be), m)] += 1
    if any(Fraction(cnt[d], m * m) != c_d(m, d) for d in divisors(m)) or sum(c_d(m, d) for d in divisors(m)) != 1:
        ok = False
check("c_d = share of the residue vectors with gcd(alpha, beta, m) = d, and sum c_d = 1 (11 odd moduli)", ok)
ok = True
for m in (9, 15, 27, 45, 49, 75):
    for a in range(m):
        for b in range(0, m, max(1, m // 9)):
            d = math.gcd(math.gcd(a, b), m)
            hist = [0] * m
            for w1 in range(m):
                for w2 in range(m):
                    hist[(a * w1 + b * w2) % m] += 1
            if any(hist[c] != (m * d if c % d == 0 else 0) for c in range(m)):
                ok = False
check("(A): w -> a w1 + b w2 takes each multiple of d = gcd(a, b, m) on m d residue vectors and no other value", ok)


def sig1(n):
    return sum(d for d in divisors(n))


ok = True
okm = True
for p, kmax in ((3, 4), (5, 3), (7, 2), (11, 2), (2, 5)):
    for k in range(1, kmax + 1):
        m = p ** k
        for e in range(m):
            v = k
            for i in range(k):
                if e % p ** (i + 1) != 0:
                    v = i
                    break
            if e == 0:
                v = k
            if v < k:
                want = (1 - Fraction(1, p * p)) * Fraction(sig1(p ** v), p ** v)
                want2 = (1 + Fraction(1, p)) * (1 - Fraction(1, p ** (v + 1)))
            else:
                want = (1 - Fraction(1, p * p)) * Fraction(sig1(p ** (k - 1)), p ** (k - 1)) + Fraction(1, p ** k)
                want2 = 1 + Fraction(1, p) - Fraction(1, p ** (k + 1))
            if L(m, e) != want or want != want2:
                ok = False
        if sum(L(m, e) for e in range(m)) != m:
            okm = False
check("L_(p^k)(e): both displayed cases (v < k and v >= k), p^k in 3..81, 5..125, 7, 49, 11, 121, 2..32", ok)
check("mean of L_(p^k) over the classes = 1", okm)
ok = True
for m in (15, 21, 45, 63, 75, 105, 225, 12, 36, 60):
    pp = []
    mm = m
    for p in prime_factors(m):
        pk = 1
        while mm % p == 0:
            mm //= p
            pk *= p
        pp.append(pk)
    for e in range(m):
        prod = Fraction(1)
        for pk in pp:
            prod *= L(pk, e % pk)
        if prod != L(m, e):
            ok = False
    if sum(L(m, e) for e in range(m)) != m:
        ok = False
check("L_m multiplicative over prime powers and of mean 1 (15, 21, 45, 63, 75, 105, 225, 12, 36, 60)", ok)


def table_square(m, u=1):
    mult = [0] * m
    for a1 in range(m):
        for a2 in range(m):
            mult[(a1 * a2) % m] += 1
    multu = [0] * m
    for t in range(m):
        multu[(u * t) % m] += mult[t]
    out = [0] * m
    for t1 in range(m):
        if mult[t1]:
            for t2 in range(m):
                out[(t1 + t2) % m] += mult[t1] * multu[t2]
    return out


ok = True
for m in range(1, 61):
    sq = table_square(m)
    if any(sq[c] != m ** 3 * L(m, c) for c in range(m)):
        ok = False
check("m^3 L_m(c) = #{(a1, a2, b1, b2) mod m : a1 a2 + b1 b2 = c} for EVERY 1 <= m <= 60 (even m included)", ok)
ok = True
for m in (9, 15, 25, 27, 49):
    for u in range(1, m):
        if math.gcd(u, m) == 1:
            if table_square(m, u) != table_square(m, 1):
                ok = False
check("same counts for a1 a2 + u b1 b2, every unit u (9, 15, 25, 27, 49)", ok)


def four_squares(m):
    sqc = [0] * m
    for a in range(m):
        sqc[(a * a) % m] += 1
    two = [0] * m
    for i in range(m):
        if sqc[i]:
            for j in range(m):
                two[(i + j) % m] += sqc[i] * sqc[j]
    four = [0] * m
    for i in range(m):
        if two[i]:
            for j in range(m):
                four[(i + j) % m] += two[i] * two[j]
    return four


ok = all(four_squares(m) == table_square(m) for m in (3, 5, 7, 9, 15, 21, 25, 27, 45, 49, 63, 75, 81, 121))
check("(F): a^2 + b^2 + c^2 + d^2 and a1 a2 + b1 b2 have the same class counts modulo 14 odd moduli", ok)
check("(F): they differ modulo 4 and 8 (the prime 2)", four_squares(4) != table_square(4)
      and four_squares(8) != table_square(8))
f4 = four_squares(16)
check("(F): four squares have 2-adic density 1 at every odd class modulo 16",
      all(f4[c] == 16 ** 3 for c in range(1, 16, 2)))
ok = True
for n in range(1, 200, 2):
    r4 = 0
    lim = int(math.isqrt(n))
    for a in range(-lim, lim + 1):
        for b in range(-lim, lim + 1):
            r = n - a * a - b * b
            if r < 0:
                continue
            for c in range(-lim, lim + 1):
                d2 = r - c * c
                if d2 < 0:
                    continue
                d = math.isqrt(d2)
                if d * d == d2:
                    r4 += 1 if d == 0 else 2
    if r4 != 8 * sig1(n):
        ok = False
check("Jacobi: r_4(n) = 8 sigma_1(n) for odd n < 200", ok)
# Remark 1
ok = True
for p in (3, 5, 7, 11):
    sf = 1 - Fraction(1, p * p) + Fraction(1, p)
    f1 = (1 + Fraction(1, p)) * (1 - Fraction(1, p ** 2))
    f2 = (1 + Fraction(1, p)) * (1 - Fraction(1, p ** 3))
    if f1 - sf != -Fraction(1, p ** 3) or f2 <= sf:
        ok = False
    for k in (1, 2, 3):
        m = p ** k
        div = [L(m, e) for e in range(0, m, p)]
        if sum(div) / len(div) != sf:
            ok = False
check("Remark 1: factor at v = 1 is the squarefree factor minus p^-3, larger from v = 2 on; the mean of L_(p^k) "
      "over the classes divisible by p is 1 - p^-2 + 1/p for k = 1, 2, 3", ok)

# ---------------------------------------------------------------- E6
print("E6 base b: region, value, content, second digit")
ok = True
okV = True
okc = True
for b in (2, 3, 4, 5, 6, 7, 9, 10):
    for h1 in range(50):
        for h2 in range(50):
            for s in range(5):
                e0 = 1 if s == 0 else 0
                Nb = xb(b, h1) * xb(b, h2) + b * ((1 << s) * y(h1) * y(h2) - e0) // 2
                a = Fraction(b * xb(b, h1), 2)
                bp = Fraction(b * (1 << s) * y(h1), 2)
                k = xb(b, h1) * (1 - Fraction(b, 2)) - Fraction(b, 2) * e0
                if a * x(h2) + bp * y(h2) + k != Nb or k > 0 or a <= 0 or bp <= 0:
                    ok = False
                if xb(b, h2) != (1 - Fraction(b, 2)) + Fraction(b, 2) * x(h2):
                    ok = False
                if s >= 1:
                    V = Nb
                    co = (xb(b, h1), b * (1 << (s - 1)) * y(h1))
                else:
                    V = 2 * Nb + b
                    co = (2 * xb(b, h1), b * y(h1))
                if V != co[0] * xb(b, h2) + co[1] * y(h2):
                    okV = False
                if Nb < b * (1 << s) * y(h1) * y(h2) // 2 - b // 2 - 1 or 2 * Nb < b * (1 << s) * y(h1) * y(h2) - b:
                    okV = False
                for m in (3, 7, 9, 11, 21, 49, 77):
                    if math.gcd(m, b) != 1:
                        continue
                    gb = math.gcd(math.gcd(xb(b, h1), y(h1)), m)
                    if math.gcd(math.gcd(co[0], co[1]), m) != gb or V % gb != 0:
                        okc = False
check("N_b = a x2 + b' y2 + k with a = b x_b(h1)/2, b' = b 2^(s-1) y1 (b y1/2 at s = 0), "
      "k = x_b(h1)(1 - b/2) - (b/2)[s = 0] <= 0", ok)
check("V = N_b (s >= 1), 2 N_b + b (s = 0) is the stated linear form in (x_b(h2), y2); "
      "N_b >= b 2^(s-1) y1 y2 - b/2", okV)
check("content of the form modulo odd m prime to b is gcd(g_b(h1), m), and type d => d | V (exact)", okc)
ok = True
for b in (3, 5, 7, 9, 15):
    for m in (3, 5, 7, 9, 11, 13, 21, 25, 27, 49):
        if math.gcd(m, 2 * b) != 1:
            continue
        img = set()
        for qq in range(m):
            for hh in range(m):
                img.add(((1 + b * (hh - qq)) % m, (2 * hh + 1) % m))
        if len(img) != m * m:
            ok = False
check("(q, h) -> (1 + b (h - q), 2h + 1) is a bijection of (Z/m)^2 for odd m prime to b", ok)
ok = True
okcount = True
for b in (3, 5, 7, 9, 11, 15):
    inv2 = pow(2, -1, b)
    for h1 in range(30):
        for h2 in range(30):
            for s in range(5):
                e0 = 1 if s == 0 else 0
                al1 = psi(h1) + psi(h2) + ((1 << s) * y(h1) * y(h2) - e0) // 2
                if s >= 1:
                    f = inv2 * x(h2) + ((1 << (s - 1)) * y(h1)) * y(h2) + (psi(h1) - inv2)
                else:
                    f = inv2 * x(h2) + (inv2 * y(h1)) * y(h2) + (psi(h1) - 1)
                if (al1 - f) % b != 0:
                    ok = False
    for cy in range(b):
        for const in range(b):
            hist = [0] * b
            for w1 in range(b):
                for w2 in range(b):
                    hist[(inv2 * w1 + cy * w2 + const) % b] += 1
            if any(hh != b for hh in hist):
                okcount = False
check("second digit: the two displayed congruences for alpha_1 modulo odd b (b = 3, 5, 7, 9, 11, 15)", ok)
check("second digit: an affine form with the unit coefficient 2^(-1) at x2 takes each class on exactly b of the "
      "b^2 residue vectors", okcount)
okp = True
okodd = True
for b in (2, 3, 5, 6, 7, 10):
    for h1 in range(20):
        for h2 in range(20):
            for s in range(4):
                Nb = xb(b, h1) * xb(b, h2) + b * ((1 << s) * y(h1) * y(h2) - (1 if s == 0 else 0)) // 2
                if Nb % b != 1 or (2 * Nb + b) % b != 2 % b:
                    okp = False
                V = Nb if s >= 1 else 2 * Nb + b
                n = V
                while n % 2 == 0:
                    n //= 2
                if math.gcd(n, 2 * b) != 1:
                    okodd = False
check("N_b = 1 (mod b) and 2 N_b + b = 2 (mod b) on 9600 triples (b = 2, 3, 5, 6, 7, 10)", okp)
check("the odd part n of the value (N_b on sheet 1, 2 N_b + b on sheet 0) is prime to 2b", okodd)

# ---------------------------------------------------------------- E7
print("E7 balanced states, early boosts")
N24 = 1 << 20
psi_tab = [0] * N24
for h in range(1, N24):
    k = h >> 1
    psi_tab[h] = psi_tab[k] + (k if h % 2 == 0 else (k >> 1) + 1)
bal3 = [h for h in range(N24) if 1 + 3 * psi_tab[h] == 2 * h + 1]
check("balanced states of base 3 below 2^20: 0, 3, 45, 117, 309, 9045, 62805, 111957, 382293",
      bal3 == [0, 3, 45, 117, 309, 9045, 62805, 111957, 382293], str(bal3))
check("base b >= 4: only h = 0 is balanced below 2^20 (b = 4..12)",
      all([h for h in range(N24) if 1 + b * psi_tab[h] == 2 * h + 1] == [0] for b in range(4, 13)))
bal2 = [h for h in range(N24) if 1 + 2 * psi_tab[h] == 2 * h + 1]
check("base 2: h = 0 and the powers of 2", bal2 == [0] + [1 << j for j in range(20)])
check("psi(3) = 2, y(3) = 7; x_b(3) = 1 + 2b is divisible by 7 iff b = 3 (mod 7) (b <= 60)",
      psi(3) == 2 and y(3) == 7 and all(((1 + 2 * b) % 7 == 0) == (b % 7 == 3) for b in range(2, 61)))
check("base 3: x_b(3) = y(3) = 7; base 10: x_b(3) = 21 = 3 y(3); base 5: h = 6 has x_b = 26 = 2 y, y = 13",
      xb(3, 3) == 7 and xb(10, 3) == 21 and psi(6) == 5 and xb(5, 6) == 26 and y(6) == 13)
first = {}
for (b, p) in ((3, 7), (10, 7), (5, 13)):
    for h in range(1, 10000):
        if xb(b, h) % p == 0 and y(h) % p == 0:
            first[(b, p)] = h
            break
check("first p-divisible states: base 3, p = 7: h = 3; base 10, p = 7: h = 3; base 5, p = 13: h = 6",
      first == {(3, 7): 3, (10, 7): 3, (5, 13): 6}, str(first))
# first states of type 25 (base 2)
f25 = next(h for h in range(1, N24) if (2 * psi_tab[h] + 1) % 25 == 0 and (2 * h + 1) % 25 == 0)
check("Remark 2: the first state with 25 | g is h = 512, x = y = 1025; h = 2 has x = y = 5",
      f25 == 512 and x(512) == 1025 and y(512) == 1025 and x(2) == 5 and y(2) == 5, str(f25))

print()
if FAIL:
    print("FAILURES: %d" % len(FAIL))
    for f in FAIL:
        print("   " + f)
    sys.exit(1)
print("ALL EXACT CHECKS PASSED")
