#!/usr/bin/env python3
"""Independent small audits (own code, exact arithmetic):
 A. signed inverse fibre F(n,s)=(2n+s,-s): U_{-s}(2n+s)=U_s(n), F^2=S_s, Jacobsthal root fibre,
    bijection (n,s)->(y,k); forward-conjugacy hostile.
 B. a_k=(10^k-7)/3: primes k=2..8, composite witnesses k=9..17, a_18 prime (deterministic MR +
    Pocklington with the stated n-1 factorization), 17|a_k iff k=9 mod 16, no prime <=13 divides.
 C. four-vertex tournaments: |Pf| of the sign matrix, switching invariance, class table, and
    relation to H (Hamiltonian paths) and c3.
 D. rational Collatz on marked triples: C'<C/3 (plus, 3|t), C'<C/4 (minus legal, 3|t);
    fixed points 1/(2^k-3); exact r-step entry for denominators 3^r.
 E. ordered carries: R(u01v)-R(u10v)=2^j 3^b, M(u01v)-M(u10v)=-2^j 3^h (all words L<=10);
    root-4 H-image -4/15; tail 266/243.
 F. Mills cubic greedy chain 2,11,1361,2521008887 (next prime above p^3).
 G. divisor balance F=S+U iff N in {p, p^3, p^2qr} for N<=20000.
 H. orbit count (10^k-7)/3 for k<=4 by brute force over words.
 I. square roots of permutations of cycle type 2^2 5 18^2 (36) and 2^2 5^2 18^2 (216) by formula,
    formula validated by brute force on all permutations of n<=7.
 J. residual language at depth 14: 734 classes, 142 accepted ports.
 K. sewing residual rank 2 (small lengths) and the 9 -> 4 certificate.
"""
from fractions import Fraction
from itertools import permutations, combinations, product
from math import gcd, isqrt, comb


def v2(m):
    c = 0
    while m % 2 == 0:
        m //= 2
        c += 1
    return c


def U(n, s):
    m = 3 * n + s
    while m % 2 == 0:
        m //= 2
    return m


# ---------------- A
def audit_A():
    ok = True
    for s in (1, -1):
        for n in range(1, 200001, 2):
            if s == -1 and n == 1:
                pass
            if U(2 * n + s, -s) != U(n, s):
                ok = False
            if v2(3 * (2 * n + s) - s) != v2(3 * n + s) + 1:
                ok = False
            if 4 * n + s != 2 * (2 * n + s) - s:
                ok = False
    print("A1 U_{-s}(2n+s)=U_s(n), exponent+1, F^2=(4n+s,s): all odd n<=200000 both signs:", ok)
    # bijection (n,s) -> (y,k) on a box
    seen = {}
    dup = False
    for s in (1, -1):
        for n in range(1, 20001, 2):
            y = U(n, s)
            k = v2(3 * n + s)
            if (y, k) in seen:
                dup = True
            seen[(y, k)] = (n, s)
    # every (y,k) with y odd, 3 !| y, small, k<=12 is hit
    miss = []
    for y in range(1, 200, 2):
        if y % 3 == 0:
            continue
        for k in range(1, 13):
            s = 1 if (2 ** k * y) % 3 == 1 else -1
            n = (2 ** k * y - s) // 3
            if n <= 20000 and seen.get((y, k)) != (n, s):
                miss.append((y, k))
    print("A2 germ map injective on box:", not dup, "; inverse formula consistent:", not miss)
    root = []
    for k in range(1, 11):
        s = 1 if (2 ** k) % 3 == 1 else -1
        root.append(((2 ** k - s) // 3, s))
    print("A3 root fibre (y=1):", root)
    print("A4 hostile: U_-(5)=%d, F(5,-)=(9,+), U_+(9)=%d; next: U_-(7)=%d vs U_+(7)=%d"
          % (U(5, -1), U(9, 1), U(7, -1), U(7, 1)))
    # shortcut form: T_-(2n+1) = 2 T_+(n) for odd n; T_+(2n-1) = 2 T_-(n) for odd n
    T = lambda n, s: n // 2 if n % 2 == 0 else (3 * n + s) // 2
    okT = all(T(2 * n + 1, -1) == 2 * T(n, 1) and T(2 * n - 1, 1) == 2 * T(n, -1) for n in range(1, 100001, 2))
    print("A5 shortcut intertwining T_-(2n+1)=2T_+(n), T_+(2n-1)=2T_-(n) (odd n<=1e5):", okT)


# ---------------- B
def is_prime_mr(n):
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]
    for p in small:
        if n % p == 0:
            return n == p
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in small:  # deterministic for n < 3.3e24
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def audit_B():
    a = lambda k: (10 ** k - 7) // 3
    prim = [(k, is_prime_mr(a(k))) for k in range(2, 19)]
    print("B1 primality k=2..18:", prim)
    wit = {9: 17, 10: 673, 11: 307, 12: 19, 13: 523, 14: 607, 15: 181, 16: 199, 17: 31}
    print("B2 witnesses divide and are proper:", all(a(k) % p == 0 and 1 < p < a(k) for k, p in wit.items()))
    n = a(18)
    fac = [2, 3, 5, 2071723, 5363222357]
    prod_ok = 2 * 3 * 5 * 2071723 * 5363222357 == n - 1
    facs_prime = all(is_prime_mr(q) for q in fac)
    lucas = pow(2, n - 1, n) == 1 and all(gcd(pow(2, (n - 1) // q, n) - 1, n) == 1 for q in fac)
    print("B3 a18 =", n, " n-1 factorization ok:", prod_ok, " factors prime:", facs_prime, " Lucas/Pocklington ok:", lucas)
    per = all((a(k) % 17 == 0) == (k % 16 == 9) for k in range(1, 400))
    print("B4 17|a_k iff k=9 mod 16 (k<400):", per)
    small = all(a(k) % p != 0 for k in range(1, 400) for p in (2, 3, 5, 7, 11, 13))
    print("B5 no prime <=13 divides a_k (k<400):", small)


# ---------------- C
def audit_C():
    V = range(4)
    pairs = list(combinations(V, 2))
    rows = []
    for mask in range(64):
        S = [[0] * 4 for _ in V]
        for i, (a, b) in enumerate(pairs):
            if (mask >> i) & 1:
                S[a][b], S[b][a] = 1, -1
            else:
                S[a][b], S[b][a] = -1, 1
        pf = S[0][1] * S[2][3] - S[0][2] * S[1][3] + S[0][3] * S[1][2]
        scores = tuple(sorted(sum(1 for j in V if S[i][j] == 1) for i in V))
        H = sum(1 for p in permutations(V) if all(S[p[i]][p[i + 1]] == 1 for i in range(3)))
        c3 = sum(1 for t in combinations(V, 3)
                 if (S[t[0]][t[1]] == S[t[1]][t[2]] == S[t[2]][t[0]]))
        # switching invariance of |Pf|
        inv = True
        for d in product((1, -1), repeat=4):
            S2 = [[d[i] * S[i][j] * d[j] for j in V] for i in V]
            pf2 = S2[0][1] * S2[2][3] - S2[0][2] * S2[1][3] + S2[0][3] * S2[1][2]
            if abs(pf2) != abs(pf):
                inv = False
        rows.append((scores, abs(pf), H, c3, inv))
    table = {}
    for r in rows:
        table.setdefault(r[:4], 0)
        table[r[:4]] += 1
    print("C1 (scores, |Pf|, H, c3): count ->", dict(sorted(table.items())))
    print("C2 |Pf| switching-invariant for all 64:", all(r[4] for r in rows))
    print("C3 |Pf| == H mod 4 for all 64:", all(r[1] % 4 == r[2] % 4 for r in rows),
          "; |Pf|=3 iff c3 odd:", all((r[1] == 3) == (r[3] % 2 == 1) for r in rows))


# ---------------- D
def audit_D():
    worst_p = Fraction(0)
    worst_m = Fraction(0)
    cnt_p = cnt_m = 0
    for t in range(3, 400, 6):  # odd multiples of 3
        for s in range(1, 800, 2):
            if gcd(s, t) != 1:
                continue
            C = Fraction(s * s + t * t, 2)
            for sg in (1, -1):
                num = 3 * s + sg * t
                if num <= 0:
                    continue
                k = v2(num)
                g = gcd(num >> k, t)
                s2, t2 = (num >> k) // g, t // g
                C2 = Fraction(s2 * s2 + t2 * t2, 2)
                r = C2 / C
                if sg == 1:
                    cnt_p += 1
                    worst_p = max(worst_p, r)
                else:
                    cnt_m += 1
                    worst_m = max(worst_m, r)
    print("D1 plus C'/C max over %d cases: %s (<1/3: %s); minus max over %d: %s (<1/4: %s)"
          % (cnt_p, float(worst_p), worst_p < Fraction(1, 3), cnt_m, float(worst_m), worst_m < Fraction(1, 4)))
    fx = []
    for k in range(3, 20):
        x = Fraction(1, 2 ** k - 3)
        y = 3 * x + 1
        e = 0
        while y.numerator % 2 == 0:
            y /= 2
            e += 1
        fx.append(y == x)
    print("D2 1/(2^k-3) fixed by U_+ for k=3..19:", all(fx))
    # exact r-step entry
    okr = True
    for r in range(0, 7):
        for s in range(1, 500, 2):
            if s % 3 == 0 and r > 0:
                continue
            x = Fraction(s, 3 ** r)
            steps = 0
            while x.denominator != 1:
                y = 3 * x + 1
                while y.numerator % 2 == 0:
                    y /= 2
                x = y
                steps += 1
            if steps != r:
                okr = False
    print("D3 denominator 3^r enters Z in exactly r plus steps (r<=6, s<500):", okr)


# ---------------- E
def carries(w):
    L = len(w)
    a = sum(w)
    R = 0
    ones_after = [0] * (L + 1)
    for j in range(L - 1, -1, -1):
        ones_after[j] = ones_after[j + 1] + w[j]
    for j, e in enumerate(w):
        if e:
            R += 2 ** j * 3 ** ones_after[j + 1]
    M = sum(e * 2 ** j * 3 ** (L - 1 - j) for j, e in enumerate(w))
    return R, M


def audit_E():
    ok = True
    n = 0
    for L in range(2, 11):
        for w in product((0, 1), repeat=L):
            for j in range(L - 1):
                if w[j] == 0 and w[j + 1] == 1:
                    w2 = list(w)
                    w2[j], w2[j + 1] = 1, 0
                    R1, M1 = carries(w)
                    R2, M2 = carries(w2)
                    b = sum(w[j + 2:])
                    h = L - j - 2
                    n += 1
                    if R1 - R2 != 2 ** j * 3 ** b or M1 - M2 != -(2 ** j) * 3 ** h:
                        ok = False
    print("E1 adjacent exchange laws over %d exchanges (L<=10):" % n, ok)
    # verify R via direct iteration on the cylinder
    okc = True
    for L in range(1, 9):
        for w in product((0, 1), repeat=L):
            R, M = carries(w)
            a = sum(w)
            # find n in cylinder
            for n0 in range(0, 2 ** L):
                x = n0
                good = True
                for e in w:
                    if x % 2 != e:
                        good = False
                        break
                    x = x // 2 if e == 0 else (3 * x + 1) // 2
                if good:
                    if 2 ** L * x != 3 ** a * n0 + R:
                        okc = False
                    break
    print("E2 2^L C^L(n)=3^a n+R(w) on cylinders (L<=8):", okc)
    x = Fraction(-3, 5)
    y = x
    for e in (1, 0):
        y = (3 * y + e) / 2
    print("E3 H-fixed point of period '10' is -3/5:", y == x, "; root-4 image (4/9)(-3/5) =", Fraction(4, 9) * x)
    print("E4 Mahler tail of 10101:", sum(Fraction(2, 3) ** (k + 1) for k, e in enumerate((1, 0, 1, 0, 1)) if e))


# ---------------- F
def audit_F():
    chain = [2]
    for _ in range(3):
        p = chain[-1]
        q = p ** 3 + 1
        while not is_prime_mr(q):
            q += 1
        assert q < (p + 1) ** 3
        chain.append(q)
    print("F1 greedy cubic Mills prefix:", chain)
    # bracket: A^(3^4) in [2521008887, 2521008888)
    lo = Fraction(1306377883863080, 10 ** 15)
    hi = Fraction(1306377883869479, 10 ** 15)
    # check lo^81 < 2521008887 and hi^81 > 2521008888 (outer bracket)
    print("F2 outer bracket: lo^81 < p4:", lo ** 81 < 2521008887, "; hi^81 > p4+1:", hi ** 81 > 2521008888)


# ---------------- G
def audit_G():
    Nmax = 20000
    spf = list(range(Nmax + 1))
    for i in range(2, isqrt(Nmax) + 1):
        if spf[i] == i:
            for j in range(i * i, Nmax + 1, i):
                if spf[j] == j:
                    spf[j] = i
    bad = []
    for N in range(2, Nmax + 1):
        divs = [d for d in range(2, N) if N % d == 0] if N < 3000 else None
        f = {}
        m = N
        while m > 1:
            p = spf[m]
            f[p] = f.get(p, 0) + 1
            m //= p
        exps = sorted(f.values())
        r = len(f)
        tau = 1
        for e in exps:
            tau *= e + 1
        F = tau - 2
        S = 2 ** r - 1 - (1 if all(e == 1 for e in exps) else 0)
        Uc = r - (1 if (r == 1 and exps == [1]) else 0)
        if divs is not None:
            F2 = len(divs)
            S2 = sum(1 for d in divs if all(d % (p * p) for p in range(2, isqrt(d) + 1)))
            U2 = sum(1 for d in divs if spf[d] == d)
            if (F, S, Uc) != (F2, S2, U2):
                bad.append(("formula", N))
        balanced = F == S + Uc
        shape = (exps == [1]) or (exps == [3]) or (exps == [1, 1, 2])
        if balanced != shape:
            bad.append(("shape", N))
    print("G1 divisor balance F=S+U iff p, p^3, p^2qr (N<=%d); formula vs direct (N<3000):" % Nmax, not bad, bad[:3])


# ---------------- H
def audit_H():
    V = [(0, 0), (1, 0), (0, 1), (1, 1)]  # 0,R,G,B
    rho = {(0, 0): (0, 0), (1, 0): (0, 1), (0, 1): (1, 1), (1, 1): (1, 0)}
    D = []
    for i in range(4):
        for j in range(i, 4):
            D.append(frozenset([V[i], V[j]]) if i != j else frozenset([V[i]]))
    # represent repeated pair {v,v} distinctly
    D = [(V[i], V[j]) for i in range(4) for j in range(i, 4)]
    canon = lambda p: tuple(sorted(p))
    act = lambda p: canon((rho[p[0]], rho[p[1]]))
    Sset = [p for p in D if p[0] != p[1]] + [((0, 0), (0, 0))]
    res = []
    for k in range(1, 5):
        words = set(product([canon(p) for p in D], repeat=k))
        for s in Sset:
            words.discard(tuple([canon(s)] * k))
        orbits = 0
        seen = set()
        for w in words:
            if w in seen:
                continue
            orb = {w}
            x = w
            for _ in range(2):
                x = tuple(act(p) for p in x)
                orb.add(x)
            seen |= orb
            orbits += 1
        res.append((k, orbits, (10 ** k - 7) // 3))
    print("H1 residual C3-orbit counts vs (10^k-7)/3:", res)


# ---------------- I
def nroots_formula(ctype):
    # ctype: dict length -> count
    from math import factorial
    def dfact(m):  # (m-1)!! for even m? here: number of perfect matchings of 2j items = (2j-1)!!
        r = 1
        for x in range(m - 1, 0, -2):
            r *= x
        return r
    total = 1
    for m, c in ctype.items():
        if m % 2 == 0:
            if c % 2:
                return 0
            total *= dfact(c) * m ** (c // 2)
        else:
            s = 0
            for j in range(0, c // 2 + 1):
                s += comb(c, 2 * j) * dfact(2 * j) * m ** j
            total *= s
    return total


def audit_I():
    ok = True
    for n in range(1, 8):
        perms = list(permutations(range(n)))
        sq = {}
        for p in perms:
            s = tuple(p[p[i]] for i in range(n))
            sq[s] = sq.get(s, 0) + 1
        for p in perms:
            seen = [False] * n
            ct = {}
            for i in range(n):
                if not seen[i]:
                    L = 0
                    j = i
                    while not seen[j]:
                        seen[j] = True
                        j = p[j]
                        L += 1
                    ct[L] = ct.get(L, 0) + 1
            if nroots_formula(ct) != sq.get(p, 0):
                ok = False
    print("I1 square-root count formula validated on all permutations n<=7:", ok)
    print("I2 roots of type 2^2 5 18^2:", nroots_formula({2: 2, 5: 1, 18: 2}),
          "; of type 2^2 5^2 18^2:", nroots_formula({2: 2, 5: 2, 18: 2}),
          "; odd-only 1,2,7 +extra 2:", nroots_formula({1: 1, 2: 2, 7: 1}),
          "; shortcut 1,3,11:", nroots_formula({1: 1, 3: 1, 11: 1}))


# ---------------- J
def audit_J(D=14):
    residual = 0
    accepted = 0
    # BFS over words, accept at first prefix with 3^ones < 2^len
    frontier = [(0, 0)]  # (length, ones)
    from collections import Counter
    cur = Counter({(0, 0): 1})
    for L in range(1, D + 1):
        nxt = Counter()
        for (l, o), c in cur.items():
            for e in (0, 1):
                o2 = o + e
                if 3 ** o2 < 2 ** L:
                    accepted += c
                else:
                    nxt[(L, o2)] += c
        cur = nxt
    residual = sum(cur.values())
    print("J1 depth %d: accepted ports %d, residual classes %d" % (D, accepted, residual))


# ---------------- K
def audit_K():
    # 9 -> 4 certificate
    def mat(w, s=1):
        A, B, Dd = 1, 0, 1
        for e in w:
            if e == 1:
                A, B = 3 * A, 3 * B + s * Dd
            Dd *= 2
        return A, B, Dd
    u = (1, 0)
    v = (1, 1, 1, 0, 1, 0, 0, 1, 0)
    Au, Bu, Du = mat(u)
    Av, Bv, Dv = mat(v)
    left = Fraction(Au * 9 + Bu, Du)
    right = Fraction(Dv * 4 - Bv, Av)
    print("K1 M_u=%s M_v=%s seam left=%s right=%s" % ((Au, Bu, Du), (Av, Bv, Dv), left, right))
    # rank of sewing residual matrix for small lengths
    import itertools
    def rank_Q(M):
        M = [list(map(Fraction, r)) for r in M]
        rk = 0
        rows, cols = len(M), len(M[0])
        for c in range(cols):
            piv = next((i for i in range(rk, rows) if M[i][c] != 0), None)
            if piv is None:
                continue
            M[rk], M[piv] = M[piv], M[rk]
            for i in range(rows):
                if i != rk and M[i][c] != 0:
                    f = M[i][c] / M[rk][c]
                    M[i] = [x - f * y for x, y in zip(M[i], M[rk])]
            rk += 1
        return rk
    ranks = set()
    for s in (1, -1):
        for (h, k) in ((1, 1), (2, 3), (3, 2), (4, 4)):
            for (n, r) in ((9, 4), (7, 11), (27, 5)):
                rowsM = []
                for uu in itertools.product((0, 1), repeat=h):
                    Au, Bu, Du = mat(uu, s)
                    row = []
                    for vv in itertools.product((0, 1), repeat=k):
                        Av, Bv, Dv = mat(vv, s)
                        row.append(Av * (Au * n + Bu) + Du * (Bv - Dv * r))
                    rowsM.append(row)
                ranks.add(rank_Q(rowsM))
    print("K2 sewing residual ranks over tested (sign, lengths, n, r):", ranks)


if __name__ == "__main__":
    audit_A()
    audit_B()
    audit_C()
    audit_D()
    audit_E()
    audit_F()
    audit_G()
    audit_H()
    audit_I()
    audit_J()
    audit_K()
