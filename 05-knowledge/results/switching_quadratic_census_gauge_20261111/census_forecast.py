#!/usr/bin/env python3
"""A000568 and A049313: level expansion, explicit uniform truncation bound, certified enclosures.

Notation: a(n) = A000568(n), b(n) = 2^C(n,2)/n!, f_n = a(n)/b(n) = sum_k L_k(n), L_k(n) = P_k(n) 4^(-kn),
where level k collects the cycle types mu (odd parts >= 3) with defect m - r = 2k (m = support, r = #cycles).
(The expansion and the polynomials P_k are Session98 / THM-4845 / Session271; this file adds the explicit bound
B_K(n) on the remainder after K levels, uniform in n and K, the ball-arithmetic enclosures, and the same for the
census of switching classes A049313.)

usage: census_forecast.py [check|bound|giant|switching|strong|heights|all]
"""
import sys, json, os
from fractions import Fraction
from math import comb, gcd, factorial

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'session270_chanpak_enumeration_20261111', 'enumeration', 'prefix.json')
A049313 = [1,1,1,2,2,6,12,79,792,19576,886288,75369960,11856006240,3467430423264,1893448825054528,
           1938818712501985736,3737086626658278741376,13606268915761294708760704,93863103860384959101157737728]  # n = 1..19 (OEIS)

def partitions(k, maxpart=None):
    if maxpart is None: maxpart = k
    if k == 0: yield (); return
    for p in range(min(k, maxpart), 0, -1):
        for rest in partitions(k - p, p): yield (p,) + rest

def falling(n, m):
    r = 1
    for i in range(m): r *= (n - i)
    return r

def types_at_level(k):
    """cycle types at level k: parts h of a partition of k -> cycles of length 2h+1. Returns (m, r, z, e)."""
    out = []
    for hs in partitions(k):
        cyc = [2*h + 1 for h in hs]; m = sum(cyc); r = len(cyc)
        z = 1
        for a in set(cyc): c = cyc.count(a); z *= a**c * factorial(c)
        e = sum((a - 1)//2 for a in cyc) + sum(gcd(cyc[i], cyc[j]) for i in range(r) for j in range(i+1, r))
        out.append((m, r, z, e))
    return out
_TYPES = {}
def level_sum(k, n, shift=0):
    """L_k(n) (exact Fraction); with shift=1 the weight of the switching-class census, 4^k L_k(n)."""
    if k not in _TYPES: _TYPES[k] = types_at_level(k)
    tot = Fraction(0)
    for (m, r, z, e) in _TYPES[k]:
        if m > n: continue
        D = (n - m)*(m - r) + comb(m, 2) - e
        tot += Fraction(falling(n, m), z) * Fraction(1, 2**D)
    return tot * (4**k if shift else 1)

def t_term(k, n):
    """t_k(n) = C(n,2)^(2k) 2^(3k^2+k-2kn) / (2k)!  >= L_k(n) for 3k <= n."""
    return Fraction(comb(n, 2)**(2*k) * 2**(3*k*k + k), factorial(2*k) * 4**(k*n))

def B_bound(K, n, shift=0):
    """Explicit bound on sum_{k>K} L_k(n) (shift=1: on sum_{k>K} 4^k L_k(n))."""
    tot = Fraction(0)
    for k in range(K + 1, n//3 + 1): tot += t_term(k, n) * (4**k if shift else 1)
    tail = Fraction(factorial(n), 2**((n*(n-1) + 2)//3))      # n! 2^(-ceil(n(n-1)/3)) : n(n-1)/3 bits lost, exponent is an integer bound
    return tot + tail * (4**((n - 1)//2) if shift else 1)

def sci(x, d=6):
    """decimal scientific string of a positive Fraction (no float underflow)"""
    if x == 0: return "0"
    from math import log10
    num, den = x.numerator, x.denominator
    e = len(str(num)) - len(str(den))
    y = Fraction(num, den) / Fraction(10)**e if e >= 0 else Fraction(num, den) * Fraction(10)**(-e)
    while y >= 10: y /= 10; e += 1
    while y < 1: y *= 10; e -= 1
    return "%.*fe%+d" % (d, float(y), e)

def exact_values():
    return [int(v) for v in json.load(open(DATA))['values']]

def burnside(n):
    """independent exact a(n) by the odd-cycle Burnside sum"""
    tot = Fraction(0)
    def odd_parts(n, maxp):
        if n == 0: yield (); return
        for p in range(min(n, maxp), 0, -1):
            if p % 2 == 1:
                for rest in odd_parts(n - p, p): yield (p,) + rest
    for lam in odd_parts(n, n):
        z = 1
        for a in set(lam): c = lam.count(a); z *= a**c * factorial(c)
        o = sum((a - 1)//2 for a in lam) + sum(gcd(lam[i], lam[j]) for i in range(len(lam)) for j in range(i+1, len(lam)))
        tot += Fraction(2**o, z)
    assert tot.denominator == 1
    return tot.numerator

def cmd_check():
    a = exact_values()
    print("control 1: independent Burnside sum equals the recorded exact values for n <= 40:", all(burnside(n) == a[n] for n in range(1, 41)))
    ok = True
    for n in range(3, 61):
        f = Fraction(a[n] * factorial(n), 2**comb(n, 2))
        s = sum(level_sum(k, n) for k in range(0, (n - 1)//2 + 1))
        ok &= (s == f)
    print("control 2: f_n = sum_k L_k(n) exactly (all levels) for 3 <= n <= 60:", ok)

def cmd_bound():
    a = exact_values()
    print("remainder after K levels, T_K(n) = f_n - sum_{k<=K} L_k(n), against the explicit bound B_K(n)")
    print("  n  K   T_K(n)        B_K(n)        B/T     leading omitted L_{K+1}(n)")
    viol = 0; cases = 0
    for n in list(range(6, 61)) + [80, 100, 150, 200, 250]:
        f = Fraction(a[n] * factorial(n), 2**comb(n, 2))
        part = Fraction(0)
        for K in range(0, 7):
            part += level_sum(K, n)
            if 2*(K+1) + 1 > n: break
            T = f - part; B = B_bound(K, n); cases += 1
            if T > B: viol += 1
            if n in (8, 12, 16, 20, 30, 50, 100, 200) and K <= 3:
                print("%4d %2d   %s  %s  %9.3g  %s" % (n, K, sci(T), sci(B), float(B/T) if T else float('inf'), sci(level_sum(K+1, n))))
    print("cases checked: %d (6 <= n <= 60 and n = 80, 100, 150, 200, 250; K <= 6); violations of T_K(n) <= B_K(n): %d" % (cases, viol))
    n0 = next(n for n in range(3, 200) if all(B_bound(0, m) < 1 for m in range(n, n + 30)))
    print("B_0(n) < 1 for all n >= %d (checked on the next 30 values)" % n0)

def cmd_giant():
    from flint import arb, fmpz, ctx
    a = exact_values()
    def enclosure(n, K, digits=30):
        """returns (E, lead) certified: a(n) = lead... x 10^E with `digits`+1 leading digits, or None"""
        nbits = max(400, int(2.2 * (2 * len(str(n)) * 3.33 + digits * 3.4 + 200)))
        ctx.prec = nbits
        N = fmpz(n); pairs = N * (N - 1) / 2 if False else fmpz(n * (n - 1) // 2)
        ln2 = arb(2).log(); ln10 = arb(10).log()
        logb = arb(pairs) * ln2 - (arb(N) + 1).lgamma()                 # ln b(n)
        # exact part C = sum_{k<=K} L_k(n) and bound on the rest, as balls
        def L_ball(k):
            tot = arb(0)
            for (m, r, z, e) in types_at_level(k):
                if m > n: continue
                # (n)_m / z * 2^(-D),  D = (n-m)(m-r) + C(m,2) - e
                ff = arb(1)
                for i in range(m): ff *= arb(fmpz(n - i))
                D = fmpz((n - m) * (m - r) + m * (m - 1) // 2 - e)
                tot += ff / arb(z) * (-(arb(D) * ln2)).exp()
            return tot
        C = arb(0)
        for k in range(1, K + 1): C += L_ball(k)
        def log_t(k):                                                    # ln t_k(n)
            kk = arb(fmpz(k))
            return 2 * kk * arb(pairs).log() + (3 * kk * kk + kk - 2 * kk * arb(N)) * ln2 - (2 * kk + 1).lgamma()
        k1 = K + 1; k2 = n // 3
        if k2 >= k1:
            mx = log_t(k1).max(log_t(k2)) if hasattr(arb, 'max') else None
            if mx is None:
                l1 = log_t(k1); l2 = log_t(k2); mx = l1 if (l1 - l2).mid() > 0 else l2
                mx = mx + abs(l1 - l2).rad()
            Bt = (arb(fmpz(k2 - k1 + 1)).log() + mx).exp()               # (number of terms) * max(endpoints), log-convexity
        else:
            Bt = arb(0)
        tail = ((arb(N) + 1).lgamma() - arb(fmpz(n * (n - 1) // 3)) * ln2).exp()
        B = Bt + tail
        lo = (logb + (1 + C).log()) / ln10
        hi = (logb + (1 + C + B).log()) / ln10
        x = lo.union(hi)                                                  # ball containing log10 a(n)
        E = x.floor().unique_fmpz()
        if E is None: return None, None, B
        mant = (ln10 * (x - arb(E))).exp()                                # in [1, 10)
        lead = (mant * arb(10)**digits).floor().unique_fmpz()
        return E, lead, B
    print("control: enclosure (K = 1) against the exact values, n = 20, 30, 50, 100, 150, 200, 250")
    for n in (20, 30, 50, 100, 150, 200, 250):
        E, lead, B = enclosure(n, 1, 30)
        s = str(a[n]); okE = (E is not None and int(E) == len(s) - 1)
        okL = "undetermined at this width" if lead is None else str(str(lead) == s[:31])
        print("  n=%d exponent ok: %s, 31 leading digits ok: %s, relative width bound 10^%.1f" % (n, okE, okL, float(B.log().mid() / arb(10).log().mid()) if B.mid() > 0 else float('-inf')))
    print("certified enclosures (Arb ball arithmetic; K = number of exact levels kept; a(n) = d.ddd... x 10^E)")
    for n, K in ((300, 1), (1000, 0), (1000, 1), (10**4, 0), (10**5, 0), (10**6, 0), (10**9, 0), (10**12, 0), (10**100, 0), (10**1000, 0)):
        E, lead, B = enclosure(n, K, 30)
        name = str(n) if n < 10**13 else "10^%d" % (len(str(n)) - 1)
        if E is None or lead is None:
            print("  n=%s K=%d: precision insufficient" % (name, K)); continue
        Es = str(E); ls = str(lead)
        logB = (B.log() / arb(10).log())
        wid = logB.mid()
        Eshow = Es if len(Es) <= 40 else "%s...%s (%d digits)" % (Es[:20], Es[-8:], len(Es))
        print("  n=%s K=%d: a(n) = %s.%s... x 10^%s ; relative width of the enclosure < 10^(%s)" % (name, K, ls[0], ls[1:], Eshow, str(wid.floor().unique_fmpz() + 1) if wid.floor().unique_fmpz() is not None else "?"))

def cmd_switching():
    print("census of switching classes of tournaments, A049313 = N_odd + N_lev (THM-479, Babai-Cameron)")
    print("  N_odd(n) = sum over odd cycle types of 2^(o - l + 1)/z ;  claim: N_odd(n) n!/2^C(n-1,2) = sum_k P_k(n) 4^(-k(n-1))")
    def odd_parts(n, maxp):
        if n == 0: yield (); return
        for p in range(min(n, maxp), 0, -1):
            if p % 2 == 1:
                for rest in odd_parts(n - p, p): yield (p,) + rest
    okid = True; rows = []
    for n in range(1, 20):
        tot = Fraction(0)
        for lam in odd_parts(n, n):
            z = 1
            for a in set(lam): c = lam.count(a); z *= a**c * factorial(c)
            o = sum((a - 1)//2 for a in lam) + sum(gcd(lam[i], lam[j]) for i in range(len(lam)) for j in range(i+1, len(lam)))
            tot += Fraction(2**(o + 1), z * 2**len(lam))
        Nodd = tot; Nlev = A049313[n-1] - Nodd          # N_odd need not be an integer when n is even
        if n >= 3:
            lhs = Fraction(Nodd * factorial(n), 2**comb(n-1, 2))
            rhs = sum(level_sum(k, n, shift=1) for k in range(0, (n-1)//2 + 1))
            okid &= (lhs == rhs)
        first = 1 + float(Fraction(64, 3) * falling(n, 3) / 4**n) if n >= 3 else 1.0
        rows.append((n, A049313[n-1], Nodd, Nlev, first))
    for (n, A, Nodd, Nlev, first) in rows:
        b = Fraction(2**comb(n-1, 2), factorial(n))
        print("  n=%2d A049313=%d N_odd=%s N_lev=%s  0<=N_lev<=2^(n^2/4): %s  A049313/b_sw=%.9f  1+(64/3)(n)_3 4^-n=%.9f" % (n, A, str(Nodd), str(Nlev), 0 <= Nlev <= 2**(n*n//4), float(A/b), first))
    print("clock-shift identity exact for 3 <= n <= 19:", okid)
    print("N_lev = 0 at every odd n <= 19:", all(r[3] == 0 for r in rows if r[0] % 2 == 1), "; N_lev >= 0 always:", all(r[3] >= 0 for r in rows))
    viol = 0
    for n in range(8, 20):
        Nodd = rows[n-1][2]; f = Fraction(Nodd * factorial(n), 2**comb(n-1, 2)); part = Fraction(0)
        for K in range(0, 4):
            part += level_sum(K, n, shift=1)
            if 2*(K+1) + 1 > n: break
            if f - part > B_bound(K, n, shift=1): viol += 1
    print("violations of the shifted bound for 8 <= n <= 19, K <= 3:", viol)

def cmd_strong():
    a = exact_values()
    # strong tournaments: A(x) = 1/(1 - S(x))
    N = 60; s = [0]*(N+1)
    for n in range(1, N+1): s[n] = a[n] - sum(s[k]*a[n-k] for k in range(1, n))
    print("strong tournaments s(n) (A051337) from A(x) = 1/(1 - S(x)); first terms:", s[1:9])
    print("  n   1 - s(n)/a(n)        n 2^(2-n)         ratio")
    for n in (10, 15, 20, 30, 40, 60):
        x = 1 - Fraction(s[n], a[n]); y = Fraction(4*n, 2**n)
        print("  %2d  %.12e  %.12e  %.9f" % (n, float(x), float(y), float(x/y)))

def cmd_heights():
    import mpmath as mp
    mp.mp.dps = 40
    print("integer-height (discrete Gaussian) moments: E[2^(M J)] for P(J=j) ~ 2^(-j^2/(2 phi)), against 2^(phi M^2/2)")
    worst = mp.mpf(0)
    for phi in (2, 4, 6, 8, 12):
        Zn = mp.nsum(lambda j: mp.mpf(2)**(-j*j/(2*mp.mpf(phi))), [-mp.inf, mp.inf])
        for M in range(-3, 9):
            mom = mp.nsum(lambda j: mp.mpf(2)**(M*j - j*j/(2*mp.mpf(phi))), [-mp.inf, mp.inf]) / Zn
            worst = max(worst, abs(mom / mp.mpf(2)**(mp.mpf(phi)*M*M/2) - 1))
    print("  largest relative deviation over phi in {2,4,6,8,12}, M = -3..8:", mp.nstr(worst, 5))
    Zt = mp.nsum(lambda j: mp.mpf(2)**(-j*(j+1)/mp.mpf(2)), [-mp.inf, mp.inf])
    worst = mp.mpf(0)
    for M in range(0, 13):
        mom = mp.nsum(lambda j: mp.mpf(2)**(M*j - j*(j+1)/mp.mpf(2)), [-mp.inf, mp.inf]) / Zt
        worst = max(worst, abs(mom / mp.mpf(2)**comb(M, 2) - 1))
    print("  triangular heights P(J=j) ~ 2^(-j(j+1)/2): E[2^(M J)] = 2^C(M,2) (labelled tournaments), M = 0..12; largest relative deviation:", mp.nstr(worst, 5))
    # a(n) from the integer-height field: E prod_p W_p^(m_p) = prod_d E[2^(M_d J_d)], M_d = #parts divisible by d
    def odd_parts(n, maxp):
        if n == 0: yield (); return
        for p in range(min(n, maxp), 0, -1):
            if p % 2 == 1:
                for rest in odd_parts(n - p, p): yield (p,) + rest
    def phi_(d): return sum(1 for i in range(1, d+1) if gcd(i, d) == 1)
    ok = True; a = exact_values()
    for n in range(1, 16):
        tot = Fraction(0)
        for lam in odd_parts(n, n):
            z = 1
            for q in set(lam): c = lam.count(q); z *= q**c * factorial(c)
            expo = Fraction(0)
            for d in range(1, n + 1, 2):
                Md = sum(1 for q in lam if q % d == 0)
                expo += Fraction(comb(Md, 2)) if d == 1 else Fraction(phi_(d) * Md * Md, 2)
            assert expo.denominator == 1
            tot += Fraction(2**expo.numerator, z)
        ok &= (tot == a[n])
    print("  a(n) = E[ Z_n(2^(H_p)) ] with H_p = sum_{d | p} J_d (exact moments), n <= 15:", ok)

if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("check", "all"): cmd_check()
    if what in ("bound", "all"): cmd_bound()
    if what in ("switching", "all"): cmd_switching()
    if what in ("strong", "all"): cmd_strong()
    if what in ("heights", "all"): cmd_heights()
    if what in ("giant", "all"): cmd_giant()
