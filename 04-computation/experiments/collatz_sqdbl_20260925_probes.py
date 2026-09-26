#!/usr/bin/env python3
"""collatz_sqdbl_20260925_probes.py -- FINITE-EXACT probes for the squares/doubles
(additive <-> multiplicative) transport foundry, session collatz-squares-doubles-20260925 (opus).

Seed (owner): "multiplication relates to the squares the way addition relates to
the doubles", i.e. the exponential transport n -> q^n sends doubles to squares,
halving to square root, and 3n+1 to q*X^3.

Probes (each with a positive control and a hostile control where one exists):
  P1  multiplicative Collatz M(X) = sqrt(X) if X is a perfect square, else rad(X)*X^3,
      run on exponent vectors: convergence iff X = m^e with m squarefree (Theorem 2).
  P2  Fermat-prime realisation: for p = 17, 257, 65537 the parity graph mod (p-1)
      is the +-sqrt graph on F_p^* under the discrete log (Lemma 1).
  P3  quadratic characters: v_2(3n+s) >= 2 iff s = chi_{-4}(n); >= 3 iff also chi_8(n) = -1;
      level-2 strategy cube = {chi_{-4}, -chi_{-4}, 3n-1, Collatz} (Lemma 3).
  P4  Jacobi law along Syracuse: (3/m_i) = (-1)^{v_i + [v_{i+1}=1]} on both sheets;
      (5/m_i) = (-1)^{v_i} for 5n+-1 (Lemma 4).
  P5  squares in the orbit: (2t+1)^2 -> 3t(t+1)+1; successor square iff Pell (Lemma 5).
  P6  real Bernstein value vs orbit product: R_L(d) = n(1 - prod(1 - 1/(3m_l))) on the
      minus sheet, R_L <= n; growth constant c = lim m_L 2^{d_L}/3^L (Proposition 6),
      with 5n+-1 controls showing full-rate divergence (c > 0 numerically).
Universe, filters and bounds are printed with each probe.
"""
import sys, math
from fractions import Fraction


def v2(x):
    return (x & -x).bit_length() - 1


def jacobi(a, n):
    """Jacobi symbol (a/n), n positive odd."""
    assert n > 0 and n % 2 == 1
    a %= n
    result = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def isqrt_exact(x):
    r = math.isqrt(x)
    return r if r * r == x else None


# ---------------------------------------------------------------- P1
def factor_small(x):
    f = {}
    d = 2
    while d * d <= x:
        while x % d == 0:
            f[d] = f.get(d, 0) + 1
            x //= d
        d += 1 if d == 2 else 2
    if x > 1:
        f[x] = f.get(x, 0) + 1
    return f


def mult_collatz_exponents(exps, sign=+1, cap=10**6, max_steps=100000):
    """Gated dynamics on exponent vectors: all even -> halve all; else e -> 3e+sign for e>0.
    Returns ('ones', steps) if the all-ones vector is reached, ('cycle', steps) if a
    state repeats, ('diverge', steps) if some exponent exceeds cap."""
    e = tuple(exps)
    seen = {e: 0}
    for step in range(1, max_steps + 1):
        if all(x == 1 for x in e):
            return ('ones', step - 1)
        if all(x % 2 == 0 for x in e):
            e = tuple(x // 2 for x in e)
        else:
            e = tuple(3 * x + sign for x in e)
        if max(e) > cap:
            return ('diverge', step)
        if e in seen:
            return ('cycle', step)
        seen[e] = step
    return ('timeout', max_steps)


def probe_p1(N=100000):
    print("== P1 multiplicative Collatz M(X) = sqrt(X) [square] else rad(X)*X^3, X <= %d ==" % N)
    print("   run on exponent vectors (gate = all exponents even); sign=+1 is 3e+1, sign=-1 is 3e-1")
    for sign in (+1, -1):
        counts = {}
        bad = []
        for X in range(2, N + 1):
            f = factor_small(X)
            exps = list(f.values())
            res, steps = mult_collatz_exponents(exps, sign)
            equal = len(set(exps)) == 1
            key = (res, equal)
            counts[key] = counts.get(key, 0) + 1
            if not equal and res != 'diverge':
                bad.append((X, res))
            if equal:
                res1, _ = mult_collatz_exponents([exps[0]], sign)
                if res1 != res:
                    bad.append((X, res, res1))
        print("   sign=%+d: outcome x (all exponents equal?) counts: %s" % (sign, sorted(counts.items())))
        print("   sign=%+d: violations of Theorem 2 prediction: %d %s" % (sign, len(bad), bad[:5]))
    for X in (32, 18, 72, 6**3, 12):
        f = factor_small(X)
        print("   X=%d exps=%s -> %s" % (X, f, mult_collatz_exponents(list(f.values()), +1)))


# ---------------------------------------------------------------- P2
def parity_graph_edges(k):
    M = 1 << k
    E = set()
    for s in range(M):
        t = (s // 2) if s % 2 == 0 else ((3 * s + 1) // 2)
        t %= (M // 2)
        E.add((s, t))
        E.add((s, t + M // 2))
    return E


def sqrt_graph_edges(p, g):
    """+-sqrt graph on F_p^*: y -> +-sqrt(y) if y is a square, else y -> +-sqrt(g*y^3)."""
    k = (p - 1).bit_length() - 1
    assert (1 << k) == p - 1
    dlog = {}
    x = 1
    for e in range(p - 1):
        dlog[x] = e
        x = x * g % p
    assert len(dlog) == p - 1, "g is not a primitive root"
    roots = {}
    for r in range(1, p):
        roots.setdefault(r * r % p, []).append(r)
    E = set()
    for y in range(1, p):
        if dlog[y] % 2 == 0:
            tgt = y
        else:
            tgt = g * pow(y, 3, p) % p
        rs = roots.get(tgt)
        assert rs is not None and len(rs) == 2, (y, tgt)
        for r in rs:
            E.add((dlog[y], dlog[r]))
    return E


def probe_p2():
    print("== P2 Fermat-prime realisation of the parity graph (Lemma 1) ==")
    for p, g in ((17, 3), (257, 3), (65537, 3)):
        k = (p - 1).bit_length() - 1
        E1 = parity_graph_edges(k)
        E2 = sqrt_graph_edges(p, g)
        same = (E1 == E2)
        print("   p=%d, k=%d, g=%d: parity graph mod 2^%d == +-sqrt graph on F_p^* under dlog_g: %s (%d edges)"
              % (p, k, g, k, same, len(E1)))
        print("      loop at 0 (y=1): %s ; loop at -1 (y=g^-1): %s" % ((0, 0) in E2, (p - 2, p - 2) in E2))
    p, g = 17, 2
    try:
        sqrt_graph_edges(p, g)
        print("   hostile g=2 mod 17: unexpectedly built")
    except AssertionError:
        print("   hostile control: g=2 mod 17 is not a primitive root -> dlog table fails, as it must")


# ---------------------------------------------------------------- P3
def chi_m4(n):
    return 1 if n % 4 == 1 else -1


def chi_8(n):
    return 1 if n % 8 in (1, 7) else -1


def probe_p3(N=1 << 20):
    print("== P3 quadratic characters and the first two bits of v_2(3n+s), odd n < %d ==" % N)
    ok2 = ok3 = True
    for n in range(1, N, 2):
        for s in (1, -1):
            v = v2(3 * n + s)
            if (v >= 2) != (s == chi_m4(n)):
                ok2 = False
            if (v >= 3) != (s == chi_m4(n) and chi_8(n) == -1):
                ok3 = False
    print("   v_2(3n+s) >= 2  <=>  s = chi_{-4}(n): %s" % ok2)
    print("   v_2(3n+s) >= 3  <=>  s = chi_{-4}(n) and chi_8(n) = -1: %s" % ok3)
    print("   level-2 strategies (sigma(1 mod 4), sigma(3 mod 4)) -> behaviour of odd n < 10^5 (steps <= 3000):")
    for s1 in (1, -1):
        for s3 in (1, -1):
            sig = {1: s1, 3: s3}
            desc = 0
            cyc = set()
            div = 0
            for n in range(3, 100001, 2):
                m = n
                steps = 0
                seen = set()
                broke = False
                while steps < 3000:
                    m = m // 2 if m % 2 == 0 else (3 * m + sig[m % 4]) // 2
                    steps += 1
                    if m < n:
                        desc += 1
                        broke = True
                        break
                    if m in seen:
                        cyc.add(min(seen))
                        broke = True
                        break
                    seen.add(m)
                    if m > 10**30:
                        div += 1
                        broke = True
                        break
                if not broke:
                    div += 1
            name = {(1, 1): 'Collatz', (-1, -1): '3n-1', (1, -1): 'chi_{-4}', (-1, 1): '-chi_{-4}'}[(s1, s3)]
            print("      %-9s descend-below-start: %d/49999, cycle minima met: %s, no descent within 3000 steps or > 10^30: %d"
                  % (name, desc, sorted(cyc)[:6], div))


# ---------------------------------------------------------------- P4
def syracuse_orbit(n, q, b, L):
    """odd iterates m_0=n, m_{i+1} = (q m_i + b)/2^{v_{i+1}}; returns (ms, vs)."""
    ms = [n]
    vs = [None]
    for _ in range(L):
        x = q * ms[-1] + b
        v = v2(x)
        ms.append(x >> v)
        vs.append(v)
    return ms, vs


def probe_p4(N=200000, L=12):
    print("== P4 Jacobi law along Syracuse orbits, odd n < %d, %d odd steps ==" % (N, L))
    for q, b in ((3, 1), (3, -1), (5, 1), (5, -1)):
        good = bad = 0
        wrong_law_bad = 0
        for n in range(1, N, 2):
            ms, vs = syracuse_orbit(n, q, b, L)
            for i in range(1, L):
                m = ms[i]
                if m % q == 0 or m <= 0:
                    continue
                lhs = jacobi(q, m)
                if q == 3:
                    pred = (-1) ** (vs[i] + (1 if vs[i + 1] == 1 else 0))
                    wrong = (-1) ** (vs[i])
                else:
                    pred = (-1) ** (vs[i])
                    wrong = (-1) ** (vs[i] + (1 if vs[i + 1] == 1 else 0))
                if lhs == pred:
                    good += 1
                else:
                    bad += 1
                if lhs != wrong:
                    wrong_law_bad += 1
        print("   %dn%+d: law holds %d, fails %d; the other law (control) fails %d" % (q, b, good, bad, wrong_law_bad))


# ---------------------------------------------------------------- P5
def pell_y(nmax):
    """(x, y) solutions of x^2 - 3y^2 = 1 with x, y > 0, in order, y <= nmax."""
    x, y = 2, 1
    out = []
    while y <= nmax:
        out.append((x, y))
        x, y = 2 * x + 3 * y, x + 2 * y
    return out


def probe_p5(S=2 * 10**6):
    print("== P5 odd squares in the Syracuse orbit (Lemma 5), odd s <= %d ==" % S)
    hits = []
    for s in range(1, S + 1, 2):
        t = (s - 1) // 2
        succ = 3 * t * (t + 1) + 1
        assert succ == (3 * s * s + 1) // 4 and (3 * s * s + 1) % 4 == 0
        u = isqrt_exact(succ)
        if u is not None:
            hits.append((s, u))
    pell = pell_y(S)
    odd_index_y = [y for j, (x, y) in enumerate(pell, start=1) if j % 2 == 1]
    print("   s with (3s^2+1)/4 a square: %s" % hits)
    print("   Pell y_j (j odd), x^2-3y^2=1: %s -> match: %s" % (odd_index_y, [s for s, u in hits] == odd_index_y))
    print("   length-3 chains (u again a hit): %s" % [(s, u) for s, u in hits if u in dict(hits)])
    N = 200000
    L = 8
    ok2 = ok3 = True
    for n in range(1, N, 2):
        ms, vs = syracuse_orbit(n, 3, 1, L)
        for i in range(1, L):
            m = ms[i]
            if (m % 8 == 1) != (vs[i + 1] == 2):
                ok2 = False
            if (m % 3 == 1) != (vs[i] % 2 == 0):
                ok3 = False
    print("   odd iterate m is a 2-adic square (m = 1 mod 8) iff v_next = 2: %s" % ok2)
    print("   odd iterate m is a 3-adic square unit (m = 1 mod 3) iff v_prev even: %s" % ok3)


# ---------------------------------------------------------------- P6
def orbit_until_cycle(n, q, b, max_odd=100000, cap=10**60):
    """odd iterates until a repeat (cycle), a cap, a nonpositive value, or a timeout."""
    ms = [n]
    vs = [0]
    seen = {n: 0}
    for i in range(max_odd):
        x = q * ms[-1] + b
        if x <= 0:
            return ms, vs, 'nonpositive'
        v = v2(x)
        m = x >> v
        ms.append(m)
        vs.append(v)
        if m in seen:
            return ms, vs, 'cycle'
        seen[m] = len(ms) - 1
        if m > cap:
            return ms, vs, 'cap'
    return ms, vs, 'timeout'


def probe_p6():
    print("== P6 real Bernstein value R(d) = sum 2^{d_l}/q^{l+1} versus the orbit product (Proposition 6) ==")
    N = 2000
    ident_ok = True
    ineq_ok = True
    cyc_R_eq_n = 0
    details = []
    for n in range(1, N + 1, 2):
        ms, vs, flag = orbit_until_cycle(n, 3, -1)
        d = 0
        R = Fraction(0)
        P = Fraction(1)
        for l in range(len(ms) - 1):
            R += Fraction(2 ** d, 3 ** (l + 1))
            P *= (1 - Fraction(1, 3 * ms[l]))
            d += vs[l + 1]
            if R != n * (1 - P):
                ident_ok = False
            if R > n:
                ineq_ok = False
        if flag == 'cycle':
            j0 = ms.index(ms[-1])
            Lc = len(ms) - 1 - j0
            Kc = sum(vs[j0 + 1:])
            dpre = sum(vs[1:j0 + 1])
            Rpre = sum(Fraction(2 ** sum(vs[1:l + 1]), 3 ** (l + 1)) for l in range(j0))
            e = 0
            inner = Fraction(0)
            for i in range(Lc):
                inner += Fraction(2 ** e, 3 ** (i + 1))
                e += vs[j0 + 1 + i]
            Rinf = Rpre + Fraction(2 ** dpre, 3 ** j0) * inner / (1 - Fraction(2 ** Kc, 3 ** Lc))
            if Rinf == n:
                cyc_R_eq_n += 1
            else:
                details.append((n, Rinf))
    print("   3n-1 on positive odd n <= %d: identity R_L = n(1 - prod_{l<L}(1 - 1/(3 m_l))) exact: %s; R_L <= n: %s"
          % (N, ident_ok, ineq_ok))
    print("   orbits entering a cycle: closed-form R(d) == n for %d of %d (mismatches: %s)" % (cyc_R_eq_n, (N + 1) // 2, details[:3]))
    ident_ok = True
    for n in range(1, N + 1, 2):
        ms, vs, flag = orbit_until_cycle(n, 3, 1)
        d = 0
        R = Fraction(0)
        P = Fraction(1)
        for l in range(len(ms) - 1):
            R += Fraction(2 ** d, 3 ** (l + 1))
            P *= (1 + Fraction(1, 3 * ms[l]))
            d += vs[l + 1]
            if R != n * (P - 1):
                ident_ok = False
    print("   3n+1 on positive odd n <= %d: identity R_L = n(prod_{l<L}(1 + 1/(3 m_l)) - 1) exact: %s" % (N, ident_ok))
    for q, b in ((5, 1), (5, -1)):
        shown = 0
        for n in range(1, 60, 2):
            ms, vs, flag = orbit_until_cycle(n, q, b, max_odd=4000, cap=10**400)
            if flag != 'cap':
                continue
            cs = []
            for L in (10, 100, 1000, len(ms) - 1):
                if L >= len(ms):
                    break
                dL = sum(vs[1:L + 1])
                cs.append((L, float(Fraction(ms[L] * 2 ** dL, q ** L))))
            recip = float(sum(Fraction(1, m) for m in ms[:len(ms) - 1]))
            print("   %dn%+d, n=%d: orbit exceeds 10^400 after %d odd steps; c_L = m_L 2^{d_L}/%d^L at L=%s ; sum 1/m_l = %.6f"
                  % (q, b, n, len(ms) - 1, q, cs, recip))
            shown += 1
            if shown >= 2:
                break
    print("   (c_L is monotone: increasing on the plus sheet, decreasing and >= 0 on the minus sheet; a positive limit = full-rate divergence)")


if __name__ == '__main__':
    which = sys.argv[1:] or ['p1', 'p2', 'p3', 'p4', 'p5', 'p6']
    for w in which:
        globals()['probe_' + w]()
        sys.stdout.flush()
