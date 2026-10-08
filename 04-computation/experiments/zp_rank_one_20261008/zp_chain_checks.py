"""Checks for THM-4610 (rank-one Haar coalescence of the base-p Collatz maps on Z_p), session opus-2026-10-08-S22.
Note: 05-knowledge/results/zp_rank_one_coalescence_20261008.md.

C_p(x) = x/p if p | x, else ((p+1)x + p - (x mod p))/p.  Pair chain u = P^k v + e, P = p+1, integer state
A = P^max(0,-k) e, digit i of v, digit j = (i + A) mod p of u (THM-4610's table).
  (A) the integer table against exact Fractions and against direct big-integer orbits (absorption = equality).
  (B) local lemmas, exhaustive on boxes of states x digits: factor law |f'| <= F|f| + (p-1)/p with the reference digit,
      direction (toward iff reference digit nonzero), flip law (1/p, 1/p, (p-2)/p when p does not divide A; no flip when it
      does), runs at k = 0 (exact multiplication), valuation runs at levels h >= 1 (w < mu_h: w-1 for every digit;
      w >= mu_h: exactly p-1 digits leave), mu_h = v_p(P^h - 1) = 1 + v_p(h) for odd p.
  (C) accessibility: the explicit descent path of statement 4 strictly decreases |e| for 2 <= |e| <= 20000, and the
      paths for e = +-1 reach (0,0); p = 3, 5, 7, 11, 13.
  (D) the level weight: kappa(theta) < 1 on (0,1) with kappa(1) = 1; least s with g(s) <= 1 at theta = 1/2 for p = 2..13
      (p = 2 reproduces THM-4581's s = 0.8966).
  (E) positive-integer cycles of C_p for starts <= 10^5 (p = 3, 5, 7, 11, 13).
  (F) dictionary facts for the other readings of Z_7, Z_11, Z_13 (note section 6): the 16 fixed points of T^4 for the
      Terras map T (3x+1) are c_w/(16 - 3^l), one per class mod 16; besides 0, -1 and the trivial cycle {1, 2} they form
      three 4-cycles with denominators 13, 7, -11 (l = 1, 2, 3), and {1, 2} is the integral member of the l = 2 family;
      the circulant tournaments on Z_7, Z_11, Z_13 with connection sets <2> = QR_7, <3> = QR_11 and <3> u 2<3> have affine
      automorphism groups of orders 21, 55, 39 (Z_p semidirect <2> or <3>).
Prints ALL CHECKS PASSED.  Runtime about 1 minute."""
import sys, math, random, time
from fractions import Fraction as Fr

FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def C(p, x):
    i = x % p
    return x // p if i == 0 else ((p + 1) * x + p - i) // p


def istep(p, k, A, i):
    P = p + 1
    j = (i + A) % p
    if i == 0 and j == 0:
        assert A % p == 0
        return k, A // p
    if k >= 0:
        if i and j:
            n = P * A + (p - j) - P ** k * (p - i)
            nk = k
        elif i == 0:
            n = P * A + p - j
            nk = k + 1
        else:
            n = (A - P ** (k - 1) * (p - i)) if k >= 1 else (P * A - (p - i))
            nk = k - 1
    else:
        h = -k
        if i and j:
            n = P * A + P ** h * (p - j) - (p - i)
            nk = k
        elif i == 0:
            n = A + P ** (h - 1) * (p - j)
            nk = k + 1
        else:
            n = P * A - (p - i)
            nk = k - 1
    assert n % p == 0
    return nk, n // p


def e_of(p, k, A):
    return Fr(A, (p + 1) ** max(0, -k))


def frac_step(p, k, e, i):
    """the chain from the defining relation u = P^k v + e, by exact rationals (independent of the table)."""
    P = p + 1
    M = Fr(P) ** k
    em = e.numerator * pow(e.denominator, -1, p) % p
    j = (i + em) % p
    if i == 0 and j == 0:
        return k, e / p
    if i and j:
        return k, (P * e + (p - j) - M * (p - i)) / p
    if i == 0:
        return k + 1, (P * e + (p - j)) / p
    return k - 1, (e - M * Fr(p - i, P)) / p


def vp(n, p):
    if n == 0:
        return 10 ** 9
    v = 0
    while n % p == 0:
        n //= p; v += 1
    return v


def section_A(rng):
    print('(A) the integer table against exact rationals and direct orbits')
    bad_t = 0; n_t = 0
    for p in (3, 5, 7, 11, 13):
        for _ in range(1500):
            k = rng.randint(-4, 4); A = rng.randint(-500, 500)
            e = e_of(p, k, A)
            for _ in range(40):
                i = rng.randrange(p)
                k1, e1 = frac_step(p, k, e, i)
                k2, A2 = istep(p, k, A, i)
                n_t += 1
                if k1 != k2 or e1 != e_of(p, k2, A2):
                    bad_t += 1; break
                k, e, A = k1, e1, A2
    check(bad_t == 0, f'table = exact rational chain on {n_t} random steps, p = 3, 5, 7, 11, 13')
    bad = 0; merges = 0; cases = 0
    for p in (3, 5, 7, 11, 13):
        W = 260
        for _ in range(240):
            k0 = rng.randint(-3, 3); A0 = rng.randint(-40, 40)
            if k0 == 0 and A0 == 0:
                A0 = 1
            h0 = max(0, -k0); P = p + 1
            # choose y so that u_0 = P^k0 y + e_0 is an integer: P^h0 | y + A0 when k0 < 0
            y = rng.randrange(p ** W)
            if k0 < 0:
                y = y - (y + A0) % (P ** h0)
                if y < 0:
                    y += P ** h0
            u = (P ** k0 * y + A0) if k0 >= 0 else (y + A0) // P ** h0
            assert (Fr(u) == Fr(P) ** k0 * y + e_of(p, k0, A0))
            k, A, v = k0, A0, y
            cases += 1
            for n in range(W - 20):
                if Fr(u) != Fr(P) ** k * v + e_of(p, k, A):
                    bad += 1; break
                if k == 0 and A == 0:
                    merges += 1
                    if u != v:
                        bad += 1
                    break
                if u == v:                      # an equality without absorption would contradict the relation
                    bad += 1; break
                k, A = istep(p, k, A, v % p)
                u, v = C(p, u), C(p, v)
    check(bad == 0, f'{cases} random starts on direct big-integer orbits (240 digits): the relation holds at every step, '
                    f'absorption is equality ({merges} merges), and no equality occurs before absorption')


def section_B():
    print('(B) local lemmas, exhaustive on boxes')
    for p in (2, 3, 5, 7, 11, 13):
        P = p + 1; c = Fr(p - 1, p)
        bad = 0; worst = Fr(0); n = 0
        K = 5; Amax = 150 if p <= 7 else 90
        for k in range(-K, K + 1):
            for A in range(-Amax, Amax + 1):
                if k == 0 and A == 0:
                    continue
                f = Fr(A, P ** k) if k >= 0 else Fr(A, P ** (-k))
                flips = {'up': 0, 'down': 0, 'stay': 0}
                for i in range(p):
                    j = (i + A) % p
                    k2, A2 = istep(p, k, A, i)
                    f2 = Fr(A2, P ** k2) if k2 >= 0 else Fr(A2, P ** (-k2))
                    if k == 0:
                        F = Fr(1, p) if (k2 != 0 or i == 0) else Fr(P, p)
                        if k2 == 0 and i != j and (i == 0 or j == 0):
                            bad += 1
                    else:
                        ref = i if k > 0 else j
                        F = Fr(1, p) if ref == 0 else Fr(P, p)
                        if k2 != k and ((abs(k2) < abs(k)) != (ref != 0)):
                            bad += 1
                    if abs(f2) - F * abs(f) > c:
                        bad += 1
                    worst = max(worst, abs(f2) - F * abs(f))
                    if k2 > k:
                        flips['up'] += 1
                    elif k2 < k:
                        flips['down'] += 1
                    else:
                        flips['stay'] += 1
                    n += 1
                if A % p:
                    if p >= 2 and (flips['up'] != 1 or flips['down'] != 1):
                        bad += 1
                elif flips['up'] or flips['down']:
                    bad += 1
        # runs at k = 0 and valuation runs
        for h in range(0, 16):
            mu = vp(P ** h - 1, p) if h else None
            if h and p > 2 and mu != 1 + vp(h, p):
                bad += 1
            for sgn in ((1,) if h == 0 else (1, -1)):
                k = sgn * h
                for A in range(-1200, 1201):
                    if A == 0 or A % p:
                        continue
                    outs = [istep(p, k, A, i) for i in range(p)]
                    if h == 0:
                        if any(k2 != 0 or (A2 * p != A and A2 * p != P * A) for k2, A2 in outs):
                            bad += 1
                        continue
                    w = vp(A, p)
                    ws = [vp(A2, p) for _, A2 in outs]
                    if any(k2 != k for k2, _ in outs):
                        bad += 1
                    if w < mu:
                        if any(x != w - 1 for x in ws):
                            bad += 1
                    elif sum(1 for x in ws if x < mu) != p - 1:
                        bad += 1
        check(bad == 0 and worst == c,
              f'p = {p}: factor law, direction and flip law on {n} state-digit pairs (worst additive = (p-1)/p exactly); '
              f'runs at k = 0 and valuation runs (mu_h = v_p(P^h - 1){" = 1 + v_p(h)" if p > 2 else ""}) for h <= 15')


def descent(p, e):
    k, A = 0, e
    if e > 0:
        k, A = istep(p, k, A, (-A) % p)
        assert k == -1
        while A % p == 0:
            k, A = istep(p, k, A, 0)
            assert k == -1
        k, A = istep(p, k, A, 0)
    else:
        k, A = istep(p, k, A, 0)
        assert k == 1
        while A % p == 0:
            k, A = istep(p, k, A, 0)
            assert k == 1
        k, A = istep(p, k, A, (-A) % p)
    assert k == 0
    return A


def section_C():
    print('(C) accessibility by explicit descent')
    for p in (3, 5, 7, 11, 13):
        bad = 0; worst = 0.0
        for e in range(-20000, 20001):
            if abs(e) < 2 or e % p == 0:
                continue
            e2 = descent(p, e)
            if not (abs(e2) < abs(e) and (e2 > 0) == (e > 0) or e2 == 0):
                bad += 1
            bound_pos = (( (p + 1) * e - 1 + p * p - p) / (p * p)) if e > 0 else None
            if e > 0 and e2 > bound_pos + 1e-9:
                bad += 1
            if abs(e) >= 1000:
                worst = max(worst, abs(e2) / abs(e))
        s1 = istep(p, *istep(p, 0, 1, 0), p - 2)
        s2 = istep(p, *istep(p, 0, -1, 1), 0)
        # full descent to (0,0) for every |e| <= 3000 following the explicit rules
        full_bad = 0
        for e0 in range(-3000, 3001):
            e = e0; steps = 0
            while e != 0 and steps < 10000:
                if e % p == 0:
                    e //= p
                elif e == 1 or e == -1:
                    e = 0
                else:
                    e = descent(p, e)
                steps += 1
            if e != 0:
                full_bad += 1
        check(bad == 0 and s1 == (0, 0) and s2 == (0, 0) and full_bad == 0,
              f'p = {p}: the descent strictly decreases |e| for 2 <= |e| <= 20000 (positive bound (Pe - 1 + p^2 - p)/p^2 holds; '
              f'worst ratio {worst:.3f} for |e| >= 1000); (0,1) -> (1,2) -> (0,0) and (0,-1) -> (-1,-2) -> (0,0); every |e| <= 3000 descends to 0')


def section_D():
    print('(D) the level weight')
    ok = True; rows = []
    for p in (2, 3, 5, 7, 11, 13):
        P = p + 1
        kap = lambda t: (1 / p) * p ** (-t) + ((p - 1) / p) * (P / p) ** t
        conv = all(kap(t / 100) < 1 for t in range(1, 100)) and abs(kap(1) - 1) < 1e-12 and abs(kap(0) - 1) < 1e-12
        th = 0.5
        g = lambda s: (1 / p) * p ** (-th) * s + (1 / p) * (P / p) ** th / s + ((p - 2) / p) * (P / p) ** th
        lo, hi = 0.01, 1.0          # g is convex in s; find the least root of g(s) = 1 below 1 by bisection
        if g(lo) <= 1:
            smin = lo
        else:
            for _ in range(100):
                mid = (lo + hi) / 2
                if g(mid) <= 1:
                    hi = mid
                else:
                    lo = mid
            smin = hi
        lam0 = max(kap(th), (2 / p) * p ** (-th) * 0.95 + ((p - 2) / p) * (P / p) ** th)
        rows.append((p, round(kap(th), 4), round(smin, 4), round(lam0, 4)))
        ok &= conv and smin < 1 and g(0.95) <= 1 + 1e-12 and lam0 < 1
    print('     (p, kappa(1/2), least s with g(s) <= 1, lambda_0 at s = 0.95):', rows)
    check(ok and abs(rows[0][2] - 0.8966) < 1e-3, 'kappa < 1 on (0,1) with kappa(0) = kappa(1) = 1; a level weight s < 1 exists; p = 2 gives s = 0.8966 (THM-4581)')


def section_E():
    print('(E) positive-integer cycles of C_p (starts <= 10^5)')
    res = {}
    for p in (3, 5, 7, 11, 13):
        term = {}                      # value -> min element of the cycle it ends in (memo for values < 10^7)
        cyc = {}
        for x0 in range(1, 100001):
            if x0 in term:
                continue
            path = []; pos = {}
            x = x0
            while x not in term and x not in pos:
                pos[x] = len(path); path.append(x); x = C(p, x)
            if x in term:
                m = term[x]
            else:                      # new cycle: path[pos[x]:]
                c = path[pos[x]:]
                m = min(c); cyc[m] = len(c)
            for z in path:
                if z < 10 ** 7:
                    term[z] = m
        res[p] = dict(sorted(cyc.items()))
    print('     cycles (min element: length):', res)
    check(res[7] == {1: 7} and res[13] == {1: 13} and res[11].get(642) == 57 and res[3].get(7) == 9,
          'p = 7, 13: only the trivial cycle {1..p}; p = 11: also a 57-cycle with minimum 642; p = 3: also a 9-cycle with minimum 7')


def section_F():
    print('(F) dictionary: period-4 points of 3x+1 and the Collatz-multiplier circulants')
    def T(x):
        return x / 2 if (x.numerator % 2 == 0) else (3 * x + 1) / 2
    pts = []
    for r in range(16):
        w = []; z = r
        for _ in range(4):
            w.append(z & 1); z = ((3 * z + 1) >> 1) if z & 1 else (z >> 1)
        l = sum(w)
        c = sum(3 ** (l - 1 - idx) * 2 ** t for idx, t in enumerate([t for t in range(4) if w[t]]))
        x = Fr(c, 16 - 3 ** l)
        # x must be a fixed point of T^4 (2-adic parity of a rational with odd denominator = parity of numerator*den^-1)
        y = x
        for _ in range(4):
            y = y / 2 if (y.numerator * y.denominator) % 2 == 0 else (3 * y + 1) / 2
        pts.append((r, l, x, y == x))
    dens = sorted(set((l, (16 - 3 ** l)) for _, l, _, _ in pts))
    fam7 = sorted(x for _, l, x, _ in pts if l == 2)
    ok = all(f for *_, f in pts) and len(set(x for _, _, x, _ in pts)) == 16
    ok &= [d for _, d in dens] == [15, 13, 7, -11, -65]
    ok &= Fr(1) in fam7 and Fr(2) in fam7 and sum(1 for x in fam7 if x.denominator == 7) == 4
    ok &= sorted(x for _, l, x, _ in pts if l == 1) == [Fr(1, 13), Fr(2, 13), Fr(4, 13), Fr(8, 13)]
    ok &= sorted(x for _, l, x, _ in pts if l == 3) == [Fr(-38, 11), Fr(-29, 11), Fr(-23, 11), Fr(-19, 11)]
    print(f'     T^4 fixed points by l: ' + '; '.join(f'l={l}: ' + ', '.join(str(x) for _, ll, x, _ in pts if ll == l) for l in range(5)))
    orders = {}
    for q, S in ((7, {1, 2, 4}), (11, {1, 3, 4, 5, 9}), (13, {1, 3, 9, 2, 6, 5})):
        tourn = all(((s_ in S) != ((-s_) % q in S)) for s_ in range(1, q))
        aut = sum(1 for a in range(1, q) for b in range(q) if {(a * s_) % q for s_ in S} == S)
        orders[q] = (tourn, aut)
    ok &= orders == {7: (True, 21), 11: (True, 55), 13: (True, 39)}
    print(f'     circulant tournaments (is a tournament, affine automorphisms): {orders}')
    check(ok, 'T^4 fixed points c_w/(16 - 3^l): denominators 15, 13, 7, -11, -65; {1, 2} lies in the l = 2 family with the 5/7 cycle; '
              'the 1/13 and -19/11 cycles; affine automorphism groups of orders 21, 55, 39')


if __name__ == '__main__':
    t0 = time.time()
    rng = random.Random(20261008)
    section_A(rng)
    section_B()
    section_C()
    section_D()
    section_E()
    section_F()
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
