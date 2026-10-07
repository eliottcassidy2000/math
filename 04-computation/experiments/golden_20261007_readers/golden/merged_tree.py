#!/usr/bin/env python3
"""The merged tree 332 -> 233 -> 425 <- 223, the number 105, pair-chain data, and base-rate tests."""
import random
from fractions import Fraction

def std_orbit(n, stop=1):
    out = [n]
    while n != stop:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        out.append(n)
    return out

def terras(x):
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2

def terras_orbit(n, steps):
    out = [n]
    for _ in range(steps):
        n = terras(n); out.append(n)
    return out

def odd_word(a, b):
    """valuation word of the odd map U from odd a to odd b (b on a's odd orbit)."""
    w = []; x = a
    while x != b:
        y = 3 * x + 1; k = 0
        while y % 2 == 0:
            y //= 2; k += 1
        w.append(k); x = y
        if len(w) > 500:
            return None
    return w

FIB = set(); a, b = 1, 2
while a < 10 ** 7:
    FIB.add(a); a, b = b, a + b
LUC = set(); a, b = 2, 1
while a < 10 ** 7:
    LUC.add(a); a, b = b, a + b
LUC.add(1)

def fib_index(n):
    a, b, k = 0, 1, 0
    while a < n:
        a, b, k = b, a + b, k + 1
    return k if a == n else None

def luc_index(n):
    a, b, k = 2, 1, 0
    while a < n and k < 100:
        a, b, k = b, a + b, k + 1
    return k if a == n else None

def tag(n):
    t = []
    if n in FIB and n > 3: t.append(f'F_{fib_index(n)}')
    if n in LUC and n > 3: t.append(f'L_{luc_index(n)}')
    if n == 121: t.append('11^2')
    if n == 242: t.append('3^5-1')
    if n == 243: t.append('3^5')
    return '/'.join(t)

def pair_chain_check(u, v, T):
    """Run THM-4581's chain from the exact relation u = 3^k v + e (k=0 start) for T Terras steps using
    the actual orbits; return list of (k, e) and check consistency with the table."""
    k, e = 0, Fraction(u - v)
    hist = [(k, e)]
    uu, vv = u, v
    for _ in range(T):
        sigma = int(e) % 2 if e.denominator == 1 else int((e * 3 ** max(0, -k)).numerator) % 2  # parity of e (2-adic unit denom)
        # parity of e as a 2-adic number with odd denominator
        num, den = e.numerator, e.denominator
        sigma = num % 2  # den odd
        beta = vv % 2
        if sigma == 0 and beta == 0:
            e = e / 2
        elif sigma == 0 and beta == 1:
            e = (3 * e + 1 - Fraction(3) ** k) / 2
        elif sigma == 1 and beta == 0:
            k, e = k + 1, (3 * e + 1) / 2
        else:
            k, e = k - 1, (e - Fraction(3) ** (k - 1)) / 2
        uu, vv = terras(uu), terras(vv)
        assert Fraction(uu) == Fraction(3) ** k * vv + e, (u, v, uu, vv, k, e)
        hist.append((k, e))
        if k == 0 and e == 0:
            break
    return hist

def main():
    print('== Standard orbits (3n+1, n/2), golden-tagged values')
    for n in (332, 233, 223, 425, 105, 27):
        o = std_orbit(n)
        tags = [(i, x, tag(x)) for i, x in enumerate(o) if tag(x)]
        print(f'{n}: length {len(o)-1}, peak {max(o)}; tagged: ' + ', '.join(f'{x}({t})@{i}' for i, x, t in tags))
    o27 = std_orbit(27); s27 = set(o27)
    for n in (332, 233, 223, 425, 105):
        o = std_orbit(n)
        j = next((i for i, x in enumerate(o) if x in s27), None)
        print(f'   {n} joins the 27-trajectory at {o[j]} after {j} steps' if o[j] != 1 and o[j] in s27 and o[j] > 8 else f'   {n} does not join 27 above 8 (first common {o[j]})')
    print('\n== Odd-map words')
    for a, b in ((83, 233), (233, 425), (223, 425), (105, 1), (233, 377), (377, 425)):
        w = odd_word(a, b)
        print(f'   {a} -> {b}: word {w} (r={len(w)}, A={sum(w)})')
    # 332 = 4*83
    print('\n== Terras times / odd counts to the meeting points')
    def tt(n, target):
        t = 0; r = 0; x = n
        while x != target:
            r += x % 2; x = terras(x); t += 1
            if t > 10000: return None
        return t, r
    for n, m in ((332, 233), (233, 425), (223, 425), (332, 425)):
        print(f'   T^t({n}) = {m}: (t, odd steps) =', tt(n, m))
    # equal-time relations at the meetings
    print('\n== Pair-chain data at the meetings (THM-4581 language)')
    for a, b, m in ((223, 233, 425), (332, 233, 233), (233, 105, None)):
        if m is None:
            continue
        ta, ra = tt(a, m); tb, rb = tt(b, m)
        if ta <= tb:
            u, d = a, tb - ta
            v = terras_orbit(b, d)[-1]
            t = ta
            ru = ra; rv = rb - sum(x % 2 for x in terras_orbit(b, d)[:-1])
        else:
            v, d = b, ta - tb
            u = terras_orbit(a, d)[-1]
            t = tb
            ru = ra - sum(x % 2 for x in terras_orbit(a, d)[:-1]); rv = rb
        k0 = rv - ru
        e0 = Fraction(u) - Fraction(3) ** k0 * v
        print(f'   {a} vs {b} meet at {m}: Terras times {ta},{tb}; clock shift {d}; equal-time pair (u,v)=({u},{v}) '
              f'merges at Terras time {t}; k0 = {k0}, e0 = u - 3^k0 v = {e0}')
        # verify with chain from (k0, e0)
        k, e = k0, e0; uu, vv = u, v
        for s in range(t):
            sigma = e.numerator % 2; beta = vv % 2
            if sigma == 0 and beta == 0: e = e / 2
            elif sigma == 0 and beta == 1: e = (3 * e + 1 - Fraction(3) ** k) / 2
            elif sigma == 1 and beta == 0: k, e = k + 1, (3 * e + 1) / 2
            else: k, e = k - 1, (e - Fraction(3) ** (k - 1)) / 2
            uu, vv = terras(uu), terras(vv)
            assert Fraction(uu) == Fraction(3) ** k * vv + e
        print(f'      chain state at time {t}: (k,e) = ({k},{e}), u_t = v_t = {uu}')
    print('\n== 2-adic proximity (common Terras parity prefix length = v_2(a-b))')
    S = [27, 31, 41, 105, 223, 233, 332, 425]
    def v2(x):
        x = abs(x); k = 0
        while x % 2 == 0: x //= 2; k += 1
        return k
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            a, b = S[i], S[j]
            if v2(a - b) >= 5:
                ta = terras_orbit(a, v2(a - b)); tb = terras_orbit(b, v2(a - b))
                r = sum(x % 2 for x in ta[:-1])
                print(f'   v2({b}-{a}) = {v2(b-a)}; T^{v2(a-b)}: {ta[-1]} vs {tb[-1]}, difference {tb[-1]-ta[-1]} = 3^{r} * {(tb[-1]-ta[-1])//3**r}')
    print('   residues mod 64:', {n: n % 64 for n in S}, ' mod 128:', {n: n % 128 for n in S})

    # ---------- base-rate tests ----------
    print('\n== Base rates (standard orbits of n in [200,500])')
    special = {'233=F13': 233, '377=F14': 377, '47=L8': 47, '322=L12': 322, '121=11^2': 121, '242=3^5-1': 242,
               '425': 425, '9232 (27-trunk peak)': 9232}
    N0, N1 = 200, 500
    orbs = {n: set(std_orbit(n)) for n in range(N0, N1 + 1)}
    for name, val in special.items():
        frac = sum(1 for n in orbs if val in orbs[n]) / len(orbs)
        print(f'   P(orbit contains {name}) = {frac:.3f}')
    # count of 'golden' values >= 29 (Fibonacci >= 34, Lucas >= 29, 121, 242) per orbit, vs shifted controls
    def golden_count(o, shift=0):
        c = 0
        for x in o:
            y = x - shift
            if y >= 29 and (y in FIB or y in LUC or y in (121, 242)):
                c += 1
        return c
    import statistics
    for shift in (0, 1, -1, 2, -2, 3):
        counts = [golden_count(orbs[n], shift) for n in orbs]
        mean = statistics.mean(counts)
        print(f'   shift {shift:+d}: mean # of (golden+shift) values >=29 per orbit = {mean:.2f}; '
              f'332 has {golden_count(orbs[332], shift)}, 233 has {golden_count(orbs[233], shift)}; '
              f'fraction of orbits with >= 6: {sum(c >= 6 for c in counts)/len(counts):.3f}')
    # Fibonacci consecutive pair test: F_k and F_{k+1} with F_{k+1} on the orbit of F_k
    print('   consecutive Fibonacci pairs with F_{k+1} on orbit(F_k):',
          [(fib_a, fib_b) for fib_a, fib_b in [(34, 55), (55, 89), (89, 144), (144, 233), (233, 377), (377, 610),
           (610, 987), (987, 1597), (1597, 2584), (2584, 4181), (4181, 6765)] if fib_b in set(std_orbit(fib_a))])
    # base rate: for a in [150,400], b = round(phi*a): P(b in orbit(a))
    phi = (1 + 5 ** 0.5) / 2
    hits = sum(1 for a in range(150, 401) if round(phi * a) in set(std_orbit(a)))
    print(f'   control: P(round(phi a) in orbit(a)), a in [150,400] = {hits/251:.3f}')
    hits2 = sum(1 for a in range(150, 401) for c in (1.5, 1.6, 1.7) if round(c * a) in set(std_orbit(a))) / (251 * 3)
    print(f'   control: P(round(c a) in orbit(a)), c in {{1.5,1.6,1.7}} = {hits2:.3f}')
    # how many n <= 1000 have orbits through 9232
    thr = sum(1 for n in range(1, 1001) if 9232 in set(std_orbit(n)))
    print(f'   n <= 1000 whose orbits pass 9232 (the 27-trunk): {thr}')
    # 105: orbit features base rate for numbers ~100
    o105 = std_orbit(105)
    print('\n== 105: odd orbit', [x for x in o105 if x % 2], ' length', len(o105) - 1)
    orbs2 = {n: set(std_orbit(n)) for n in range(80, 131)}
    for val in (101, 19, 29, 11, 17, 13, 76):
        print(f'   P(orbit of n in [80,130] contains {val}) = {sum(1 for n in orbs2 if val in orbs2[n])/len(orbs2):.3f}')

if __name__ == '__main__':
    main()
