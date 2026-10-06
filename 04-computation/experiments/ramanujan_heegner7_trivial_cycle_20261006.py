#!/usr/bin/env python3
"""Ramanujan's constant, the split primes {1,2,4} mod 7, and the trivial Collatz cycle's codes.

Session opus-2026-10-06-S17.  Note: 05-knowledge/results/ramanujan_heegner7_trivial_cycle_20261006.md
Builds on collatz_paley_bridge_20261001.md (the trivial cycle's real code {1,2,4}/7 and 2-adic code -{1,2,4}/7).

Sections:
 1. Heegner table: j(tau_d), e^(pi sqrt d) vs -j + 744, splitting of 2 and 3, parity of j (also the 4 non-maximal orders).
 2. Gross-Zagier (1985) check: primes dividing j(tau_d) are non-split in Q(sqrt(-d)) and Q(sqrt(-3)).
 3. Q(sqrt(-7)) from the trivial cycle: decomposition group <2> = {1,2,4}; the Gauss period g = z + z^2 + z^4
    satisfies g^2 + g + 2 = 0, g + 1 = tau_7 (a Heegner point) and j(g) = -15^3; Paley spectrum; F_8 trace;
    2-adic roots of x^2 + x + 2 and the 2-adic code Phi(1) = -1/7 = (1 + 2g)^(-2).
 4. Visibility: a Heegner field lies in a decomposition field of 2 iff 2 splits in it, i.e. only for d = 7;
    class numbers of squarefree d = 7 mod 8; the norm lemma (smallest split prime = (d+1)/4).
 5. Gauss periods of the codes of every integer cycle, in both codings: exact degree, CM point or not; Propositions 4-5:
    census of all doubling orbits mod odd m < 3001 (CM periods, class number one, which are codes of T-cycles).
 6. Ordinary elliptic curves over F_2 (all 32 Weierstrass tuples) and E_7 = 49a1 over F_(2^k).
 7. E_7 over F_p: a_p = 0 iff p = 3,5,6 mod 7; 4p = a_p^2 + 7B^2 on split p.
 8. Klein quartic x^3y + y^3z + z^3x = 0: over F_p (cubic-twist formula) and over F_(2^n) (L(T) = 1 + 5T^3 + 8T^6).
 9. Supersingular congruences for j(tau_163) and j(tau_7).
10. E8 chain from {1,2,4}; gamma_2 = E4/eta^8 as a q-series (248, 744 = 3*248, 196884); gamma_2 at Heegner points.
11. Weber's f(sqrt(-7))^24 = 2^12: e^(pi sqrt 7) = 2^12 - 24 - 276 e^(-pi sqrt 7) - ...; j(sqrt(-7)) = 255^3.
12. Ramanujan's tau mod 7: tau(n) = n sigma_9(n) mod 7; tau(p) = 0 mod 7 iff p = 3,5,6 mod 7 (or 7).
13. Chudnovsky's constants: (640320^3 + 1728)/163 = m^2, 12 * 545140134 = 163 m; Gross-Zagier for j - 1728.
14. Landau-Ramanujan constant of Q(sqrt(-7)).
15. Numerology probes.

Reproduce: python3 04-computation/experiments/ramanujan_heegner7_trivial_cycle_20261006.py   (a few minutes)
"""
import math, cmath
from collections import Counter
import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 60
HEEGNER = [1, 2, 3, 7, 11, 19, 43, 67, 163]
JKNOWN = {1: 1728, 2: 8000, 3: 0, 7: -3375, 11: -32768, 19: -884736, 43: -884736000,
          67: -147197952000, 163: -262537412640768000}
FAIL = []


def check(label, ok):
    if not ok:
        FAIL.append(label)
    return ok


def tau_of(d):
    return mp.mpc(0, mp.sqrt(d)) if d in (1, 2) else mp.mpc(mp.mpf(1) / 2, mp.sqrt(d) / 2)


def kron(D, p):
    """Kronecker symbol (D/p) for a prime p."""
    if p == 2:
        if D % 2 == 0:
            return 0
        return 1 if D % 8 in (1, 7) else -1
    return sp.jacobi_symbol(D % p, p) if D % p else 0


def disc(d):
    return -4 * d if d in (1, 2) else -d


def series_mul(a, b, n=None):
    n = min(len(a), len(b)) if n is None else n
    return [sum(a[i] * b[k - i] for i in range(k + 1)) for k in range(n)]


def gamma2(tau, terms=None):
    """Weber gamma_2 = E4/eta^8 with eta = exp(2 pi i tau/24) prod(1 - q^n) (correct branch)."""
    q = mp.exp(2j * mp.pi * tau)
    E4 = 1 + 240 * mp.nsum(lambda n: n ** 3 * q ** n / (1 - q ** n), [1, mp.inf])
    eta = mp.exp(2j * mp.pi * tau / 24) * mp.nprod(lambda n: 1 - q ** n, [1, mp.inf])
    return E4 / eta ** 8


print('== 1. Heegner table')
for d in HEEGNER:
    t = tau_of(d)
    j = mp.kleinj(t) * 1728
    jr = int(mp.nint(mp.re(j)))
    check(f'j({d})', jr == JKNOWN[d])
    e = mp.e ** (mp.pi * mp.sqrt(d))
    D = disc(d)
    cr = round(abs(jr) ** (1 / 3)) * (1 if jr >= 0 else -1) if jr else 0
    print(f'   d={d:3d}  j={jr:>22d}  cube root {cr:>8d}'
          f'  e^(pi sqrt d) - (-j + 744) = {mp.nstr(e - (-jr + 744), 8):>14s}'
          f'  2: {["inert", "ramified", "split"][kron(D, 2) + 1]:8s} 3: {["inert", "ramified", "split"][kron(D, 3) + 1]:8s}'
          f'  j odd: {jr % 2 == 1}')
print('   (for d = 1, 2 the point is i sqrt d and q > 0, so the comparison there is with j - 744, not -j + 744)')
NONMAX = {-12: (mp.mpc(0, mp.sqrt(3)), 54000), -16: (mp.mpc(0, 2), 287496),
          -27: (mp.mpc(mp.mpf(1) / 2, 3 * mp.sqrt(3) / 2), -12288000), -28: (mp.mpc(0, mp.sqrt(7)), 16581375)}
for D, (t, jk) in NONMAX.items():
    jr = int(mp.nint(mp.re(mp.kleinj(t) * 1728)))
    check(f'j({D})', jr == jk)
    print(f'   non-maximal order D = {D}: j = {jr} = {sp.factorint(jr) if jr > 0 else "-" + str(sp.factorint(-jr))}, odd: {jr % 2 == 1}')
odd13 = sorted([disc(d) for d in HEEGNER if JKNOWN[d] % 2] + [D for D, (t, jk) in NONMAX.items() if jk % 2])
print('   among all 13 class-number-one discriminants, j is odd exactly for D =', odd13, '(the two orders of Q(sqrt(-7)))', check('odd13', odd13 == [-28, -7]))

print('\n== 2. Gross-Zagier: primes dividing j(tau_d) are non-split in Q(sqrt(-d)) and in Q(sqrt(-3))')
for d in HEEGNER:
    jr = JKNOWN[d]
    if jr == 0:
        print(f'   d={d}: j = 0 (tau is the cube root of unity point)')
        continue
    fac = sp.factorint(abs(jr))
    rows, ok = [], True
    for p in fac:
        ok &= kron(disc(d), p) != 1 and kron(-3, p) != 1
        rows.append(f'{p}^{fac[p]}')
    check(f'GZ {d}', ok)
    print(f'   d={d:3d}: |j| = {" * ".join(rows):28s} all non-split in both fields: {ok};  largest prime {max(fac)} <= 3|D|/4 = {3 * abs(disc(d)) / 4}')
print('   d = 7: the only primes <= 21/4 that are non-split in Q(sqrt(-7)) are 3 and 5 (2 splits), and j(tau_7) = -3^3 5^3.')

print('\n== 3. Q(sqrt(-7)) as the decomposition field of 2 for the trivial cycle')
z = cmath.exp(2j * math.pi / 7)
g7 = z + z ** 2 + z ** 4
print('   doubling orbit of 1 mod 7:', sorted({pow(2, k, 7) for k in range(3)}), ' = quadratic residues mod 7:', sorted({x * x % 7 for x in range(1, 7)}),
      '  ord_7(2) =', sp.n_order(2, 7), '= (7 - 1)/2')
print(f'   Gauss period g = z + z^2 + z^4 = {g7:.12f};  (-1 + sqrt(-7))/2 = {(-1 + cmath.sqrt(-7)) / 2:.12f};  |g^2 + g + 2| = {abs(g7 ** 2 + g7 + 2):.2e};  g*conj(g) = {abs(g7) ** 2:.12f}')
zm = mp.exp(2j * mp.pi / 7)
gm = zm + zm ** 2 + zm ** 4
jg = mp.kleinj(gm) * 1728
print('   Im g = sqrt(7)/2 > 0:', check('Im g', abs(mp.im(gm) - mp.sqrt(7) / 2) < mp.mpf(10) ** -50),
      ';  g + 1 = tau_7 = (1 + sqrt(-7))/2:', check('g+1', abs(gm + 1 - tau_of(7)) < mp.mpf(10) ** -50),
      ';  j(g) =', mp.nstr(jg, 25), check('j(g)', abs(jg + 3375) < mp.mpf(10) ** -40))
A7 = np.array([[1 if (j - i) % 7 in (1, 2, 4) else 0 for j in range(7)] for i in range(7)])
ev = sorted(np.linalg.eigvals(A7), key=lambda x: (round(x.real, 6), round(x.imag, 6)))
print('   Paley P7 = Cay(Z7,{1,2,4}) eigenvalues:', [complex(round(x.real, 9), round(x.imag, 9)) for x in ev], ' = 3, g (x3), conj g (x3)')


def f8mul(x, y):
    r = 0
    for i in range(3):
        if (y >> i) & 1:
            r ^= x << i
    for i in (4, 3):
        if (r >> i) & 1:
            r ^= 0b1011 << (i - 3)
    return r


def f8pow(x, n):
    r = 1
    for _ in range(n):
        r = f8mul(r, x)
    return r


tr = {i: (f8pow(2, i) ^ f8pow(2, 2 * i % 7) ^ f8pow(2, 4 * i % 7)) for i in range(1, 7)}
print('   trace Tr_{F8/F2}(a^i) for a^3+a+1=0:', tr, ' -> trace-zero exponents', sorted(i for i, v in tr.items() if v == 0))
M40 = 2 ** 40
sols = []
for r0 in (0, 1):
    r = r0
    for k in range(1, 40):
        for c in (r, r + 2 ** k):
            if (c * c + c + 2) % 2 ** (k + 1) == 0:
                r = c
                break
        else:
            raise AssertionError('Hensel lift failed')
    sols.append(r)
for r in sols:
    v = (r & -r).bit_length() - 1 if r else 40
    print(f'   2-adic root of x^2+x+2: g = {r} mod 2^40, v_2(g) = {v};  (1 + 2g)^2 = -7 mod 2^40: {check("hensel", (1 + 2 * r) ** 2 % M40 == (-7) % M40)};'
          f'  Frobenius 1 + g = -conj(g) is a 2-adic {"unit" if (1 + r) % 2 else "non-unit"}')
# 2-adic code of 1 under T: parity sequence 1,0,0,1,0,0,...
x, phi = 1, 0
for i in range(40):
    phi |= (x % 2) << i
    x = x // 2 if x % 2 == 0 else 3 * x + 1
print('   Phi(1) = sum parity_i 2^i = -1/7 mod 2^40:', check('Phi(1)', phi == (-pow(7, -1, M40)) % M40),
      ';  (1 + 2g)^2 Phi(1) = 1 mod 2^40 for both roots g:', check('Phi2', all((1 + 2 * r) ** 2 * phi % M40 == 1 for r in sols)))
print('   real numerators {1,2,4} and 2-adic numerators -{1,2,4} = {3,5,6} are the two cosets of <2> in (Z/7)^*,',
      'i.e. the two primes above 2; complex conjugation x -> -x swaps them:', sorted((-x) % 7 for x in (1, 2, 4)))

print('\n== 4. Visibility: which Heegner fields lie inside a decomposition field of 2?')
print('   Q(sqrt D) lies in Q(zeta_m)^<2> (m odd) iff |D| divides m and chi_D(2) = +1 (2 splits).')
for d in HEEGNER:
    D = disc(d)
    vis = D % 2 != 0 and kron(D, 2) == 1
    print(f'   D = {D:5d}: chi_D(2) = {kron(D, 2):+d};  odd conductor: {D % 2 != 0};  visible to binary codes: {vis}')
check('visible only -7', [d for d in HEEGNER if disc(d) % 2 and kron(disc(d), 2) == 1] == [7])


def classno(D):
    """Number of reduced primitive forms of discriminant D < 0."""
    h, a = 0, 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a) == 0:
                c = (b * b - D) // (4 * a)
                if c >= a and not (b < 0 and a == c) and math.gcd(math.gcd(a, abs(b)), c) == 1:
                    h += 1
        a += 1
    return h


print('   sanity: h(-3), h(-4), h(-7), h(-23), h(-31), h(-163) =', [classno(D) for D in (-3, -4, -7, -23, -31, -163)])
h1 = [d for d in range(3, 20001) if d % 8 == 7 and sp.factorint(d) and max(sp.factorint(d).values()) == 1 and classno(-d) == 1]
print('   squarefree d = 7 mod 8 (2 splits in Q(sqrt(-d))), d <= 20000, with class number 1:', h1)
check('h1 only 7', h1 == [7])
print('   norm lemma: if p splits in Q(sqrt(-d)) with h = 1 then p = x^2 + xy + ((d+1)/4) y^2 with y != 0, so p >= (d+1)/4:')
for d in [7, 11, 19, 43, 67, 163]:
    ps = next(p for p in sp.primerange(2, 10 ** 4) if kron(-d, p) == 1)
    print(f'     d = {d:3d}: smallest split prime {ps:3d};  (d+1)/4 = {(d + 1) // 4:3d};  equal: {check("norm lemma", ps == (d + 1) // 4)}')

print('\n== 5. Gauss periods of the codes of every integer cycle (both codings)')


def T_ns(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


def T_sc(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def cycle_of(f, n0):
    seq = [n0]
    while True:
        n = f(seq[-1])
        if n == n0:
            return seq
        seq.append(n)


def n_distinct(vals, dec):
    return len({(round(v.real, dec), round(v.imag, dec)) for v in vals})


def prime_1_mod(m, start=10 ** 18):
    k = start // m + 1
    while not sp.isprime(k * m + 1):
        k += 1
    return k * m + 1


def exact_degree(orbit, m):
    """Exact degree of sum_{c in orbit} zeta_m^c (orbit stable under x2): reduce Z[zeta_m] -> F_l (l = 1 mod m prime),
    count distinct images of the conjugates sigma_a, a over U/<2>.  Distinct images => distinct conjugates, so the count is a
    lower bound for the degree; the number of cosets is an upper bound.  Equality = exact degree."""
    l = prime_1_mod(m)
    qs = list(sp.factorint(m))
    a0 = 2
    while True:
        w = pow(a0, (l - 1) // m, l)
        if all(pow(w, m // q, l) != 1 for q in qs):
            break
        a0 += 1
    H = {pow(2, i, m) for i in range(sp.n_order(2, m))}
    seen, reps = set(), []
    for a in range(1, m):
        if math.gcd(a, m) == 1 and a not in seen:
            reps.append(a)
            seen.update(a * h % m for h in H)
    vals = {sum(pow(w, a * c % m, l) for c in orbit) % l for a in reps}
    return len(vals), len(reps)


for name, f, n0 in [('T', T_ns, 1), ('T', T_ns, -1), ('T', T_ns, -5), ('T', T_ns, -17),
                    ('T_1', T_sc, 1), ('T_1', T_sc, -1), ('T_1', T_sc, -5), ('T_1', T_sc, -17)]:
    cyc = cycle_of(f, n0)
    w = [n % 2 for n in cyc]
    L = len(w)
    m = 2 ** L - 1
    real = sorted({sum(w[(k + i) % L] << (L - 1 - i) for i in range(L)) % m for k in range(L)})
    two = sorted({(-sum(w[(k + i) % L] << i for i in range(L))) % m for k in range(L)})
    if m == 1:
        print(f'   {name:3s} cycle {cyc}: L = 1, m = 1, codes are integers; Gauss period 1 (rational)')
        continue
    units = np.array([a for a in range(1, m) if math.gcd(a, m) == 1], dtype=np.int64)
    out = []
    for lab, R in (('real', real), ('2-adic', two)):
        Rv = np.array(R, dtype=np.int64)
        G = np.exp(2j * np.pi * ((units[:, None] * Rv[None, :]) % m) / m).sum(axis=1)
        G0 = G[0]                                    # a = 1
        deg5, deg8 = n_distinct(G, 5), n_distinct(G, 8)
        kind = ('imaginary quadratic (a CM point after conj/translation)' if deg5 == 2 and abs(G0.imag) > 1e-6 else
                'rational' if deg5 == 1 else
                f'degree {deg5} (not quadratic: j there is transcendental by Schneider if Im > 0)')
        lo, hi = exact_degree(R, m)
        check(f'exact degree {n0} {lab}', lo == hi == deg5)
        out.append(f'{lab}: period {G0.real:+.6f}{G0.imag:+.6f}i, conjugates {deg5} (={deg8} at 1e-8; exact: {lo} distinct mod a prime l = 1 mod m, of {hi} cosets) -> {kind}')
    print(f'   {name:3s} cycle from {n0} (L = {L}, word {"".join(map(str, w))}), m = {m} = {sp.factorint(m)}')
    for s in out:
        print('      ' + s)
# control: a rational cycle can carry a CM point of higher class number.  Under T, 1/5 -> 8/5 -> 4/5 -> 2/5 -> 1/5
# has word 1000 (L = 4), real code orbit {1,2,4,8}/15; <2> has index 2 in (Z/15)^*.
z15 = mp.exp(2j * mp.pi / 15)
g15 = sum(z15 ** k for k in (1, 2, 4, 8))
print('   control, rational cycle 1/5 of T (word 1000): Gauss period of {1,2,4,8}/15 =', mp.nstr(g15, 15), '= (1 + sqrt(-15))/2:',
      check('g15', abs(g15 - (1 + mp.sqrt(-15)) / 2) < mp.mpf(10) ** -50), '(the unit sum is mu(15) = +1, not -1);  h(-15) =', classno(-15))
cmrows = []
for p in sp.primerange(3, 400):
    if p % 8 == 7 and sp.n_order(2, p) == (p - 1) // 2:
        cmrows.append((p, (p - 1) // 2, classno(-p)))
print('   primes p = 7 mod 8 with ord_p(2) = (p-1)/2 (a doubling orbit of size (p-1)/2 whose Gauss period is the CM point (-1+sqrt(-p))/2),')
print('   as (p, orbit size, h(-p)):', cmrows, ';  class number 1 only at p = 7:', check('cm census', [r for r in cmrows if r[2] == 1] == [(7, 3, 1)]))
# Propositions 4-5: census of all unit cosets of <2> mod m (odd m < 3001) = real codes of all rational cycles of T_1 (lowest terms).
# A coset is the code of a T-cycle iff its word has no "11" iff no element c has c/m > 3/4.
MC = 3001
cm_rows, t_codes, bad6 = [], [], []
for m in range(3, MC, 2):
    U = [a for a in range(1, m) if math.gcd(a, m) == 1]
    H = [pow(2, i, m) for i in range(sp.n_order(2, m))]
    lab = {}
    cos = []
    for a in U:
        if a not in lab:
            cs = sorted({a * h % m for h in H})
            for c in cs:
                lab[c] = len(cos)
            cos.append(cs)
    if len(cos) == 1:
        continue
    zz = np.exp(2j * np.pi * np.arange(m) / m)
    ids = np.array([lab[a] for a in U])
    Ua = np.array(U)
    per = np.bincount(ids, weights=zz[Ua].real) + 1j * np.bincount(ids, weights=zz[Ua].imag)
    keys = {(round(v.real, 7), round(v.imag, 7)) for v in per}
    if len(keys) == 2 and abs(per[0].imag) > 1e-6:
        sqf = max(sp.factorint(m).values()) == 1
        D = int(round(((per[0] - per[1]) ** 2).real))
        mu = int(sp.mobius(m))
        ok6 = sqf and len(cos) == 2 and abs(per[0] + per[1] - mu) < 1e-8 and D % 8 == 1 and max(sp.factorint(-D).values()) == 1
        if not ok6:
            bad6.append(m)
        hD = classno(D)
        cm_rows.append((m, D, hD))
        for cs, pv in zip(cos, per):
            if not any(4 * c > 3 * m for c in cs):
                t_codes.append((m, cs, D, hD))
print(f'   Proposition 4 census, odd m < {MC}: {len(cm_rows)} denominators with an imaginary quadratic Gauss period;',
      'all with m squarefree, index 2, periods (mu(m) +- sqrt D)/2, D fundamental = 1 mod 8:', check('prop6', not bad6), bad6[:5])
h1m = [r[0] for r in cm_rows if r[2] == 1]
pred = [7] + [7 * p for p in sp.primerange(3, MC // 7 + 1) if p != 7 and sp.n_order(2, p) == p - 1 and p % 3 != 1 and 7 * p < MC]
print('   class-number-one cases (all D = -7):', sorted({r[1] for r in cm_rows if r[2] == 1}), ';  denominators', h1m,
      ';  = predicted {7} u {7p : 2 primitive mod p, p != 1 mod 3}:', check('h1 m list', h1m == pred), f'({len(pred) - 1} values 7p)')
print('   Proposition 5: CM cosets that are codes of T-cycles (no "11"): ', [(m, cs, D, h) for m, cs, D, h in t_codes],
      check('prop7', [(m, D) for m, cs, D, h in t_codes] == [(7, -7), (15, -15)]))
# the T_1 cycle of 1/11 carries tau_7 itself
from fractions import Fraction as Fr
x = Fr(1, 11)
cyc11, word11 = [], []
while True:
    cyc11.append(x)
    word11.append(x.numerator % 2)
    x = (3 * x + 1) / 2 if x.numerator % 2 else x / 2
    if x == Fr(1, 11):
        break
z21 = mp.exp(2j * mp.pi / 21)
g21 = sum(z21 ** k for k in (1, 2, 4, 8, 16, 11))
print('   T_1 cycle', [str(c) for c in cyc11], 'word', ''.join(map(str, word11)), ': code orbit {1,2,4,8,16,11}/21, Gauss period',
      mp.nstr(g21, 15), '= tau_7:', check('g21', abs(g21 - tau_of(7)) < mp.mpf(10) ** -50))

print('\n== 6. Ordinary elliptic curves over F_2, and E_7 = 49a1 over F_(2^k)')


def disc_weier(a1, a2, a3, a4, a6):
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    return -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6


rows = []
for t in range(32):
    a1, a2, a3, a4, a6 = [(t >> i) & 1 for i in range(5)]
    if disc_weier(a1, a2, a3, a4, a6) % 2 == 0:
        continue
    npts = 1 + sum(1 for xx in (0, 1) for yy in (0, 1)
                   if (yy * yy + a1 * xx * yy + a3 * yy - (xx ** 3 + a2 * xx * xx + a4 * xx + a6)) % 2 == 0)
    a = 3 - npts
    rows.append((a1, a, a * a - 8))
ordinary = [r for r in rows if r[1] % 2]
print(f'   {len(rows)} smooth Weierstrass tuples over F_2; ordinary (trace odd): {len(ordinary)}, all with a1 = 1: {all(r[0] == 1 for r in ordinary)};'
      f'  traces {sorted(Counter(r[1] for r in ordinary).items())};  Frobenius discriminants a^2 - 8: {sorted(set(r[2] for r in ordinary))}')
check('ordinary F2', ordinary and all(r[2] == -7 for r in ordinary) and all(r[0] == 1 for r in ordinary))
print('   supersingular (trace even):', sorted(Counter(r[1] for r in rows if r[1] % 2 == 0).items()),
      ' -> every ordinary curve over F_2 has Frobenius (+-1 +- sqrt(-7))/2 in Q(sqrt(-7)); x^2 - x + 2 has roots -g, -conj(g)')
POLYS = {1: 0b11, 2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011, 7: 0b10000011, 8: 0b100011101}


def gf_table(k):
    poly, q = POLYS[k], 1 << k
    tab = [[0] * q for _ in range(q)]
    for x0 in range(q):
        for y0 in range(q):
            r, xx, yy = 0, x0, y0
            while yy:
                if yy & 1:
                    r ^= xx
                yy >>= 1
                xx <<= 1
                if (xx >> k) & 1:
                    xx ^= poly
            tab[x0][y0] = r
    return tab


TABS = {k: gf_table(k) for k in range(1, 9)}


def count_E7_gf(k):
    """#E(F_(2^k)) for 49a1 mod 2: y^2 + xy = x^3 + x^2 + 1."""
    mul, q = TABS[k], 1 << k
    n = 1
    for xx in range(q):
        x2 = mul[xx][xx]
        rhs = mul[x2][xx] ^ x2 ^ 1
        for yy in range(q):
            if mul[yy][yy] ^ mul[xx][yy] == rhs:
                n += 1
    return n


tk = [2, 1]
for k in range(2, 37):
    tk.append(tk[-1] - 2 * tk[-2])
direct = {k: count_E7_gf(k) for k in range(1, 9)}
rec = {k: 2 ** k + 1 - tk[k] for k in range(1, 37)}
print('   #E_7(F_(2^k)), k = 1..8, direct count:', direct, '  matches 2^k + 1 - t_k (t_k = t_(k-1) - 2 t_(k-2)):',
      check('E7 F2k', all(direct[k] == rec[k] for k in direct)))
print('   k <= 36 with 7 | #E_7(F_(2^k)):', [k for k in rec if rec[k] % 7 == 0], '  = multiples of 3:',
      check('7 | #E iff 3 | k', [k for k in rec if rec[k] % 7 == 0] == list(range(3, 37, 3))))
print('   (Frobenius pi = 1 + g acts on the 7-isogeny kernel E[sqrt(-7)] = Z/7 as multiplication by 4 = 2^(-1), an element of order 3 of <2>)')


def ok_mul(x, y, mod):
    # (a + b g)(c + d g) with g^2 = -g - 2
    return ((x[0] * y[0] - 2 * x[1] * y[1]) % mod, (x[0] * y[1] + x[1] * y[0] - x[1] * y[1]) % mod)


e_, ordpi = (1, 1), 1
while e_ != (1, 0):
    e_ = ok_mul(e_, (1, 1), 7)
    ordpi += 1
v7_21 = sp.multiplicity(7, rec[21])
print('   order of pi = 1 + g in (O_K/7)^*:', ordpi, check('ord21', ordpi == 21), ' -> E[7] is rational over F_(2^k) iff 21 | k;  v_7 #E_7(F_(2^21)) =', v7_21, check('v7', v7_21 == 3))

print('\n== 7. CM curve E_7 = 49a1: y^2 + xy = x^3 - x^2 - 2x - 1 (j = -3375)')


def count_49a1(p):
    n = 1
    for xx in range(p):
        for yy in range(p):
            if (yy * yy + xx * yy - (xx ** 3 - xx * xx - 2 * xx - 1)) % p == 0:
                n += 1
    return n


ap7 = {p: p + 1 - count_49a1(p) for p in sp.primerange(2, 400) if p != 7}
bad = [(p, a) for p, a in ap7.items() if (a == 0) != (p % 7 in (3, 5, 6))]
print('   a_p for p < 60:', {p: a for p, a in ap7.items() if p < 60})
print('   a_p = 0 exactly for p = 3,5,6 mod 7 (p < 400, p != 7): violations', bad)
check('ap zero', not bad)
cm = [(p, a) for p, a in ap7.items() if p % 7 in (1, 2, 4) and not any(4 * p - a * a == 7 * B * B for B in range(0, math.isqrt(4 * p) + 1))]
print('   4p = a_p^2 + 7B^2 on split p < 400: violations', cm)
check('4p', not cm)

print('\n== 8. Klein quartic x^3 y + y^3 z + z^3 x = 0')


def klein_count(p):
    pts = 0
    for (xx, yy, zz) in [(1, a_, b_) for a_ in range(p) for b_ in range(p)] + [(0, 1, b_) for b_ in range(p)] + [(0, 0, 1)]:
        if (xx ** 3 * yy + yy ** 3 * zz + zz ** 3 * xx) % p == 0:
            pts += 1
    return pts


kq = []
for p in sp.primerange(2, 110):
    if p == 7:
        continue
    tw = 3 if p % 7 in (1, 6) else 0     # 1 + chi(p) + chi(p)^2 for the cubic character chi mod 7
    kq.append((p, klein_count(p), p + 1 - 3 * ap7[p], p + 1 - tw * ap7[p]))
print('   (p, #X(F_p), p+1-3a_p [untwisted], p+1-a_p(1+chi+chi^2) [cubic twist]):', kq)
print('   untwisted formula holds for all p:', all(a == b for _, a, b, c in kq),
      '   cubic-twist formula holds for all p < 110:', check('klein twist', all(a == c for _, a, b, c in kq)),
      '   #X = p+1 exactly unless p = 1 mod 7:', check('klein p+1', all(a == p + 1 for p, a, b, c in kq if p % 7 != 1)))


def klein_count_gf(n):
    mul, q = TABS[n], 1 << n
    cube = [mul[mul[t][t]][t] for t in range(q)]
    pts = 2                               # (0:1:0) and (0:0:1)
    for yy in range(q):
        for zz in range(q):
            # x = 1: y + y^3 z + z^3
            if yy ^ mul[cube[yy]][zz] ^ cube[zz] == 0:
                pts += 1
    return pts


kg = {n: klein_count_gf(n) for n in range(1, 7)}
pred = {}
u_sum = {0: 2, 1: -5}                    # power sums of the roots of u^2 + 5u + 8 (u = pi^3)
for k in range(2, 3):
    u_sum[k] = -5 * u_sum[k - 1] - 8 * u_sum[k - 2]
for n in range(1, 7):
    S = 0 if n % 3 else 3 * u_sum[n // 3]
    pred[n] = 2 ** n + 1 - S
print('   #X(F_(2^n)), n = 1..6, direct:', kg, ';  from L(T) = 1 + 5T^3 + 8T^6:', pred, check('klein F2n', kg == pred))
print('   Serre bound over F_8 for genus 3: 8 + 1 + 3*floor(2 sqrt 8) =', 9 + 3 * math.floor(2 * math.sqrt(8)), '-> the Klein quartic is maximal over F_8')

print('\n== 9. Supersingular congruences')


def ap_from_j(j, p):
    j %= p
    if j == 0:
        A, B = 0, 1
    elif j == 1728 % p:
        A, B = 1, 0
    else:
        k = (j * pow((1728 - j) % p, -1, p)) % p
        A, B = (3 * k) % p, (2 * k) % p
    sq = Counter((y * y) % p for y in range(p))
    return p + 1 - (1 + sum(sq[(xx ** 3 + A * xx + B) % p] for xx in range(p)))


def supersingular_set(p):
    return sorted(j for j in range(p) if ap_from_j(j, p) % p == 0)


ok163, ok7 = True, True
for p in sp.primerange(5, 72):
    ss = supersingular_set(p)
    j163, j7 = JKNOWN[163] % p, JKNOWN[7] % p
    ok163 &= (j163 in ss) == (kron(-163, p) != 1)
    if p != 7:
        ok7 &= (j7 in ss) == (kron(-7, p) != 1)
    print(f'   p={p:2d} supersingular j: {ss};  j(163) mod p = {j163} ({"ss" if j163 in ss else "ord"}, 163 {"inert" if kron(-163, p) == -1 else "split" if kron(-163, p) == 1 else "ram"});'
          f'  j(7) mod p = {j7} ({"ss" if j7 in ss else "ord"}, 7 {"inert" if kron(-7, p) == -1 else "split" if kron(-7, p) == 1 else "ram"})')
print('   Deuring pattern (supersingular iff non-split) holds for j(163) at 5 <= p < 72:', check('deuring163', ok163), ';  for j(7) (p != 7):', check('deuring7', ok7))
print('   first prime split in Q(sqrt(-163)):', next(p for p in sp.primerange(2, 1000) if kron(-163, p) == 1), '  (all primes < 41 inert: Euler n^2+n+41 / Rabinowitsch)')
print('   ord_163(2) =', sp.n_order(2, 163), ' (2 is a primitive root mod 163: a doubling orbit mod 163 is all of (Z/163)^*, decomposition field Q)')

print('\n== 10. E8 chain from {1,2,4}; gamma_2 = E4/eta^8')
gpoly = [1, 1, 0, 1]   # 1 + x + x^3, zeros a, a^2, a^4
codes = set()
for mm in range(16):
    c = [0] * 7
    for i in range(4):
        if (mm >> i) & 1:
            for k, gk in enumerate(gpoly):
                c[i + k] ^= gk
    codes.add(tuple(c))
wts = Counter(sum(c) for c in codes)
ext = {c + (sum(c) % 2,) for c in codes}
print('   Hamming [7,4,3] weight distribution:', dict(sorted(wts.items())), '; extended [8,4,4]:', dict(sorted(Counter(sum(c) for c in ext).items())))
roots = 16 + sum(1 for c in ext if sum(c) == 4) * 16
print('   E8 = Construction A (x in Z^8, x mod 2 in H8, norm x.x/2): norm-2 vectors = 16 + 14*16 =', roots, check('240', roots == 240))
NQ = 6
E4s = [1] + [240 * sum(dd ** 3 for dd in sp.divisors(n)) for n in range(1, NQ + 1)]
inv8 = [1] + [0] * NQ
for n in range(1, NQ + 1):
    for _ in range(8):
        for i in range(n, NQ + 1):
            inv8[i] += inv8[i - n]
g2s = series_mul(E4s, inv8)
jq = series_mul(series_mul(g2s, g2s), g2s)
print('   q^(1/3) gamma_2 = (theta_E8)(q^(1/3)/eta^8) = 1 + 248q + ... :', g2s[:4], ';  248 = 240 + 8:', 240 + 8)
print('   q j = (q^(1/3) gamma_2)^3 =', jq[:4], ';  744 = 3*248:', check('744', jq[1] == 744 == 3 * 248), ' 196884:', check('196884', jq[2] == 196884))
for lab, t, expect in [('(3 + i sqrt 163)/2', tau_of(163) + 1, -640320), ('(1 + i sqrt 163)/2', tau_of(163), None),
                       ('(3 + i sqrt 7)/2', tau_of(7) + 1, -15), ('i sqrt 7', mp.mpc(0, mp.sqrt(7)), 255)]:
    v = gamma2(t)
    tag = f'(expected {expect})' if expect is not None else '(= 640320 e^(-i pi/3): gamma_2(tau + 1) = e^(-2 pi i/3) gamma_2(tau))'
    if expect is not None:
        check(f'gamma2 {lab}', abs(v - expect) < mp.mpf(10) ** -30)
    print(f'   gamma_2({lab}) = {mp.nstr(v, 22)}  {tag}')
print('   640320 =', sp.factorint(640320), '; Ramanujan constant e^(pi sqrt 163) =', mp.nstr(mp.e ** (mp.pi * mp.sqrt(163)), 40))

print('\n== 11. Weber: f(sqrt(-7))^24 = 2^12 and the d = 7 near-integers')
x7 = mp.e ** (-mp.pi * mp.sqrt(7))
f24 = mp.nprod(lambda n: (1 + x7 ** (2 * n - 1)) ** 24, [1, mp.inf]) / x7
print('   f(sqrt(-7))^24 = e^(pi sqrt 7) prod (1 + e^(-(2n-1) pi sqrt 7))^24 =', mp.nstr(f24, 40), check('f24', abs(f24 - 4096) < mp.mpf(10) ** -40))
NX = 6
ser = [1] + [0] * NX                 # prod (1 + x^(2n-1))^24 as a series in x
for n in range(1, NX + 1):
    e = 2 * n - 1
    if e > NX:
        break
    for _ in range(24):
        for i in range(NX, e - 1, -1):
            ser[i] += ser[i - e]
print('   x f^24 = prod (1 + x^(2n-1))^24 = ', ser[:5], ' -> e^(pi sqrt 7) = 2^12 - 24 - 276 x - 2048 x^2 - ..., x = e^(-pi sqrt 7)')
ep7 = 1 / x7
print('   e^(pi sqrt 7) =', mp.nstr(ep7, 25), ';  minus 4072 =', mp.nstr(ep7 - 4072, 12),
      ';  minus (4096 - 24 - 276x - 2048x^2) =', mp.nstr(ep7 - (4096 - 24 - 276 * x7 - 2048 * x7 ** 2), 6))
j_sqrt7 = mp.kleinj(mp.mpc(0, mp.sqrt(7))) * 1728
print('   j(sqrt(-7)) =', mp.nstr(j_sqrt7, 25), '= 255^3 =', 255 ** 3, check('j sqrt -7', abs(j_sqrt7 - 255 ** 3) < mp.mpf(10) ** -30),
      ';  (4096 - 16)^3/4096 =', (4096 - 16) ** 3 // 4096, ' (Weber j = (f^24 - 16)^3/f^24)')
e2 = mp.e ** (2 * mp.pi * mp.sqrt(7))
print('   e^(2 pi sqrt 7) = e^(pi sqrt 28) =', mp.nstr(e2, 20), ';  minus (255^3 - 744) =', mp.nstr(e2 - (255 ** 3 - 744), 10))
Xr = sp.symbols('X')
print('   roots of (X - 16)^3 = 255^3 X:', [sp.nsimplify(r) for r in sp.solve((Xr - 16) ** 3 - 255 ** 3 * Xr, Xr)], '-> the only positive root is 4096 = 2^12')
print('   gamma_2: -15 = -(2^4 - 1) at (3 + sqrt(-7))/2, 255 = 2^8 - 1 at sqrt(-7);  Ramanujan G_7 = 2^(-1/4) f(sqrt(-7)) =', mp.nstr(mp.root(f24, 24) / mp.root(2, 4), 20), '= 2^(1/4) =', mp.nstr(mp.root(2, 4), 20))
k7 = (mp.jtheta(2, 0, x7) / mp.jtheta(3, 0, x7)) ** 2           # singular modulus, nome e^(-pi sqrt 7)
print('   singular modulus k_7 = theta_2^2/theta_3^2 =', mp.nstr(k7, 25), '= (3 - sqrt 7)/(4 sqrt 2):', check('k7', abs(k7 - (3 - mp.sqrt(7)) / (4 * mp.sqrt(2))) < mp.mpf(10) ** -50),
      ';  4 k^2 k\'^2 =', mp.nstr(4 * k7 ** 2 * (1 - k7 ** 2), 20), '= 1/64 = G_7^(-24)')
Sram = mp.nsum(lambda n: (42 * n + 5) * (mp.rf(mp.mpf(1) / 2, n) / mp.factorial(n)) ** 3 / mp.mpf(64) ** n, [0, mp.inf])
print('   Ramanujan (1914): sum (42k+5) ((1/2)_k/k!)^3 / 64^k =', mp.nstr(Sram, 45), ';  16/pi =', mp.nstr(16 / mp.pi, 45), check('ram series', abs(Sram - 16 / mp.pi) < mp.mpf(10) ** -50))


def weber_w(m):
    xm = mp.e ** (-mp.pi * mp.sqrt(m))
    return mp.root(mp.nprod(lambda n: (1 + xm ** (2 * n - 1)) ** 24, [1, mp.inf]) / xm, 24) / mp.sqrt(2)


print('   f(sqrt(-m))/sqrt 2 for m = 7 mod 8 (2 splits in Q(sqrt(-m))): minimal polynomial found by an integer-relation search')
units_all = True
for m in (7, 15, 23, 31, 39, 47, 55, 71):
    hm = classno(-4 * m)
    with mp.workdps(200):
        wv = weber_w(m)
        found, var = None, 'w'
        for deg in range(1, 3 * hm + 1):
            found = mp.findpoly(wv, deg, maxcoeff=10 ** 4)
            if found:
                break
        if not found:                       # PSLQ can miss a degree-12 relation; search in u = w^3 instead
            var = 'w^3'
            for deg in range(1, 3 * hm + 1):
                found = mp.findpoly(wv ** 3, deg, maxcoeff=10 ** 4)
                if found:
                    break
        resid = abs(mp.polyval(found, wv if var == 'w' else wv ** 3)) if found else None
        wv_s = mp.nstr(wv, 18)
    unit = bool(found) and abs(found[0]) == 1 and abs(found[-1]) == 1 and resid < mp.mpf(10) ** -150
    units_all &= unit
    print(f'     m = {m:2d}: h(-4m) = {hm}, w = {wv_s}, minimal polynomial of {var} (high to low) {found if found else "not found"}, residual {mp.nstr(resid, 3) if found else "-"}, unit: {unit if found else "undecided"}')
print('   all eight are units:', check('weber units', units_all))

print('\n== 12. Ramanujan tau mod 7 (and mod 49 on the non-residues)')
NT, MOD = 6000, 49
P = np.zeros(NT + 1, dtype=np.int64)
for k in range(-90, 91):
    e = k * (3 * k - 1) // 2
    if 0 <= e <= NT:
        P[e] += -1 if k % 2 else 1


def cmul(a, b):
    return np.convolve(a % MOD, b % MOD)[:NT + 1] % MOD


P2 = cmul(P, P)
P4 = cmul(P2, P2)
P8 = cmul(P4, P4)
P16 = cmul(P8, P8)
P24 = cmul(P16, P8)
tau_mod = {n: int(P24[n - 1]) for n in range(1, NT + 1)}
known = {1: 1, 2: -24, 3: 252, 4: -1472, 5: 4830, 6: -6048, 7: -16744, 8: 84480, 9: -113643, 10: -115920}
print('   tau(1..10) mod 49 agree with the known values:', check('tau known', all(tau_mod[n] == known[n] % MOD for n in known)))
sig9 = [0] * (NT + 1)
for dd in range(1, NT + 1):
    pw = pow(dd, 9, MOD)
    for mm in range(dd, NT + 1, dd):
        sig9[mm] = (sig9[mm] + pw) % MOD
c7 = all((tau_mod[n] - n * sig9[n]) % 7 == 0 for n in range(1, NT + 1))
c49 = all((tau_mod[n] - n * sig9[n]) % 49 == 0 for n in range(1, NT + 1) if n % 7 in (3, 5, 6))
c49_all = sum(1 for n in range(1, NT + 1) if (tau_mod[n] - n * sig9[n]) % 49)
print(f'   tau(n) = n sigma_9(n) mod 7 for all n <= {NT}:', check('tau mod 7', c7),
      f';  mod 49 for all n = 3,5,6 mod 7:', check('tau mod 49', c49), f'  (mod 49 fails for {c49_all} other n <= {NT})')
zero_set = sorted({p % 7 for p in sp.primerange(2, NT + 1) if tau_mod[p] % 7 == 0})
nonzero_set = sorted({p % 7 for p in sp.primerange(2, NT + 1) if tau_mod[p] % 7 != 0})
print('   residues mod 7 of primes p <= %d with tau(p) = 0 mod 7:' % NT, zero_set, ';  with tau(p) != 0 mod 7:', nonzero_set,
      check('tau split', zero_set == [0, 3, 5, 6] and nonzero_set == [1, 2, 4]))
print('   tau(p) = 2p mod 7 for p = 1,2,4 mod 7:', check('tau 2p', all((tau_mod[p] - 2 * p) % 7 == 0 for p in sp.primerange(2, NT + 1) if p % 7 in (1, 2, 4))),
      ';  mechanism: tau(p) = p + p^4 = p(1 + p^3) mod 7 and p^3 = (p/7) mod 7 (Euler, exponent 3 = (7-1)/2)')

print('\n== 13. Chudnovsky: 1/pi = 12 sum (-1)^k (6k)! (13591409 + 545140134 k) / ((3k)! k!^3 640320^(3k+3/2))')
C = 640320
num = C ** 3 + 1728
m2, rem = divmod(num, 163)
mroot = math.isqrt(m2)
print('   1728 - j(tau_163) = 640320^3 + 1728 = 163 * m^2 with m =', mroot, sp.factorint(mroot), check('chud square', rem == 0 and mroot * mroot == m2))
Bc, Ac = 545140134, 13591409
print('   545140134 =', sp.factorint(Bc), ';  12 * 545140134 = 163 m:', check('12B', 12 * Bc == 163 * mroot), ';  13591409 =', sp.factorint(Ac))
for p in sp.factorint(mroot):
    print(f'     p = {p:3d}: in Q(sqrt(-163)) {"inert" if kron(-163, p) == -1 else "split" if kron(-163, p) == 1 else "ram"}, in Q(i) {"inert" if kron(-4, p) == -1 else "split" if kron(-4, p) == 1 else "ram"}')
check('GZ 1728', all(kron(-163, p) != 1 and kron(-4, p) != 1 for p in sp.factorint(mroot)))
gzn = {y: 163 - y * y for y in range(0, 13)}
print('   Gross-Zagier norms (d1 d2 - x^2)/4 = 163 - y^2 (d2 = 4, x = 2y):', {y: sp.factorint(v) for y, v in gzn.items()})
print('   every prime of m divides some 163 - y^2:', check('GZ norms', all(any(v % p == 0 for v in gzn.values()) for p in sp.factorint(mroot))),
      ';  127 = 163 - 6^2;  y with 7 | 163 - y^2:', [y for y in gzn if gzn[y] % 7 == 0], '(163 = 2 = 3^2 mod 7 is a square: 163 mod 7 lies in {1,2,4})')
for K in range(1, 4):
    s = mp.fsum((-1) ** k * mp.factorial(6 * k) * (Ac + Bc * k) / (mp.factorial(3 * k) * mp.factorial(k) ** 3 * mp.mpf(C) ** (3 * k + mp.mpf(3) / 2)) for k in range(K))
    print(f'     {K} term(s): |12 S - 1/pi| = {mp.nstr(abs(12 * s - 1 / mp.pi), 5)}')
print('   j(tau_163) = 1728 mod 7^2:', JKNOWN[163] % 49 == 1728 % 49, '(7 | m; 7 is inert in Q(sqrt(-163)) because -163 = 5 mod 7 is a non-residue)')

print('\n== 14. Landau-Ramanujan constant of Q(sqrt(-7))')
L1 = mp.pi / mp.sqrt(7)              # L(1, chi_{-7}) = 2 pi h / (w sqrt 7) with h = 1, w = 2
prod = mp.mpf(1)
for qq in sp.primerange(2, 2_000_000):
    if qq % 7 in (3, 5, 6):
        prod *= 1 / (1 - mp.mpf(qq) ** -2)
C7 = mp.sqrt(L1 * mp.mpf(7) / 6 * prod) / mp.sqrt(mp.pi)
prod4 = mp.mpf(1)
for qq in sp.primerange(3, 2_000_000):
    if qq % 4 == 3:
        prod4 *= 1 / (1 - mp.mpf(qq) ** -2)
C4 = mp.sqrt(prod4) / mp.sqrt(2)
print('   C_7 = (1/sqrt(pi)) * sqrt( L(1,chi_-7) * 7/6 * prod_{q=3,5,6 mod 7} (1-q^-2)^-1 ) =', mp.nstr(C7, 12), ' (Euler product truncated at 2*10^6)')
print('   control: Landau-Ramanujan K (sums of two squares) by the same formula =', mp.nstr(C4, 12), ' (known 0.764223653589...)')
X = 10 ** 7
spf = list(range(X + 1))
for i in range(2, int(X ** 0.5) + 1):
    if spf[i] == i:
        for k in range(i * i, X + 1, i):
            if spf[k] == k:
                spf[k] = i
cnt = 0
checkpoints = {10 ** 4: None, 10 ** 5: None, 10 ** 6: None, 10 ** 7: None}
for n in range(1, X + 1):
    mm = n
    good = True
    while mm > 1:
        p = spf[mm]
        e = 0
        while mm % p == 0:
            mm //= p
            e += 1
        if p % 7 in (3, 5, 6) and e % 2:
            good = False
            break
    cnt += good
    if n in checkpoints:
        checkpoints[n] = cnt
for n, c in checkpoints.items():
    print(f'   B_7({n:.0e}) = {c};  B * sqrt(log x)/x = {c * math.sqrt(math.log(n)) / n:.6f}')
print('   (B_7(x) counts n <= x of the form x^2 + xy + 2y^2; convergence to C_7 is slow, with a 1/log x correction, as for sums of two squares)')

print('\n== 15. Numerology probes')
vals = {}
for A in range(0, 60):
    for l in range(0, 40):
        v = abs(2 ** A - 3 ** l)
        if v <= 200 and v not in vals:
            vals[v] = (A, l)
print('   Heegner numbers of the form |2^A - 3^l| (A, l):', {d: vals.get(d) for d in HEEGNER},
      ';  share of 1..20 of that form:', sum(1 for v in range(1, 21) if v in vals), '/ 20')
print('   163 mod 7 =', 163 % 7, '(a residue: 163 splits in Q(sqrt(-7)), 163 = 10^2 + 7*3^2 =', 10 ** 2 + 7 * 9, ');',
      ' 7 in Q(sqrt(-163)):', 'inert' if kron(-163, 7) == -1 else 'split', ' (reciprocity twins: (-7/163) = (163/7), (-163/7) = (7/163))')
print('   640320 mod 7 =', 640320 % 7, '; j(163) mod 7 =', JKNOWN[163] % 7, '= 1728 mod 7 =', 1728 % 7, '(the only supersingular j mod 7)')
print('   2^18 - 1 (code denominator of the -17 cycle under T) =', sp.factorint(2 ** 18 - 1), ';  2^4 - 1 = 15, 2^8 - 1 = 255 (gamma_2 values at the two d = -7, -28 points)')
sq23 = [d for d in range(1, 200) if sp.factorint(d) and max(sp.factorint(d).values()) == 1 and kron(-d if d % 4 == 3 else -4 * d, 2) == 1 and kron(-d if d % 4 == 3 else -4 * d, 3) == 1]
print('   smallest squarefree d with both 2 and 3 split in Q(sqrt(-d)):', sq23[:4], ';  h(-23) =', classno(-23))

print('\nFAILED CHECKS:', FAIL if FAIL else 'none')
print('ALL CHECKS PASSED' if not FAIL else 'SOME CHECKS FAILED')
