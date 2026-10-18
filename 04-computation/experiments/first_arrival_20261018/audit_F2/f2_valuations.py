#!/usr/bin/env python3
"""Audit F2, checks (1)-(3): the primes above 19 of F, valuations, sheet moves -- independent of nineteen_places.py.

F_361 is realised here as F_19[s]/(s^2 - 2) (2 is a non-residue mod 19; the session used t^2 = 8).
Ring homomorphisms Z[zeta, j, w] -> F_361 are found by brute force over all 361 field elements; primes above 19 are
their Frobenius orbits.  Since the tower order is 19-maximal (f2_disc.py) and 19 is unramified, for x in F with
19^k || (coordinate denominator), y = 19^k x is 19-integral and v_P(x) = -k iff y does not vanish at P (else >= 1 - k).
Exact valuations that one reduction cannot decide are settled by a second element (19/u, 19 u), and by PARI
(f2_pari_check.gp) independently.
"""
import sys, os, itertools
from collections import Counter
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f2_field import E, ONE, J, W, EXP, N, build, roots_of_unity, norm_cc, conj

P = 19
NR = 2                                     # s^2 = 2


def fmul(x, y): return ((x[0] * y[0] + NR * x[1] * y[1]) % P, (x[0] * y[1] + x[1] * y[0]) % P)
def fadd(x, y): return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)
def fneg(x): return ((-x[0]) % P, (-x[1]) % P)
def fpow(x, k):
    r = (1, 0)
    while k:
        if k & 1: r = fmul(r, x)
        x = fmul(x, x); k >>= 1
    return r
def frob(x): return fpow(x, P)


ELS = [(a, b) for a in range(P) for b in range(P)]
assert all(fmul((0, 1), (0, 1)) == (2, 0) for _ in [0]) and pow(2, (P - 1) // 2, P) == P - 1, '2 is a non-residue mod 19'
# order-15 elements (primitive 15th roots of unity)
Z15 = [x for x in ELS if x != (0, 0) and fpow(x, 15) == (1, 0) and fpow(x, 5) != (1, 0) and fpow(x, 3) != (1, 0)]
JR = [x for x in ELS if fmul(x, x) == ((-11) % P, 0)]
WR = [x for x in ELS if fadd(fadd(fpow(x, 3), fmul((2, 0), fpow(x, 2))), fadd(fmul(((-4) % P, 0), x), ((-4) % P, 0))) == (0, 0)]


def hom_value(x, hm):
    """Image of x (with 19-free coordinate denominators) under zeta, j, w -> hm."""
    z, j, w = hm
    assert x.d % P, 'coordinate denominator divisible by 19'
    dinv = pow(x.d, -1, P)
    s = (0, 0)
    for k, v in enumerate(x.n):
        if v:
            a, b, c = EXP[k]
            s = fadd(s, fmul(((v * dinv) % P, 0), fmul(fpow(z, a), fmul(fpow(j, b), fpow(w, c)))))
    return s


def k19(x):
    D, k = x.d, 0
    while D % P == 0: D //= P; k += 1
    return k


def val(x, hm):
    """(exact valuation or None, lower bound)."""
    k = k19(x)
    y = x.scale(P ** k)
    if hom_value(y, hm) != (0, 0):
        return -k, -k
    return None, 1 - k


def vstr(x, hm):
    v, lb = val(x, hm)
    return str(v) if v is not None else '>=%d' % lb


def main():
    print('F_361 = F_19[s]/(s^2-2): primitive 15th roots %d, sqrt(-11) %s, roots of w^3+2w^2-4w-4: %s' % (len(Z15), JR, WR))
    homs = [(z, j, w) for z in Z15 for j in JR for w in WR]
    # sanity: every defining relation vanishes
    for z, j, w in homs:
        zz = [fpow(z, k) for k in range(9)]
        phi15 = [1, -1, 0, 1, -1, 1, 0, -1, 1]                 # x^8 - x^7 + x^5 - x^4 + x^3 - x + 1, low to high
        s = (0, 0)
        for k, cf in enumerate(phi15): s = fadd(s, fmul(((cf) % P, 0), zz[k]))
        assert s == (0, 0)
    primes, seen = [], set()
    for hm in homs:
        if hm in seen: continue
        orb = [hm]; nxt = tuple(frob(x) for x in hm)
        while nxt != hm:
            orb.append(nxt); nxt = tuple(frob(x) for x in nxt)
        seen.update(orb); primes.append((hm, len(orb)))
    print('ring homomorphisms Z[zeta,j,w] -> F_361: %d ; Frobenius orbits (= primes above 19): %d ; orbit sizes %s'
          % (len(homs), len(primes), sorted(Counter(l for _, l in primes).items())))

    def typ(hm):
        z, j, w = hm; hr = fadd(w, (3, 0))
        if hr == j: return '+j'
        if hr == fneg(j): return '-j'
        if hr == (7, 0): return '7'
        return '?'
    # type is Frobenius invariant
    for hm, _ in primes:
        assert typ(hm) == typ(tuple(frob(x) for x in hm))
    print('types:', sorted(Counter(typ(hm) for hm, _ in primes).items()))

    O = build()
    b, c, u, u2, a = O['b'], O['c'], O['u'], O['u2'], O['a']
    h = O['h']
    # ---- table of valuations
    rows = []
    for hm, _ in primes:
        rows.append((typ(hm), vstr(b, hm), vstr(c, hm), vstr(u, hm), vstr(u2, hm), vstr(a, hm), vstr(b / c, hm)))
    print('per-prime valuations (type, v(b), v(c), v(u), v(u2), v(a), v(b/c)), grouped:')
    for k, n in sorted(Counter(rows).items()):
        print('   ', k, 'x', n)
    # ---- exact valuation of u at the +j primes via 19/u and at the -j primes via 19 u
    ui19 = u.inv().scale(19); u19 = u.scale(19)
    print('19/u has 19-free denominators:', k19(ui19) == 0, '; 19 u has 19-free denominators:', k19(u19) == 0)
    for t in ('+j', '-j', '7'):
        ps = [hm for hm, _ in primes if typ(hm) == t]
        print('  type %2s: 19/u nonzero at all: %s ; 19 u nonzero at all: %s ; (h-j) nonzero: %s ; (h+j) nonzero: %s'
              % (t, all(hom_value(ui19, hm) != (0, 0) for hm in ps), all(hom_value(u19, hm) != (0, 0) for hm in ps),
                 all(hom_value(h - J, hm) != (0, 0) for hm in ps), all(hom_value(h + J, hm) != (0, 0) for hm in ps)))
    print('  => v(u) = +1 exactly at the +j primes iff 19/u is nonzero there; v(u) = -1 at -j iff 19u nonzero there.')
    # ---- 19-integrality / units of the named sets
    M, C, A = O['M'], O['C'], O['A']
    roots = roots_of_unity()
    allp = [hm for hm, _ in primes]
    def unit_everywhere(x): return k19(x) == 0 and all(hom_value(x, hm) != (0, 0) for hm in allp)
    print('points of C 19-integral: %s ; points of A 19-integral: %s ; spindle points 19-integral: %s'
          % (all(k19(x) == 0 for x in C), all(k19(x) == 0 for x in A), all(k19(x) == 0 for x in M)))
    print('a = conj(d) unit at all 24: %s ; u2 unit at all 24: %s ; 30 roots of unity units at all 24: %s'
          % (unit_everywhere(a), unit_everywhere(u2), all(unit_everywhere(r) for r in roots)))
    Dl = []
    for x in M:
        for y in M:
            dd = y - x
            if dd not in Dl: Dl.append(dd)
    nz = [x for x in Dl if not x.is_zero()]
    print('|M - M| = %d (including 0); nonzero differences %d ; all units at all 24 primes: %s'
          % (len(Dl), len(nz), all(unit_everywhere(x) for x in nz)))
    spin_pairs = [(i, k) for i in range(7) for k in range(i)]
    print('the 21 spindle pairs m != m\': all differences units at all 24 primes:',
          all(unit_everywhere(M[i] - M[k]) for i, k in spin_pairs))
    # ---- b and c: simple poles everywhere; b/c unit
    print('b: 19 || den: %s, 19 b nonzero at all 24: %s ; c: 19 || den: %s, 19 c nonzero at all 24: %s'
          % (k19(b) == 1, all(hom_value(b.scale(19), hm) != (0, 0) for hm in allp),
             k19(c) == 1, all(hom_value(c.scale(19), hm) != (0, 0) for hm in allp)))
    r = b / c
    print('b/c and c/b 19-integral: %s %s ; b/c unit at all 24: %s' % (k19(r) == 0, k19(r.inv()) == 0, unit_everywhere(r)))
    # ---- sheet moves (check 3)
    P7 = [hm for hm in allp if typ(hm) == '7']
    zero = E.zero()
    bD = {i: b * x for i, x in enumerate(Dl)}
    cD = {i: c * x for i, x in enumerate(Dl)}
    dist7 = Counter(); integral_all = 0; per_type_dist = {t: Counter() for t in ('+j', '-j', '7')}
    kmax = 0
    for i2 in range(len(Dl)):
        for i3 in range(len(Dl)):
            if Dl[i2].is_zero() and Dl[i3].is_zero(): continue
            y = bD[i2] + cD[i3]
            k = k19(y); kmax = max(kmax, k)
            assert k <= 1
            yy = y.scale(19 ** k)
            integ = {hm: (k == 0 or hom_value(yy, hm) == (0, 0)) for hm in allp}
            dist7[sum(integ[hm] for hm in P7)] += 1
            for t in per_type_dist:
                per_type_dist[t][sum(integ[hm] for hm in allp if typ(hm) == t)] += 1
            if all(integ.values()): integral_all += 1
    tot = sum(dist7.values())
    print('sheet moves: %d nonzero pairs (dm2, dm3) in (M-M)^2 ; max 19-power in denominators %d' % (tot, kmax))
    print('  number of h=7 primes at which b dm2 + c dm3 is integral -> count:', sorted(dist7.items()))
    print('  same distribution at the +j / -j primes:', sorted(per_type_dist['+j'].items()), sorted(per_type_dist['-j'].items()))
    print('  pairs integral at all 24 primes (i.e. 19-integral):', integral_all)
    # ---- single-factor statement: b dm has a pole at all 24 primes for every nonzero dm
    print('single factor: every nonzero b*dm, c*dm has a pole at all 24 primes:',
          all(k19(b * x) == 1 and all(hom_value((b * x).scale(19), hm) != (0, 0) for hm in allp) for x in nz) and
          all(k19(c * x) == 1 and all(hom_value((c * x).scale(19), hm) != (0, 0) for hm in allp) for x in nz))
    # direction families: v(u) = 0 at the h = 7 primes (h - j and h + j both nonzero there, printed above) and u2 is a
    # 19-unit, so v(zeta^l u^k u2^m q) = v(q) >= 0 there for every k, m in Z and 19-integral q (multiplicativity).


if __name__ == '__main__':
    main()
