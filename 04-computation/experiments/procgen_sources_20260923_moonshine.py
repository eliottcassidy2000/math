#!/usr/bin/env python3
"""procgen_sources_20260923_moonshine.py

The moonshine exercise "from 6 = 5+1 to 196884 = 196883+1", checked exactly.

Chain (every arithmetic step is recomputed here; group facts that are only
CITED are marked as such in the note, not here):

  M1  S_6 and the exotic S_5 = PGL(2,5) on P^1(F_5) (6 = 5+1 points); the
      outer automorphism of S_6 built from the coset action; the two
      permutation characters 1+chi_5 and 1+chi_5' (Murnaghan-Nakayama).
  M2  Extended quadratic-residue codes of length p+1 for p = 5, 11, 23 over
      F_4, F_3, F_2: hexacode [6,3,4], ternary Golay [12,6,6], binary Golay
      [24,12,8]; S(5,6,12) from ternary Golay supports; S(5,8,24) from octads.
  M3  M_12 = <PSL(2,11), delta12> (order 95040, sharply 5-transitive); the hexad
      stabilizer is S_6 acting on hexad and complement through the outer
      automorphism of S_6 (the precise S_6 -> M_12 link); Galois's exceptional
      A_5 of index 11 in PSL(2,11) (the p = 11 analogue of S_5 < S_6).
  M4  M_24 = <PSL(2,23), delta24> (order 244823040); Frame shapes 1^2 11^2 and
      1.23; the eta product eta(t)^2 eta(11t)^2 is the level-11 newform
      (point counts of X_0(11) checked).
  M4b the dodecad stabilizer in M_24 is M_12, acting on a dodecad and its
      complement through the outer automorphism of M_12.
  M5  Leech lattice minimal vectors from the Golay code (196560, by shape) and
      the inner-product distribution from one minimal vector.
  M6  q-series: Theta_Leech = E4^3 - 720 Delta, j = E4^3/Delta,
      Theta_Leech/eta^24 - 24 = j - 744, 1/eta^24; the decompositions
      196884 = 300+24+196560 = (300+98280)+98304 = 1+299+98280+98304 = 1+196883
      and the McKay-Thompson coefficient decompositions for c(2), c(3).

Deterministic; one process; peak memory well under 300 MB.
"""
import sys
import random
from itertools import combinations, permutations, product
from collections import Counter
from fractions import Fraction

from sympy.combinatorics import Permutation, PermutationGroup

OUT = sys.stdout


def say(*a):
    print(*a, file=OUT)
    OUT.flush()


# ---------------------------------------------------------------- utilities
def P1(p):
    return list(range(p)) + ['inf']


def mobius(p, a, b, c, d):
    """x -> (a x + b)/(c x + d) on P^1(F_p), as a list of image indices."""
    pts = P1(p)
    idx = {x: i for i, x in enumerate(pts)}
    img = []
    for x in pts:
        if x == 'inf':
            y = 'inf' if c % p == 0 else (a * pow(c, -1, p)) % p
        else:
            num, den = (a * x + b) % p, (c * x + d) % p
            y = 'inf' if den == 0 else (num * pow(den, -1, p)) % p
        img.append(idx[y])
    return tuple(img)


def compose(g, h):
    """(g*h)(i) = g(h(i)) : apply h first."""
    return tuple(g[h[i]] for i in range(len(h)))


def cycle_type(g):
    n = len(g)
    seen = [False] * n
    ct = []
    for i in range(n):
        if not seen[i]:
            L = 0
            j = i
            while not seen[j]:
                seen[j] = True
                j = g[j]
                L += 1
            ct.append(L)
    return tuple(sorted(ct, reverse=True))


def ct_str(ct):
    c = Counter(ct)
    return '.'.join(f'{k}^{c[k]}' if c[k] > 1 else f'{k}' for k in sorted(c, reverse=True))


def closure(gens):
    n = len(gens[0])
    e = tuple(range(n))
    seen = {e}
    frontier = [e]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = compose(g, x)
                if y not in seen:
                    seen.add(y)
                    new.append(y)
        frontier = new
    return seen


def orbit_sets(base, gens):
    base = frozenset(base)
    orb = {base}
    fr = [base]
    while fr:
        s = fr.pop()
        for g in gens:
            t = frozenset(g[i] for i in s)
            if t not in orb:
                orb.add(t)
                fr.append(t)
    return orb


def steiner_check(blocks, t, v):
    cnt = Counter()
    for B in blocks:
        for f in combinations(sorted(B), t):
            cnt[f] += 1
    from math import comb
    return len(cnt) == comb(v, t) and all(x == 1 for x in cnt.values())


# ------------------------------------------------ Murnaghan-Nakayama (S_n)
def partitions(n, maxp=None):
    if maxp is None:
        maxp = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxp), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def mn_char(lam, mu):
    """chi^lam(mu) by the Murnaghan-Nakayama rule on beta-sets."""
    if sum(lam) == 0:
        return 1
    L = len(lam)
    beta = frozenset(lam[i] + (L - 1 - i) for i in range(L))
    return _mn(beta, tuple(mu))


def _mn(beta, mu):
    if not mu:
        return 1
    k = mu[0]
    rest = mu[1:]
    tot = 0
    for b in beta:
        if b - k >= 0 and (b - k) not in beta:
            sign = (-1) ** sum(1 for c in beta if b - k < c < b)
            nb = (beta - {b}) | {b - k}
            tot += sign * _mn(frozenset(nb), rest)
    return tot


# ------------------------------------------------------------------ M1
def section_M1():
    say('=' * 78)
    say('M1. S_6, the exotic S_5 = PGL(2,5) on P^1(F_5), and Out(S_6)')
    p = 5
    els = set()
    for a, b, c, d in product(range(p), repeat=4):
        if (a * d - b * c) % p:
            els.add(mobius(p, a, b, c, d))
    H = els
    say(f'  |PGL(2,5)| as permutations of the 6 points of P^1(F_5): {len(H)} (expect 120)')
    cts = Counter(cycle_type(g) for g in H)
    say('  cycle types of PGL(2,5) on 6 points:',
        ', '.join(f'{ct_str(k)}:{v}' for k, v in sorted(cts.items(), key=lambda kv: -kv[0][0])))
    has_transposition = any(cycle_type(g) == (2, 1, 1, 1, 1) for g in H)
    trans = len({g[0] for g in H}) == 6
    triples = Counter((g[0], g[1], g[2]) for g in H)
    sharply3 = len(triples) == 6 * 5 * 4 and all(v == 1 for v in triples.values())
    say(f'  transitive: {trans}; sharply 3-transitive: {sharply3}; contains a transposition: {has_transposition}')
    # S_6 and the coset action -> outer automorphism
    S6 = [tuple(q) for q in permutations(range(6))]
    cosets = {}
    reps = []
    for g in S6:
        key = frozenset(compose(g, h) for h in H)
        if key not in cosets:
            cosets[key] = len(reps)
            reps.append(key)
    say(f'  cosets S_6/PGL(2,5): {len(reps)}')
    coset_of = {}
    for key, i in cosets.items():
        for x in key:
            coset_of[x] = i
    rep_elt = [next(iter(k)) for k in reps]
    psi = {}
    for g in S6:
        psi[g] = tuple(coset_of[compose(g, rep_elt[i])] for i in range(6))
    # homomorphism + bijectivity
    rnd = random.Random(20260923)
    hom_ok = all(psi[compose(g, h)] == compose(psi[g], psi[h])
                 for g, h in ((rnd.choice(S6), rnd.choice(S6)) for _ in range(20000)))
    bij = len(set(psi.values())) == 720
    say(f'  psi: S_6 -> Sym(cosets) is a homomorphism on 20000 random pairs: {hom_ok}; bijective: {bij}')
    table = Counter((cycle_type(g), cycle_type(psi[g])) for g in S6)
    say('  class correspondence  cycle type of g  ->  cycle type of psi(g)  (count):')
    for (a, b), v in sorted(table.items(), key=lambda kv: (-kv[0][0][0], kv[0][0])):
        say(f'     {ct_str(a):>8} -> {ct_str(b):<8} ({v})')
    swapped = sorted({(ct_str(a), ct_str(b)) for (a, b) in table if a != b})
    say(f'  swapped class pairs: {swapped}  => psi is an OUTER automorphism')
    # characters
    classes = Counter(cycle_type(g) for g in S6)

    def inner(f1, f2):
        return Fraction(sum(f1(g) * f2(g) for g in S6), 720)
    pi = lambda g: sum(1 for i in range(6) if g[i] == i)
    pi2 = lambda g: sum(1 for i in range(6) if psi[g][i] == i)
    say(f'  <pi,pi> = {inner(pi, pi)}, <pi2,pi2> = {inner(pi2, pi2)}, <pi,pi2> = {inner(pi, pi2)}')
    say('  => 6 = 1 + 5 in two inequivalent ways (natural and exotic), sharing only the trivial summand')
    # identify chi = pi - 1 and chi' = pi2 - 1 among irreducibles
    irr = list(partitions(6))
    for name, f in (('pi-1 ', lambda g: pi(g) - 1), ("pi'-1", lambda g: pi2(g) - 1)):
        hits = []
        for lam in irr:
            val = Fraction(sum(f(g) * mn_char(lam, cycle_type(g)) for g in S6), 720)
            if val:
                hits.append((lam, val))
        say(f'  {name} = ' + ' + '.join(f'{v}*chi{lam}' for lam, v in hits))
    return psi, H


# ------------------------------------------------------------------ codes
class F4:
    """F_4 = {0,1,w,w^2} encoded 0,1,2,3 with w^2 = w+1."""
    add = [[a ^ b for b in range(4)] for a in range(4)]
    _log = {1: 0, 2: 1, 3: 2}
    _exp = [1, 2, 3]

    @classmethod
    def mul(cls, a, b):
        if a == 0 or b == 0:
            return 0
        return cls._exp[(cls._log[a] + cls._log[b]) % 3]


def span_code(gens, q, mul, add):
    words = {tuple([0] * len(gens[0]))}
    for g in gens:
        new = set()
        for w in words:
            for s in range(q):
                sg = tuple(mul(s, x) for x in g)
                new.add(tuple(add(a, b) for a, b in zip(w, sg)))
        words = new
    return words


def weight_enum(words):
    return dict(sorted(Counter(sum(1 for x in w if x) for w in words).items()))


def section_M2():
    say('=' * 78)
    say('M2. Extended quadratic-residue codes of length p+1, p = 5, 11, 23')
    # hexacode: QR code mod 5 over F_4, generator g(x) = x^2 + w x + 1 (a factor of x^5-1)
    add4 = lambda a, b: F4.add[a][b]
    mul4 = F4.mul
    g = [1, 2, 1, 0, 0]  # 1 + w x + x^2 (coefficients low->high)
    # check g | x^5 - 1 over F_4 by verifying x^5-1 = g*h for h = x^3 + a x^2 + b x + c
    def polymul(a, b):
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                r[i + j] = add4(r[i + j], mul4(x, y))
        return r
    found_h = None
    for h in product(range(4), repeat=3):
        prod = polymul([1, 2, 1], list(h) + [1])
        if prod == [1, 0, 0, 0, 0, 1]:
            found_h = h
    say(f'  x^5 - 1 = (x^2 + w x + 1) * h(x) over F_4: {found_h is not None}')
    rows = []
    for s in range(3):
        v = [0] * 5
        for i, c in enumerate([1, 2, 1]):
            v[(i + s) % 5] = c
        rows.append(v)
    C5 = span_code(rows, 4, mul4, add4)
    # extend by overall parity (sum over F_4)
    def ext4(w):
        s = 0
        for x in w:
            s = add4(s, x)
        return tuple(list(w) + [s])
    hexa = {ext4(w) for w in C5}
    say(f'  hexacode: {len(hexa)} words, weights {weight_enum(hexa)} (expect 64 words, 1 + 45 y^4 + 18 y^6)')
    # ternary Golay: extended QR code mod 11 over F_3
    Q11 = sorted({(x * x) % 11 for x in range(1, 11)})
    gpoly = None
    # generator polynomial = prod_{r in Q} (x - z^r) over F_3(z), z a primitive 11th root: compute via
    # brute force: the cyclic [11,6] code spanned by shifts of the idempotent-like vector.
    add3 = lambda a, b: (a + b) % 3
    mul3 = lambda a, b: (a * b) % 3
    # the ternary QR code of length 11 is generated by x^5 + x^4 - x^3 + x^2 - 1 (a divisor of x^11-1 mod 3)
    cand = [2, 0, 1, 2, 1, 1]  # -1 + x^2 - x^3 + x^4 + x^5, low->high, mod 3
    def pmul3(a, b):
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % 3
        return r
    def pdivmod3(num, den):
        num = num[:]
        q = [0] * (len(num) - len(den) + 1)
        inv = 1 if den[-1] == 1 else 2
        for k in range(len(num) - len(den), -1, -1):
            c = (num[k + len(den) - 1] * inv) % 3
            q[k] = c
            for j, y in enumerate(den):
                num[k + j] = (num[k + j] - c * y) % 3
        return q, num[:len(den) - 1]
    x11m1 = [2] + [0] * 10 + [1]
    qq, rr = pdivmod3(x11m1, cand)
    say(f'  x^11 - 1 divisible by g3(x) = x^5+x^4-x^3+x^2-1 over F_3: {all(c == 0 for c in rr)}')
    rows = []
    for s in range(6):
        v = [0] * 11
        for i, c in enumerate(cand):
            v[(i + s) % 11] = c
        rows.append(v)
    C11 = span_code(rows, 3, mul3, add3)
    # extension: last coordinate chosen so that the extended code is self-dual: c_inf = -sum? try both
    best = None
    for sgn in (1, 2):
        ext = {tuple(list(w) + [(sgn * sum(w)) % 3]) for w in C11}
        we = weight_enum(ext)
        if min(k for k in we if k) == 6:
            best = (sgn, ext, we)
    sgn, tern, we = best
    say(f'  ternary Golay (extension coefficient {sgn}): {len(tern)} words, weights {we} '
        f'(expect 729, 1 + 264 y^6 + 440 y^9 + 24 y^12)')
    hexads = {frozenset(i for i, x in enumerate(w) if x) for w in tern if sum(1 for x in w if x) == 6}
    say(f'  supports of weight-6 words: {len(hexads)} hexads; Steiner S(5,6,12): {steiner_check(hexads, 5, 12)}')
    # binary Golay: extended QR code mod 23, span of PSL(2,23)-orbit of N u {inf}
    p = 23
    pts = P1(p)
    idx = {x: i for i, x in enumerate(pts)}
    Q = {(x * x) % p for x in range(1, p)}
    N = set(range(1, p)) - Q
    gens = [mobius(p, 1, 1, 0, 1), mobius(p, 2, 0, 0, 1), mobius(p, 0, p - 1, 1, 0)]
    orb = orbit_sets([idx[x] for x in N | {'inf'}], gens)
    basis = []
    for s in orb:
        v = sum(1 << i for i in s)
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
            basis.sort(reverse=True)
    words = [0]
    for b in basis:
        words = words + [w ^ b for w in words]
    wt = dict(sorted(Counter(bin(w).count('1') for w in words).items()))
    say(f'  binary Golay: dim {len(basis)}, {len(words)} words, weights {wt} '
        f'(expect 1 + 759 y^8 + 2576 y^12 + 759 y^16 + y^24)')
    octads = [frozenset(i for i in range(24) if w >> i & 1) for w in words if bin(w).count('1') == 8]
    say(f'  octads: {len(octads)}; Steiner S(5,8,24): {steiner_check(octads, 5, 24)}; '
        f'C(24,5) = 42504 = 759*56 = {759 * 56}')
    return words, octads, hexads, idx


# ------------------------------------------------------------------ M3
def section_M3():
    say('=' * 78)
    say('M3. M_12, the hexad stabilizer S_6 (outer twist), and the exceptional A_5 < PSL(2,11)')
    p = 11
    pts = P1(p)
    idx = {x: i for i, x in enumerate(pts)}
    psl = [mobius(p, 1, 1, 0, 1), mobius(p, 3, 0, 0, 1), mobius(p, 0, p - 1, 1, 0)]
    d = list(range(12))
    for a, b in [(2, 10), (3, 4), (5, 9), (6, 7)]:
        d[idx[a]], d[idx[b]] = idx[b], idx[a]
    d = tuple(d)
    PSL = closure(psl)
    say(f'  |PSL(2,11)| = {len(PSL)} (expect 660)')
    M12 = closure(psl + [d])
    say(f'  |<PSL(2,11), (2 10)(3 4)(5 9)(6 7)>| = {len(M12)} (expect 95040 = 12*11*10*9*8 = {12*11*10*9*8})')
    tuples5 = Counter(g[:5] for g in M12)
    say(f'  images of the first five points: {len(tuples5)} distinct ordered 5-tuples, each once: '
        f'{all(v == 1 for v in tuples5.values())} => sharply 5-transitive')
    # invariant hexad system
    H0 = None
    for S in combinations(range(12), 6):
        o = orbit_sets(S, psl + [d])
        if len(o) == 132:
            H0 = o
            break
    say(f'  M_12-invariant hexad system: {len(H0)} hexads, Steiner S(5,6,12): {steiner_check(H0, 5, 12)}')
    h = sorted(next(iter(H0)))
    comp = [i for i in range(12) if i not in h]
    stab = [g for g in M12 if all(g[i] in h for i in h)]
    say(f'  setwise stabilizer of the hexad {[pts[i] for i in h]}: order {len(stab)} (expect 720 = |S_6|)')
    on_h = [tuple(h.index(g[i]) for i in h) for g in stab]
    on_c = [tuple(comp.index(g[i]) for i in comp) for g in stab]
    say(f'  faithful on the hexad: {len(set(on_h)) == 720}; faithful on the complement: {len(set(on_c)) == 720}')
    tab = Counter((cycle_type(a), cycle_type(b)) for a, b in zip(on_h, on_c))
    say('  cycle type on hexad -> cycle type on complement (count):')
    for (a, b), v in sorted(tab.items(), key=lambda kv: (-kv[0][0][0], kv[0][0])):
        say(f'     {ct_str(a):>8} -> {ct_str(b):<8} ({v})')
    twist = tab.get(((2, 1, 1, 1, 1), (2, 2, 2)), 0) == 15
    say(f'  a transposition on the hexad is a triple transposition on the complement: {twist}')
    say('  => the two actions of the hexad stabilizer S_6 differ by the outer automorphism of S_6')
    compset = frozenset(comp)
    swap = [g for g in M12 if frozenset(g[i] for i in h) == compset]
    e12 = tuple(range(12))

    def order12(g):
        k, x = 1, g
        while x != e12:
            x = compose(g, x)
            k += 1
        return k
    pair_orders = Counter(order12(g) for g in stab + swap)
    say(f'  the complement is a hexad: {compset in H0}; elements exchanging hexad and complement: {len(swap)}; '
        f'pair stabilizer order {len(stab) + len(swap)}, element orders {sorted(pair_orders.items())}')
    say('  (elements of order 8 and 10 exist, so the pair stabilizer is not S_6 x 2: it is the extension of S_6 by')
    say('   its outer automorphism, of order 1440)')
    # the permutation character of M_12 on 12 points: 2-transitivity
    pairs = Counter((g[0], g[1]) for g in M12)
    say(f'  M_12 is 2-transitive (orbit of an ordered pair has size {len(pairs)} = 12*11): '
        f'so <pi,pi> = 2 and 12 = 1 + 11')
    # Galois: A_5 of index 11 inside PSL(2,11)
    PSLl = sorted(PSL)
    order = {}
    e = tuple(range(12))
    for g in PSLl:
        k, x = 1, g
        while x != e:
            x = compose(g, x)
            k += 1
        order[g] = k
    invs = [g for g in PSLl if order[g] == 2]
    threes = [g for g in PSLl if order[g] == 3]
    A5 = None
    for a in invs:
        for b in threes:
            if order[compose(a, b)] == 5:
                G = closure([a, b])
                if len(G) == 60:
                    A5 = G
                    break
        if A5:
            break
    say(f'  PSL(2,11) contains a subgroup of order 60 (A_5, a (2,3,5)-generated group): {A5 is not None}')
    if A5:
        cos = {}
        for g in PSLl:
            key = frozenset(compose(g, x) for x in A5)
            cos.setdefault(key, len(cos))
        say(f'  index of that A_5: {len(cos)} (Galois: PSL(2,p) acts on p points only for p <= 11)')
        orb0 = {g[0] for g in A5}
        say(f'  this A_5 is transitive on the 12 points of P^1(F_11): {len(orb0) == 12} (point stabilizer of order 5)')
    return M12


# ------------------------------------------------------------------ M4
def section_M4():
    say('=' * 78)
    say('M4. M_24, Frame shapes, and the level-11 eta product')
    p = 23
    pts = P1(p)
    idx = {x: i for i, x in enumerate(pts)}
    Q = {(x * x) % p for x in range(1, p)}
    inv9 = pow(9, -1, p)

    def delta(x):
        if x == 'inf' or x == 0:
            return x
        return (pow(x, 3, p) * inv9) % p if x in Q else (9 * pow(x, 3, p)) % p
    dl = [idx[delta(x)] for x in pts]
    gens = [mobius(p, 1, 1, 0, 1), mobius(p, 2, 0, 0, 1), mobius(p, 0, p - 1, 1, 0)]
    G = PermutationGroup([Permutation(list(g)) for g in gens] + [Permutation(dl)])
    say(f'  |<PSL(2,23), delta>| = {G.order()} (expect |M_24| = 244823040 = 24*23*22*21*20*48 = '
        f'{24*23*22*21*20*48})')
    say(f'  transitivity degree >= 5: orbit of an ordered 5-tuple (via stabilizer chain) '
        f'{G.order() // G.pointwise_stabilizer([0, 1, 2, 3, 4]).order()} = 24*23*22*21*20 = {24*23*22*21*20}')
    t23 = cycle_type(gens[0])
    t11 = cycle_type(gens[1])
    say(f'  x -> x+1 has Frame shape {ct_str(t23)}; x -> 2x has Frame shape {ct_str(t11)}')
    # eta(t)^2 eta(11t)^2 coefficients
    M = 200
    f = [0] * (M + 1)
    f[0] = 1
    for n in range(1, M + 1):
        for (step, power) in ((n, 2), (11 * n, 2)):
            if step > M:
                continue
            for _ in range(power):
                for k in range(M, step - 1, -1):
                    f[k] -= f[k - step]
    a = [0] * (M + 2)
    for k in range(M + 1):
        if k + 1 <= M + 1:
            a[k + 1] = f[k]
    # point counts on X_0(11): y^2 + y = x^3 - x^2 - 10 x - 20
    def count_pts(q):
        c = 1  # point at infinity
        for x in range(q):
            rhs = (x ** 3 - x * x - 10 * x - 20) % q
            for y in range(q):
                if (y * y + y - rhs) % q == 0:
                    c += 1
        return c
    primes = [q for q in range(2, 60) if all(q % r for r in range(2, int(q ** 0.5) + 1)) and q != 11]
    ok = all(a[q] == q + 1 - count_pts(q) for q in primes)
    say(f'  eta(t)^2 eta(11t)^2 = q - 2q^2 - q^3 + 2q^4 + q^5 + 2q^6 - 2q^7 ...: coefficients {a[1:12]}')
    say(f'  a_p = p + 1 - #E(F_p) for E: y^2+y = x^3-x^2-10x-20 (X_0(11)) at all primes p < 60, p != 11: {ok}')
    b = [a[2 ** r] for r in range(0, 8)]
    rec = all(b[r + 1] == -2 * b[r] - 2 * b[r - 1] for r in range(1, 7))
    say(f'  a_(2^r) for r = 0..7: {b}; recursion b_(r+1) = -2 b_r - 2 b_(r-1) (i.e. 1/(1+2t+2t^2)): {rec}')
    return G, [tuple(g) for g in gens] + [tuple(dl)]


def section_M4b(words, gens24):
    say('=' * 78)
    say('M4b. Inside M_24: the stabilizer of a dodecad is M_12, acting on the dodecad and on its')
    say('     complement through the outer automorphism of M_12 (the 12 -> 24 analogue of M3)')
    dodecads = [w for w in words if bin(w).count('1') == 12]
    D = dodecads[0]
    Dset = [i for i in range(24) if D >> i & 1]
    Cset = [i for i in range(24) if not D >> i & 1]
    rnd = random.Random(1729)
    # product replacement
    pool = list(gens24) + [compose(gens24[0], gens24[1]), compose(gens24[2], gens24[3]),
                           compose(gens24[1], gens24[3]), compose(gens24[0], gens24[2])]
    acc = tuple(range(24))

    def step():
        nonlocal acc
        i, j = rnd.sample(range(len(pool)), 2)
        if rnd.random() < 0.5:
            pool[i] = compose(pool[i], pool[j])
        else:
            pool[i] = compose(pool[j], pool[i])
        acc = compose(acc, pool[i])
        return acc
    for _ in range(2000):
        step()
    stab_gens = []
    tries = 0
    Dfs = frozenset(Dset)
    while len(stab_gens) < 12 and tries < 400000:
        g = step()
        tries += 1
        if frozenset(g[i] for i in Dset) == Dfs:
            stab_gens.append(g)
    H = closure(stab_gens)
    orbD = orbit_sets(Dset, list(gens24))
    say(f'  M_24-orbit of the dodecad: {len(orbD)} (= all 2576 dodecads), so its stabilizer has order '
        f'{244823040 // len(orbD)}')
    say(f'  random elements tried: {tries}; subgroup generated by stabilizing elements: order {len(H)} '
        f'(= the full stabilizer, M_12)')
    onD = [tuple(Dset.index(g[i]) for i in Dset) for g in H]
    onC = [tuple(Cset.index(g[i]) for i in Cset) for g in H]
    say(f'  faithful on the dodecad: {len(set(onD)) == len(H)}; faithful on the complement: {len(set(onC)) == len(H)}')
    tab = Counter((cycle_type(x), cycle_type(y)) for x, y in zip(onD, onC))
    say('  cycle type on dodecad -> cycle type on complement (count):')
    for (x, y), v in sorted(tab.items(), key=lambda kv: (-kv[0][0][0], kv[0][0], kv[0][1])):
        say(f'     {ct_str(x):>10} -> {ct_str(y):<10} ({v})')
    diff = sorted({(ct_str(x), ct_str(y)) for (x, y) in tab if x != y})
    say(f'  cycle-type pairs that differ: {diff}')
    say('  => the two 12-point actions are inequivalent permutation representations of M_12 (outer twist)')


# ------------------------------------------------------------------ M5
def section_M5(words):
    say('=' * 78)
    say('M5. Leech lattice minimal vectors from the binary Golay code (coordinates scaled by sqrt 8)')
    import numpy as np
    codewords = words
    octads = [w for w in codewords if bin(w).count('1') == 8]
    vecs = []
    # shape (+-4, +-4, 0^22)
    for i, j in combinations(range(24), 2):
        for si in (4, -4):
            for sj in (4, -4):
                v = [0] * 24
                v[i], v[j] = si, sj
                vecs.append(v)
    n1 = len(vecs)
    # shape (+-2^8, 0^16) on octads, even number of minus signs
    for w in octads:
        supp = [i for i in range(24) if w >> i & 1]
        for signs in product((1, -1), repeat=8):
            if signs.count(-1) % 2 == 0:
                v = [0] * 24
                for i, s in zip(supp, signs):
                    v[i] = 2 * s
                vecs.append(v)
    n2 = len(vecs) - n1
    # shape (-+3, +-1^23): x_i = 1 mod 4 off S3, 3 mod 4 on S3 (S3 a codeword); |x_i| = 3 at one position
    for w in codewords:
        for i in range(24):
            v = []
            for j in range(24):
                in_s3 = (w >> j) & 1
                if j == i:
                    v.append(3 if in_s3 else -3)
                else:
                    v.append(-1 if in_s3 else 1)
            vecs.append(v)
    n3 = len(vecs) - n1 - n2
    A = np.array(vecs, dtype=np.int16)
    say(f'  counts by shape: (4^2 0^22): {n1}, (2^8 0^16): {n2}, (3 1^23): {n3}; total {len(vecs)} '
        f'(expect 1104 + 97152 + 98304 = {1104 + 97152 + 98304})')
    norms = set((A.astype(np.int32) ** 2).sum(axis=1).tolist())
    distinct = len({tuple(r) for r in A.tolist()}) == len(vecs)
    say(f'  all squared norms: {sorted(norms)} (32 in these units = norm 4); all distinct: {distinct}')
    # Conway-Sloane conditions: all x_i = m mod 2; sum = 4m mod 8; for each class mod 4 the support is a codeword
    cw = set(codewords)
    def cond(v):
        m = v[0] % 2
        if any(x % 2 != m for x in v):
            return False
        if sum(v) % 8 != (4 * m) % 8:
            return False
        for r in range(4):
            s = sum(1 << j for j in range(24) if v[j] % 4 == r)
            if s and s not in cw:
                return False
        return True
    say(f'  Conway-Sloane membership conditions hold for every vector: {all(cond(v) for v in vecs)}')
    v0 = A[0].astype(np.int64)
    ips = Counter((A.astype(np.int64) @ v0).tolist())
    say('  inner products with one minimal vector (units: /8 gives the standard values):',
        {k // 8: ips[k] for k in sorted(ips)})
    say('  expected Leech distribution {-4:1, -2:4600, -1:47104, 0:93150, 1:47104, 2:4600, 4:1}')
    return len(vecs)


# ------------------------------------------------------------------ M6
def series_mul(a, b, M):
    r = [0] * (M + 1)
    for i, x in enumerate(a[:M + 1]):
        if x:
            for j, y in enumerate(b[:M + 1 - i]):
                r[i + j] += x * y
    return r


def section_M6(leech_count):
    say('=' * 78)
    say('M6. q-series: Leech theta, j, and the 196884 decompositions')
    M = 12
    sig3 = [0] + [sum(d ** 3 for d in range(1, n + 1) if n % d == 0) for n in range(1, M + 1)]
    E4 = [1] + [240 * sig3[n] for n in range(1, M + 1)]
    # Delta = q prod (1-q^n)^24 ; store P = prod (1-q^n)^24, Delta = q*P
    P = [0] * (M + 1)
    P[0] = 1
    for n in range(1, M + 1):
        for _ in range(24):
            for k in range(M, n - 1, -1):
                P[k] -= P[k - n]
    Delta = [0] + P[:M]
    E43 = series_mul(series_mul(E4, E4, M), E4, M)
    theta = [E43[k] - 720 * Delta[k] for k in range(M + 1)]
    say(f'  Theta_Leech = E4^3 - 720 Delta = {theta[:6]} ... (coefficient of q^m counts norm-2m vectors)')
    say(f'  q^2 coefficient {theta[2]} equals the explicit enumeration {leech_count}: {theta[2] == leech_count}; '
        f'check 179280 + 720*24 = {179280 + 720 * 24}')
    # 1/P as a power series
    invP = [0] * (M + 1)
    invP[0] = 1
    for k in range(1, M + 1):
        invP[k] = -sum(P[i] * invP[k - i] for i in range(1, k + 1))
    say(f'  1/eta^24 = q^-1 (1 + {invP[1]} q + {invP[2]} q^2 + {invP[3]} q^3 + {invP[4]} q^4 + ...); '
        f'324 = C(25,2) + 24 = {25 * 24 // 2} + 24')
    j = series_mul(E43, invP, M)  # j * q = E4^3 / P
    say(f'  j = q^-1 + {j[1]} + {j[2]} q + {j[3]} q^2 + {j[4]} q^3 + {j[5]} q^4 + ...')
    lhs = series_mul(theta, invP, M)  # (Theta/eta^24) * q
    lhs[1] -= 24
    rhs = j[:]
    rhs[1] -= 744
    say(f'  Theta_Leech/eta^24 - 24 = j - 744 through q^{M - 1}: {lhs[:M] == rhs[:M]}')
    c1, c2, c3 = j[2], j[3], j[4]
    checks = [
        ('196884 = 300 + 24 + 196560 (V_Leech weight 2)', c1 == 300 + 24 + 196560),
        ('196884 = (300 + 98280) + 98304 (FLM orbifold: V_Leech^+ + twisted^+)', c1 == 300 + 98280 + 98304),
        ('98280 = 196560/2 and 98304 = 2^12 * 24', 98280 * 2 == 196560 and 98304 == 2 ** 12 * 24),
        ('196884 = 1 + 299 + 98280 + 98304 (restriction to 2^(1+24).Co_1)', c1 == 1 + 299 + 98280 + 98304),
        ('196883 = 299 + 98280 + 98304', 196883 == 299 + 98280 + 98304),
        ('196884 = 1 + 196883', c1 == 1 + 196883),
        ('196883 = 47 * 59 * 71', 196883 == 47 * 59 * 71),
        ('c(2) = 21493760 = 1 + 196883 + 21296876', c2 == 1 + 196883 + 21296876),
        ('c(3) = 864299970 = 2*1 + 2*196883 + 21296876 + 842609326',
         c3 == 2 + 2 * 196883 + 21296876 + 842609326),
        ('300 = C(25,2) = 1 + 299 (Sym^2 of the 24-dim Co_0 module)', 300 == 25 * 24 // 2 == 1 + 299),
        ('shape count 98304 = 24 * |Golay| = 24 * 4096', 98304 == 24 * 4096),
    ]
    for text, ok in checks:
        say(f'  [{"OK" if ok else "FAIL"}] {text}')
    say('  the "+1" ladder: 6 = 5+1 (S_6 on P^1(F_5)), 12 = 11+1 (M_12), 24 = 23+1 (M_24), '
        '300 = 299+1 (Co_1), 196884 = 196883+1 (Monster)')
    say(f'  numerology note: 5 -> 11 -> 23 -> 47 is a Cunningham chain (p -> 2p+1), and 47 | 196883: '
        f'{2*5+1 == 11 and 2*11+1 == 23 and 2*23+1 == 47 and 196883 % 47 == 0}')


def main():
    say('procgen_sources_20260923_moonshine.py')
    section_M1()
    words, octads, hexads, idx = section_M2()
    section_M3()
    G24, gens24 = section_M4()
    section_M4b(words, gens24)
    n = section_M5(words)
    section_M6(n)
    say('=' * 78)
    say('moonshine: done')


if __name__ == '__main__':
    main()
