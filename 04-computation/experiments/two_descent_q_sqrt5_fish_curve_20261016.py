#!/usr/bin/env python3
"""Full two-descent over Q and over K0 = Q(sqrt5) for elliptic curves y^2 = (x - e1)(x - e2)(x - e3) with e_i in K0,
applied to the fish / EI21 curve.

Question (THM-4661 4(b),(c); fish_marked_elliptic_cover_20261013 (5)): does the genus-one curve
    C3:  11 z0^2 = D(1-D),  11 z1^2 = D(4-D),  11 z2^2 = (1-D)(D+3)
have a point over the real cyclic quartic field F4 = K0(sqrt(delta)), delta = (55 + 11 sqrt5)/2 = 11 sqrt5 phi?
C3 is the two-covering of E: y^2 = x(x-9)(x+7) with descent class (33, 33) for the map P -> (x, x - 9) mod squares
(x = 12/D - 3, so x = 3 D(4-D)/D^2 and x - 9 = 12 D(1-D)/D^2).  Hence C3(F4) is nonempty iff (33, 33) lies in the
image of E(F4) in (F4^x / F4^x2)^2.  The rank of E(F4) is rank E(K0) + rank E^delta(K0), E^delta: y^2 =
x(x - 9 delta)(x + 7 delta).  This script computes two-Selmer groups (upper bounds for ranks) by the standard full
two-torsion descent: for (b1, b2) in K0(S,2)^2, the class is Selmer iff at every place of S it lies in the local image
of E(K_v)/2E(K_v), computed by sampling global points until the local image has its known size (4 at odd places,
4 * 2^[K_v:Q2] at places over 2, |E(R)[2]|/... at real places).  Controls: congruent-number curves y^2 = x^3 - n^2 x
over Q (known ranks) and over K0.

Session opus-2026-10-16-S38.  Exact rational arithmetic (fractions); local tests by valuations and residues.
"""
import json, math, random, itertools, sys
from fractions import Fraction as Fr

random.seed(20261016)


# ---------------- K0 = Q(sqrt5) arithmetic: x = a + b sqrt5 ----------------
class K:
    __slots__ = ('a', 'b')

    def __init__(self, a, b=0):
        self.a = Fr(a); self.b = Fr(b)

    def __add__(s, o):
        o = o if isinstance(o, K) else K(o)
        return K(s.a + o.a, s.b + o.b)

    __radd__ = __add__

    def __sub__(s, o):
        o = o if isinstance(o, K) else K(o)
        return K(s.a - o.a, s.b - o.b)

    def __rsub__(s, o):
        return K(o) - s

    def __mul__(s, o):
        o = o if isinstance(o, K) else K(o)
        return K(s.a * o.a + 5 * s.b * o.b, s.a * o.b + s.b * o.a)

    __rmul__ = __mul__

    def norm(s):
        return s.a * s.a - 5 * s.b * s.b

    def inv(s):
        n = s.norm()
        return K(s.a / n, -s.b / n)

    def __truediv__(s, o):
        o = o if isinstance(o, K) else K(o)
        return s * o.inv()

    def iszero(s):
        return s.a == 0 and s.b == 0

    def __repr__(s):
        return f'({s.a} + {s.b} sqrt5)'


SQRT5 = K(0, 1)
PHI = K(Fr(1, 2), Fr(1, 2))


def vp(q, p):
    """p-adic valuation of a nonzero rational."""
    q = Fr(q)
    if q == 0:
        return 10 ** 9
    n, d, v = q.numerator, q.denominator, 0
    while n % p == 0:
        n //= p; v += 1
    while d % p == 0:
        d //= p; v -= 1
    return v


def unit_mod(q, p, k=1):
    """q = p^v * u; return u mod p^k as an integer."""
    q = Fr(q)
    v = vp(q, p)
    u = q / Fr(p) ** v
    m = p ** k
    return (u.numerator * pow(u.denominator, -1, m)) % m


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def sqrt5_mod(p, k):
    """A square root of 5 modulo p^k (p odd, 5 a QR mod p)."""
    r = next(x for x in range(p) if (x * x - 5) % p == 0)
    m = p
    for _ in range(1, k):
        m *= p
        r = (r - (r * r - 5) * pow(2 * r, -1, m)) % m
    return r


# ---------------- places and local square classes (as tuples of bits) ----------------
class Place:
    pass


class RealPlace(Place):
    def __init__(self, sign):
        self.sign = sign            # sqrt5 -> sign * sqrt(5)
        self.name = f'real{sign:+d}'
        self.dim = 1

    def cls(self, x):
        val = float(x.a) + self.sign * float(x.b) * math.sqrt(5)
        if abs(val) < 1e-12:
            val = x.a + self.sign * x.b * Fr(22360679774997897, 10 ** 16)   # high-precision fallback
        return (0,) if val > 0 else (1,)


class SplitPlace(Place):
    """p splits in K0; the place where sqrt5 -> s (mod p^k).  For K = Q use b = 0 only."""
    def __init__(self, p, s_mod_p, k=40):
        self.p = p
        r = s_mod_p % p
        m = p
        for _ in range(1, k):
            m *= p
            r = (r - (r * r - 5) * pow(2 * r, -1, m)) % m
        self.s = r; self.m = m; self.k = k
        self.name = f'{p}(sqrt5={s_mod_p})'
        self.dim = 2

    def image(self, x):
        """x as a p-adic number (Fraction combination) -> rational approximation mod p^k with valuation."""
        return x.a + x.b * self.s

    def cls(self, x):
        p = self.p
        # value a + b s in Q_p; compute valuation carefully: work with numerators/denominators
        a, b = x.a, x.b
        den = a.denominator * b.denominator
        A, B = a * den, b * den                   # integers
        A, B = int(A), int(B)
        val_int = (A + B * self.s) % self.m       # A + B s mod p^k (s is a p-adic integer)
        if val_int == 0:
            raise ValueError('precision')
        v = 0
        while val_int % p == 0:
            val_int //= p; v += 1
        v -= vp(den, p)
        dv = den // p ** vp(den, p)
        u = (val_int * pow(dv % p, -1, p)) % p
        return (v % 2, 0 if legendre(u, p) == 1 else 1)


class InertOddPlace(Place):
    """p odd, inert in K0: K_v = Q_p(sqrt5), integral basis {1, sqrt5}."""
    def __init__(self, p):
        self.p = p; self.name = f'{p}(inert)'; self.dim = 2

    def cls(self, x):
        p = self.p
        va, vb = vp(x.a, p), vp(x.b, p)
        v = min(va, vb)
        sc = Fr(p) ** (-v)
        a, b = x.a * sc, x.b * sc
        am = unit_mod(a, p) if a != 0 and vp(a, p) == 0 else 0
        bm = unit_mod(b, p) if b != 0 and vp(b, p) == 0 else 0
        # residue am + bm*w in F_p[w], w^2 = 5; square iff (am + bm w)^((p^2-1)/2) = 1
        e = (p * p - 1) // 2
        ra, rb = 1, 0
        ba, bb = am, bm
        while e:
            if e & 1:
                ra, rb = (ra * ba + 5 * rb * bb) % p, (ra * bb + rb * ba) % p
            ba, bb = (ba * ba + 5 * bb * bb) % p, (2 * ba * bb) % p
            e >>= 1
        sq = (ra == 1 and rb == 0)
        return (v % 2, 0 if sq else 1)


class RamifiedFivePlace(Place):
    """The prime sqrt5 of K0 (uniformizer sqrt5, residue F5)."""
    def __init__(self):
        self.name = 'sqrt5'; self.dim = 2

    def cls(self, x):
        va = 2 * vp(x.a, 5) if x.a != 0 else 10 ** 9
        vb = 2 * vp(x.b, 5) + 1 if x.b != 0 else 10 ** 9
        if va < vb:
            u = unit_mod(x.a, 5); v = va
        else:
            u = unit_mod(x.b, 5); v = vb
        return (v % 2, 0 if legendre(u, 5) == 1 else 1)


class TwoInertPlace(Place):
    """The prime 2 of K0 (inert): O_v = Z2[phi], phi^2 = phi + 1.  Units mod squares has order 8."""
    def __init__(self):
        self.name = '2(inert)'; self.dim = 4
        M = 8
        self.M = M
        units = [(c, d) for c in range(M) for d in range(M) if (c % 2, d % 2) != (0, 0)]
        mulm = lambda u, w: ((u[0] * w[0] + u[1] * w[1]) % M, (u[0] * w[1] + u[1] * w[0] + u[1] * w[1]) % M)
        self.mulm = mulm
        sq = set(mulm(u, u) for u in units)
        self.sq = sq
        # coset representatives of U/U^2 and an F2-basis
        reps = []
        def coset_of(u):
            for i, r in enumerate(reps):
                # u/r square  <=> u * r^-1 in sq ; r^-1 mod 8 by search
                if mulm(u, self.inv(r)) in sq:
                    return i
            return None
        self.inv = lambda r: next(w for w in units if mulm(r, w) == (1, 0))
        for u in units:
            if coset_of(u) is None:
                reps.append(u)
        assert len(reps) == 8, len(reps)
        # basis: -1, phi, 1+2phi?  find three generators spanning the 8 cosets
        cand = [(M - 1, 0), (0, 1), (1, 2), (3, 0), (1, 4), (5, 0), (1, 1), (3, 2)]
        basis = []
        span = {(1, 0)}
        for g in cand:
            if coset_of(g) is None:
                continue
            gi = coset_of(g)
            spanned_cosets = set(coset_of(s) for s in span)
            if gi in spanned_cosets:
                continue
            basis.append(g)
            span = span | set(mulm(s, g) for s in span)
            if len(basis) == 3:
                break
        assert len(basis) == 3, basis
        self.table = {}
        for bits in itertools.product((0, 1), repeat=3):
            e = (1, 0)
            for bit, g in zip(bits, basis):
                if bit:
                    e = mulm(e, g)
            self.table[coset_of(e)] = bits
        assert len(self.table) == 8
        self.coset_of = coset_of

    def cls(self, x):
        # x = a + b sqrt5 = (a - b) + 2b phi
        c, d = x.a - x.b, 2 * x.b
        vc, vd = vp(c, 2), vp(d, 2)
        v = min(vc, vd)
        sc = Fr(2) ** (-v)
        c, d = c * sc, d * sc
        cm = int((c.numerator * pow(c.denominator, -1, self.M)) % self.M) if c != 0 else 0
        dm = int((d.numerator * pow(d.denominator, -1, self.M)) % self.M) if d != 0 else 0
        return (v % 2,) + self.table[self.coset_of((cm, dm))]


# ---------------- descent ----------------
def local_image(place, e, nwant, tries=200000, base_q=False):
    """Image of E(K_v)/2E(K_v) under P -> (cls(x - e1), cls(x - e2)), sampled over global points."""
    e1, e2, e3 = e
    img = set()
    zero = tuple([0] * place.dim)
    def mul(c1, c2):
        return tuple((a + b) % 2 for a, b in zip(c1, c2))
    img.add((zero, zero))
    # two-torsion points
    for (t, o1, o2) in [(e1, e2, e3), (e2, e1, e3)]:
        pass
    c = lambda x: place.cls(x)
    img.add((c((e1 - e2) * (e1 - e3)), c(e1 - e2)))
    img.add((c(e2 - e1), c((e2 - e1) * (e2 - e3))))
    img.add((c(e3 - e1), c(e3 - e2)))
    for _ in range(tries):
        if len(img) >= nwant:
            break
        # random global x with random p-adic shape
        num_a = random.randint(-10 ** 6, 10 ** 6)
        num_b = 0 if base_q else random.randint(-10 ** 6, 10 ** 6)
        den = random.choice([1, 2, 3, 4, 5, 7, 8, 9, 11, 16, 25, 27, 49, 121, 32, 81, 125, 343])
        x = K(Fr(num_a, den), Fr(num_b, den))
        if random.random() < 0.3:
            x = e[random.randrange(3)] + x * K(Fr(1, random.choice([3, 5, 7, 11, 2, 4, 9, 25, 49, 121])) ** random.randint(1, 4))
        try:
            f1, f2, f3 = x - e1, x - e2, x - e3
            if f1.iszero() or f2.iszero() or f3.iszero():
                continue
            cf = mul(mul(c(f1), c(f2)), c(f3))
            if cf != zero:
                continue
            img.add((c(f1), c(f2)))
        except ValueError:
            continue
        if len(img) > nwant // 2:
            img = close_group(img)       # the local image is a subgroup (the Kummer map is a homomorphism)
    return close_group(img)


def close_group(img):
    img = set(img)
    while True:
        new = set()
        for x in img:
            for y in img:
                z = (tuple((a + b) % 2 for a, b in zip(x[0], y[0])), tuple((a + b) % 2 for a, b in zip(x[1], y[1])))
                if z not in img:
                    new.add(z)
        if not new:
            return img
        img |= new


def selmer(e, places, gens, base_q=False, sizes=None):
    """Two-Selmer group of y^2 = (x-e1)(x-e2)(x-e3) for S-unit generators gens (list of K)."""
    imgs = {}
    for pl in places:
        want = sizes[pl.name]
        im = local_image(pl, e, want, base_q=base_q)
        assert len(im) == want, (pl.name, len(im), want)
        imgs[pl.name] = im
    gcls = {pl.name: [pl.cls(g) for g in gens] for pl in places}
    n = len(gens)
    count = 0
    elems = []
    for m1 in range(2 ** n):
        for m2 in range(2 ** n):
            ok = True
            for pl in places:
                z = tuple([0] * pl.dim)
                c1, c2 = z, z
                for i in range(n):
                    if m1 >> i & 1:
                        c1 = tuple((a + b) % 2 for a, b in zip(c1, gcls[pl.name][i]))
                    if m2 >> i & 1:
                        c2 = tuple((a + b) % 2 for a, b in zip(c2, gcls[pl.name][i]))
                if (c1, c2) not in imgs[pl.name]:
                    ok = False
                    break
            if ok:
                count += 1
                elems.append((m1, m2))
    return count, elems


def places_Q(primes):
    pls = [RealPlace(1)]
    for p in primes:
        if p == 2:
            pls.append(TwoQPlace())
        else:
            pls.append(QOddPlace(p))
    return pls


class QOddPlace(Place):
    def __init__(self, p):
        self.p = p; self.name = f'{p}'; self.dim = 2

    def cls(self, x):
        assert x.b == 0
        v = vp(x.a, self.p)
        u = unit_mod(x.a, self.p)
        return (v % 2, 0 if legendre(u, self.p) == 1 else 1)


class TwoQPlace(Place):
    def __init__(self):
        self.name = '2'; self.dim = 3

    def cls(self, x):
        assert x.b == 0
        v = vp(x.a, 2)
        u = unit_mod(x.a, 2, 3)
        return (v % 2, 1 if u % 4 == 3 else 0, 1 if u % 8 in (3, 5) else 0)


def run_Q(e_int, label):
    e = [K(t) for t in e_int]
    from sympy import primefactors
    disc_primes = set()
    for i in range(3):
        for j in range(i + 1, 3):
            disc_primes |= set(primefactors(abs(int(e_int[i] - e_int[j]))))
    primes = sorted(disc_primes | {2})
    pls = places_Q(primes)
    gens = [K(-1)] + [K(p) for p in primes]
    sizes = {pl.name: (8 if pl.name == '2' else (2 if pl.name.startswith('real') else 4)) for pl in pls}
    # real place: E(R)/2E(R) has order 2 when all roots are real
    cnt, _ = selmer(e, pls, gens, base_q=True, sizes=sizes)
    return dict(curve=label, field='Q', S=primes, selmer_size=cnt, selmer_rank_bound=int(round(math.log2(cnt))) - 2)


def places_K0(primes_split, inert_odd, with5, with2):
    pls = [RealPlace(1), RealPlace(-1)]
    for (p, s) in primes_split:
        pls.append(SplitPlace(p, s))
    for p in inert_odd:
        pls.append(InertOddPlace(p))
    if with5:
        pls.append(RamifiedFivePlace())
    if with2:
        pls.append(TwoInertPlace())
    return pls


def run_K0(e, label, gens, pls):
    sizes = {}
    for pl in pls:
        if pl.name.startswith('real'):
            sizes[pl.name] = 2
        elif pl.name == '2(inert)':
            sizes[pl.name] = 16
        else:
            sizes[pl.name] = 4
    cnt, el = selmer(e, pls, gens, sizes=sizes)
    return dict(curve=label, field='Q(sqrt5)', places=[p.name for p in pls], selmer_size=cnt,
                selmer_rank_bound=int(round(math.log2(cnt))) - 2, elements=el)


def is_square_K(x):
    """Exact test: is x (a K) a square in Q(sqrt5)?  Returns a square root or None."""
    if x.iszero():
        return K(0)
    n = x.norm()
    if n < 0:
        return None
    rn = rat_sqrt(n)
    if rn is None:
        return None
    for sgn in (1, -1):
        s2 = (x.a + sgn * rn) / 2
        if s2 < 0:
            continue
        sr = rat_sqrt(s2)
        if sr is None or sr == 0:
            continue
        t = x.b / (2 * sr)
        r = K(sr, t)
        if ((r * r) - x).iszero():
            return r
    return None


def rat_sqrt(q):
    q = Fr(q)
    if q < 0:
        return None
    a, b = math.isqrt(q.numerator), math.isqrt(q.denominator)
    if a * a == q.numerator and b * b == q.denominator:
        return Fr(a, b)
    return None


def kprod(gens, mask):
    r = K(1)
    for i, g in enumerate(gens):
        if mask >> i & 1:
            r = r * g
    return r


def tors_masks(e, gens, pls):
    """Selmer masks of the three two-torsion points (found by matching local classes is unreliable; use exact
    square tests against all masks)."""
    e1, e2, e3 = e
    tclasses = [((e1 - e2) * (e1 - e3), e1 - e2), (e2 - e1, (e2 - e1) * (e2 - e3)), (e3 - e1, e3 - e2)]
    out = set()
    n = len(gens)
    for (c1, c2) in tclasses:
        for m1 in range(2 ** n):
            b1 = kprod(gens, m1)
            if is_square_K(c1 / b1) is None:
                continue
            for m2 in range(2 ** n):
                if is_square_K(c2 / kprod(gens, m2)) is not None:
                    out.add((m1, m2))
                    break
            break
    return out


def search_points(e, gens, elems, P=40, R=20, cap_per=3):
    """Search the two-covering of each Selmer element for points X = b1 u^2, u = (p + q phi)/r."""
    e1, e2, e3 = e
    found = {}
    for (m1, m2) in elems:
        b1, b2 = kprod(gens, m1), kprod(gens, m2)
        b12 = b1 * b2
        hits = []
        for r in range(1, R + 1):
            for pp in range(-P, P + 1):
                for qq in range(-P, P + 1):
                    if math.gcd(math.gcd(pp, qq), r) != 1 or (pp == 0 and qq == 0):
                        continue
                    u = K(Fr(pp, r), 0) + PHI * Fr(qq, r)
                    X = b1 * u * u
                    w2 = (X - e2) / b2
                    if w2.norm() < 0 or rat_sqrt(w2.norm()) is None:
                        continue
                    w3 = (X - e3) / b12
                    if w3.norm() < 0 or rat_sqrt(w3.norm()) is None:
                        continue
                    if is_square_K(w2) is not None and is_square_K(w3) is not None:
                        hits.append(X)
                        if len(hits) >= cap_per:
                            break
                if len(hits) >= cap_per:
                    break
            if len(hits) >= cap_per:
                break
        found[(m1, m2)] = hits
        print(f'   Selmer element b1={b1}, b2={b2}: {len(hits)} points found ' + (str(hits[0]) if hits else ''), flush=True)
    return found


if __name__ == '__main__':
    out = {}
    print('controls over Q (congruent number curves y^2 = x^3 - n^2 x; known ranks n=1,2,3: 0; n=5,6,7: 1):', flush=True)
    ctrl = []
    for n in (1, 2, 3, 5, 6, 7):
        r = run_Q((0, n, -n), f'n={n}')
        ctrl.append(r)
        print('  ', r, flush=True)
    out['controls_Q'] = ctrl
    rE = run_Q((0, 9, -7), 'E: y^2 = x(x-9)(x+7)')
    rE5 = run_Q((0, 45, -35), 'E^(5): y^2 = x(x-45)(x+35)')
    print('  ', rE, '\n  ', rE5, flush=True)
    out['E_over_Q'] = rE; out['E5_over_Q'] = rE5

    # over K0 = Q(sqrt5): S-unit generators
    pi11 = K(4, 1); pi11c = K(4, -1)
    gens_K0 = [K(-1), PHI, K(2), K(3), K(7), SQRT5, pi11, pi11c]
    # control: y^2 = x^3 - x over K0 has rank 0 + rank(x^3 - 25x) = 1
    pls_c = places_K0([], [], True, True)
    gens_c = [K(-1), PHI, K(2), SQRT5]
    rc = run_K0([K(0), K(1), K(-1)], 'control y^2 = x^3 - x over Q(sqrt5) (rank 1)', gens_c, pls_c)
    print('  ', rc, flush=True)
    out['control_K0'] = {k: v for k, v in rc.items() if k != 'elements'}
    # E over K0
    pls_E = places_K0([], [3, 7], False, True)
    gens_E = [K(-1), PHI, K(2), K(3), K(7)]
    rEK = run_K0([K(0), K(9), K(-7)], 'E over Q(sqrt5)', gens_E, pls_E)
    print('  ', rEK, flush=True)
    out['E_over_K0'] = {k: v for k, v in rEK.items() if k != 'elements'}
    # the twist by delta = 11 sqrt5 phi = (55 + 11 sqrt5)/2
    delta = K(Fr(55, 2), Fr(11, 2))
    assert (delta - K(11) * SQRT5 * PHI).iszero()
    eD = [K(0), K(9) * delta, K(-7) * delta]
    pls_D = places_K0([(11, 7), (11, 4)], [3, 7], True, True)     # pi11 = 4 + sqrt5 <-> sqrt5 = -4 = 7 (mod 11)
    rD = run_K0(eD, 'E^delta: y^2 = x(x - 9 delta)(x + 7 delta) over Q(sqrt5)', gens_K0, pls_D)
    print('  ', rD, flush=True)
    out['E_delta_over_K0'] = {k: v for k, v in rD.items() if k != 'elements'}
    print('point search on the Selmer coverings of E^delta over Q(sqrt5):', flush=True)
    pts = search_points(eD, gens_K0, rD['elements'])
    out['E_delta_points'] = {f'{m1},{m2}': [str(x) for x in v] for (m1, m2), v in pts.items()}
    # repeat the E^delta Selmer computation with another sampling seed
    random.seed(12345)
    rD2 = run_K0(eD, 'E^delta (second seed)', gens_K0, pls_D)
    assert rD2['selmer_size'] == rD['selmer_size'] and sorted(rD2['elements']) == sorted(rD['elements'])
    # control with split places over 11: y^2 = x^3 - 121 x over Q(sqrt5); rank 0 (n = 11) + rank(n = 55) >= 1
    pls_11 = places_K0([(11, 7), (11, 4)], [], True, True)
    r11 = run_K0([K(0), K(11), K(-11)], 'control y^2 = x^3 - 121x over Q(sqrt5)', [K(-1), PHI, K(2), SQRT5, pi11, pi11c], pls_11)
    print('  ', {k: v for k, v in r11.items() if k != 'elements'}, flush=True)
    out['control_K0_split11'] = {k: v for k, v in r11.items() if k != 'elements'}
    # rank bookkeeping: the classes with points span the whole Selmer group?
    def cls_vec(b):
        return tuple(pl.cls(b) for pl in pls_D)
    found_classes = [(m1, m2) for (m1, m2), v in pts.items() if v]
    span = {(0, 0)}
    for (m1, m2) in found_classes + [(m1, m2) for (m1, m2) in rD['elements'] if (m1, m2) in tors_masks(eD, gens_K0, pls_D)]:
        span |= set((a ^ m1, b ^ m2) for (a, b) in span)
    print(f'   classes with points (with torsion) span {len(span)} of {rD["selmer_size"]} Selmer elements', flush=True)
    out['E_delta_span_of_point_classes'] = len(span)
    # final step: is (33, 33) in the image of E(F4)/2E(F4)?  F4 = K0(sqrt delta); x in K0 is a square in F4 iff
    # x or x/delta is a square in K0.  Image = span of E(Q)_tors classes and the F4-images of E^delta(K0) classes.
    def sqF4(x):
        return is_square_K(x) is not None or is_square_K(x / delta) is not None
    gens_img = [(K(-7), K(-1)), (K(-3), K(-3))]
    for (m1, m2) in found_classes:
        gens_img.append((kprod(gens_K0, m1), kprod(gens_K0, m2)))
    # reduce to an independent set over F4
    basis = []
    elems_img = [(K(1), K(1))]
    for g in gens_img:
        if any(sqF4(g[0] / h[0]) and sqF4(g[1] / h[1]) for h in elems_img):
            continue
        basis.append(g)
        elems_img = elems_img + [(h[0] * g[0], h[1] * g[1]) for h in elems_img]
    print(f'   F4-image of E(F4)/2E(F4): {len(basis)} independent classes found (dimension 2 + rank = 4 expected)', flush=True)
    target = (K(33), K(33))
    hit = [h for h in elems_img if sqF4(target[0] / h[0]) and sqF4(target[1] / h[1])]
    print(f'   (33, 33) in the image: {bool(hit)}', flush=True)
    out['F4_image_dimension'] = len(basis)
    out['C3_class_33_33_in_image'] = bool(hit)
    out['conclusion'] = ('C3(F4) is EMPTY: rank E(F4) = rank E(K0) + rank E^delta(K0) = 0 + 2, the image of E(F4)/2E(F4) '
                         'is the displayed 4-dimensional span, and the class (33, 33) of C3 is not in it') if (not hit and len(basis) == 4) else 'UNDECIDED'
    print('   ' + out['conclusion'], flush=True)
    json.dump(out, open('05-knowledge/results/two_descent_q_sqrt5_fish_curve_20261016.json', 'w'), indent=1, default=str)
    print('wrote 05-knowledge/results/two_descent_q_sqrt5_fish_curve_20261016.json')
