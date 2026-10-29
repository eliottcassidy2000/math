"""Audit of THM-4845 (B), self-converse tournaments (A002785).  Independent of gf_gaussian.py.

A. Labelled brute force: for every sigma in S_n (n <= 6) and one sigma per cycle type (n = 7, 8), count labelled
   tournaments T with sigma.T = T^op directly (bit vectors over pairs); compare with the claimed condition
   (at most one fixed point, every other cycle = 2 mod 4) and the claimed count 2^((m^T G m - m_1)/2).
B. Unlabelled brute force: every class from gentourng (n <= 9), canonical forms of T and T^op by nauty labelg.
C. The cycle-type (twisted Burnside) sum, n <= 30, against the OEIS A002785 data field (fetched 2026-10-29).
D. The Gaussian form E[(1 + W_1 x) exp(sum_{j = 2 mod 4} e^{X_j} x^j / j)] expanded as a product of exponential
   series, with the moment rule taken from the EXPLICIT realisation X_j = sqrt(ln 2) sum_{d | j} sqrt(phi(d)) xi_d,
   i.e. Cov(X_a, X_b) = ln2 * sum_{d | a, d | b} phi(d) (no gcd used), and E[exp(<m, X>)] = exp(Var/2).
"""
import itertools
import math
import subprocess
from fractions import Fraction
import numpy as np

OEIS_A002785 = [1, 1, 2, 2, 8, 12, 88, 176, 2752, 8784, 279968, 1492288, 95458560, 872687552, 111698291584,
                1787154671104, 457509297625088, 13013584213369088, 6662951988432581120,
                341143107490935724032, 349330527429800077778944, 32519496073514216703585280,
                66599927718689841715666229248]          # n = 1..23, OEIS data field, read 2026-10-29


def cycle_type(perm):
    n = len(perm)
    seen = [False] * n
    out = []
    for i in range(n):
        if not seen[i]:
            L = 0
            j = i
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                L += 1
            out.append(L)
    return tuple(sorted(out))


def claimed_valid(ct):
    return ct.count(1) <= 1 and all(c % 4 == 2 for c in ct if c != 1)


def claimed_count(ct):
    m = {}
    for c in ct:
        m[c] = m.get(c, 0) + 1
    q = sum(math.gcd(a, b) * m[a] * m[b] for a in m for b in m)
    e2 = q - m.get(1, 0)
    assert e2 % 2 == 0
    return 2 ** (e2 // 2)


def count_fixed(n, perm, Tbits, pairs, pidx):
    perm_e = np.array([pidx[tuple(sorted((perm[i], perm[j])))] for (i, j) in pairs], dtype=np.int64)
    flip = np.array([1 if perm[i] > perm[j] else 0 for (i, j) in pairs], dtype=np.uint8)
    # sigma.T has, at pair perm_e[e], the bit Tbits[e] ^ flip[e]; T^op has 1 - Tbits[perm_e[e]]
    lhs = Tbits ^ flip[None, :]
    rhs = 1 - Tbits[:, perm_e]
    return int(np.all(lhs == rhs, axis=1).sum())


print('A. labelled brute force')
for n in range(1, 9):
    pairs = list(itertools.combinations(range(n), 2))
    pidx = {e: k for k, e in enumerate(pairs)}
    m = len(pairs)
    if n == 8:
        # 2^28 tournaments are too many to hold; check only via the pair-orbit structure below (part A')
        break
    Tbits = ((np.arange(2 ** m)[:, None] >> np.arange(m)[None, :]) & 1).astype(np.uint8)
    total = 0
    bad = 0
    if n <= 6:
        for perm in itertools.permutations(range(n)):
            c = count_fixed(n, perm, Tbits, pairs, pidx)
            ct = cycle_type(perm)
            exp = claimed_count(ct) if claimed_valid(ct) else 0
            bad += (c != exp)
            total += c
        sc = Fraction(total, math.factorial(n))
    else:
        # one representative per cycle type, weighted by the class size n!/z
        sc = Fraction(0)
        def parts(n, maxp=None):
            if maxp is None:
                maxp = n
            if n == 0:
                yield ()
                return
            for p_ in range(min(n, maxp), 0, -1):
                for rest in parts(n - p_, p_):
                    yield (p_,) + rest
        for lam in parts(n):
            perm = list(range(n))
            pos = 0
            for L in lam:
                cyc = list(range(pos, pos + L))
                for k in range(L):
                    perm[cyc[k]] = cyc[(k + 1) % L]
                pos += L
            c = count_fixed(n, perm, Tbits, pairs, pidx)
            ct = cycle_type(perm)
            exp = claimed_count(ct) if claimed_valid(ct) else 0
            bad += (c != exp)
            mult = {}
            for L in lam:
                mult[L] = mult.get(L, 0) + 1
            z = 1
            for L, k in mult.items():
                z *= L ** k * math.factorial(k)
            sc += Fraction(c, z)
    print('   n=%d: SC(n) = %s ; sigmas (or cycle types) violating the claimed condition/count: %d' % (n, sc, bad))
    del Tbits

print("A'. pair-orbit check of the condition and count for n <= 14 (orbits computed by brute force)")
def orbit_check(n):
    bad = 0
    def parts(n, maxp=None):
        if maxp is None:
            maxp = n
        if n == 0:
            yield ()
            return
        for p_ in range(min(n, maxp), 0, -1):
            for rest in parts(n - p_, p_):
                yield (p_,) + rest
    for lam in parts(n):
        perm = list(range(n))
        pos = 0
        for L in lam:
            cyc = list(range(pos, pos + L))
            for k in range(L):
                perm[cyc[k]] = cyc[(k + 1) % L]
            pos += L
        # orbits of sigma on ordered pairs; an orientation is consistent iff following an arc around the orbit
        # with a flip at every step returns to itself
        seen = set()
        ok = True
        norb = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (i, j) in seen:
                    continue
                norb += 1
                # walk: state = ordered arc (x, y) meaning x -> y; sigma maps arc x->y to the arc sigma(y)->sigma(x)
                x, y = i, j
                steps = 0
                while True:
                    seen.add((min(x, y), max(x, y)))
                    x, y = perm[y], perm[x]
                    steps += 1
                    if {x, y} == {i, j}:
                        break
                if (x, y) != (i, j):
                    ok = False
        exp = claimed_count(lam) if claimed_valid(lam) else 0
        got = 2 ** norb if ok else 0
        bad += (got != exp)
    return bad
for n in range(1, 15):
    print('   n=%d: cycle types violating the claim: %d' % (n, orbit_check(n)))


print('B. unlabelled brute force (gentourng -z, nauty labelg canonical forms of T and T^op)')
def d6_decode(s):
    assert s[0] == '&'
    n = ord(s[1]) - 63
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        bits += [(v >> (5 - k)) & 1 for k in range(6)]
    A = np.array(bits[:n * n], dtype=np.int64).reshape(n, n)
    return A


def d6_encode(A):
    n = A.shape[0]
    bits = list(A.reshape(-1))
    while len(bits) % 6:
        bits.append(0)
    out = '&' + chr(n + 63)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = 2 * v + int(b)
        out += chr(v + 63)
    return out


def canon(lines):
    res = subprocess.run(['labelg', '-q'], input='\n'.join(lines) + '\n', capture_output=True, text=True,
                         check=True)
    return [ln for ln in res.stdout.split('\n') if ln.strip()]


for n in range(2, 10):
    raw = subprocess.run(['gentourng', '-q', '-z', str(n)], capture_output=True, text=True, check=True).stdout
    T = [ln.strip() for ln in raw.split('\n') if ln.strip()]
    mats = [d6_decode(s) for s in T]
    for A in mats[:50]:
        assert ((A + A.T) == (1 - np.eye(n, dtype=np.int64))).all()   # really a tournament
    Top = [d6_encode(A.T) for A in mats]
    assert d6_encode(mats[0]) == T[0]
    cT, cTop = canon(T), canon(Top)
    assert len(cT) == len(T) == len(cTop)
    distinct = len(set(cT))
    sc = sum(1 for a, b in zip(cT, cTop) if a == b)
    print('   n=%d: classes=%d (distinct canonical forms %d), converse closed: %s, self-converse = %d'
          % (n, len(T), distinct, set(cTop) == set(cT), sc))


print('C. twisted Burnside sum over cycle types vs OEIS A002785')
def partitions_mult(n):
    """Yield multiplicity dicts of partitions of n."""
    def rec(n, maxp):
        if n == 0:
            yield {}
            return
        for p_ in range(min(n, maxp), 0, -1):
            for rest in rec(n - p_, p_):
                d = dict(rest)
                d[p_] = d.get(p_, 0) + 1
                yield d
    yield from rec(n, n)


def sc_burnside(n):
    tot = Fraction(0)
    for m in partitions_mult(n):
        if m.get(1, 0) > 1 or any(p_ % 4 != 2 for p_ in m if p_ != 1):
            continue
        q = sum(math.gcd(a, b) * m[a] * m[b] for a in m for b in m)
        e2 = q - m.get(1, 0)
        z = 1
        for a, k in m.items():
            z *= a ** k * math.factorial(k)
        tot += Fraction(2 ** (e2 // 2), z)
    assert tot.denominator == 1
    return int(tot)


vals = [sc_burnside(n) for n in range(1, 31)]
print('   n=1..16 :', vals[:16])
print('   agrees with OEIS data n=1..23:', vals[:23] == OEIS_A002785,
      '; first disagreement:', next((k + 1 for k in range(23) if vals[k] != OEIS_A002785[k]), None))
print('   n=24..30:', vals[23:])


print('D. Gaussian form from the explicit divisor realisation')
def phi(d):
    return sum(1 for k in range(1, d + 1) if math.gcd(k, d) == 1)


def cov(a, b):          # in units of ln 2: sum over common divisors d of phi(d)
    return sum(phi(d) for d in range(1, min(a, b) + 1) if a % d == 0 and b % d == 0)


NMAX = 23
parts_J = [1] + [j for j in range(2, NMAX + 1) if j % 4 == 2]
# series: dict degree -> dict monomial(tuple of multiplicities over parts_J) -> Fraction
series = {0: {tuple([0] * len(parts_J)): Fraction(1)}}
for idx, j in enumerate(parts_J):
    new = {}
    if j == 1:
        factor = [(0, Fraction(1)), (1, Fraction(1))]          # (1 + W_1 x): exponents 0, 1, coefficient 1
    else:
        factor = [(k, Fraction(1, j ** k * math.factorial(k))) for k in range(0, NMAX // j + 1)]
    for deg, mons in series.items():
        for k, cf in factor:
            d2 = deg + j * k
            if d2 > NMAX:
                continue
            bucket = new.setdefault(d2, {})
            for mon, c in mons.items():
                mon2 = list(mon)
                mon2[idx] += k
                mon2 = tuple(mon2)
                bucket[mon2] = bucket.get(mon2, Fraction(0)) + c * cf
    series = new
C = [[cov(a, b) for b in parts_J] for a in parts_J]
assert all(C[x][y] == math.gcd(parts_J[x], parts_J[y]) for x in range(len(parts_J)) for y in range(len(parts_J)))
gauss = []
for n in range(0, NMAX + 1):
    tot = Fraction(0)
    for mon, c in series.get(n, {}).items():
        quad = sum(mon[x] * mon[y] * C[x][y] for x in range(len(parts_J)) for y in range(len(parts_J)))
        # E[W_1^{m1} prod e^{m_j X_j}] = 2^{-m1/2} * 2^{quad/2}
        e2 = quad - mon[0]
        assert e2 % 2 == 0, 'half-integral exponent'
        tot += c * 2 ** (e2 // 2)
    assert tot.denominator == 1
    gauss.append(int(tot))
print('   coefficients n=0..23:', gauss)
print('   n=1..23 agree with OEIS A002785:', gauss[1:] == OEIS_A002785, '; constant term', gauss[0])
