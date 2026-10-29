"""HYP-9446 (1)-(3) data.

(1) Pre-registered before running (2026-10-29): every LOCKED tournament cycle is a cycle of the rotational
    tournament R_q.  Locked: for some appended letter v and some prime p with ord_p(2) >= 3, p | D(wv),
    p not dividing 6 (2^v - 3), and every letter = v (mod ord_p(2)), so that every marking cancels p (THM-4849 (E)).
    Tested over every tournament cycle (every necklace containing a connection set) for primes q <= 23 and every
    v <= 2 * ord_p(2) + q, which covers every residue class of v mod ord_p(2).
(2) Interval-locking primes: the fraction of primes q < 6000 whose R_q interval cycle has no locked prime under (1),
    i.e. gcd(2^m - 1, 3^m - 1) = 1 with m = (q+1)/2 (a prime-restricted Ailon-Rudnick count).
(3) Reduced denominators of the Paley cycle and of its converse, for q = 3 (mod 4) prime <= 47.
"""
import itertools
from math import gcd
from sympy import factorint, n_order, primerange, divisors


def C_of(w):
    s, out, k = 0, 0, len(w)
    for j, v in enumerate(w):
        out += 3 ** (k - 1 - j) * 2 ** s
        s += v
    return out


def gaps(S, q):
    S = sorted(S)
    return tuple((S[(i + 1) % len(S)] - S[i]) % q or q for i in range(len(S)))


def canon(S, q):
    return min(tuple(sorted((x - t) % q for x in S)) for t in S)


def connection_sets(q):
    k = (q - 1) // 2
    for choice in itertools.product((0, 1), repeat=k):
        yield frozenset(i if c == 0 else q - i for i, c in zip(range(1, k + 1), choice))


print('(1) locked tournament cycles (pre-registered: all are R_q cycles)')
for q in primerange(5, 24):
    k = (q - 1) // 2
    I = frozenset(range(1, k + 1))
    R = {canon(frozenset((m * x) % q for x in I), q) for m in range(1, q)}
    neck = {}
    for S in connection_sets(q):
        neck.setdefault(canon(S, q), S)
    found = []
    for c, S in neck.items():
        w = gaps(sorted(c), q)
        diffs = [abs(a - b) for a in w for b in w]
        d = 0
        for x in diffs:
            d = gcd(d, x)
        if d == 0:
            d = 10 ** 9   # all letters equal: impossible here since k does not divide q
        # candidate primes: ord_p(2) divides d and ord_p(2) >= 3; d is at most q
        cands = set()
        for o in divisors(d) if d < 10 ** 9 else []:
            if o >= 3:
                cands.update(p for p in factorint(2 ** o - 1) if n_order(2, p) == o)
        hits = []
        for p in sorted(cands):
            o = n_order(2, p)
            for v in range(1, 2 * o + q + 1):
                if (v - w[0]) % o:
                    continue
                if (2 ** v - 3) % p == 0:
                    continue
                Dn = 2 ** (q + v) - 3 ** (k + 1)
                if Dn % p == 0:
                    D = 2 ** q - 3 ** k
                    rots = [w[i:] + w[:i] for i in range(k)]
                    xi = pow((2 ** v - 3) % p, -1, p)
                    assert all(C_of(r) * pow(D, -1, p) % p == xi for r in rots)
                    hits.append((p, v % o))
                    break
        if hits:
            found.append((w, c in R, sorted(set(hits))))
    nonR = [f for f in found if not f[1]]
    print('  q=%2d: %d tournament cycles, %d locked, of which non-R_q: %d  %s'
          % (q, len(neck), len(found), len(nonR), [(f[0], f[2]) for f in found][:6]))

print('\n(2) primes q < 6000 with gcd(2^m - 1, 3^m - 1) = 1, m = (q+1)/2 (no interval-locked prime)')
tot = none = 0
first = []
for q in primerange(3, 6000):
    m = (q + 1) // 2
    tot += 1
    if gcd(2 ** m - 1, 3 ** m - 1) == 1:
        none += 1
        if len(first) < 15:
            first.append(q)
print('  %d of %d primes; first: %s' % (none, tot, first))

print('\n(3) Paley cycle and converse: reduced denominators')
for q in primerange(7, 48):
    if q % 4 != 3:
        continue
    k = (q - 1) // 2
    D = 2 ** q - 3 ** k
    QR = frozenset((x * x) % q for x in range(1, q))
    NQR = frozenset(range(1, q)) - QR
    out = []
    for name, S in (('QR', QR), ('NQR', NQR)):
        w = gaps(sorted(S), q)
        g = gcd(C_of(w), D)
        out.append('%s: D/%d' % (name, g))
    print('  q=%2d  D_q = %d = %s ; %s' % (q, D, factorint(D), ' ; '.join(out)))
