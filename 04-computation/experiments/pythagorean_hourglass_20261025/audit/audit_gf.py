# Audit helper (independent of the lane's code): finite fields GF(p^r) by log tables, the quadratic character, and the
# Paley tournament.  Elements are encoded as integers 0..q-1 via base-p digits (digit i = coefficient of x^i).
# The irreducible polynomial is searched from the LAST tail in lexicographic order (the lane searches from the first),
# verified by Rabin's test, and a primitive element is found by brute force.  Audit 2026-10-25.
import itertools
import numpy as np


def _polymulmod(a, b, f, p):
    r = len(f) - 1
    c = [0] * (2 * r)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    c[i + j] = (c[i + j] + x * y) % p
    for k in range(len(c) - 1, r - 1, -1):
        t = c[k]
        if t:
            for i in range(r + 1):
                c[k - r + i] = (c[k - r + i] - t * f[i]) % p
    return c[:r]


def _polypowmod(a, e, f, p):
    r = len(f) - 1
    res = [1] + [0] * (r - 1)
    while e:
        if e & 1:
            res = _polymulmod(res, a, f, p)
        a = _polymulmod(a, a, f, p)
        e >>= 1
    return res


def _polygcd_is_one(a, f, p):
    """gcd(a, f) == 1 over F_p (a, f coefficient lists low->high)"""
    def trim(u):
        u = [x % p for x in u]
        while u and u[-1] == 0:
            u.pop()
        return u
    a, b = trim(f), trim(a)
    while b:
        # a mod b
        inv = pow(b[-1], p - 2, p)
        while len(a) >= len(b) and a:
            t = a[-1] * inv % p
            s = len(a) - len(b)
            for i in range(len(b)):
                a[s + i] = (a[s + i] - t * b[i]) % p
            a = trim(a)
        a, b = b, a
    return len(a) == 1


def irreducible_rabin(f, p):
    r = len(f) - 1
    x = [0, 1] + [0] * (r - 2) if r >= 2 else [0]
    # x^(p^r) == x mod f
    y = x[:]
    for _ in range(r):
        y = _polypowmod(y, p, f, p)
    if y != x:
        return False
    for d in [d for d in range(2, r + 1) if r % d == 0 and all(d % e for e in range(2, d))]:
        z = x[:]
        for _ in range(r // d):
            z = _polypowmod(z, p, f, p)
        diff = [(zi - xi) % p for zi, xi in zip(z, x)]
        if not _polygcd_is_one(diff, f, p):
            return False
    return True


class GF:
    def __init__(self, p, r):
        self.p, self.r, self.q = p, r, p ** r
        q = self.q
        self.digits = np.array([[(x // p ** i) % p for i in range(r)] for x in range(q)], dtype=np.int64)
        self.pw = p ** np.arange(r, dtype=np.int64)
        if r == 1:
            self.f = None
            mul = lambda a, b: [(a[0] * b[0]) % p]
        else:
            f = None
            for tail in itertools.product(range(p - 1, -1, -1), repeat=r):   # reverse order: differs from the lane
                cand = list(tail) + [1]
                if cand[0] == 0:
                    continue
                if irreducible_rabin(cand, p):
                    f = cand
                    break
            self.f = f
            mul = lambda a, b: _polymulmod(a, b, f, p)
        enc = lambda c: int(sum(int(ci) * p ** i for i, ci in enumerate(c)))
        # primitive element
        order_primes = [s for s in range(2, q) if (q - 1) % s == 0 and all(s % t for t in range(2, int(s ** 0.5) + 1))]
        for g in range(2 if q > 2 else 1, q):
            gd = list(self.digits[g])
            ok = True
            for s in order_primes:
                e = (q - 1) // s
                # g^e
                res = [1] + [0] * (r - 1)
                base = gd[:]
                ee = e
                while ee:
                    if ee & 1:
                        res = mul(res, base)
                    base = mul(base, base)
                    ee >>= 1
                if enc(res) == 1:
                    ok = False
                    break
            if ok:
                self.g = g
                break
        self.exp = np.zeros(q - 1, dtype=np.int64)
        self.log = np.full(q, -1, dtype=np.int64)
        cur = [1] + [0] * (r - 1)
        gd = list(self.digits[self.g])
        for i in range(q - 1):
            e_ = enc(cur)
            self.exp[i] = e_
            self.log[e_] = i
            cur = mul(cur, gd)
        assert (self.log[1:] >= 0).all(), 'not primitive'
        self.chi = np.zeros(q, dtype=np.int64)
        self.chi[1:] = np.where(self.log[1:] % 2 == 0, 1, -1)

    def sub(self, y, x):
        """y - x elementwise (arrays of encodings, broadcast)"""
        d = (self.digits[np.asarray(y)] - self.digits[np.asarray(x)]) % self.p
        return (d * self.pw).sum(axis=-1)

    def add(self, y, x):
        d = (self.digits[np.asarray(y)] + self.digits[np.asarray(x)]) % self.p
        return (d * self.pw).sum(axis=-1)

    def mul(self, a, b):
        a = np.asarray(a); b = np.asarray(b)
        out = np.zeros(np.broadcast(a, b).shape, dtype=np.int64)
        nz = (a != 0) & (b != 0)
        la = np.broadcast_to(self.log[a], out.shape); lb = np.broadcast_to(self.log[b], out.shape)
        out[nz] = self.exp[(la[nz] + lb[nz]) % (self.q - 1)]
        return out

    def neg(self, x):
        d = (-self.digits[np.asarray(x)]) % self.p
        return (d * self.pw).sum(axis=-1)

    def inv(self, x):
        x = np.asarray(x)
        return self.exp[(-self.log[x]) % (self.q - 1)]

    def chi_of(self, z):
        return self.chi[np.asarray(z)]


def prime_field(p):
    return GF(p, 1)
