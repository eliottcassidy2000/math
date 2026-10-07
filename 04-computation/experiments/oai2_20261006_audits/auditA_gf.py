"""Minimal GF(p^n) arithmetic (audit A, independent of the repo code)."""
import itertools

def polymod(a, m, p):
    # a, m: lists of coefficients low->high; m monic
    a = a[:]
    dm = len(m) - 1
    while len(a) - 1 >= dm and any(a):
        if a[-1] % p == 0:
            a.pop(); continue
        c = a[-1] % p
        shift = len(a) - 1 - dm
        for i in range(len(m)):
            a[shift + i] = (a[shift + i] - c * m[i]) % p
        a.pop()
    while len(a) < dm: a.append(0)
    return [x % p for x in a[:dm]]

def polymul(a, b, p):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % p
    return r

def is_irreducible(m, p):
    # brute force: no factor of degree <= n/2 (small n only)
    n = len(m) - 1
    for d in range(1, n // 2 + 1):
        for tail in itertools.product(range(p), repeat=d):
            f = list(tail) + [1]
            # remainder of m mod f
            if not any(polymod(m, f, p)):
                return False
    return True

class GF:
    def __init__(self, p, n):
        self.p, self.n, self.q = p, n, p ** n
        if n == 1:
            self.m = [0, 1]
        else:
            for tail in itertools.product(range(p), repeat=n):
                m = list(tail) + [1]
                if m[0] != 0 and is_irreducible(m, p):
                    self.m = m; break
        q = self.q
        self.vec = [[(x // p ** k) % p for k in range(n)] for x in range(q)]
        idx = lambda v: sum(c * p ** k for k, c in enumerate(v))
        self.add = [[idx([(a + b) % p for a, b in zip(self.vec[x], self.vec[y])]) for y in range(q)] for x in range(q)]
        if n == 1:
            self.mul = [[(x * y) % p for y in range(q)] for x in range(q)]
        else:
            self.mul = [[idx(polymod(polymul(self.vec[x], self.vec[y], p), self.m, p)) for y in range(q)] for x in range(q)]
        self.neg = [next(y for y in range(q) if self.add[x][y] == 0) for x in range(q)]
        # sanity: field (every nonzero invertible)
        for x in range(1, q):
            assert any(self.mul[x][y] == 1 for y in range(1, q)), "not a field"

    def pw(self, x, e):
        r = 1
        b = x
        while e:
            if e & 1: r = self.mul[r][b]
            b = self.mul[b][b]; e >>= 1
        return r
