# Exact arithmetic in a multiquadratic field L = Q(sqrt a_1,...,sqrt a_k) inside C (sqrt of negative = i*sqrt|a|).
from fractions import Fraction as Fr
import math, itertools, random

class MQ:
    def __init__(self, avec):
        self.a = list(avec); self.k = len(avec)
        self.n = 1 << self.k
        # numeric values of basis elements e_S
        self.num = []
        for S in range(self.n):
            z = complex(1, 0)
            for i in range(self.k):
                if S >> i & 1:
                    ai = self.a[i]
                    z *= (complex(0, math.sqrt(-ai)) if ai < 0 else math.sqrt(ai))
            self.num.append(z)
        # multiplication table e_S e_T = coef * e_{S^T}
        self.mt = {}
        for S in range(self.n):
            for T in range(self.n):
                c = 1
                for i in range(self.k):
                    if (S >> i & 1) and (T >> i & 1):
                        c *= self.a[i]
                self.mt[(S, T)] = c
        self.csign = []
        for S in range(self.n):
            s = 1
            for i in range(self.k):
                if (S >> i & 1) and self.a[i] < 0: s = -s
            self.csign.append(s)
    def elt(self, d):  # dict S->Fraction
        v = [Fr(0)]*self.n
        for S, c in d.items(): v[S] += Fr(c)
        return tuple(v)
    def add(self, x, y): return tuple(a+b for a, b in zip(x, y))
    def sub(self, x, y): return tuple(a-b for a, b in zip(x, y))
    def mul(self, x, y):
        r = [Fr(0)]*self.n
        for S, a in enumerate(x):
            if a == 0: continue
            for T, b in enumerate(y):
                if b == 0: continue
                r[S ^ T] += a*b*self.mt[(S, T)]
        return tuple(r)
    def conj(self, x): return tuple(s*a for s, a in zip(self.csign, x))
    def norm2(self, x): return self.mul(x, self.conj(x))   # |x|^2 as element (lies in Q if x in L ... in L+ generally)
    def is_one(self, x): return x[0] == 1 and all(c == 0 for c in x[1:])
    def inv(self, x):
        # x^{-1} via multiplying by all Galois conjugates: N = prod_{g} g(x)
        prod = self.elt({0: 1})
        for g in range(1, self.n):
            gx = tuple(a*(-1 if bin(S & g).count('1') % 2 else 1) for S, a in enumerate(x))
            prod = self.mul(prod, gx)
        N = self.mul(x, prod)
        assert all(c == 0 for c in N[1:]) and N[0] != 0
        return tuple(c / N[0] for c in prod)
    def numval(self, x): return sum(float(c)*z for c, z in zip(x, self.num))
    def one(self): return self.elt({0: 1})
