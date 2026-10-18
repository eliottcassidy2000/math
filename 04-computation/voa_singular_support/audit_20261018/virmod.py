"""Independent implementation (audit) of the universal Virasoro vacuum module V^c.

PBW monomial: tuple (k1 >= k2 >= ... >= kr >= 2)  <->  L_{-k1} ... L_{-kr}|0>.
Coefficients live in a field given by `one` (flint.fmpq(1) for exact, flint.nmod(1,P) mod P).
Commutator: [L_m, L_n] = (m-n) L_{m+n} + c/12 (m^3-m) delta_{m+n,0}.
"""
import sys
import flint

sys.setrecursionlimit(100000)


def partitions_min2(D, maxpart=None):
    """Non-increasing tuples of integers >= 2 summing to D."""
    if maxpart is None:
        maxpart = D
    if D == 0:
        return [()]
    out = []
    for k in range(min(D, maxpart), 1, -1):
        for rest in partitions_min2(D - k, k):
            out.append((k,) + rest)
    return out


class VirVac:
    def __init__(self, c, one):
        self.one = one
        self.c = c * one
        self.cache = {}

    def act(self, n, mono):
        """L_n applied to the PBW monomial `mono`; returns dict mono -> coeff."""
        key = (n, mono)
        r = self.cache.get(key)
        if r is not None:
            return r
        one = self.one
        if n == 0:
            res = {mono: sum(mono) * one} if sum(mono) != 0 else {}
        elif not mono:
            res = {(-n,): one} if n <= -2 else {}
        else:
            k1 = mono[0]
            rest = mono[1:]
            if n <= -2 and -n >= k1:
                res = {(-n,) + mono: one}
            else:
                res = {}
                # L_n L_{-k1} rest = L_{-k1}(L_n rest) + [L_n, L_{-k1}] rest
                for m2, c2 in self.act(n, rest).items():
                    for m3, c3 in self.act(-k1, m2).items():
                        v = res.get(m3, 0) + c2 * c3
                        if v == 0:
                            res.pop(m3, None)
                        else:
                            res[m3] = v
                coef = n + k1
                if coef != 0:
                    for m3, c3 in self.act(n - k1, rest).items():
                        v = res.get(m3, 0) + coef * c3
                        if v == 0:
                            res.pop(m3, None)
                        else:
                            res[m3] = v
                if n == k1:
                    cc = self.c * (n ** 3 - n) / 12
                    if cc != 0:
                        v = res.get(rest, 0) + cc
                        if v == 0:
                            res.pop(rest, None)
                        else:
                            res[rest] = v
        self.cache[key] = res
        return res

    def actv(self, n, vec):
        out = {}
        for mono, c in vec.items():
            for m2, c2 in self.act(n, mono).items():
                v = out.get(m2, 0) + c * c2
                if v == 0:
                    out.pop(m2, None)
                else:
                    out[m2] = v
        return out


def central_charge(p, q):
    return flint.fmpq(1) - flint.fmpq(6 * (p - q) ** 2, p * q)


# ---------- commutative side: C[x_0, x_1, ...], x_j <-> L_{-j-2} -------------
# A polynomial monomial is stored as the SAME tuple (k1>=...>=kr) with x_{k-2} per part.

def poly_mul_mono(m1, m2):
    return tuple(sorted(m1 + m2, reverse=True))


def poly_mul(p1, p2):
    out = {}
    for a, ca in p1.items():
        for b, cb in p2.items():
            m = poly_mul_mono(a, b)
            v = out.get(m, 0) + ca * cb
            if v == 0:
                out.pop(m, None)
            else:
                out[m] = v
    return out


def poly_d(p):
    """Derivation d x_j = (j+1) x_{j+1};  in part notation k=j+2: d(k) = (k-1)*(k+1)."""
    out = {}
    for m, c in p.items():
        for i, k in enumerate(m):
            if i > 0 and m[i - 1] == k:
                continue  # handle repeated parts once with multiplicity
            mult = m.count(k)
            lst = list(m)
            lst.remove(k)
            lst.append(k + 1)
            nm = tuple(sorted(lst, reverse=True))
            v = out.get(nm, 0) + c * mult * (k - 1)
            if v == 0:
                out.pop(nm, None)
            else:
                out[nm] = v
    return out


def symbol(vec, r):
    """Length-r component of a V^c vector, read as a polynomial (same keys)."""
    return {m: c for m, c in vec.items() if len(m) == r}
