"""tw.py  -  independent engine for the audit of kind-pasteur-2026-11-11-S24 (written from scratch; the author's
scripts are not imported).

H^2(X, d - dW ^) = O(X) / ( Im(d/dx - W_x) + Im(d/dy - W_y) )   in the frame dx ^ dy, for

  RingA : C[x, y]                       basis x^i y^j                                   key (i, j)
  RingB : C[x, 1/x, y]                  basis x^i y^j, i in Z                           key (i, j)
  RingC : C[x, y, 1/q], q = 1 + x y     TORUS-WEIGHT basis  m(e, n) = x^e q^n (e >= 0),  y^(-e) q^n (e < 0),  n in Z

RingC's basis is not the author's normal form (x^i y^j w^k with k = 0 or i j = 0): it diagonalises the torus action
(x, y) -> (s x, y / s), under which q is invariant.  Elements are dicts {key: exact Python int}; all ring arithmetic
is exact over Z.  Only the linear algebra is done modulo a prime (two primes are used, both different from the
author's 2147483629).

Truncation.  A weight function wt on basis monomials is fixed (several are used).  M(N') = span of the two generators
applied to all basis monomials of weight <= N'.  The echelon form has pivots at the monomials that are largest for the
order (wt, key), so the rows whose pivot has weight <= N form a basis of M(N') cap V_N, V_N = span of monomials of
weight <= N, and
        c(N, N') = dim V_N - #{pivots of weight <= N} = dim V_N / (M(N') cap V_N)  >=  rank of V_N in H^2.
"""
import heapq
from math import comb

P1 = 2147483647          # 2^31 - 1
P2 = 1000000007


# ------------------------------------------------------------------ rings (exact integer arithmetic)
class RingA:
    name = "C[x,y]"

    def __init__(self, wx=1, wy=1):
        self.wx, self.wy = wx, wy

    def wt(self, m):
        return self.wx * m[0] + self.wy * m[1]

    def basis(self, N):
        return [(i, j) for i in range(N // self.wx + 1) for j in range((N - self.wx * i) // self.wy + 1)]

    def mulm(self, a, b):
        return {(a[0] + b[0], a[1] + b[1]): 1}

    def dxm(self, m):
        return {(m[0] - 1, m[1]): m[0]} if m[0] else {}

    def dym(self, m):
        return {(m[0], m[1] - 1): m[1]} if m[1] else {}


class RingB(RingA):
    name = "C[x,1/x,y]"

    def wt(self, m):
        return self.wx * abs(m[0]) + self.wy * m[1]

    def basis(self, N):
        out = []
        for j in range(N // self.wy + 1):
            r = (N - self.wy * j) // self.wx
            out += [(i, j) for i in range(-r, r + 1)]
        return out


class RingC:
    """C[x, y, 1/(1+xy)] in the torus-weight basis m(e, n)."""
    name = "C[x,y,1/(1+xy)]"

    def __init__(self, we=1, wpos=2, wneg=1):
        self.we, self.wpos, self.wneg = we, wpos, wneg

    def wt(self, m):
        e, n = m
        return self.we * abs(e) + (self.wpos * n if n >= 0 else -self.wneg * n)

    def basis(self, N):
        out = []
        for e in range(-(N // self.we), N // self.we + 1):
            r = N - self.we * abs(e)
            out += [(e, n) for n in range(-(r // self.wneg), r // self.wpos + 1)]
        return out

    def mulm(self, a, b):
        (e1, n1), (e2, n2) = a, b
        if e1 * e2 >= 0:
            return {(e1 + e2, n1 + n2): 1}
        k = min(abs(e1), abs(e2))                       # x^k y^k = (q - 1)^k
        return {(e1 + e2, n1 + n2 + j): comb(k, j) * (-1) ** (k - j) for j in range(k + 1)}

    def dxm(self, m):
        e, n = m
        if e >= 1:
            out = {(e - 1, n): e + n, (e - 1, n - 1): -n}
        else:                                           # d/dx q^n = n y q^(n-1);  d/dx (y^k q^n) = n y^(k+1) q^(n-1)
            out = {(e - 1, n - 1): n}
        return {k: c for k, c in out.items() if c}

    def dym(self, m):
        e, n = m
        if e <= -1:
            out = {(e + 1, n): -e + n, (e + 1, n - 1): -n}
        else:
            out = {(e + 1, n - 1): n}
        return {k: c for k, c in out.items() if c}


# ------------------------------------------------------------------ element arithmetic
def add(a, b, s=1):
    out = dict(a)
    for m, c in b.items():
        v = out.get(m, 0) + s * c
        if v:
            out[m] = v
        else:
            out.pop(m, None)
    return out


def scal(a, s):
    return {m: c * s for m, c in a.items()} if s else {}


def mul(R, a, b):
    out = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            cc = ca * cb
            for m, c in R.mulm(ma, mb).items():
                v = out.get(m, 0) + cc * c
                if v:
                    out[m] = v
                else:
                    out.pop(m, None)
    return out


def power(R, a, k, one):
    r = dict(one)
    for _ in range(k):
        r = mul(R, r, a)
    return r


def deriv(R, a, which):
    out = {}
    f = R.dxm if which == "x" else R.dym
    for m, c in a.items():
        for mm, cc in f(m).items():
            v = out.get(mm, 0) + c * cc
            if v:
                out[mm] = v
            else:
                out.pop(mm, None)
    return out


# ------------------------------------------------------------------ linear algebra modulo a prime
class Echelon:
    def __init__(self, R, p):
        self.R, self.p = R, p
        self.piv = {}            # sortkey -> (row dict {sortkey: coeff mod p}, normalised to leading coefficient 1)
        self.keyof = {}          # monomial -> sortkey
        self.monof = {}

    def sk(self, m):
        k = self.keyof.get(m)
        if k is None:
            # order: weight first, then the monomial itself; packed into one integer
            k = (self.R.wt(m) << 40) | ((m[0] + (1 << 18)) << 20) | (m[1] + (1 << 18))
            self.keyof[m] = k
            self.monof[k] = m
        return k

    def torow(self, el):
        p = self.p
        row = {}
        for m, c in el.items():
            c %= p
            if c:
                row[self.sk(m)] = c
        return row

    def _reduce(self, row, full):
        """eliminate pivots from the top; if not full, stop at the first non-pivot leading term.
        returns (row, lead) if not full, else the fully reduced row."""
        p, piv = self.p, self.piv
        heap = [-k for k in row]
        heapq.heapify(heap)
        out = {}
        while heap:
            k = -heapq.heappop(heap)
            c = row.get(k)
            if not c:
                row.pop(k, None)                    # pivot rows only contain smaller keys: k cannot come back
                continue
            pr = piv.get(k)
            if pr is None:
                if not full:
                    return row, k
                out[k] = c
                del row[k]
                continue
            for kk, cc in pr.items():
                if kk == k:
                    continue
                old = row.get(kk)
                v = ((old or 0) - c * cc) % p
                if v:
                    row[kk] = v
                    if old is None:
                        heapq.heappush(heap, -kk)
                elif old is not None:
                    row[kk] = 0                     # keep the key (already in the heap); zero is skipped on pop
            del row[k]
        if full:
            return out
        return {}, None

    def insert(self, el):
        row, lead = self._reduce(self.torow(el), False)
        if lead is None:
            return False
        inv = pow(row[lead], self.p - 2, self.p)
        self.piv[lead] = {k: c * inv % self.p for k, c in row.items() if c}
        return True

    def nf(self, el):
        """normal form: dict {monomial: coeff mod p} supported on non-pivot monomials."""
        out = self._reduce(self.torow(el), True)
        return {self.monof[k]: c for k, c in out.items() if c}

    def npiv_upto(self, N):
        return sum(1 for k in self.piv if (k >> 40) <= N)

    def nonpivots_upto(self, N):
        return sorted(m for m in self.R.basis(N) if self.sk(m) not in self.piv)


def build(R, W, Np, p):
    """echelon form of M(N')."""
    Wx, Wy = deriv(R, W, "x"), deriv(R, W, "y")
    E = Echelon(R, p)
    mons = sorted(R.basis(Np), key=lambda m: (R.wt(m), m))
    for m in mons:
        a = {m: 1}
        E.insert(add(deriv(R, a, "x"), mul(R, Wx, a), -1))
        E.insert(add(deriv(R, a, "y"), mul(R, Wy, a), -1))
    return E


def codim(R, E, N):
    return len(R.basis(N)) - E.npiv_upto(N)


def maxwt(R, el):
    return max(R.wt(m) for m in el) if el else 0


def table(R, W, Ns, Nps, p, label):
    print(label)
    E = None
    for Np in Nps:
        E = build(R, W, Np, p)
        print("     N' = %3d : %s   (%d pivots)" % (Np, "  ".join("c(%d)=%d" % (N, codim(R, E, N)) for N in Ns if N <= Np), len(E.piv)), flush=True)
    return E


# ------------------------------------------------------------------ self-test of the ring arithmetic against sympy
def selftest():
    import sympy as sp
    import random
    x, y = sp.symbols("x y")
    q = 1 + x * y
    random.seed(7)
    R = RingC()

    def to_sym(el):
        return sum(c * (x ** e if e >= 0 else y ** (-e)) * q ** n for (e, n), c in el.items())

    def rnd():
        return {(random.randint(-4, 4), random.randint(-3, 3)): random.randint(-5, 5) or 1 for _ in range(4)}

    ok = True
    for _ in range(40):
        a, b = rnd(), rnd()
        ok = ok and sp.simplify(to_sym(mul(R, a, b)) - to_sym(a) * to_sym(b)) == 0
        ok = ok and sp.simplify(to_sym(deriv(R, a, "x")) - sp.diff(to_sym(a), x)) == 0
        ok = ok and sp.simplify(to_sym(deriv(R, a, "y")) - sp.diff(to_sym(a), y)) == 0
    # uniqueness of the representation: a nonzero dict is a nonzero function (spot check through evaluation)
    for _ in range(40):
        a = rnd()
        val = to_sym(a).subs({x: sp.Rational(3, 7), y: sp.Rational(-5, 11)})
        val2 = to_sym(a).subs({x: sp.Rational(2, 5), y: sp.Rational(9, 4)})
        ok = ok and (val != 0 or val2 != 0)
    print("RingC self-test (products, d/dx, d/dy against sympy on 40 random pairs):", "OK" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    selftest()
