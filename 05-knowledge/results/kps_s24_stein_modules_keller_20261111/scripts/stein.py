"""stein.py - exact (mod p) computation of twisted de Rham cohomology H^2(X, d - dW ^) for three coordinate rings:

  ring "A": C[x, y]                         (X = C^2)
  ring "B": C[x, 1/x, y]                    (X = C^* x C)
  ring "C": C[x, y, w] / (w (1 + x y) - 1)  (X = C^2 minus the hyperbola x y = -1)

In each case dx ^ dy trivialises the 2-forms, so
      H^2 = O(X) / M,   M = Im(d/dx - W_x) + Im(d/dy - W_y).
Elements are dicts {(i, j, k): coefficient mod P}; monomial x^i y^j w^k (k = 0 in rings A, B; i in Z in ring B).
Normal form in ring C: k = 0, or i j = 0 (rule x y w = 1 - w).

Truncation: V_N = span of basis monomials of weight <= N (weight |i| + j + k).  M_(N') = span of the generators
applied to basis monomials of weight <= N'.  c(N, N') = dim V_N - dim (M_(N') cap V_N) is computed by elimination
with pivots in the monomials of largest weight; it decreases in N' to the rank of V_N in H^2."""
import sys
from functools import lru_cache

P = 2147483629          # a prime below 2^31


class Ring:
    def __init__(self, kind):
        self.kind = kind

    # ----- monomial normal form -----
    @lru_cache(maxsize=None)
    def _nf(self, i, j, k):
        if self.kind != "C" or k == 0 or i == 0 or j == 0:
            return (((i, j, k), 1),)
        # x^i y^j w^k = x^(i-1) y^(j-1) (w^(k-1) - w^k)
        out = {}
        for (m, c) in self._nf(i - 1, j - 1, k - 1):
            out[m] = (out.get(m, 0) + c) % P
        for (m, c) in self._nf(i - 1, j - 1, k):
            out[m] = (out.get(m, 0) - c) % P
        return tuple((m, c) for m, c in out.items() if c)

    def mono(self, i, j, k=0, c=1):
        return {m: (c * cc) % P for (m, cc) in self._nf(i, j, k)}

    @staticmethod
    def add(a, b, s=1):
        out = dict(a)
        for m, c in b.items():
            v = (out.get(m, 0) + s * c) % P
            if v:
                out[m] = v
            else:
                out.pop(m, None)
        return out

    def mul(self, a, b):
        out = {}
        for (i, j, k), c in a.items():
            for (i2, j2, k2), c2 in b.items():
                cc = c * c2 % P
                for (m, c3) in self._nf(i + i2, j + j2, k + k2):
                    v = (out.get(m, 0) + cc * c3) % P
                    if v:
                        out[m] = v
                    else:
                        out.pop(m, None)
        return out

    def scal(self, a, s):
        s %= P
        return {m: c * s % P for m, c in a.items() if c * s % P}

    def dx(self, a):
        out = {}
        for (i, j, k), c in a.items():
            terms = []
            if i:
                terms.append(((i - 1, j, k), c * i))
            if k:
                terms.append(((i, j + 1, k + 1), -c * k))        # d/dx w = -y w^2
            for (m0, c0) in terms:
                for (m, c3) in self._nf(*m0):
                    v = (out.get(m, 0) + c0 * c3) % P
                    if v:
                        out[m] = v
                    else:
                        out.pop(m, None)
        return out

    def dy(self, a):
        out = {}
        for (i, j, k), c in a.items():
            terms = []
            if j:
                terms.append(((i, j - 1, k), c * j))
            if k:
                terms.append(((i + 1, j, k + 1), -c * k))        # d/dy w = -x w^2
            for (m0, c0) in terms:
                for (m, c3) in self._nf(*m0):
                    v = (out.get(m, 0) + c0 * c3) % P
                    if v:
                        out[m] = v
                    else:
                        out.pop(m, None)
        return out

    def power(self, a, n):
        r = self.mono(0, 0)
        for _ in range(n):
            r = self.mul(r, a)
        return r

    @staticmethod
    def weight(m):
        return abs(m[0]) + m[1] + m[2]

    def basis(self, N):
        out = []
        if self.kind == "A":
            for i in range(N + 1):
                for j in range(N + 1 - i):
                    out.append((i, j, 0))
        elif self.kind == "B":
            for i in range(-N, N + 1):
                for j in range(N + 1 - abs(i)):
                    out.append((i, j, 0))
        else:
            for i in range(N + 1):
                for j in range(N + 1 - i):
                    out.append((i, j, 0))
            for k in range(1, N + 1):
                for i in range(N + 1 - k):
                    out.append((i, 0, k))
                for j in range(1, N + 1 - k):
                    out.append((0, j, k))
        return out


def key(m):
    return (abs(m[0]) + m[1] + m[2], m)


class Echelon:
    """rows reduced so that each has a distinct leading monomial (largest key)."""
    def __init__(self):
        self.piv = {}

    def reduce(self, row):
        row = dict(row)
        while row:
            lead = max(row, key=key)
            if lead not in self.piv:
                return row, lead
            prow = self.piv[lead]
            f = row[lead] * pow(prow[lead], P - 2, P) % P
            for m, c in prow.items():
                v = (row.get(m, 0) - f * c) % P
                if v:
                    row[m] = v
                else:
                    row.pop(m, None)
        return row, None

    def insert(self, row):
        row, lead = self.reduce(row)
        if lead is not None:
            self.piv[lead] = row
            return True
        return False

    def normal_form(self, row):
        """full reduction (also non-leading terms) against the pivots."""
        row = dict(row)
        out = {}
        while row:
            lead = max(row, key=key)
            c = row.pop(lead)
            if lead in self.piv:
                prow = self.piv[lead]
                f = c * pow(prow[lead], P - 2, P) % P
                for m, cc in prow.items():
                    if m == lead:
                        continue
                    v = (row.get(m, 0) - f * cc) % P
                    if v:
                        row[m] = v
                    else:
                        row.pop(m, None)
            else:
                out[lead] = c
        return out


def twisted(R, W, Nprime):
    """echelon form of M_(N') = span{(d/dx - W_x) m, (d/dy - W_y) m : m basis monomial of weight <= N'}."""
    Wx, Wy = R.dx(W), R.dy(W)
    E = Echelon()
    for m in R.basis(Nprime):
        a = R.mono(*m)
        E.insert(R.add(R.dx(a), R.mul(Wx, a), -1))
        E.insert(R.add(R.dy(a), R.mul(Wy, a), -1))
    return E


def codim(R, E, N):
    """dim V_N - number of pivots of weight <= N."""
    return len(R.basis(N)) - sum(1 for m in E.piv if Ring.weight(m) <= N)


def table(R, W, Ns, Nps, label):
    print(label)
    for Np in Nps:
        E = twisted(R, W, Np)
        print("   N' = %2d:" % Np, " ".join("c(%d)=%d" % (N, codim(R, E, N)) for N in Ns if N <= Np))
        sys.stdout.flush()
    return E


if __name__ == "__main__":
    A = Ring("A")
    x, y, one = A.mono(1, 0), A.mono(0, 1), A.mono(0, 0)
    # textbook controls (expected: 1, 4, 1, 1, 1)
    table(A, A.mul(x, y), (2, 4, 6), (6, 10, 14), "W = x y  (Morse; expected 1)")
    table(A, A.add(A.power(x, 3), A.power(y, 3)), (3, 5, 7), (8, 12, 16), "W = x^3 + y^3  (expected (3-1)^2 = 4)")
    table(A, A.add(x, A.mul(A.mul(x, x), y)), (2, 4, 6), (8, 12, 16, 20), "W = x + x^2 y  (Broughton: no critical point, fibre C^*; expected 1)")
    table(A, A.mul(A.mul(x, x), y), (2, 4, 6), (8, 12, 16), "W = x^2 y  (non-isolated critical line, fibre C^*; expected 1)")
    f = A.add(x, A.mul(y, y))                 # automorphism (x + y^2, y)
    g = A.add(y, A.power(f, 3))               # then (f, y + f^3)
    table(A, A.mul(f, g), (2, 4, 6), (10, 14, 18, 22), "W = f g, (f, g) = (x + y^2, y + (x + y^2)^3)  (automorphism; expected 1)")
