# AUDIT a3 (independent auditor, 2026-10-21): replay the group-ring proof of THM-4746 Theorem 1 numerically.
# For each odd prime p <= 31: K = F_p (eps=+1) or F_{p^2} (eps=-1); g0 with g0^(p-1) = eps; write
# tau = (2^C(m,2) mod p)_m as sum_{a != 0} c_a E_{g0 a} (Vandermonde solve on m = 0..p-2) and CHECK it on m <= 4p;
# t = sum c_a z^a in K[z]/(z^p - 1) (z^a <-> [g0 a]); invert t by the nilpotent series t^-1 = sum_j (1-t)^j;
# CHECK: augmentation(t) = 1, t*t^-1 = 1, (z-1)^p = 0, and s(m) = -sum_{a != 0} d_a (g0 a)^m for 1 <= m <= 4p
# against the binomial recurrence; report the m = 0 value sum_a d_a (should be 1, i.e. e - sigma at m = 0).
from math import comb

def primes(n):
    return [q for q in range(3, n + 1) if all(q % d for d in range(2, int(q ** .5) + 1))]

class Fq:
    """F_p or F_{p^2} = F_p[r]/(r^2 - nr)"""
    def __init__(self, p, quad):
        self.p = p; self.quad = quad
        if quad:
            self.nr = next(a for a in range(2, p) if pow(a, (p - 1) // 2, p) == p - 1)
    def mul(self, x, y):
        p = self.p
        if not self.quad:
            return (x[0] * y[0] % p, 0)
        return ((x[0] * y[0] + self.nr * x[1] * y[1]) % p, (x[0] * y[1] + x[1] * y[0]) % p)
    def add(self, x, y):
        return ((x[0] + y[0]) % self.p, (x[1] + y[1]) % self.p)
    def neg(self, x):
        return ((-x[0]) % self.p, (-x[1]) % self.p)
    def pow(self, x, e):
        r = (1, 0)
        while e:
            if e & 1: r = self.mul(r, x)
            x = self.mul(x, x); e >>= 1
        return r
    def inv(self, x):
        q = self.p ** (2 if self.quad else 1)
        return self.pow(x, q - 2)
    def elements(self):
        p = self.p
        if self.quad:
            return [(a, b) for a in range(p) for b in range(p)]
        return [(a, 0) for a in range(p)]

def smod(M, p):
    T = [pow(2, j * (j - 1) // 2, p) for j in range(M + 1)]
    s = [0] * (M + 1)
    for m in range(1, M + 1):
        s[m] = (T[m] - sum(comb(m, k) * s[k] * T[m - k] for k in range(1, m))) % p
    return s

def solve(F, A, rhs):
    n = len(A); A = [row[:] + [rhs[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = next(r for r in range(col, n) if A[r][col] != (0, 0))
        A[col], A[piv] = A[piv], A[col]
        iv = F.inv(A[col][col])
        A[col] = [F.mul(iv, x) for x in A[col]]
        for r in range(n):
            if r != col and A[r][col] != (0, 0):
                f = A[r][col]
                A[r] = [F.add(A[r][j], F.neg(F.mul(f, A[col][j]))) for j in range(n + 1)]
    return [A[i][n] for i in range(n)]

def polmul(F, a, b, p):
    r = [(0, 0)] * p
    for i, x in enumerate(a):
        if x == (0, 0): continue
        for j, y in enumerate(b):
            if y == (0, 0): continue
            r[(i + j) % p] = F.add(r[(i + j) % p], F.mul(x, y))
    return r

allok = True
for p in primes(31):
    eps = 1 if pow(2, (p - 1) // 2, p) == 1 else -1
    F = Fq(p, eps == -1)
    epsel = (1, 0) if eps == 1 else (p - 1, 0)
    g0 = next(x for x in F.elements() if x != (0, 0) and F.pow(x, p - 1) == epsel)
    nodes = [F.mul(g0, (a, 0)) for a in range(1, p)]
    M = 4 * p + 2
    tau = [(pow(2, m * (m - 1) // 2, p), 0) for m in range(M + 1)]
    A = [[F.pow(x, m) for x in nodes] for m in range(p - 1)]
    cvec = solve(F, A, tau[:p - 1])
    # check representation on all m <= M
    rep_ok = all(tau[m] == _sum for m in range(M + 1)
                 for _sum in [__import__('functools').reduce(F.add, [F.mul(cvec[i], F.pow(nodes[i], m)) for i in range(p - 1)], (0, 0))])
    t = [(0, 0)] + cvec  # coefficient of z^a, a = 0..p-1
    aug = __import__('functools').reduce(F.add, t, (0, 0))
    one = [(1, 0)] + [(0, 0)] * (p - 1)
    nil = [F.add(one[i], F.neg(t[i])) for i in range(p)]  # 1 - t
    zm1 = [(p - 1, 0), (1, 0)] + [(0, 0)] * (p - 2)
    pw = one
    for _ in range(p): pw = polmul(F, pw, zm1, p)
    nilpotent_ok = all(x == (0, 0) for x in pw)
    inv = one[:]; term = one[:]
    for _ in range(1, p + 1):
        term = polmul(F, term, nil, p); inv = [F.add(inv[i], term[i]) for i in range(p)]
    unit_ok = polmul(F, t, inv, p) == one
    s = smod(M, p)
    val_ok = True
    for m in range(1, M + 1):
        v = (0, 0)
        for a in range(1, p):
            v = F.add(v, F.mul(inv[a], F.pow(nodes[a - 1], m)))
        if F.neg(v) != (s[m], 0):
            val_ok = False
    m0 = __import__('functools').reduce(F.add, inv, (0, 0))
    ok = rep_ok and aug == (1, 0) and nilpotent_ok and unit_ok and val_ok and m0 == (1, 0)
    allok &= ok
    print(f"p={p} eps={eps:+d} K={'F_p^2' if F.quad else 'F_p'} g0={g0}: tau=iota(t) {rep_ok}, aug(t)=1 {aug == (1, 0)}, "
          f"(z-1)^p=0 {nilpotent_ok}, t invertible {unit_ok}, s(m)=-sum d_a (g0 a)^m for 1<=m<={M}: {val_ok}, "
          f"sum_x d_x (m=0) = {m0}, d_0 = {inv[0]}")
print("ALL GROUP-RING CHECKS PASSED" if allok else "SOME CHECK FAILED")
