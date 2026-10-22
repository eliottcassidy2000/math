"""The Arnold cat map in modular data.  Session opus-2026-10-22-S55.

In SL2(Z), with s = [[0,-1],[1,0]] and t = [[1,1],[0,1]], the cat map A = [[2,1],[1,1]] equals the commutator
t s t^-1 s^-1 and the square of the Fibonacci matrix Q = [[1,1],[1,0]].  For modular data (S, T) of a modular tensor
category the image C = T S T^-1 S^-1 does not depend on the scalar normalisations of S and T (central charge, Gauss sum
phase): it is a canonical unitary matrix.  Its order is compared with the Pisano bounds
    ord_N(A)  (= pi(N)/2 for N >= 3; N = order of the normalised twists, projective conductor; C^(ord_N A) is scalar),
    ord_N'(A) (N' = order of T including exp(-2 pi i c/24), the linear conductor; ord C divides it).
Families: SU(2)_k, SU(3)_k, Virasoro minimal models M(p,p'), Fibonacci, pointed categories from lattice discriminant
forms (definite and indefinite Gram matrices), quantum doubles D(Z_N) and D(S_3), Deligne products.
Every datum is checked (after audit_S55): S symmetric and unitary, S^2 a permutation matrix, (S T)^3 = Theta S^2 with
Theta = exp(2 pi i c/8), and Verlinde integrality for rank <= 100.
Usage: python catmap_mtc.py"""
import sys, math, itertools, cmath
from fractions import Fraction
import numpy as np


# ---------------------------------------------------------------- SL2(Z/N) side
def pisano(N):
    if N == 1:
        return 1
    a, b, k = 0, 1, 0
    while True:
        a, b = b, (a + b) % N
        k += 1
        if a == 0 and b == 1:
            return k


def mat_order(M, N, proj=False, cap=10 ** 6):
    """order of the 2x2 integer matrix M in SL2(Z/N) (proj: modulo +-I)."""
    X = [[1, 0], [0, 1]]
    for k in range(1, cap):
        X = [[(X[0][0] * M[0][0] + X[0][1] * M[1][0]) % N, (X[0][0] * M[0][1] + X[0][1] * M[1][1]) % N],
             [(X[1][0] * M[0][0] + X[1][1] * M[1][0]) % N, (X[1][0] * M[0][1] + X[1][1] * M[1][1]) % N]]
        if X == [[1 % N, 0], [0, 1 % N]] or (proj and X == [[(-1) % N, 0], [0, (-1) % N]]):
            return k
    raise ValueError


CAT = [[2, 1], [1, 1]]


# ---------------------------------------------------------------- modular data
def denominators(hs):
    return math.lcm(*[Fraction(h).limit_denominator(10 ** 6).denominator for h in hs])


class MD:
    def __init__(self, name, S, h, c):
        self.name, self.S, self.h, self.c = name, np.array(S, dtype=complex), [Fraction(x) for x in h], Fraction(c)
        self.T = np.diag([cmath.exp(2j * math.pi * float(x)) for x in self.h])
        r = len(self.h)
        U = self.S @ self.S.conj().T
        assert np.allclose(U, np.eye(r), atol=1e-9), (name, 'S not unitary')
        assert np.allclose(self.S, self.S.T, atol=1e-9), (name, 'S not symmetric')
        self.Nproj = denominators(self.h)
        self.Nlin = denominators([x - self.c / 24 for x in self.h])

    def catmap(self):
        T, S = self.T, self.S
        return T @ S @ np.linalg.inv(T) @ np.linalg.inv(S)

    def relations(self, verlinde_rank=100):
        """modular relations: S^2 a permutation, (S T)^3 = Theta S^2 with Theta = exp(2 pi i c/8), Verlinde fusion
        coefficients non-negative integers (rank <= verlinde_rank; None above it)."""
        S, T, r = self.S, self.T, len(self.h)
        S2 = S @ S
        perm = bool(np.allclose(np.abs(S2), np.round(np.abs(S2)), atol=1e-8) and np.allclose(np.abs(S2).sum(axis=1), 1)
                    and np.allclose(np.abs(S2).sum(axis=0), 1) and np.allclose(S2.imag, 0, atol=1e-8))
        ST3 = np.linalg.matrix_power(S @ T, 3)
        theta = cmath.exp(2j * math.pi * float(self.c) / 8)
        theta_ok = bool(np.allclose(ST3, theta * S2, atol=1e-7))
        ver = None
        if r <= verlinde_rank:
            Nabc = np.einsum('ax,bx,cx->abc', S, S, S.conj() / S[0][None, :], optimize=True)
            ver = bool(np.allclose(Nabc, np.round(Nabc.real), atol=1e-6) and (np.round(Nabc.real) >= 0).all())
        return perm, theta_ok, ver


def root_order(z, cap):
    a = (cmath.phase(z) / (2 * math.pi)) % 1.0
    f = Fraction(a).limit_denominator(cap)
    assert abs(float(f) - a) < 1e-7 and abs(abs(z) - 1) < 1e-7, (z, f)
    return f.denominator


def order_of(C, cap):
    ev = np.linalg.eigvals(C)
    m = math.lcm(*[root_order(z, cap) for z in ev])
    assert np.allclose(np.linalg.matrix_power(C, m), np.eye(len(C)), atol=1e-6)
    return m


def su2(k):
    n = k + 2
    S = [[math.sqrt(2 / n) * math.sin(math.pi * (a + 1) * (b + 1) / n) for b in range(k + 1)] for a in range(k + 1)]
    h = [Fraction(a * (a + 2), 4 * n) for a in range(k + 1)]
    return MD(f'SU(2)_{k}', S, h, Fraction(3 * k, k + 2))


def su3(k):
    n = k + 3
    wts = [(a, b) for a in range(k + 1) for b in range(k + 1 - a)]
    F = np.array([[2, 1], [1, 2]]) / 3.0                    # quadratic form matrix in Dynkin labels

    def ip(x, y):
        return float(np.array(x) @ F @ np.array(y))

    def weyl(v):
        out = []
        cur = [(tuple(v), 1)]
        seen = {tuple(v): 1}
        while cur:
            nxt = []
            for w, sg in cur:
                for refl in (lambda a, b: (-a, a + b), lambda a, b: (a + b, -b)):
                    u = refl(*w)
                    if u not in seen:
                        seen[u] = -sg; nxt.append((u, -sg))
            cur = nxt
        return list(seen.items())
    rho = (1, 1)
    S = np.zeros((len(wts), len(wts)), dtype=complex)
    for i, l in enumerate(wts):
        lr = (l[0] + 1, l[1] + 1)
        orb = weyl(lr)
        for j, m in enumerate(wts):
            mr = (m[0] + 1, m[1] + 1)
            S[i, j] = sum(sg * cmath.exp(-2j * math.pi * ip(w, mr) / n) for w, sg in orb)
    S = S / np.linalg.norm(S[0])
    ph = S[0, 0] / abs(S[0, 0])
    S = S / ph
    h = [Fraction(int(round(3 * ip(l, (l[0] + 2, l[1] + 2)))), 6 * n) for l in wts]   # <l, l+2rho>/(2n)
    return MD(f'SU(3)_{k}', S, h, Fraction(8 * k, k + 3))


def minimal(p, pp):
    """Virasoro minimal model M(p, p'), gcd = 1, 2 <= p < p'."""
    reps = []
    for r in range(1, pp):
        for s in range(1, p):
            if (pp - r, p - s) not in reps:
                reps.append((r, s))
    h = [Fraction((p * r - pp * s) ** 2 - (p - pp) ** 2, 4 * p * pp) for r, s in reps]
    S = [[2 * math.sqrt(2 / (p * pp)) * (-1) ** (1 + s * rr + r * ss) * math.sin(math.pi * p * r * rr / pp)
          * math.sin(math.pi * pp * s * ss / p) for (rr, ss) in reps] for (r, s) in reps]
    c = Fraction(1) - Fraction(6 * (p - pp) ** 2, p * pp)
    return MD(f'M({p},{pp})', S, h, c)


def fibonacci():
    phi = (1 + 5 ** 0.5) / 2
    D = math.sqrt(2 + phi)
    return MD('Fib = (G2)_1', [[1 / D, phi / D], [phi / D, -1 / D]], [0, Fraction(2, 5)], Fraction(14, 5))


def smith_discriminant(Gram):
    """discriminant group L*/L of an even lattice with integer Gram matrix (any signature), as a list of coset
    representatives in L* coordinates (x = Gram^-1 * integer vector), with q(x) = x^T Gram x / 2 mod 1."""
    G = np.array(Gram, dtype=object)
    n = len(G)
    import sympy as sp
    M = sp.Matrix(G.tolist())
    Minv = M.inv()
    det = abs(int(M.det()))
    # enumerate L*/L: vectors Minv*v for v in Z^n / M Z^n -- brute force over a box (small determinants only)
    reps, seen = [], set()
    B = det
    for v in itertools.product(range(B), repeat=n):
        x = Minv * sp.Matrix(v)
        key = tuple(sp.nsimplify(c - sp.floor(c)) for c in (x))
        # canonical class: reduce x modulo Z^n (L* / L in dual coordinates: x ~ x + integer vectors? no: L = Z^n)
        if key in seen:
            continue
        seen.add(key); reps.append(sp.Matrix(key))
        if len(reps) == det:
            break
    q = [Fraction(str(sp.nsimplify((x.T * M * x)[0] / 2))) % 1 for x in reps]
    bil = [[Fraction(str(sp.nsimplify((x.T * M * y)[0]))) % 1 for y in reps] for x in reps]
    return reps, q, bil


def pointed(name, Gram):
    import sympy as sp
    reps, q, bil = smith_discriminant(Gram)
    n = len(reps)
    S = [[cmath.exp(-2j * math.pi * float(bil[i][j])) / math.sqrt(n) for j in range(n)] for i in range(n)]
    sig = sum(1 if e > 0 else -1 for e in np.linalg.eigvalsh(np.array(Gram, dtype=float)))
    gm = sum(cmath.exp(2j * math.pi * float(x)) for x in q) / math.sqrt(n)        # Gauss-Milgram: exp(2 pi i sig/8)
    assert abs(gm - cmath.exp(2j * math.pi * sig / 8)) < 1e-9, (name, gm, sig)
    md = MD(name, S, q, Fraction(sig))
    md.signature = sig
    return md


def double_ZN(N):
    """D(Z_N): anyons (a, b) in Z_N^2, h = ab/N, S = exp(-2 pi i (a b' + a' b)/N)/N."""
    an = [(a, b) for a in range(N) for b in range(N)]
    S = [[cmath.exp(-2j * math.pi * (a * bb + aa * b) / N) / N for (aa, bb) in an] for (a, b) in an]
    return MD(f'D(Z_{N})', S, [Fraction(a * b, N) for a, b in an], 0)


def double_S3():
    """D(S_3) modular data (Coste-Gannon-Ruelle conventions): 8 anyons."""
    w = cmath.exp(2j * math.pi / 3)
    # labels: A (e,1), B (e,sign), C (e,2-dim), D (3-cycle class, chi=1), E (3-cycle, w), F (3-cycle, w^2),
    #         G (transposition class, +), H (transposition class, -)
    S = (1 / 6) * np.array([
        [1, 1, 2, 2, 2, 2, 3, 3],
        [1, 1, 2, 2, 2, 2, -3, -3],
        [2, 2, 4, -2, -2, -2, 0, 0],
        [2, 2, -2, 4, -2, -2, 0, 0],
        [2, 2, -2, -2, -2, 4, 0, 0],
        [2, 2, -2, -2, 4, -2, 0, 0],
        [3, -3, 0, 0, 0, 0, 3, -3],
        [3, -3, 0, 0, 0, 0, -3, 3]], dtype=complex)
    h = [0, 0, 0, 0, Fraction(1, 3), Fraction(2, 3), 0, Fraction(1, 2)]
    return MD('D(S_3)', S, h, 0)


def deligne(a, b):
    S = np.kron(a.S, b.S)
    h = [x + y for x in a.h for y in b.h]
    md = MD(f'{a.name} x {b.name}', S, h, a.c + b.c)
    return md


def proj_order(C, cap):
    """smallest m with C^m a scalar matrix."""
    X = np.eye(len(C), dtype=complex)
    for m in range(1, cap + 1):
        X = X @ C
        if np.allclose(X, X[0, 0] * np.eye(len(C)), atol=1e-7):
            return m
    raise ValueError


def report(md, out):
    C = md.catmap()
    o = order_of(C, cap=10 ** 5)
    po = proj_order(C, o)
    Np, Nl = md.Nproj, md.Nlin
    mN = mat_order(CAT, Np) if Np > 1 else 1               # = pi(N)/2 for N >= 3
    mL = mat_order(CAT, Nl) if Nl > 1 else 1
    sat = 'SATURATES' if po == mN else f'deficient by {mN // po}' if mN % po == 0 else 'NOT A DIVISOR'
    perm, theta_ok, ver = md.relations()
    print(f'{md.name:30s} rank {len(md.h):4d}  N {Np:4d}  pi(N) {pisano(Np):5d}  ord_SL2(Z/N) A {mN:4d}  proj ord C {po:4d}  '
          f'ord C {o:4d}  Nlin {Nl:5d}  ord_SL2(Z/Nlin) A {mL:4d}  [{sat}; ord C | ord_Nlin A: {mL % o == 0}]  '
          f'[S^2 perm {perm}; (ST)^3 = e(c/8) S^2 {theta_ok}; Verlinde {ver}]', flush=True)
    out.append(dict(name=md.name, rank=len(md.h), N=Np, pi=pisano(Np), mN=mN, po=po, o=o, Nl=Nl, mL=mL,
                    S2_perm=perm, ST3_theta=theta_ok, verlinde=ver))


E8 = [[2, -1, 0, 0, 0, 0, 0, 0], [-1, 2, -1, 0, 0, 0, 0, 0], [0, -1, 2, -1, 0, 0, 0, -1], [0, 0, -1, 2, -1, 0, 0, 0],
      [0, 0, 0, -1, 2, -1, 0, 0], [0, 0, 0, 0, -1, 2, -1, 0], [0, 0, 0, 0, 0, -1, 2, 0], [0, 0, -1, 0, 0, 0, 0, 2]]
LATTICES = [('A1 (definite)', [[2]]), ('A1 + U (Lorentzian)', [[2, 0, 0], [0, 0, 1], [0, 1, 0]]),
            ('A1(-1) (negative definite)', [[-2]]), ('A2 (definite)', [[2, -1], [-1, 2]]),
            ('A2 + U (Lorentzian)', [[2, -1, 0, 0], [-1, 2, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]),
            ('A4 (definite)', [[2, -1, 0, 0], [-1, 2, -1, 0], [0, -1, 2, -1], [0, 0, -1, 2]]),
            ('D4 (definite)', [[2, -1, 0, 0], [-1, 2, -1, -1], [0, -1, 2, 0], [0, -1, 0, 2]]),
            ('golden [[2,1],[1,-2]] (indef.)', [[2, 1], [1, -2]]),
            ('E8 (definite, unimodular)', E8), ('U = II_1,1 (Lorentzian)', [[0, 1], [1, 0]])]


if __name__ == '__main__':
    out = []
    print('cat map A = t s t^-1 s^-1 =', np.array([[1, 1], [0, 1]]) @ np.array([[0, -1], [1, 0]]) @ np.linalg.inv(np.array([[1, 1], [0, 1]])) @ np.linalg.inv(np.array([[0, -1], [1, 0]])))
    print('Pisano periods pi(N), N = 1..30:', [pisano(N) for N in range(1, 31)])
    print('ord of A in SL2(Z/N), N = 2..30:', [mat_order(CAT, N) for N in range(2, 31)])
    print('pi(N) = 2 ord_N(A) for N >= 3:', all(pisano(N) == 2 * mat_order(CAT, N) for N in range(3, 200)))
    print('no N has pi(N) = 15 (N <= 10^4):', all(pisano(N) != 15 for N in range(1, 10001)),
          '; least N with pi(N) = 120:', min(N for N in range(1, 1000) if pisano(N) == 120),
          '; pi(4), pi(6), pi(30) =', pisano(4), pisano(6), pisano(30))
    print()
    fams = [fibonacci()] + [minimal(p, pp) for p in range(2, 7) for pp in range(p + 1, 13) if math.gcd(p, pp) == 1 and p * pp <= 66]
    for md in fams:
        report(md, out)
    for k in range(1, 61):
        report(su2(k), out)
    for k in range(1, 11):
        report(su3(k), out)
    for name, G in LATTICES:
        report(pointed(name, G), out)
    for N in range(2, 31):
        report(double_ZN(N), out)
    report(double_S3(), out)
    report(deligne(fibonacci(), deligne(double_ZN(2), double_ZN(3))), out)
    import json
    json.dump(out, open('catmap_census.json', 'w'), indent=1)
    sat = [r for r in out if r['po'] == r['mN']]
    print()
    print(f'{len(out)} modular data; projective order of C divides ord_SL2(Z/N)(A) = pi(N)/2 in all: '
          f'{all(r["mN"] % r["po"] == 0 for r in out)}; ord C divides ord_SL2(Z/Nlin)(A) in all: {all(r["mL"] % r["o"] == 0 for r in out)}; '
          f'saturating: {len(sat)} ({sum(r["N"] >= 3 for r in sat)} with N >= 3)')
    print(f'modular relations hold in all: S^2 permutation {all(r["S2_perm"] for r in out)}; (ST)^3 = e(c/8) S^2 '
          f'{all(r["ST3_theta"] for r in out)}; Verlinde integrality {all(r["verlinde"] for r in out if r["verlinde"] is not None)} '
          f'({sum(r["verlinde"] is not None for r in out)} data of rank <= 100)')
