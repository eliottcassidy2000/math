"""Weyl-vector cosets of the irreducible root lattices and the Coxeter Frame shape (session opus-2026-10-22-S53).

For an irreducible simply-laced root lattice Q of rank l, Coxeter number h, Weyl vector rho, Coxeter element c with
char. polynomial det(x - c) = prod_a (x^a - 1)^(m_a) (Frame shape pi = prod a^(m_a)), check

  (I)  sum_{x in Q} q^{(h/2)|x + rho/h|^2} = eta(h tau)^l eta_pi(tau),   eta_pi = prod_a eta(a tau)^(m_a);
       equivalently theta_{Q + rho/h}(tau) = eta(tau)^l eta_pi(tau/h)  (classical: principal specialisation of the
       basic representation, Lepowsky / Kac-Kazhdan-Lepowsky-Wilson; for A_(t-1) these are t-cores, GKS 1990);
  (II) every coset Q + w (w a minuscule weight) has the same Weyl-vector-shifted theta series as Q (the special
       automorphisms of the extended Dynkin diagram fix the alcove barycentre rho/h);
  (III) pi is anti-self-dual: a -> h/a sends pi to pi^(-1), and prod_a a^(m_a) = det(1 - c) = |Q^*/Q|.

Pure Python (fractions); runs in seconds. Usage: python principal_cores.py [TERMS=40]"""
import sys
from fractions import Fraction as Fr
import itertools

TERMS = int(sys.argv[1]) if len(sys.argv) > 1 else 40
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


def cartan(t, r):
    G = [[0] * r for _ in range(r)]
    for i in range(r):
        G[i][i] = 2
    if t == 'A':
        edges = [(i, i + 1) for i in range(r - 1)]
    elif t == 'D':
        edges = [(i, i + 1) for i in range(r - 2)] + [(r - 3, r - 1)]
    else:   # E_r: chain 0-1-...-(r-2), node r-1 attached to node 2
        edges = [(i, i + 1) for i in range(r - 2)] + [(2, r - 1)]
    for i, j in edges:
        G[i][j] = G[j][i] = -1
    return G


def frame(t, r):
    if t == 'A':
        return {r + 1: 1, 1: -1}
    if t == 'D':
        out = {}
        for a, m in ((2 * r - 2, 1), (2, 1), (r - 1, -1), (1, -1)):
            out[a] = out.get(a, 0) + m
        return {a: m for a, m in out.items() if m}
    return {6: {12: 1, 2: 1, 3: 1, 6: -1, 4: -1, 1: -1},
            7: {18: 1, 3: 1, 2: 1, 9: -1, 6: -1, 1: -1},
            8: {30: 1, 5: 1, 3: 1, 2: 1, 15: -1, 10: -1, 6: -1, 1: -1}}[r]


def coxeter_number(t, r):
    return {'A': r + 1, 'D': 2 * r - 2, 'E': {6: 12, 7: 18, 8: 30}.get(r)}[t]


def inv(M):
    n = len(M)
    A = [[Fr(M[i][j]) for j in range(n)] + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]


def det(M):
    n = len(M)
    A = [[Fr(x) for x in row] for row in M]
    d = Fr(1)
    for c in range(n):
        p = next((i for i in range(c, n) if A[i][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            A[c], A[p] = A[p], A[c]; d = -d
        d *= A[c][c]
        for i in range(c + 1, n):
            f = A[i][c] / A[c][c]
            A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return d


def charpoly_coxeter_matches(t, r):
    """Coxeter element c = s_1 ... s_r on the root lattice; compare det(x - c) with the Frame shape."""
    import sympy as sp
    G = cartan(t, r)
    x = sp.symbols('x')
    M = sp.eye(r)
    for i in range(r):          # reflection s_i in root coordinates: alpha_j -> alpha_j - G[i][j] alpha_i
        S = sp.eye(r)
        for j in range(r):
            S[i, j] -= G[i][j]
        M = M * S
    cp = sp.expand((x * sp.eye(r) - M).det())
    fp = sp.Integer(1)
    for a, m in frame(t, r).items():
        fp *= (x ** a - 1) ** m
    return sp.simplify(cp - sp.cancel(fp)) == 0


def coset_counts(G, h, shift, B):
    """counts[n] = #{a in Z^r : (h/2)|a + shift + rho/h|^2 - const = n}, shift in root coords (Fractions).
    Returns dict n -> count for exponents (exact rationals) below B, relative to the minimum over Q + rho/h."""
    r = len(G)
    Gi = inv(G)
    rho = [sum(Gi[i][j] for j in range(r)) for i in range(r)]       # (rho, alpha_i) = 1
    cen = [-(s + p / h) for s, p in zip(shift, rho)]               # minimise |a - cen|^2
    rho2 = sum(rho[i] * G[i][j] * rho[j] for i in range(r) for j in range(r))
    base = Fr(h, 2) * rho2 / (h * h)                               # value at a = 0, shift = 0
    # Cholesky-free bound search: brute force over a box sized by the smallest eigen-bound is wasteful; use
    # recursive Fincke-Pohst in exact rationals with LDL^T
    L = [[Fr(0)] * r for _ in range(r)]
    D = [Fr(0)] * r
    for j in range(r):
        D[j] = G[j][j] - sum(L[j][k] ** 2 * D[k] for k in range(j))
        for i in range(j + 1, r):
            L[i][j] = (G[i][j] - sum(L[i][k] * L[j][k] * D[k] for k in range(j))) / D[j]
    # q(a) = sum_j D[j] (y_j)^2 with y_j = (a_j - cen_j) + sum_{i>j} L[i][j] (a_i - cen_i)
    bound = Fr(2, h) * (B + base) + 1          # |a - cen|^2 <= bound (generous)
    out = {}
    a = [0] * r

    def rec(j, partial):
        if j < 0:
            val = Fr(h, 2) * sum((a[i] - cen[i]) * G[i][k] * (a[k] - cen[k]) for i in range(r) for k in range(r))
            e = val - base
            if e < B:
                out[e] = out.get(e, 0) + 1
            return
        off = sum(L[i][j] * (a[i] - cen[i]) for i in range(j + 1, r))
        c0 = cen[j] - off
        rem = bound - partial
        if rem < 0:
            return
        import math
        w = math.sqrt(float(rem / D[j])) + 1e-9
        lo, hi = math.ceil(float(c0) - w), math.floor(float(c0) + w)
        for v in range(lo, hi + 1):
            y = v - c0
            p2 = partial + D[j] * y * y
            if p2 <= bound:
                a[j] = v
                rec(j - 1, p2)
        a[j] = 0

    rec(r - 1, Fr(0))
    return out


def eta_series(fr, l, h, L):
    """coefficients of prod_n (1-q^{hn})^l prod_a prod_n (1-q^{an})^{m_a}, q^0..q^{L-1}."""
    s = [0] * L; s[0] = 1
    terms = [(h, l)] + list(fr.items())
    for a, m in terms:
        for k in range(1, L // a + 1):
            for _ in range(abs(m)):
                if m > 0:
                    s = [s[t] - (s[t - a * k] if t >= a * k else 0) for t in range(L)]
                else:
                    for t in range(a * k, L):
                        s[t] += s[t - a * k]
    return s


TYPES = [('A', r) for r in range(1, 9)] + [('A', 12)] + [('D', r) for r in range(4, 11)] + [('E', 6), ('E', 7), ('E', 8)]
for t, r in TYPES:
    G = cartan(t, r)
    h = coxeter_number(t, r)
    fr = frame(t, r)
    # (III)
    dual = {}
    for a, m in fr.items():
        assert h % a == 0
        dual[h // a] = dual.get(h // a, 0) + m
    antiself = all(dual.get(a, 0) == -m for a, m in fr.items()) and len(dual) == len(fr)
    detQ = det(G)
    prod = Fr(1)
    for a, m in fr.items():
        prod *= Fr(a) ** m
    cp_ok = charpoly_coxeter_matches(t, r)
    check(antiself and prod == detQ and sum(a * m for a, m in fr.items()) == r and cp_ok,
          f'{t}{r}: Frame shape {fr} = char. poly of the Coxeter element; anti-self-dual; prod a^m = det = {detQ}')
    # (I)
    B = TERMS
    cnt = coset_counts(G, h, [Fr(0)] * r, B)
    integral = all(e.denominator == 1 for e in cnt)
    pred = eta_series(fr, r, h, B)
    got = [cnt.get(Fr(n), 0) for n in range(B)]
    check(integral and got == pred, f'{t}{r} (h={h}): sum_Q q^((h/2)|x+rho/h|^2) = eta(h tau)^{r} eta_pi(tau) through q^{B - 1}: {got[:12]}...')
    # (II): cosets by minuscule weights (columns of the inverse Cartan matrix for the special nodes)
    Gi = inv(G)
    special = list(range(r)) if t == 'A' else [0, r - 2, r - 1] if t == 'D' else {6: [0, 4], 7: [5], 8: []}[r]
    for j in special:
        w = [Gi[i][j] for i in range(r)]
        cw = coset_counts(G, h, w, B)
        same = all(e.denominator == 1 for e in cw) and [cw.get(Fr(n), 0) for n in range(B)] == got
        check(same, f'{t}{r}: coset Q + omega_{j + 1} has the same rho/h-shifted theta series')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
