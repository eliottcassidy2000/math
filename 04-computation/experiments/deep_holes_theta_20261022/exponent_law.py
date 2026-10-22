"""The exponent law behind the shallow multiples of the Weyl vector (session opus-2026-10-22-S53; THM-4759 (5)).

R irreducible simply-laced, rank l, Coxeter number h, exponents m_1..m_l. For a divisor e of h put
J_e(m) = r (e - r) with r = m mod e.

(A) Exponent inequality:  sum_i J_e(m_i) >= l e (e - 1) / 6, with equality iff e = 1 or e = h/p for a bad prime p of R
    with p || h.  (Proved by hand for A_n and D_n in THM-4759; checked here exactly for A_1..A_200, D_4..D_200,
    E_6, E_7, E_8.)
(B) Twisted-module bound: c^(h/e) has eigenvalues exp(2 pi i m_i / e), so (Kostant) every coset Q + w + rho/e, w in P/Q,
    is the weight lattice of a c^(h/e)-twisted V_Q-module, of conformal weight >= the vacuum anomaly
    sum_i J_e(m_i) / (4 e^2). Hence min_{x in P} |x + rho/e|^2 >= sum_i J_e(m_i) / (2 e^2) >= (l/12)(1 - 1/e).
    Checked here by exact enumeration of every coset for A_1..A_8, D_4..D_10, E_6, E_7, E_8, with the equality cases.
Pure Python (fractions). Usage: python exponent_law.py"""
import math
from fractions import Fraction as Fr

FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


def exponents(t, r):
    if t == 'A':
        return list(range(1, r + 1))
    if t == 'D':
        return list(range(1, 2 * r - 2, 2)) + [r - 1]
    return {6: [1, 4, 5, 7, 8, 11], 7: [1, 5, 7, 9, 11, 13, 17], 8: [1, 7, 11, 13, 17, 19, 23, 29]}[r]


def cox(t, r):
    return {'A': r + 1, 'D': 2 * r - 2, 'E': {6: 12, 7: 18, 8: 30}.get(r)}[t]


def badp(t, r):
    return [] if t == 'A' else [2] if t == 'D' else ([2, 3] if r in (6, 7) else [2, 3, 5])


# (A) exponent inequality
cnt = 0
okA = True
eq_cases = {}
for t, rng in (('A', range(1, 201)), ('D', range(4, 201)), ('E', (6, 7, 8))):
    for r in rng:
        h, E, l = cox(t, r), exponents(t, r), r
        assert len(E) == l and sum(E) == l * h // 2
        for e in range(1, h + 1):
            if h % e:
                continue
            S = sum((m % e) * (e - m % e) for m in E)
            lhs6, rhs6 = 6 * S, l * e * (e - 1)
            pred_eq = e == 1 or any(e == h // p and math.gcd(p, h // p) == 1 for p in badp(t, r))
            okA &= lhs6 >= rhs6 and ((lhs6 == rhs6) == pred_eq)
            if lhs6 == rhs6 and e > 1:
                eq_cases.setdefault(t, set()).add((r, e))
            cnt += 1
check(okA, f'(A) sum_i J_e(m_i) >= l e (e-1)/6, equality iff e = 1 or e = h/p (p bad, p || h): {cnt} pairs (R, e); '
           f'equality cases e > 1: E: {sorted(eq_cases.get("E", []))}, D: n even (e = n - 1), A: none')
# the strange formula as the case e = h: sum m (h - m) = l h (h + 1) / 6
okS = all(sum(m * (cox(t, r) - m) for m in exponents(t, r)) * 6 == r * cox(t, r) * (cox(t, r) + 1)
          for t, rng in (('A', range(1, 60)), ('D', range(4, 60)), ('E', (6, 7, 8))) for r in rng)
check(okS, 'e = h: sum_i m_i (h - m_i) = l h (h + 1)/6 (equivalent to the Freudenthal-de Vries strange formula)')


# (B) exact minima over the weight lattice
def cartan(t, r):
    G = [[0] * r for _ in range(r)]
    for i in range(r):
        G[i][i] = 2
    if t == 'A':
        edges = [(i, i + 1) for i in range(r - 1)]
    elif t == 'D':
        edges = [(i, i + 1) for i in range(r - 2)] + [(r - 3, r - 1)]
    else:
        edges = [(i, i + 1) for i in range(r - 2)] + [(2, r - 1)]
    for i, j in edges:
        G[i][j] = G[j][i] = -1
    return G


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


def min_norm(G, cen):
    """min over a in Z^r of (a - cen)^T G (a - cen), exact; recursive Fincke-Pohst with LDL^T."""
    r = len(G)
    L = [[Fr(0)] * r for _ in range(r)]
    D = [Fr(0)] * r
    for j in range(r):
        D[j] = Fr(G[j][j]) - sum((L[j][k] ** 2 * D[k] for k in range(j)), Fr(0))
        for i in range(j + 1, r):
            L[i][j] = (Fr(G[i][j]) - sum((L[i][k] * L[j][k] * D[k] for k in range(j)), Fr(0))) / D[j]
    best = [None]
    # initial bound: rounding
    a0 = [round(c) for c in cen]
    best[0] = sum((a0[i] - cen[i]) * G[i][k] * (a0[k] - cen[k]) for i in range(r) for k in range(r))
    a = [0] * r

    def rec(j, partial):
        if j < 0:
            val = partial
            if val < best[0]:
                best[0] = val
            return
        off = sum(L[i][j] * (a[i] - cen[i]) for i in range(j + 1, r))
        c0 = cen[j] - off
        rem = best[0] - partial
        if rem < 0:
            return
        w = math.sqrt(float(rem / D[j])) + 1e-9
        for v in range(math.ceil(float(c0) - w), math.floor(float(c0) + w) + 1):
            y = v - c0
            p2 = partial + D[j] * y * y
            if p2 <= best[0]:
                a[j] = v
                rec(j - 1, p2)
        a[j] = 0

    rec(r - 1, Fr(0))
    return best[0]


okB = True
lines = []
for t, r in [('A', r) for r in range(1, 9)] + [('D', r) for r in range(4, 11)] + [('E', 6), ('E', 7), ('E', 8)]:
    G = cartan(t, r)
    Gi = inv(G)
    h, E = cox(t, r), exponents(t, r)
    rho = [sum(Gi[i][j] for j in range(r)) for i in range(r)]            # root coordinates of rho
    special = list(range(r)) if t == 'A' else [0, r - 2, r - 1] if t == 'D' else {6: [0, 4], 7: [5], 8: []}[r]
    weights = [[Fr(0)] * r] + [[Gi[i][j] for i in range(r)] for j in special]   # minuscule coset representatives
    eqs = []
    for e in [e for e in range(1, h + 1) if h % e == 0]:
        S = sum((m % e) * (e - m % e) for m in E)
        bound_tw = Fr(S, 2 * e * e)
        bound_l = Fr(r, 12) * (1 - Fr(1, e))
        mn = min(min_norm(G, [-(w[i] + rho[i] / e) for i in range(r)]) for w in weights)
        okB &= mn >= bound_tw >= bound_l
        if bound_tw == bound_l:
            eqs.append((e, mn == bound_tw))
    lines.append(f'{t}{r}: e with sum J = l e(e-1)/6: {[e for e, _ in eqs]} (min over P attains it: {all(x for _, x in eqs)})')
    okB &= all(x for _, x in eqs)
check(okB, '(B) min_{x in P} |x + rho/e|^2 >= sum J/(2e^2) >= (l/12)(1 - 1/e) for every e | h, with the minimum equal to '
           '(l/12)(1 - 1/e) in every equality case: ' + '; '.join(lines))
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
