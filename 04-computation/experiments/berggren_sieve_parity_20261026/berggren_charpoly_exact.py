# Exact characteristic polynomials of 3T = P_U + P_A + P_D and 3T_chi = P_U - P_A + P_D on P1(F_q), factored over Q.
import sys
from sympy import Matrix, symbols, factor_list, primerange, Poly
INF = 'inf'
def act(q, t, g):
    if g == 'D': return INF if t == INF else (t + 2) % q
    if t == INF: return 2 % q
    if t == 0: return INF
    inv = pow(t, -1, q)
    return (2 - inv) % q if g == 'U' else (2 + inv) % q
x = symbols('x')
for q in primerange(3, int(sys.argv[1]) if len(sys.argv) > 1 else 30):
    pts = list(range(q)) + [INF]; idx = {p: i for i, p in enumerate(pts)}; N = len(pts)
    for name, sg in (('3T', (1, 1, 1)), ('3T_chi', (1, -1, 1))):
        M = [[0] * N for _ in range(N)]
        for p in pts:
            for g, s in zip('UAD', sg): M[idx[p]][idx[act(q, p, g)]] += s
        cp = Matrix(M).charpoly(x).as_expr()
        fl = factor_list(cp)
        desc = ' * '.join(f"({f})^{e}" if e > 1 else f"({f})" for f, e in fl[1])
        # minimal polynomial degree of nilpotent part: rank of M^k
        Mm = Matrix(M); ranks = []; P = Mm
        for k in range(1, 6):
            ranks.append(P.rank()); P = P * Mm
        print(f"q = {q}, {name}: charpoly = {desc};  ranks of M^1..M^5: {ranks}")
