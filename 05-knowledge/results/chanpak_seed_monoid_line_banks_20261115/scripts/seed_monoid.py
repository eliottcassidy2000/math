# The seed monoid.  A seed is (block C, head pattern X).  Over R = Z[e]/(e^2 - 1) (e marks an odd number of path
# vertices) put   M(C, X) = [[h, b, c], [0, G, a], [0, 0, h]]   with
#   h = H(C),  G = G_even + G_odd e,  a = sum_{s notin X} lam[s],  b = sum_{s in X} kap[s],
#   c = sum_{u in X, v notin X} W(u, v)  (all path parities kept; the even part of c is cut_C(X)).
# Claim: M(C1 -> C2, X1 u X2) = M(C1, X1) M(C2, X2), and H(C + x) = h + 2 * (even part of c).
# A tape vertex is the one-vertex seed: non-contact M0 = [[1,0,0],[0,1+e,e],[0,0,1]], contact M1 = [[1,e,0],[0,1+e,0],[0,0,1]].
import random
from master_formula_check import block_stats, brute
class R:
    __slots__ = ('ev', 'od')
    def __init__(s, ev=0, od=0): s.ev, s.od = ev, od
    def __add__(s, t): return R(s.ev + t.ev, s.od + t.od)
    def __mul__(s, t): return R(s.ev * t.ev + s.od * t.od, s.ev * t.od + s.od * t.ev)
    def __eq__(s, t): return s.ev == t.ev and s.od == t.od
    def __repr__(s): return "%d+%de" % (s.ev, s.od)
def mmul(A, B): return [[sum((A[i][k] * B[k][j] for k in range(3)), R()) for j in range(3)] for i in range(3)]
def seed_matrix(C, m, X):
    h, W, G, lam, kap = block_stats(C, m); full = (1 << m) - 1
    a = R(sum(lam[s][0] for s in range(m) if not X >> s & 1), sum(lam[s][1] for s in range(m) if not X >> s & 1))
    b = R(sum(kap[s][0] for s in range(m) if X >> s & 1), sum(kap[s][1] for s in range(m) if X >> s & 1))
    c = R(sum(W[u][v][0] for u in range(m) for v in range(m) if u != v and (X >> u & 1) and not (X >> v & 1)),
          sum(W[u][v][1] for u in range(m) for v in range(m) if u != v and (X >> u & 1) and not (X >> v & 1)))
    return [[R(h[full]), b, c], [R(), R(G[0], G[1]), a], [R(), R(), R(h[full])]]
M0 = [[R(1), R(), R()], [R(), R(1, 1), R(0, 1)], [R(), R(), R(1)]]
M1 = [[R(1), R(0, 1), R()], [R(), R(1, 1), R()], [R(), R(), R(1)]]
if __name__ == '__main__':
    rng = random.Random(9); bad = 0; tests = 0
    for trial in range(40):
        m1 = rng.randint(1, 4); m2 = rng.randint(1, 4); m = m1 + m2
        C = [[0] * m for _ in range(m)]
        for u in range(m):
            for v in range(u + 1, m):
                if u < m1 <= v: C[u][v] = 1
                elif rng.random() < 0.5: C[u][v] = 1
                else: C[v][u] = 1
        C1 = [row[:m1] for row in C[:m1]]; C2 = [row[m1:] for row in C[m1:]]
        X1 = rng.randrange(1 << m1); X2 = rng.randrange(1 << m2)
        P = mmul(seed_matrix(C1, m1, X1), seed_matrix(C2, m2, X2)); Mfull = seed_matrix(C, m, X1 | (X2 << m1))
        tests += 1; bad += (P != Mfull)
    print("M(C1 -> C2) = M(C1) M(C2): %d random pairs, failures %d" % (tests, bad))
    # tape (letters top to bottom) -> block -> tape, against brute force
    bad2 = 0
    for trial in range(30):
        m = rng.randint(1, 4); C = [[0] * m for _ in range(m)]
        for u in range(m):
            for v in range(u + 1, m):
                if rng.random() < 0.5: C[u][v] = 1
                else: C[v][u] = 1
        i = rng.randint(0, 3); j = rng.randint(0, 3); X = rng.randrange(1 << m)
        U = {hh for hh in range(1, i + 1) if rng.random() < 0.5}; L = {d for d in range(1, j + 1) if rng.random() < 0.5}
        P = [[R(1), R(), R()], [R(), R(1), R()], [R(), R(), R(1)]]
        for hh in range(i, 0, -1): P = mmul(P, M1 if hh in U else M0)
        P = mmul(P, seed_matrix(C, m, X))
        for d in range(1, j + 1): P = mmul(P, M1 if d in L else M0)
        H = P[0][0].ev + 2 * P[0][2].ev
        bad2 += (H != brute(C, m, i, j, U, X, L))
    print("H = h + 2 * even part of the corner of the matrix product: 30 random chains against brute force, failures %d" % bad2)
    # the two characters e = +1, -1 on tape letters: (b, c) -> (2b, c + b) / (2b + 1, c)   and   c^- = number of '10' factors
    def chars(word):
        P = [[R(1), R(), R()], [R(), R(1), R()], [R(), R(), R(1)]]
        for ch in word: P = mmul(P, M1 if ch == '1' else M0)
        b, c = P[0][1], P[0][2]; return (b.ev + b.od, c.ev + c.od, b.ev - b.od, c.ev - c.od, 1 + 2 * c.ev)
    for wd in ("10", "100", "1010", "1001100", "11101001000"):
        bp, cp, bm, cm, H = chars(wd); print("word %s: (b+, c+) = (%d, %d), (b-, c-) = (%d, %d), H = 1 + c+ + c- = %d" % (wd, bp, cp, bm, cm, H))
