# (1) composition law for ordinal sums of blocks; (2) how many odd numbers the explicit two-tape bank reaches.
import random
from master_formula_check import block_stats
rng = random.Random(3); bad = 0
for trial in range(40):
    m1 = rng.randint(1, 4); m2 = rng.randint(1, 4); m = m1 + m2
    C = [[0] * m for _ in range(m)]
    for a in range(m):
        for b in range(a + 1, m):
            if a < m1 <= b: C[a][b] = 1                      # C1 -> C2
            elif rng.random() < 0.5: C[a][b] = 1
            else: C[b][a] = 1
    C1 = [row[:m1] for row in C[:m1]]; C2 = [row[m1:] for row in C[m1:]]
    h, W, G, lam, kap = block_stats(C, m); h1, W1, G1, lam1, kap1 = block_stats(C1, m1); h2, W2, G2, lam2, kap2 = block_stats(C2, m2)
    X1 = rng.randrange(1 << m1); X2 = rng.randrange(1 << m2); X = X1 | (X2 << m1)
    cut = lambda W, X, n: sum(W[a][b][0] for a in range(n) for b in range(n) if a != b and (X >> a & 1) and not (X >> b & 1))
    be = sum(kap1[s][0] for s in range(m1) if X1 >> s & 1); bo = sum(kap1[s][1] for s in range(m1) if X1 >> s & 1)
    ae = sum(lam2[t][0] for t in range(m2) if not X2 >> t & 1); ao = sum(lam2[t][1] for t in range(m2) if not X2 >> t & 1)
    pred = h2[-1] * cut(W1, X1, m1) + h1[-1] * cut(W2, X2, m2) + be * ae + bo * ao
    bad += (cut(W, X, m) != pred)
    bad += (h[-1] != h1[-1] * h2[-1]) + (G[0] != G1[0] * G2[0] + G1[1] * G2[1]) + (G[1] != G1[0] * G2[1] + G1[1] * G2[0])
print("composition law  cut(C1 -> C2) = h2 cut1 + h1 cut2 + b_e a_e + b_o a_o,  h = h1 h2,  G = G1 (x) G2 by parity: failures", bad)
def J(n):
    b = [i for i in range(n.bit_length()) if n >> i & 1]
    return sum(1 if y - x == 1 else 2 ** (y - x - 2) for i, x in enumerate(b) for y in b[i + 1:])
Jc = {}
def Jm(n):
    v = Jc.get(n)
    if v is None: v = Jc[n] = J(n)
    return v
for X in (10 ** 5, 10 ** 6, 10 ** 7):
    two = set(); one = set(); pres = 0
    A = 0
    while True:
        # K is increasing in B for fixed A up to the energy term; scan B while the main term is below ~ 2X
        if 174 * A - 90 * A > 2 * X and A > 0: break
        B = 0; any_hit = False
        while 432 * A * B + 303 * B - 90 * B <= 2 * X + 200:
            K = 432 * A * B + 174 * A + 303 * B + 45 - 30 * (A & 1) - 15 * (B & 1) - 90 * (Jm(A) + Jm(B))
            for d in (90, 92, 94):
                v = K + d
                if 0 < v <= X:
                    two.add(v); pres += 1
                    if B == 0: one.add(v)
            B += 1
        A += 1
        if 84 * A > 2 * X + 200: break
    print("X = %d: two-tape bank reaches %d distinct odd numbers (%.2f%% of the odd numbers; %d presentations); with the upper tape only: %d (%.2f%%)" % (X, len(two), 200.0 * len(two) / X, pres, len(one), 200.0 * len(one) / X))
