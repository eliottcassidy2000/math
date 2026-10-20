"""Redei for 2-coloured tournaments: P_T(x) = sum over Hamiltonian paths of x^(number of red arcs). Redei: P_T(1) is odd.
What is P_T mod 2 as a polynomial over F_2?  Random 2-coloured tournaments, n <= 8; also the sign-twisted version
Q_T = P_T(-1) = (#even-red paths) - (#odd-red paths)."""
import random, sys
from collections import Counter
from functools import lru_cache
random.seed(5)
def ham_paths_poly(n, orient, colour):
    @lru_cache(maxsize=None)
    def f(mask, last):
        if mask == (1 << n) - 1: return {0: 1}
        tot = {}
        for v in range(n):
            if not (mask >> v) & 1 and orient[last][v]:
                c = colour[last][v]
                for d, coef in f(mask | (1 << v), v).items():
                    tot[d + c] = tot.get(d + c, 0) + coef
        return tot
    P = {}
    for s in range(n):
        for d, coef in f(1 << s, s).items(): P[d] = P.get(d, 0) + coef
    return [P.get(d, 0) for d in range(n)]
NTR = int(sys.argv[1]) if len(sys.argv) > 1 else 300
for n in range(3, 9):
    pats = Counter(); qvals = Counter()
    for _ in range(NTR):
        orient = [[False]*n for _ in range(n)]; colour = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5: orient[i][j] = True
                else: orient[j][i] = True
                c = random.randrange(2); colour[i][j] = colour[j][i] = c
        P = ham_paths_poly(n, orient, colour)
        assert sum(P) % 2 == 1
        pats[tuple(c % 2 for c in P)] += 1
        qvals[sum(c*(-1)**d for d, c in enumerate(P)) % 4] += 1
    oddc = Counter(sum(p) for p in pats.elements())
    print(f"n={n}: {len(pats)} distinct P mod 2 patterns; most common {pats.most_common(3)}; #odd coefficients {dict(sorted(oddc.items()))}; P(-1) mod 4: {dict(sorted(qvals.items()))}")
