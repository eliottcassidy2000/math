# Graph-side twin of the tournament TT5 formula: the 4-clique count of the Paley GRAPH G_p (p = 1 mod 4), fitted exactly
# as K4(G_p) = p(p-1) [alpha p^2 + beta p + gamma + delta a_p(eta(4z)^6)] / c and verified on all p = 1 mod 4 < 700.
# (Classical: Evans-Pulham-Sheehan 1981; this script re-derives the constants rather than quoting them.)
import sys, itertools
from fractions import Fraction
import numpy as np
sys.path.insert(0, '.')
from tt5_modular_search import eta_product, solve
from sympy import primerange
a3 = eta_product([(4, 6)], 700)
def k4(p):
    sq = {(x * x) % p for x in range(1, p)}; Q = sorted(sq)
    N01 = np.array([x for x in Q if (x - 1) % p in sq])
    D = (N01[None, :] - N01[:, None]) % p
    through = int(np.isin(D, Q).sum()) // 2           # K4 containing the edge {0,1}
    return through * (p * (p - 1) // 4) // 6
for p in (13, 17, 29, 37):
    sq = {(x * x) % p for x in range(1, p)}
    brute = sum(1 for S in itertools.combinations(range(p), 4) if all((b - a) % p in sq for a, b in itertools.combinations(S, 2)))
    assert brute == k4(p), (p, brute, k4(p))
rows = [(p, k4(p)) for p in primerange(5, 700) if p % 4 == 1]
train = rows[3:7]
A_ = [[Fraction(p) ** 2 * p * (p - 1), Fraction(p) * p * (p - 1), Fraction(p * (p - 1)), Fraction(p * (p - 1) * a3[p])] for p, t in train]
sol = solve(A_, [Fraction(t) for p, t in train])
ok = all(sol[0] * p ** 3 * (p - 1) + sol[1] * p ** 2 * (p - 1) + sol[2] * p * (p - 1) + sol[3] * p * (p - 1) * a3[p] == t for p, t in rows)
den = 1
for s in sol: den = max(den, s.denominator)
print(f"K4(Paley graph G_p) = p(p-1)[({sol[0] * den}) p^2 + ({sol[1] * den}) p + ({sol[2] * den}) + ({sol[3] * den}) a_p(eta(4z)^6)]/{den}"
      f"  -- exact on all {len(rows)} primes p = 1 mod 4 below 700: {ok}")
