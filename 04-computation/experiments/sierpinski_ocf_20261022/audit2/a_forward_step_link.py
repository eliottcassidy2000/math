# Audit2 (A), remark: the forest identity is a substitution in the forward-step polynomial F(T,x) = sum_sigma x^fwd
# (THM-4756):  sum_sigma prod_i (1 + t S_i) = sum_sigma (1+t)^fwd (1-t)^bwd, hence
#     sum_L (n-L)! N_L t^L = (1-t)^(n-1) F(T, (1+t)/(1-t)),
# and t = 1 gives 2^(n-1) H.  Checked exactly (sympy) on random tournaments n <= 7 against direct forest enumeration.
import random, itertools, sys
import sympy as sp
from math import factorial
sys.path.insert(0, '.')
from a_forest_identity import forest_sums, tour_from_bits
t, x = sp.symbols('t x')
random.seed(5)
ok = True
for n in (3, 4, 5, 6, 7):
    for r in range(3):
        m = n * (n - 1) // 2
        A = tour_from_bits(n, [random.randint(0, 1) for _ in range(m)])
        S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]
        N, _, _ = forest_sums(S, n)
        F = 0
        for p in itertools.permutations(range(n)):
            F += x ** sum(A[p[i]][p[i + 1]] for i in range(n - 1))
        lhs = sp.expand(sum(factorial(n - L) * N[L] * t ** L for L in range(n)))
        rhs = sp.expand(sp.cancel((1 - t) ** (n - 1) * F.subs(x, (1 + t) / (1 - t))))
        good = sp.simplify(lhs - rhs) == 0
        ok &= good
        print(f"n={n} #{r}: {'ok' if good else 'MISMATCH'}")
print('ALL OK' if ok else 'FAIL')
