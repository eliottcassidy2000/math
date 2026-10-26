# Exact formula for the transitive 5-subtournaments of the Paley tournament:
#   TT5(P_p) = A(p) + B(p) a_p(f8),  f8 = eta(2z)^4 eta(4z)^4 (weight 4, level 8),
# with A of degree 5 and B of degree <= 2, solved exactly on 9 primes and verified exactly on all primes p = 3 (mod 4)
# below 1500 (paley_tt5_data.out).  mac-mini-2026-10-25-pyth.
import re, sys
from fractions import Fraction
sys.path.insert(0, '.')
from tt5_modular_search import eta_product, solve
rows = [(int(m.group(1)), int(m.group(2))) for m in (re.match(r'p=(\d+) TT5=(\d+)', l) for l in open('paley_tt5_data.out')) if m]
M = max(p for p, t in rows) + 2
a = eta_product([(2, 4), (4, 4)], M)
for degB in (0, 1, 2):
    nunk = 6 + degB + 1
    train = rows[4:4 + nunk]
    A_ = [[Fraction(p) ** k for k in range(6)] + [Fraction(p) ** k * a[p] for k in range(degB + 1)] for p, t in train]
    b_ = [Fraction(t) for p, t in train]
    sol = solve(A_, b_)
    if sol is None: print('singular', degB); continue
    def F(p): return sum(sol[k] * Fraction(p) ** k for k in range(6)) + sum(sol[6 + k] * Fraction(p) ** k for k in range(degB + 1)) * a[p]
    ok = all(F(p) == t for p, t in rows)
    print(f"deg B = {degB}: exact on all {len(rows)} primes: {ok}")
    if ok:
        print("  A(p) =", ' + '.join(f"({sol[k]}) p^{k}" for k in range(6)))
        print("  B(p) =", ' + '.join(f"({sol[6 + k]}) p^{k}" for k in range(degB + 1)))
        # factor A and B with sympy
        import sympy as sp
        x = sp.symbols('p')
        Ap = sum(sp.Rational(sol[k].numerator, sol[k].denominator) * x ** k for k in range(6))
        Bp = sum(sp.Rational(sol[6 + k].numerator, sol[6 + k].denominator) * x ** k for k in range(degB + 1))
        print("  A(p) factored:", sp.factor(Ap)); print("  B(p) factored:", sp.factor(Bp))
        break
