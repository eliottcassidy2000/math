# Exact search: T_{K5-e}(p) = c0 + c1 p + c2 p^2 + d * eps(p) * a_p(f) over all odd primes p in the data, with f an
# eta-product newform (weights 2, 3, 4; weight 2 also times p) and eps a quadratic twist (1, chi_-4, chi_8, chi_-8,
# chi_-3, chi_12).  Also checks T_K4 = a_p(eta(4z)^6) and T_K5 = (p - 4) T_K4.  mac-mini-2026-10-25-pyth.
import re, sys
from fractions import Fraction
sys.path.insert(0, '.')
from tt5_modular_search import eta_product, forms as eta_forms, solve
rows = []
for l in open('paley_graph_sums.out'):
    m = re.match(r'p=(\d+) mod4=\d K4=(-?\d+) K5=(-?\d+) K5-e=(-?\d+) K5-2K2=(-?\d+)', l)
    if m: rows.append(tuple(int(v) for v in m.groups()))
M = max(r[0] for r in rows) + 2
f4 = eta_product([(4, 6)], M)
print("T_K4 = a_p(eta(4z)^6) for all p:", all(r[1] == f4[r[0]] for r in rows),
      "; T_K5 = (p - 4) T_K4:", all(r[2] == (r[0] - 4) * r[1] for r in rows),
      "; T_(K5-e) = T_(K5-2K2):", all(r[3] == r[4] for r in rows))
extra = {'w4 N5': [(1, 4), (5, 4)], 'w4 N6': [(1, 2), (2, 2), (3, 2), (6, 2)], 'w4 N8': [(2, 4), (4, 4)], 'w4 N9': [(3, 8)]}
allforms = dict(eta_forms); allforms.update(extra)
coeffs = {name: eta_product(spec, M) for name, spec in allforms.items()}
def kron(D, p):
    from sympy import jacobi_symbol
    return jacobi_symbol(D % p, p) if D % p else 0
twists = {'1': lambda p: 1, 'chi_-4': lambda p: kron(-4, p), 'chi_8': lambda p: kron(8, p), 'chi_-8': lambda p: kron(-8, p),
          'chi_-3': lambda p: kron(-3, p), 'chi_12': lambda p: kron(12, p)}
data = [(r[0], r[3]) for r in rows if r[0] > 5]
hits = 0
for name, a in coeffs.items():
    for tw, eps in twists.items():
        for j in (0, 1):
            train = data[2:6]
            A = [[Fraction(1), Fraction(p), Fraction(p) ** 2, Fraction(eps(p) * a[p] * p ** j)] for p, t in train]
            b = [Fraction(t) for p, t in train]
            sol = solve(A, b)
            if sol is None: continue
            ok = all(sol[0] + sol[1] * p + sol[2] * p * p + sol[3] * eps(p) * a[p] * p ** j == t for p, t in data)
            if ok:
                hits += 1
                print(f"EXACT: T_(K5-e)(p) = {sol[0]} + {sol[1]} p + {sol[2]} p^2 + ({sol[3]}) {tw}(p) p^{j} a_p[{name}]")
print("exact identities:", hits)
