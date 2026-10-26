# Search: TT6(P_p) = A(p) [deg 6] + B(p) a_p(f8) [deg <= 2] + C(p) g(p) [deg <= 1], g a CM Hecke-character trace of
# weight 3 or 5 (Q(sqrt-2), Q(sqrt-7), Q(sqrt-3)) or an eta product, solved exactly and verified on all primes.
import re, sys, math
from fractions import Fraction
sys.path.insert(0, '.')
from tt5_modular_search import eta_product, solve
rows = [(int(m.group(1)), int(m.group(2))) for m in (re.match(r'p=(\d+) TT6=(\d+)', l) for l in open('paley_tt6_data.out')) if m]
M = max(p for p, t in rows) + 2
f8 = eta_product([(2, 4), (4, 4)], M)
def cm_trace(p, D, k):
    """a_p = pi^(k-1) + conj for p split in Q(sqrt -D) (class number one), else 0; returns all candidate values for the
    unit ambiguity (one for D = 2, 7; up to units for D = 3)"""
    sols = []
    if D in (2, 1):
        for x in range(0, math.isqrt(p) + 1):
            y2 = (p - x * x)
            if y2 % D == 0 and y2 >= 0 and math.isqrt(y2 // D) ** 2 == y2 // D:
                y = math.isqrt(y2 // D)
                z = complex(x, y * math.sqrt(D))
                sols.append(round(2 * (z ** (k - 1)).real))
    elif D in (3, 7):
        # p = (a^2 + D b^2)/4
        for a in range(0, 2 * math.isqrt(p) + 2):
            r = 4 * p - a * a
            if r >= 0 and r % D == 0 and math.isqrt(r // D) ** 2 == r // D:
                b = math.isqrt(r // D)
                z = complex(a / 2, b * math.sqrt(D) / 2)
                sols.append(round(2 * (z ** (k - 1)).real))
    return sols
cands = {}
for D in (2, 3, 7):
    for k in (3, 5):
        vals = {}
        for p, t in rows:
            s = cm_trace(p, D, k)
            vals[p] = s[0] if s else 0       # first representation; sign/unit ambiguity handled by the d coefficient for D = 2, 7
        cands[f'CM Q(sqrt-{D}) weight {k}'] = vals
eta5 = {'eta(z)^4eta(2z)^2eta(4z)^4': [(1, 4), (2, 2), (4, 4)]}
for name, spec in eta5.items():
    c = eta_product(spec, M); cands[name] = {p: c[p] for p, t in rows}
from sympy import jacobi_symbol
twists = {'1': lambda p: 1, 'chi8': lambda p: jacobi_symbol(8, p), 'chi-8': lambda p: jacobi_symbol(-8 % p, p),
          'chi-3': lambda p: jacobi_symbol(-3 % p, p), 'chi-7': lambda p: jacobi_symbol(-7 % p, p), 'chi12': lambda p: jacobi_symbol(12, p)}
hits = 0
for name, g in cands.items():
    for tw, eps in twists.items():
        for degC in (0, 1):
            nun = 7 + 3 + degC + 1
            train = rows[3:3 + nun]
            A_ = [[Fraction(p) ** k for k in range(7)] + [Fraction(p) ** k * f8[p] for k in range(3)] + [Fraction(p) ** k * eps(p) * g[p] for k in range(degC + 1)] for p, t in train]
            sol = solve(A_, [Fraction(t) for p, t in train])
            if sol is None: continue
            def F(p): return (sum(sol[k] * Fraction(p) ** k for k in range(7)) + sum(sol[7 + k] * Fraction(p) ** k for k in range(3)) * f8[p]
                              + sum(sol[10 + k] * Fraction(p) ** k for k in range(degC + 1)) * eps(p) * g[p])
            nok = sum(1 for p, t in rows if F(p) == t)
            if nok == len(rows):
                hits += 1
                print(f"EXACT ({len(rows)} primes): g = {name}, twist {tw}, deg C = {degC}: A = {[str(s) for s in sol[:7]]}, B = {[str(s) for s in sol[7:10]]}, C = {[str(s) for s in sol[10:]]}")
print("hits:", hits)
