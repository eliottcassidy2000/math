# Search for an exact identity TT5(P_p) = sum_(k<=5) c_k p^k + d p^j a_f(p) over primes p = 3 (mod 4), a_f the p-th
# coefficient of an eta-product newform (weights 2, 3, 4).  Unknowns c_0..c_5, d solved exactly (Fractions) from 7
# primes; the identity is accepted only if it holds exactly at all other primes.  mac-mini-2026-10-25-pyth.
import re, sys
from fractions import Fraction
rows = [(int(m.group(1)), int(m.group(2))) for m in (re.match(r'p=(\d+) TT5=(\d+)', l) for l in open('paley_tt5_data.out')) if m]
M = max(p for p, t in rows) + 2

def eta_product(spec, M):
    """spec: list of (scale, power); returns coefficients of prod eta(scale z)^power as a q-series starting at q^shift"""
    shift = sum(s * e for s, e in spec)
    assert shift % 24 == 0
    shift //= 24
    c = [0] * (M + 1); c[0] = 1
    for s, e in spec:
        # multiply by prod_n (1 - q^(s n))^e
        for _ in range(abs(e)):
            for n in range(1, M // s + 1):
                step = s * n
                if e > 0:
                    for k in range(M, step - 1, -1): c[k] -= c[k - step]
                else:
                    for k in range(step, M + 1): c[k] += c[k - step]
    out = [0] * (M + 1)
    for k in range(M + 1 - shift): out[k + shift] = c[k]
    return out

forms = {
    'w2 N11 eta(z)^2eta(11z)^2': [(1, 2), (11, 2)],
    'w2 N14': [(1, 1), (2, 1), (7, 1), (14, 1)],
    'w2 N15': [(1, 1), (3, 1), (5, 1), (15, 1)],
    'w2 N20': [(2, 2), (10, 2)],
    'w2 N24': [(2, 1), (4, 1), (6, 1), (12, 1)],
    'w2 N27 (CM -3)': [(3, 2), (9, 2)],
    'w2 N32 (CM -4)': [(4, 2), (8, 2)],
    'w2 N36 (CM -3)': [(6, 4)],
    'w3 N7 (CM -7)': [(1, 3), (7, 3)],
    'w3 N8 (CM -2)': [(1, 2), (2, 1), (4, 1), (8, 2)],
    'w3 N12 (CM -3)': [(2, 3), (6, 3)],
    'w3 N16 (CM -4)': [(4, 6)],
    'w4 N5': [(1, 4), (5, 4)],
    'w4 N6': [(1, 2), (2, 2), (3, 2), (6, 2)],
    'w4 N8': [(2, 4), (4, 4)],
    'w4 N9 (CM -3)': [(3, 8)],
}
coeffs = {name: eta_product(spec, M) for name, spec in forms.items()}

def solve(A, b):
    n = len(A); M_ = [row[:] + [bb] for row, bb in zip(A, b)]
    for i in range(n):
        piv = next((r for r in range(i, n) if M_[r][i] != 0), None)
        if piv is None: return None
        M_[i], M_[piv] = M_[piv], M_[i]
        for r in range(n):
            if r != i and M_[r][i] != 0:
                f = M_[r][i] / M_[i][i]
                M_[r] = [x - f * y for x, y in zip(M_[r], M_[i])]
    return [M_[i][n] / M_[i][i] for i in range(n)]

hits = []
train = rows[3:10]          # avoid the smallest primes in the solve
for name, a in coeffs.items():
    for j in range(0, 4):
        A = [[Fraction(p) ** k for k in range(6)] + [Fraction(p) ** j * a[p]] for p, t in train]
        b = [Fraction(t) for p, t in train]
        sol = solve(A, b)
        if sol is None: continue
        ok = all(sum(sol[k] * Fraction(p) ** k for k in range(6)) + sol[6] * Fraction(p) ** j * a[p] == t for p, t in rows)
        if ok:
            hits.append((name, j, sol))
            print(f"EXACT: TT5(P_p) = {' + '.join(f'({sol[k]}) p^{k}' for k in range(6))} + ({sol[6]}) p^{j} a(p)   [{name}]")
print(f"searched {len(coeffs)} forms x 4 powers; exact identities found: {len(hits)}")
