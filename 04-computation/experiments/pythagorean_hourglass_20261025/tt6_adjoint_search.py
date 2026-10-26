# Six-vertex Paley census: is TT6(P_p) (p = 3 mod 4) an exact combination of polynomials times
#   a_p(f8)      f8 = eta(2z)^4 eta(4z)^4                (weight 4, level 8: the five-vertex Apery form),
#   |b(p)|^2     b = LMFDB 32.3.c.a, the weight-3 level-32 chi_-4 newform (non-CM) of HYP-9276 (K6-P6 adjoint),
#   c6(p)        the weight-6 level-8 newform 8.6.a.a,
# with the polynomial part allowed to depend on (2/p)?  Exact rational column-span test on all primes 19..751.
# mac-mini-2026-10-25-pyth.
import re, sys
from fractions import Fraction
from sympy import jacobi_symbol

rows = [(int(m.group(1)), int(m.group(2))) for m in (re.match(r'p=(\d+) TT6=(\d+)', l) for l in open('paley_tt6_data.out')) if m]
M = max(p for p, t in rows) * 3 + 10

def mul(a, b):
    c = [0] * (M + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(0, M + 1 - i):
                if b[j]: c[i + j] += x * b[j]
    return c

def eta_product(spec):
    shift = sum(s * e for s, e in spec); assert shift % 24 == 0; shift //= 24
    c = [0] * (M + 1); c[0] = 1
    for s, e in spec:
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

def theta(scale):
    c = [0] * (M + 1); n = 0
    while scale * n * n <= M:
        c[scale * n * n] += 1 if n == 0 else 2; n += 1
    return c

f8 = eta_product([(2, 4), (4, 4)])
# b = eta(4z)^2 eta(8z)^2 (theta(4z)^2 + i (theta(2z)^2 - theta(4z)^2)); at p = 3 mod 4 only the i-part survives
e48 = eta_product([(4, 2), (8, 2)])
t2, t4 = theta(2), theta(4)
t2s, t4s = mul(t2, t2), mul(t4, t4)
re_b = mul(e48, t4s); im_b = mul(e48, [x - y for x, y in zip(t2s, t4s)])
assert (re_b[1], im_b[1], re_b[3], im_b[3], re_b[5], im_b[5], re_b[7], im_b[7]) == (1, 0, 0, 4, 2, 0, 0, -8), 'b expansion'

# weight-6 level-8: S_6(Gamma0(8)) = f8 * M_2(Gamma0(8)), M_2 spanned by E2(z) - d E2(dz), d = 2, 4, 8
def sigma(n): return sum(d for d in range(1, n + 1) if n % d == 0)
E2 = [1] + [-24 * sigma(n) for n in range(1, M + 1)]
def E2d(d):
    c = [0] * (M + 1)
    for n in range(0, M // d + 1): c[d * n] = E2[n]
    return c
basis6 = [mul(f8, [x - d * y for x, y in zip(E2, E2d(d))]) for d in (2, 4, 8)]
# old forms: eta(2z)^12 and eta(4z)^12 (a_3 = -12); find the T_3-eigenvector with another eigenvalue
import sympy as sp
K = 12
B = sp.Matrix([[f[n] for n in range(1, K + 1)] for f in basis6])          # 3 x K
def T3(f): return [f[3 * n] + (3 ** 5 * f[n // 3] if n % 3 == 0 else 0) for n in range(0, M // 3 + 1)]
TB = sp.Matrix([[T3(f)[n] for n in range(1, K + 1)] for f in basis6])
# solve TB = Mt * B for the 3x3 matrix Mt (row space)
Mt = (TB * B.T) * (B * B.T).inv()
assert Mt * B == TB
ev = Mt.T.eigenvects()
new = None
for val, mult, vecs in ev:
    if val != -12:
        assert mult == 1; v = vecs[0]
        cand = [sum(v[i] * basis6[i][n] for i in range(3)) for n in range(M + 1)]
        c1 = cand[1]; c6 = [sp.nsimplify(x / c1) for x in cand]
        new = (val, c6)
print("weight-6 level-8 newform: T_3 eigenvalue", new[0], "q-expansion", [int(new[1][n]) for n in range(1, 12)])
c6 = [int(x) for x in new[1][:M // 3 + 1]]
# Hecke check of c6 at p = 5, 7: a(p^2) = a(p)^2 - p^5, multiplicativity a(15) = a(3) a(5)
assert c6[25] == c6[5] ** 2 - 5 ** 5 and c6[15] == c6[3] * c6[5] and c6[21] == c6[3] * c6[7]
from sympy import isprime
assert all(f8[p * p] == f8[p] ** 2 - p ** 3 for p in (3, 5, 7, 11, 13))

ps = [p for p, t in rows]; ts = [Fraction(t) for p, t in rows]
lam = {p: im_b[p] ** 2 for p in ps}
chi8 = {p: int(jacobi_symbol(2, p)) for p in ps}
cols_all = {}
for k in range(7): cols_all[f'p^{k}'] = {p: p ** k for p in ps}
for k in range(7): cols_all[f'(2/p) p^{k}'] = {p: chi8[p] * p ** k for p in ps}
for k in range(4): cols_all[f'p^{k} a(f8)'] = {p: p ** k * f8[p] for p in ps}
for k in range(4): cols_all[f'p^{k} |b|^2'] = {p: p ** k * lam[p] for p in ps}
for k in range(3): cols_all[f'p^{k} c6'] = {p: p ** k * c6[p] for p in ps}
for k in range(3): cols_all[f'(2/p) p^{k} a(f8)'] = {p: chi8[p] * p ** k * f8[p] for p in ps}
for k in range(3): cols_all[f'(2/p) p^{k} |b|^2'] = {p: chi8[p] * p ** k * lam[p] for p in ps}
cols_all['a(f8)^2'] = {p: f8[p] ** 2 for p in ps}

def in_span(names):
    A = [[Fraction(cols_all[nm][p]) for nm in names] + [t] for p, t in zip(ps, ts)]
    m, n = len(A), len(names)
    r = 0; piv = []
    for c in range(n + 1):
        pr = next((i for i in range(r, m) if A[i][c] != 0), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]
        inv = 1 / A[r][c]; A[r] = [x * inv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    if n in piv: return None          # target not in span
    sol = {names[c]: A[i][n] for i, c in enumerate(piv)}
    return sol, len(piv)

tests = {
 'poly6 + f8(deg2)': ['p^%d' % k for k in range(7)] + ['p^%d a(f8)' % k for k in range(3)],
 'poly6 + f8(deg2) + |b|^2(deg1)': ['p^%d' % k for k in range(7)] + ['p^%d a(f8)' % k for k in range(3)] + ['p^%d |b|^2' % k for k in range(2)],
 'poly6 + f8(deg3) + |b|^2(deg2)': ['p^%d' % k for k in range(7)] + ['p^%d a(f8)' % k for k in range(4)] + ['p^%d |b|^2' % k for k in range(3)],
 'poly6 + f8(deg2) + |b|^2(deg1) + c6(deg1)': ['p^%d' % k for k in range(7)] + ['p^%d a(f8)' % k for k in range(3)] + ['p^%d |b|^2' % k for k in range(2)] + ['p^%d c6' % k for k in range(2)],
 'poly6+-(2/p) + f8 + |b|^2 + c6': ['p^%d' % k for k in range(7)] + ['(2/p) p^%d' % k for k in range(7)] + ['p^%d a(f8)' % k for k in range(3)] + ['p^%d |b|^2' % k for k in range(2)] + ['p^%d c6' % k for k in range(2)],
 'everything': list(cols_all),
}
for name, names in tests.items():
    res = in_span(names)
    print(f"{name:45s} ({len(names)} columns, {len(ps)} primes): ", end='')
    if res is None: print("NOT in span")
    else:
        sol, rk = res
        print(f"EXACT (rank {rk}); nonzero:", {k: str(v) for k, v in sol.items() if v != 0})
