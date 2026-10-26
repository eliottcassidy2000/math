# Does the seven-vertex Paley census carry new arithmetic?  THM-4780's lemmas give
#   2^21 TT7 = A7(q) + B7(q) q(q-1) S_W4 + C7(q) q(q-1) S_(K6-P6) + sum over the 37 genuine 7-vertex sums,
# deg A7 <= 7, deg B7 <= 2, deg C7 <= 1.  Test exact membership of TT7 in the span of q^i, q(q-1)q^i a_q (i <= 2),
# q(q-1)q^i |b(q)|^2 (i <= 1); if it fails, try extra candidate columns.  mac-mini-2026-10-25-pyth.
import re, sys
from fractions import Fraction
rows = [(int(m.group(1)), int(m.group(2))) for m in (re.match(r'p=(\d+) TT7=(\d+)', l) for l in open(sys.argv[1] if len(sys.argv) > 1 else 'paley_tt7_data.out')) if m]
rows = [(3, 0), (7, 0), (11, 0)] + rows
P = max(p for p, t in rows)
M = P + 10
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
        for _ in range(e):
            for n in range(1, M // s + 1):
                for k in range(M, s * n - 1, -1): c[k] -= c[k - s * n]
    out = [0] * (M + 1)
    for k in range(M + 1 - shift): out[k + shift] = c[k]
    return out
def theta(sc):
    c = [0] * (M + 1); n = 0
    while sc * n * n <= M: c[sc * n * n] += 1 if n == 0 else 2; n += 1
    return c
f8 = eta_product([(2, 4), (4, 4)])
e48 = eta_product([(4, 2), (8, 2)]); t2, t4 = theta(2), theta(4)
imb = mul(e48, [x - y for x, y in zip(mul(t2, t2), mul(t4, t4))])
assert (f8[3], f8[7], imb[3], imb[7]) == (-4, 24, 4, -8)
e12 = eta_product([(2, 12)])                        # 4.6.a.a
cm3 = eta_product([(1, 2), (2, 1), (4, 1), (8, 2)])  # 8.3.d.a, CM by Q(sqrt-2)
def sigma(n): return sum(d for d in range(1, n + 1) if n % d == 0)
E2 = [1] + [-24 * sigma(n) for n in range(1, M + 1)]
def E2d(d):
    c = [0] * (M + 1)
    for n in range(0, M // d + 1): c[d * n] = E2[n]
    return c
# weight-6 level-8 newform 8.6.a.a: f8 * (E2(z) - 2E2(2z)) etc.; the T_3 = 20 eigenvector, normalised
import sympy as sp
basis6 = [mul(f8, [x - d * y for x, y in zip(E2, E2d(d))]) for d in (2, 4, 8)]
K = 12
Bm = sp.Matrix([[f[n] for n in range(1, K + 1)] for f in basis6])
TB = sp.Matrix([[(f[3 * n] if 3 * n <= M else 0) + (3 ** 5 * f[n // 3] if n % 3 == 0 else 0) for n in range(1, K + 1)] for f in basis6])
Mt = (TB * Bm.T) * (Bm * Bm.T).inv()
vec = [v for val, mult, vs in Mt.T.eigenvects() if val == 20 for v in vs][0]
c6 = [sum(vec[i] * basis6[i][n] for i in range(3)) for n in range(M + 1)]
c6 = [int(sp.nsimplify(x / c6[1])) for x in c6]
assert c6[3] == 20 and c6[5] == -74
cols = {}
for i in range(8): cols[f'q^{i}'] = lambda q, i=i: q ** i
for i in range(3): cols[f'q(q-1)q^{i} a'] = lambda q, i=i: q * (q - 1) * q ** i * f8[q]
for i in range(2): cols[f'q(q-1)q^{i} |b|^2'] = lambda q, i=i: q * (q - 1) * q ** i * imb[q] ** 2
extra = {'q(q-1) c6': lambda q: q * (q - 1) * c6[q], 'q(q-1) eta(2z)^12': lambda q: q * (q - 1) * e12[q],
         'q(q-1) a^2': lambda q: q * (q - 1) * f8[q] ** 2, 'q(q-1) cm3(-2)^2': lambda q: q * (q - 1) * cm3[q] ** 2,
         'q(q-1) q cm3(-2)': lambda q: q * (q - 1) * q * cm3[q], 'q(q-1) a |b|^2/q': lambda q: (q - 1) * f8[q] * imb[q] ** 2}
def in_span(names, cf):
    A = [[Fraction(cf[nm](p)) for nm in names] + [Fraction(t)] for p, t in rows]
    m, n = len(A), len(names); r = 0; piv = []
    for c in range(n + 1):
        pr = next((i for i in range(r, m) if A[i][c] != 0), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]; inv = 1 / A[r][c]; A[r] = [x * inv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    if n in piv: return None, len(piv) - 1
    return {names[c]: A[i][n] for i, c in enumerate(piv)}, len(piv)
base = list(cols)
sol, rk = in_span(base, cols)
print(f"{len(rows)} values of q (max {P}); base span: {len(base)} columns, rank {rk}:", 'EXACT ' + str({k: str(v) for k, v in sol.items() if v}) if sol else 'NOT in span')
if not sol:
    allc = dict(cols); allc.update(extra)
    for nm in extra:
        s2, rk2 = in_span(base + [nm], allc)
        print(f"  + {nm}: {'EXACT ' + str({k: str(v) for k, v in s2.items() if v}) if s2 else 'NOT in span'}")
    s3, rk3 = in_span(base + list(extra), allc)
    print(f"  + all extras ({len(base) + len(extra)} columns): {'EXACT ' + str({k: str(v) for k, v in s3.items() if v}) if s3 else 'NOT in span'}")
