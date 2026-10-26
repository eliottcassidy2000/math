# Isolate the genuine seven-vertex part of the Paley census.  With the audit's exact reduced form
#   2^21 TT7 = A7 + B7 T_W4 + C7 T_(K6-P6) + G,   G = sum over the 37 genuine dense sums m_H T_H,
#   A7 = q(q-1)(q-3)(q^4 - 87q^3 + 2785q^2 - 38827q + 198842), B7 = 28q^2 - 1288q + 13238, C7 = 1912 - 64q,
# and T_W4 = q(q-1)(-4 - a_q) (HYP-9274 (B)), T_(K6-P6) = q(q-1)(|b|^2 - q^2 - 7) (THM-4653 + HYP-9276, FINITE-EXACT
# far beyond these q), compute g(q) = G/(q(q-1)) and test it against polynomials plus candidate newform orbits.
# mac-mini-2026-10-25-pyth.
import re, sys, math
from fractions import Fraction
sys.argv = [sys.argv[0]]
exec(open('tt7_candidate_fit.py').read().split("print(f\"{len(ps)} values")[0])
A7 = lambda q: q * (q - 1) * (q - 3) * (q ** 4 - 87 * q ** 3 + 2785 * q ** 2 - 38827 * q + 198842)
B7 = lambda q: 28 * q * q - 1288 * q + 13238
C7 = lambda q: 1912 - 64 * q
g = {}
for p in ps:
    TW4 = p * (p - 1) * (-4 - F8[p]); TP6 = p * (p - 1) * (imb[p] ** 2 - p * p - 7)
    G = 2 ** 21 * target[p] - A7(p) - B7(p) * TW4 - C7(p) * TP6
    assert G % (p * (p - 1)) == 0, p
    g[p] = G // (p * (p - 1))
print("g(q) = G/(q(q-1)) at the first values:", [(p, g[p]) for p in ps[:8]])
print("audit check G(7) = 25536, G(31) = 861522240:", g[7] * 42 == 25536, g[31] * 31 * 30 == 861522240)
print("g(q)/q^(5/2):", [round(g[p] / p ** 2.5, 3) for p in ps[:40]])
# divisibility patterns
for d in (2, 4, 8, 16, 32, 64, 3, 5, 7):
    print(f"  g = 0 mod {d} at {sum(1 for p in ps if g[p] % d == 0)}/{len(ps)} q")
def span_g(cols):
    names = list(cols)
    A = [[Fraction(cols[nm][p]) for nm in names] + [Fraction(g[p])] for p in ps]
    m, n = len(A), len(names); r = 0; piv = []
    for c in range(n + 1):
        pr = next((i for i in range(r, m) if A[i][c] != 0), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]; inv = 1 / A[r][c]; A[r] = [x * inv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    if n in piv: return None
    return {names[c]: A[i][n] for i, c in enumerate(piv)}
poly = {f'q^{i}': {p: p ** i for p in ps} for i in range(4)}
extra_base = {'a': {p: F8[p] for p in ps}, 'q a': {p: p * F8[p] for p in ps}, '|b|^2': {p: imb[p] ** 2 for p in ps}}
def orbit_cols_g(lab):
    k = lab_k(lab); out = {}
    for j in range(0, (6 - k) // 2 + 1):
        out[f'{lab}:q^{j}Tr'] = {p: p ** j * forms[lab][p][0] for p in ps}
    if k <= 3:
        for j in range(0, (7 - 2 * k) // 2 + 1):
            out[f'{lab}:q^{j}Nr'] = {p: p ** j * forms[lab][p][1] for p in ps}
    return out
base = dict(poly); base.update(extra_base)
print("g in span(poly deg<=3, a, qa, |b|^2):", 'EXACT' if span_g(base) else 'NOT')
hits = 0
for lab in forms:
    c = dict(base); c.update(orbit_cols_g(lab))
    s = span_g(c)
    if s: hits += 1; print("EXACT with", lab, {k: str(v) for k, v in s.items() if v})
print("single-orbit hits:", hits)
for group in [(6,), (5,), (4,), (3,), (2,), (5, 6), (4, 6), (3, 6), (3, 5), (2, 6)]:
    c = dict(base)
    for lab in forms:
        if lab_k(lab) in group: c.update(orbit_cols_g(lab))
    if len(c) >= len(ps) - 3:
        print(f"weights {group}: {len(c)} columns, too many"); continue
    s = span_g(c)
    print(f"weights {group} jointly ({len(c)} columns):", ('EXACT ' + str({k: str(v) for k, v in s.items() if v})) if s else 'NOT')
