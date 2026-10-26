# Exact span tests for the seven-vertex Paley census against every Galois orbit of 2-power-level newforms
# (candidate_forms.out from candidate_forms.gp), on top of the THM-4780 base span
#   q^i (i <= 7), q(q-1) q^i a_q (i <= 2), q(q-1) q^i |b(q)|^2 (i <= 1).
# Per orbit: q(q-1) q^j Tr(a_p) for j <= floor((6-k)/2), and for k <= 3 also q(q-1) q^j Tr|a_p|^2, j <= floor((7-2k)/2).
# mac-mini-2026-10-25-pyth.
import re, sys, itertools
from fractions import Fraction
rows = []
for fn in ('paley_tt7_data.out', 'paley_tt7_data2.out'):
    for l in open(fn):
        m = re.match(r'p=(\d+) TT7=(\d+)', l)
        if m: rows.append((int(m.group(1)), int(m.group(2))))
rows = sorted(set([(3, 0), (7, 0), (11, 0)] + rows))
ps = [p for p, t in rows]; target = {p: t for p, t in rows}
forms = {}
for l in open('candidate_forms.out'):
    lab, p, tr, nr = l.split(); p = int(p)
    forms.setdefault(lab, {})[p] = (int(tr), int(nr))
def lab_k(lab): return int(re.match(r'N\d+k(\d+)', lab).group(1))
f8 = {p: forms['N8k4c1o1d1'][p][0] for p in ps} if 'N8k4c1o1d1' in forms else None
# a_q(f8), |b|^2 from their own expansions (as in THM-4780's scripts)
sys.argv = [sys.argv[0], 'paley_tt7_data.out']
M = max(ps) + 10
def mul(a, b):
    c = [0] * (M + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(0, M + 1 - i):
                if b[j]: c[i + j] += x * b[j]
    return c
def eta_product(spec):
    shift = sum(s * e for s, e in spec) // 24
    c = [0] * (M + 1); c[0] = 1
    for s, e in spec:
        for _ in range(e):
            for n in range(1, M // s + 1):
                for k in range(M, s * n - 1, -1): c[k] -= c[k - s * n]
    return [0] * shift + c[:M + 1 - shift]
def theta(sc):
    c = [0] * (M + 1); n = 0
    while sc * n * n <= M: c[sc * n * n] += 1 if n == 0 else 2; n += 1
    return c
F8 = eta_product([(2, 4), (4, 4)])
e48 = eta_product([(4, 2), (8, 2)]); t2, t4 = theta(2), theta(4)
imb = mul(e48, [x - y for x, y in zip(mul(t2, t2), mul(t4, t4))])
base = {}
for i in range(8): base[f'q^{i}'] = {p: p ** i for p in ps}
for i in range(3): base[f'q(q-1)q^{i}a'] = {p: p * (p - 1) * p ** i * F8[p] for p in ps}
for i in range(2): base[f'q(q-1)q^{i}|b|^2'] = {p: p * (p - 1) * p ** i * imb[p] ** 2 for p in ps}

def in_span(cols):
    names = list(cols)
    A = [[Fraction(cols[nm][p]) for nm in names] + [Fraction(target[p])] for p in ps]
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

def orbit_cols(lab):
    k = lab_k(lab); out = {}
    for j in range(0, (6 - k) // 2 + 1):
        out[f'{lab}:q^{j}Tr'] = {p: p * (p - 1) * p ** j * forms[lab][p][0] for p in ps}
    if k <= 3:
        for j in range(0, (7 - 2 * k) // 2 + 1):
            out[f'{lab}:q^{j}Nr'] = {p: p * (p - 1) * p ** j * forms[lab][p][1] for p in ps}
    return out
print(f"{len(ps)} values of q, max {max(ps)}; {len(forms)} candidate orbits")
print("base span alone:", 'EXACT' if in_span(base) else 'NOT in span')
hits = 0
for lab in forms:
    cols = dict(base); cols.update(orbit_cols(lab))
    sol = in_span(cols)
    if sol:
        hits += 1
        print("EXACT with", lab, {k: str(v) for k, v in sol.items() if v})
print("single-orbit hits:", hits)
# all orbits of one weight together (when the column count allows)
for k in (2, 3, 4, 5, 6):
    cols = dict(base)
    for lab in forms:
        if lab_k(lab) == k: cols.update(orbit_cols(lab))
    if len(cols) >= len(ps) - 5:
        print(f"weight {k}: {len(cols)} columns, too many for {len(ps)} values"); continue
    sol = in_span(cols)
    print(f"weight {k} all orbits ({len(cols)} columns):", 'EXACT ' + str({kk: str(v) for kk, v in sol.items() if v}) if sol else 'NOT in span')
