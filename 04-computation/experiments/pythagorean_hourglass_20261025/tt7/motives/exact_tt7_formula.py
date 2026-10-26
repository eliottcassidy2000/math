# Exact O(q^2) formula for the seven-vertex Paley census.  Each genuine sum S_H is fitted exactly as
#   S_H = (combination of the base B = 1, q, q^2, q^3, a84, q a84, a86, a46, |b|^2) + (combination of Y02, Y05, Y09, Y18)
# over all primes where S_H is known (O(q^5)/O(q^4) data, 7..459); then g(q) = sum m_H S_H and
#   2^21 TT7 = A7 + B7 T_W4 + C7 T_(K6-P6) + q(q-1) g(q)
# is checked against the 72 directly counted values TT7(P_q), q <= 811.  mac-mini-2026-10-25-pyth.
import json, os, re, subprocess, tempfile
from fractions import Fraction
G = json.load(open('genuine37.json'))
S = {}
for l in open('sums_7_199.jsonl'):
    d = json.loads(l)
    for i, v in enumerate(d['S']): S.setdefault(i, {})[d['q']] = v
for fn in [f for f in os.listdir('.') if f.startswith('fast_') and f.endswith('.jsonl')]:
    for l in open(fn):
        d = json.loads(l)
        for i, v in d['S'].items(): S.setdefault(int(i), {})[d['q']] = v
Yrows = {}
for fn in sorted(f for f in os.listdir('.') if f.startswith('cheapY_') and f.endswith('.txt')):
    for l in open(fn):
        r = list(map(int, l.split())); Yrows[r[0]] = {'Y02': r[1], 'Y05': r[3], 'Y09': r[5], 'Y18': r[6]}
P = max(Yrows) + 2
gpf = os.path.join(tempfile.gettempdir(), 'forms_ex.gp')
open(gpf, 'w').write(f'default(parisize,1000000000);\nf=mfeigenbasis(mfinit([8,4,1],0))[1];g=mfeigenbasis(mfinit([8,6,1],0))[1];h=mfeigenbasis(mfinit([4,6,1],0))[1];b=mfeigenbasis(mfinit([32,3,-4],0))[1];\n'
                     f'print(mfcoefs(f,{P}));print(mfcoefs(g,{P}));print(mfcoefs(h,{P}));c=mfcoefs(b,{P});print(vector({P}+1,n,norml2(subst(lift(c[n]),y,I))));\nquit;\n')
out = subprocess.run(['gp', '-q', gpf], capture_output=True, text=True).stdout.strip().splitlines()
f84, f86, f46, B2 = (eval(l.replace('I', '1j')) for l in out[-4:])
B2 = [int(round(abs(x))) for x in B2]
def basecols(q): return {'1': 1, 'q': q, 'q^2': q * q, 'q^3': q ** 3, 'a84': f84[q], 'q a84': q * f84[q], 'a86': f86[q], 'a46': f46[q], '|b|^2': B2[q]}
names = ['1', 'q', 'q^2', 'q^3', 'a84', 'q a84', 'a86', 'a46', '|b|^2', 'Y02', 'Y05', 'Y09', 'Y18']
def col(q): c = basecols(q); c.update(Yrows[q]); return c
def solve(target, qs):
    A = [[Fraction(col(q)[n]) for n in names] + [Fraction(target[q])] for q in qs]
    m, n = len(A), len(names); r = 0; piv = []
    for c in range(n + 1):
        pr = next((i for i in range(r, m) if A[i][c] != 0), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]; inv = 1 / A[r][c]; A[r] = [x * inv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                fct = A[i][c]; A[i] = [x - fct * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    if n in piv: return None
    return {names[c]: A[i][n] for i, c in enumerate(piv)}
genuine = [i for i, H in enumerate(G) if not H['all_even']]
gform = {n: Fraction(0) for n in names}
ok = True
for i in genuine:
    qs = sorted(q for q in S[i] if q in Yrows)
    sol = solve(S[i], qs)
    if sol is None:
        print(f"H{i:02d}: NOT in base + Y span ({len(qs)} primes)"); ok = False; continue
    for n, v in sol.items(): gform[n] += G[i]['m'] * v
    nz = {k: str(v) for k, v in sol.items() if v}
    print(f"H{i:02d} (m = {G[i]['m']:+d}, {len(qs)} primes): " + ' + '.join(f'({v})*{k}' for k, v in nz.items()))
print("\ng(q) = sum_H m_H S_H = " + ' + '.join(f'({v})*{k}' for k, v in gform.items() if v))
# check against the counted TT7 values
TT7 = {}
for fn in ('../paley_tt7_data.out', '../paley_tt7_data2.out'):
    for l in open(fn):
        m = re.match(r'p=(\d+) TT7=(\d+)', l)
        if m: TT7[int(m.group(1))] = int(m.group(2))
TT7.update({3: 0, 7: 0, 11: 0})
A7 = lambda q: q * (q - 1) * (q - 3) * (q ** 4 - 87 * q ** 3 + 2785 * q ** 2 - 38827 * q + 198842)
B7 = lambda q: 28 * q * q - 1288 * q + 13238
C7 = lambda q: 1912 - 64 * q
good = bad = 0
for q, t in sorted(TT7.items()):
    if q not in Yrows or q < 7: continue
    gq = sum(gform[n] * col(q)[n] for n in names)
    TW4 = q * (q - 1) * (-4 - f84[q]); TP6 = q * (q - 1) * (B2[q] - q * q - 7)
    rhs = A7(q) + B7(q) * TW4 + C7(q) * TP6 + q * (q - 1) * gq
    if rhs == 2 ** 21 * t: good += 1
    else: bad += 1
print(f"O(q^2) formula for 2^21 TT7 against direct counts: {good} agree, {bad} disagree (q in 7..{max(TT7)})")
json.dump({k: str(v) for k, v in gform.items()}, open('g_formula.json', 'w'))
