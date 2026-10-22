#!/usr/bin/env python3
"""AUDIT pass 2, items 2 and 3: evaluate the packing claims of THM-4751 (7), results note section 13 and MISTAKE-660
from the own alpha tables (a13_alpha_ie_n.json, method Z; cross-checked against method E, a13b_alpha_explicit_n.out,
for n <= 13), and compare with the session's alpha_decomp_*.out.
Output: a13c_alpha_claims.out"""
import os, re, json, glob
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__))
SES = os.path.join(HERE, '..')
out = open(os.path.join(HERE, 'a13c_alpha_claims.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
def load(n):
    fn = os.path.join(HERE, f'a13_alpha_ie_{n}.json')
    if not os.path.exists(fn): return None
    return {tuple(json.loads(k)): v for k, v in json.load(open(fn)).items()}
def explicit(n):
    fn = os.path.join(HERE, f'a13b_alpha_explicit_{n}.out'); d = {}
    if not os.path.exists(fn): return d
    for line in open(fn):
        m = re.match(r'n=\d+ S=\[([\d, ]+)\].*alpha = \[([\d, ]+)\]', line)
        if m: d[tuple(int(x) for x in m.group(1).split(','))] = [int(x) for x in m.group(2).split(',')]
    return d
def session():
    d = {}
    for f in glob.glob(os.path.join(SES, 'alpha_decomp_*.out')):
        cur = None
        for line in open(f):
            m = re.match(r'n=(\d+) S=\{([\d,]+)\}', line)
            if m: cur = (int(m.group(1)), tuple(int(x) for x in m.group(2).split(','))); d[cur] = {}; continue
            m = re.match(r'\s+alpha_(\d+)=(\d+)', line)
            if m and cur: d[cur][int(m.group(1))] = int(m.group(2))
    return d
ses = session()
flat_named = {9: (1, 2, 3, 5), 11: (1, 3, 4, 5, 9), 13: (1, 2, 3, 5, 6, 9), 17: (1, 2, 3, 5, 7, 8, 11, 13),
              19: (1, 4, 5, 6, 7, 9, 11, 16, 17), 21: (1, 2, 3, 7, 10, 12, 13, 15, 16, 17)}
for n in range(7, 22, 2):
    T = load(n)
    if T is None: continue
    U = [a for a in range(1, n) if gcd(a, n) == 1]
    canon = lambda S: min(tuple(sorted((a * s) % n for s in S)) for a in U)
    m = (n - 1) // 2; R = canon(range(1, m + 1)); K = n // 3
    Emin = min(v['E'] for v in T.values()); EminS = [S for S, v in T.items() if v['E'] == Emin]
    ex = explicit(n)
    agreeE = [S for S in ex if S in T and ex[S] == T[S]['alpha']]
    ses_n = {canon(S): v for (nn, S), v in ses.items() if nn == n}
    agreeS = [S for S in ses_n if S in T and all(T[S]['alpha'][k] == ses_n[S][k] for k in ses_n[S])]
    say(f"n={n}: {len(T)} classes in own table; all H checks passed (see a13_alpha_ie_{n}.out); explicit method agrees on {len(agreeE)}/{len(ex)}; "
        f"session alpha_decomp agrees on {len(agreeS)}/{len(ses_n)} classes printed there")
    tag = lambda S: 'R' if S == R else ('Emin' if S in EminS else str(list(S)))
    full = len(T) > 3
    if full:
        for k in range(1, K + 1):
            best = max(v['alpha'][k] for v in T.values()); arg = [S for S, v in T.items() if v['alpha'][k] == best]
            say(f"   alpha_{k}: max {best} at {[tag(S) for S in arg]} (unique: {len(arg) == 1})")
    fl = canon(flat_named[n]) if n in flat_named else None
    if fl in T and R in T:
        a, b = T[R]['alpha'], T[fl]['alpha']
        d12 = 2 * (a[1] - b[1]) + 4 * (a[2] - b[2]); dk = sum((a[k] - b[k]) << k for k in range(3, K + 1))
        say(f"   split vs flat {list(fl)}: 2da1 + 4da2 = {d12:,}; sum_(k>=3) 2^k dak = {dk:,}; dH = {d12 + dk:,}; "
            f"alpha_2(R) = {a[2]:,} {'<' if a[2] < b[2] else '>'} alpha_2(flat) = {b[2]:,}; per-k sign of dak (k>=3): {['+' if a[k] > b[k] else '-' for k in range(3, K + 1)]}")
        if n == 13:
            say(f"   n=13: 8*da3 = {8 * (a[3] - b[3]):,}, 16*da4 = {16 * (a[4] - b[4]):,}; without alpha_4: dH = {d12 + 8 * (a[3] - b[3]):,}")
    if n == 15 and R in T:
        for S in EminS:
            say(f"   n=15 alpha_2(R) = {T[R]['alpha'][2]:,} vs energy minimiser {list(S)}: {T[S]['alpha'][2]:,}")
out.close()
