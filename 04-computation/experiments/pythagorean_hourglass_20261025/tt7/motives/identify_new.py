# Identify the new seven-vertex motives.  Load the survivor sums (O(q^5) run for q <= 199, O(q^4) run beyond), take
# the genuine sums NOT explained by the base span B = {1, q, q^2, q^3, a84, q a84, a86, a46, |b|^2}, and test which
# candidate newform orbits (tt7/candidate_forms.out: 2-power levels <= 256, weights 2-6, characters 1, chi_8,
# chi_-4, chi_-8; Tr a_p and Tr |a_p|^2) fill the new directions:
#   newdim(C) = rank(B + C + NEW) - rank(B + C); a set C explains everything when newdim(C) = 0.
# Columns per orbit of weight k: q^j Tr for j <= floor((6-k)/2); for k <= 3 also q^j Nr for j <= floor((7-2k)/2).
# mac-mini-2026-10-25-pyth.
import json, sys, re, os, subprocess, tempfile
from sympy import Matrix
G = json.load(open('genuine37.json'))
S = {}
for l in open('sums_7_199.jsonl'):
    d = json.loads(l)
    for i, v in enumerate(d['S']): S.setdefault(i, {})[d['q']] = v
for fn in [f for f in os.listdir('.') if f.startswith('fast_') and f.endswith('.jsonl')]:
    for l in open(fn):
        d = json.loads(l)
        for i, v in d['S'].items(): S.setdefault(int(i), {})[d['q']] = v
NEW = [2, 4, 5, 6, 9, 10, 16, 17, 18, 20, 21, 24, 28, 30, 31]          # H36 (no degree-3 vertex) handled separately
qs = sorted(set.intersection(*[set(S[i]) for i in NEW]))
forms = {}
for l in open('../candidate_forms.out'):
    lab, p, tr, nr = l.split(); forms.setdefault(lab, {})[int(p)] = (int(tr), int(nr))
def lab_k(lab): return int(re.match(r'N\d+k(\d+)', lab).group(1))
P = max(qs) + 2
gpf = os.path.join(tempfile.gettempdir(), 'forms_id.gp')
open(gpf, 'w').write(f'default(parisize,400000000);\nf=mfeigenbasis(mfinit([8,4,1],0))[1];g=mfeigenbasis(mfinit([8,6,1],0))[1];h=mfeigenbasis(mfinit([4,6,1],0))[1];\n'
                     f'print(mfcoefs(f,{P}));print(mfcoefs(g,{P}));print(mfcoefs(h,{P}));\nquit;\n')
f84, f86, f46 = (eval(l) for l in subprocess.run(['gp', '-q', gpf], capture_output=True, text=True).stdout.strip().splitlines()[-3:])
b2 = {p: forms['N32k3c-4o1d2'][p][1] // 2 for p in qs}   # 32.3.c.a orbit (d = 2: b and its conjugate): Nr = 2|b|^2
assert b2[19] == 144 and b2[23] == 1600
base = [[1] * len(qs), qs, [q * q for q in qs], [q ** 3 for q in qs], [f84[q] for q in qs], [q * f84[q] for q in qs],
        [f86[q] for q in qs], [f46[q] for q in qs], [b2[q] for q in qs]]
newcols = [[S[i][q] for q in qs] for i in NEW]
def rank(cols): return Matrix(cols).T.rank() if cols else 0
rB = rank(base); rBN = rank(base + newcols)
print(f"{len(qs)} primes {qs[0]}..{qs[-1]}; rank(B) = {rB}; new directions = {rBN - rB}")
def orbit_cols(lab):
    k = lab_k(lab); out = []
    for j in range(0, (6 - k) // 2 + 1): out.append([q ** j * forms[lab][q][0] for q in qs])
    if k <= 3:
        for j in range(0, (7 - 2 * k) // 2 + 1): out.append([q ** j * forms[lab][q][1] for q in qs])
    return out
results = []
for lab in forms:
    if lab in ('N8k4c1o1d1',): continue
    C = orbit_cols(lab)
    if not any(any(v) for v in C): continue
    rBC = rank(base + C); d = rank(base + C + newcols) - rBC
    results.append((d, lab, rBC - rB))
results.sort()
print("orbits that reduce the new dimension (new dim after adding the orbit, label, columns added):")
for d, lab, c in results:
    if d < rBN - rB: print(f"   {d}  {lab}  (+{c})")
print("best:", results[:5])
