# Do tensor products / symmetric powers of 2-power-level newforms fill the four new seven-vertex directions?
# Candidates (traces at p = 3 mod 4): products Tr_f Tr_g q^j of two orbits with motivic weight w = k_f + k_g - 2 <= 5
# (j <= (5 - w)/2), and Tr_f^2, Tr_f^3 for weight-2 orbits.  Same rank test as identify_new.py.  mac-mini-2026-10-25-pyth.
import json, re, os, itertools
from sympy import Matrix
exec(open('identify_new.py').read().split('results = []')[0].replace('print(f"{len(qs)} primes', 'pass  # print(f"{len(qs)} primes'))
real = {}
for lab in forms:
    k = lab_k(lab)
    tr = [forms[lab][q][0] for q in qs]
    if any(tr) and lab != 'N32k3c-4o1d2': real[lab] = (k, tr)
real['f84'] = (4, [f84[q] for q in qs]); real['f86'] = (6, [f86[q] for q in qs]); real['f46'] = (6, [f46[q] for q in qs])
real['h32'] = (3, None)
print(f"{len(qs)} primes; {len(real)} orbits with nonzero traces at p = 3 mod 4; new dimension {rBN - rB}")
cands = []
labs = [l for l in real if l != 'h32']
for f, g in itertools.combinations_with_replacement(labs, 2):
    kf, tf = real[f]; kg, tg = real[g]
    w = kf + kg - 2
    if w > 5: continue
    for j in range(0, (5 - w) // 2 + 1):
        cands.append((f'{f}*{g}*q^{j}', [a * b * q ** j for a, b, q in zip(tf, tg, qs)]))
for f in labs:
    kf, tf = real[f]
    if kf == 2:
        for j in range(0, 2): cands.append((f'{f}^3*q^{j}', [a ** 3 * q ** j for a, q in zip(tf, qs)]))
print(len(cands), "product candidates")
hits = []
for name, col in cands:
    if not any(col): continue
    rBC = rank(base + [col]); d = rank(base + [col] + newcols) - rBC
    if d < rBN - rB: hits.append((d, name))
print("candidates that reduce the new dimension:", sorted(hits)[:20] if hits else 'none')
