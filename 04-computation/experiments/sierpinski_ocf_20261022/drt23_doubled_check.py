# Out-of-sample test at n = 23: the doubly regular tournament obtained by normalising the doubled Paley skew-Hadamard
# matrix of order 24 (drt_c1_test.py) has TT5 = 1430, far outside the annealed classes (1466..1494).  Exact H (hk64),
# exact forest numbers (forest_numbers.c, CRT), and comparison with
#   - the empirical TT5 law of the 12 annealed non-Paley classes (HYP-9367: H = 1.575375e16 + 4.199788e9 TT5),
#   - the four-loop strict and dressed predictions (HYP-9371),
#   - the exact anatomy by forest length.
# Output: drt23_doubled_check.out.  mac-mini-2026-10-22-loop4.
import sys, json, subprocess, itertools
from fractions import Fraction
from math import factorial, comb
import numpy as np
sys.path.insert(0, '.')
import loop_engine4 as LE4
from drt_switch import d6, S_to_code
from drt_c1_test import bordered, double, normalisations, paley_prime, census5
from drt23_forest_exact import interpolate_mod, crt, PRIMES

n = 23
cls = json.load(open('drt23_classes.json'))
ex = {r['class']: r for r in json.load(open('drt23_forest_exact.json'))}
P = [c for c in cls if c['label'] == 'P_23'][0]
SP = paley_prime(23)
fam = normalisations(double(bordered(paley_prime(11))))
canon = subprocess.run(['labelg', '-q'], input='\n'.join(d6(D) for D in [SP] + fam) + '\n', capture_output=True, text=True).stdout.split()
known = {c['canon'] for c in cls}
D = [M for M, cf in zip(fam, canon[1:]) if cf != canon[0]][0]
cfD = [cf for cf in canon[1:] if cf != canon[0]][0]
code = S_to_code(D)
print('doubled-Paley DRT(23): new nauty class (not among the 13 annealed/Paley classes):', cfD not in known)
H = int(subprocess.run(['./hk64'], input=code + '\n', capture_output=True, text=True).stdout.split()[0])
TT5 = census5(D)['(0, 1, 2, 3, 4)|1']
print(f"H = {H}, TT5 = {TT5}, H/H(P_23) - 1 = {float(Fraction(H, P['H']) - 1):+.6e}")

# empirical TT5 law from the 12 annealed classes, refitted from drt23_classes.json
nonP = [c for c in cls if c['label'] != 'P_23']
A = np.array([[1.0, c['TT5']] for c in nonP]); y = np.array([float(c['H']) for c in nonP])
(a, b), *_ = np.linalg.lstsq(A, y, rcond=None)
pred = a + b * TT5
print(f"TT5 law of the annealed classes (H = {a:.6e} + {b:.6e} TT5) predicts H = {pred:.6e}: relative error {pred / H - 1:+.3e} (in-sample max residual 3.7e-7)")

# exact forest numbers
vals = {}
for p in PRIMES:
    r = subprocess.run(['./forest_numbers', str(p)], input=code + '\n', capture_output=True, text=True, check=True)
    vals[p] = {int(u): int(v) for u, v in (l.split() for l in r.stdout.split('\n') if l.strip())}
co = [interpolate_mod([vals[p][x] for x in range(n)], p) for p in PRIMES]
f = [crt([co[j][i] for j in range(3)], PRIMES)[0] for i in range(n)]
assert f[0] == H and sum(f) == factorial(n) and f == f[::-1]
N = []
for L in range(n):
    s = sum(fk * sum(comb(k, a_) * comb(n - 1 - k, L - a_) * (-1) ** (L - a_) for a_ in range(max(0, L - (n - 1 - k)), min(k, L) + 1)) for k, fk in enumerate(f))
    assert s % factorial(n - L) == 0; N.append(s // factorial(n - L))
NP = [int(x) for x in ex[P['class']]['N']]
dTT5 = TT5 - P['TT5']
print(f"Delta N_L for L = 2, 4, 6: {[N[L] - NP[L] for L in (2, 4, 6)]}; Delta N_8 = {N[8] - NP[8]} = {Fraction(N[8] - NP[8], dTT5)} x Delta TT5")
g = [Fraction(n ** L, factorial(n) // factorial(n - L)) for L in range(n)]
kapP = sum(Fraction(NP[L], n ** L) * g[L] for L in range(n))
perL = [float(Fraction(N[L] - NP[L], n ** L) * g[L] / kapP) for L in range(n)]
exact = float(Fraction(H, P['H']) - 1)
print("per forest length g(L) Delta phi_L/kappa: " + ' '.join(f"L{L}:{perL[L]:+.3e}" for L in range(8, n, 2)) + f"  (sum {sum(perL):+.4e})")
print("   fractions: " + ' '.join(f"L{L}:{perL[L] / exact:.3f}" for L in range(8, 16, 2)))
# four-loop predictions: Delta X_3 is 1.0061e-7 Delta TT5 (strict) and 2.9972e-7 Delta TT5 (dressed) per HYP-9371
fl = json.load(open('drt23_four_loop.json'))
rows = [r for r in fl['rows'] if r['label'] != 'P_23']
ss = np.polyfit([r['TT5'] for r in rows], [r['strict'] for r in rows], 1)[0]
sd = np.polyfit([r['TT5'] for r in rows], [r['dressed'] for r in rows], 1)[0]
print(f"four-loop strict {ss * dTT5:+.4e} (exact/strict {exact / (ss * dTT5):.3f}), dressed {sd * dTT5:+.4e} (dressed/exact {sd * dTT5 / exact:.3f}); exact {exact:+.4e}")
print(f"exact relative to Paley per TT5: {exact / dTT5:.4e} (annealed-class slope 2.6648e-07, line through Paley would need {exact / dTT5:.4e})")
json.dump({'code': code, 'canon': cfD, 'H': str(H), 'TT5': TT5, 'f': [str(x) for x in f], 'N': [str(x) for x in N]},
          open('drt23_doubled.json', 'w'), indent=1)
