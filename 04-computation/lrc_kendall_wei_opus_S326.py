# opus-2026-07-16-S326 -- HYP-7030 / THM-894: LRC'S KENDALL-WEI.
# The resonance matrix M[i][j] = rho(x_i, x_j) (the exact sawtooth pair law,
# diag 2/13 = single-comb density) over the cluster's OWN speeds; the Perron
# pair (lambda, w) = the self-referential resonance classifier.
# Battery: named families. Questions:
#   (Q1) does w-uniformity/concentration classify species (AP uniform-ish,
#        wells concentrated on binding speeds)?
#   (Q2) is lambda spectral beyond the first moment (compare with row-sum
#        mean = saw-adjacent)?
#   (Q3) which speeds carry the Perron mass in the deep families -- the
#        binding speeds? (the classifier's self-referential selection)
from fractions import Fraction
from math import gcd
import numpy as np, itertools

def T_of(a, b):
    s = a + b
    T, m = 0, 0
    while s - 13*m > 0:
        T += (min(2*a, s - 13*m)) * (1 if m == 0 else 2)
        m += 1
    return T

def rho(x1, x2):
    if x1 == x2: return Fraction(2, 13)
    g = gcd(x1, x2)
    a, b = x1//g, x2//g
    if a > b: a, b = b, a
    return Fraction(T_of(a, b), 13*a*b)

FAMS = {
    'AP {1..13} (tight)': list(range(1, 14)),
    '2AP-1 (loose twin-E)': list(range(1, 26, 2)),
    'GAP {1..7}u{10..16}': list(range(1, 8)) + list(range(10, 17)),
    'near-AP {21..32}': list(range(21, 33)),
    'floor {1..13}-{6}': [v for v in range(1, 14) if v != 6],
    'mediant S_13 {1..11,13,36}': list(range(1, 12)) + [13, 36],
    'row (13,9) [2U0+x,y]': sorted([2*u for u in [1,2,3,5,7,8,9,10,11,12]] + [13, 9]),
    'row (17,13) [2U0+x,y]': sorted([2*u for u in [1,2,3,5,7,8,9,10,11,12]] + [17, 13]),
}
print(f"{'family':28s} {'lam':8s} {'rowmean':8s} {'lam/rm':7s} {'w_max/w_min':11s} "
      f"{'top-2 Perron speeds':22s}")
for name, S in FAMS.items():
    n = len(S)
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            M[i][j] = float(rho(S[i], S[j]))
    ev, V = np.linalg.eigh(M)
    lam = ev[-1]
    w = np.abs(V[:, -1]); w = w / w.sum()
    rowmean = M.sum() / n   # average row sum = first-moment proxy
    conc = w.max() / w.min()
    top = sorted(zip(w, S), reverse=True)[:2]
    print(f"{name:28s} {lam:8.4f} {rowmean:8.4f} {lam/rowmean:7.4f} {conc:11.2f} "
          f"{str([(sp, round(float(ww),3)) for ww, sp in top]):22s}")

print()
print("Q1/Q3 detail -- the Perron profile of the floor family (binding structure):")
S = [v for v in range(1, 14) if v != 6]
n = len(S)
M = np.zeros((n, n))
for i in range(n):
    for j in range(n): M[i][j] = float(rho(S[i], S[j]))
ev, V = np.linalg.eigh(M)
w = np.abs(V[:, -1]); w = w / w.sum()
for sp, ww in zip(S, w): print(f"   speed {sp:3d}: w = {ww:.4f}")
print()
print("self-referential stability check (power iteration = the recursion into")
print("its own second moment): ||M w - lam w|| for the AP:")
S = list(range(1, 14)); n = 13
M = np.zeros((n, n))
for i in range(n):
    for j in range(n): M[i][j] = float(rho(S[i], S[j]))
ev, V = np.linalg.eigh(M)
lam, w = ev[-1], V[:, -1]
print(f"   residual = {np.linalg.norm(M @ w - lam * w):.2e} (fixed point exact)")
print(f"   AP Perron profile (should be near-uniform if tight = maximally")
print(f"   mutually resonant): min = {np.abs(w).min()/np.abs(w).sum():.4f}, "
      f"max = {np.abs(w).max()/np.abs(w).sum():.4f} (uniform = {1/13:.4f})")
