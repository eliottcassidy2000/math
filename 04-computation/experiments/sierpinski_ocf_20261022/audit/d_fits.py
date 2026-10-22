# Audit of HYP-9364 (a)-(d): normalised counts kappa = H 2^(n-1)/n! and polynomial extrapolations in 1/n.
# H values: n <= 27 recomputed by the auditor (b_ie_circulant, d_ie_small.out); n = 29, 31 (R) and 31 (P) from
# repo files, re-checked by d_ie_large.out when it finishes; R_33 from the session's R33.out (not used by HYP-9364).
from fractions import Fraction
from math import factorial
from mpmath import mp, mpf, matrix, lu_solve, cosh, cos, e
mp.dps = 40
HR = {3: 3, 5: 15, 7: 175, 9: 3267, 11: 93027, 13: 3711175, 15: 198464295, 17: 13689269499,
      19: 1184212824763, 21: 125547534942879, 23: 16011537490557279, 25: 2418453569285650675,
      27: 426990256001863046739, 29: 87133267870433699742519, 31: 20351875726012821643802839,
      33: 5394851076280218488551446699}
HP = {3: 3, 7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349, 27: 418425945688463432031,
      31: 19876729650849384218015865}
kap = lambda n, h: mpf(Fraction(h * 2 ** (n - 1), factorial(n)).numerator) / mpf(Fraction(h * 2 ** (n - 1), factorial(n)).denominator)
KR = {n: kap(n, h) for n, h in HR.items()}
KP = {n: kap(n, h) for n, h in HP.items()}
L_rot = cosh(1) / cos(1)
def fit0(pts):
    """polynomial through (1/n, y) of degree len-1, value at 1/n = 0"""
    k = len(pts)
    A = matrix(k, k); b = matrix(k, 1)
    for i, (n, y) in enumerate(pts):
        for j in range(k): A[i, j] = mpf(1) / mpf(n) ** j
        b[i] = y
    return lu_solve(A, b)[0]
print("kappa(R_n):", ", ".join(f"{n}:{mp.nstr(KR[n], 7)}" for n in sorted(KR)))
print("kappa(P_q):", ", ".join(f"{n}:{mp.nstr(KP[n], 7)}" for n in sorted(KP)))
print("targets: cosh1/cos1 =", mp.nstr(L_rot, 10), " e =", mp.nstr(e, 10), " ratio =", mp.nstr(L_rot / e, 10))
print("\n(b) rotational fits in 1/n (value at 0):")
for ns in [(27, 29, 31), (25, 27, 29, 31), (23, 25, 27, 29, 31), (21, 23, 25, 27, 29, 31), (29, 31, 33), (27, 29, 31, 33), (25, 27, 29, 31, 33), (23, 25, 27, 29, 31, 33), (21,23,25,27,29,31,33), (25, 27), (27, 29), (29, 31), (31, 33)]:
    print(f"  n={ns}: {mp.nstr(fit0([(n, KR[n]) for n in ns]), 8)}")
print("\n(b') successive increments of kappa(R_n) and n^2 * increment (a 1/n law gives n^2*diff/2 ~ const):")
ns = sorted(KR)
for a, b2 in zip(ns, ns[1:]):
    d = KR[b2] - KR[a]
    print(f"  {a}->{b2}: diff {mp.nstr(d, 6)}, n_mid^2*diff/2 = {mp.nstr(d * ((a + b2) / 2) ** 2 / 2, 6)}")
print("\n(c) Paley fits in 1/q:")
for qs in [(19, 23, 31), (19, 23, 27, 31), (23, 27, 31), (11, 19, 23, 27, 31), (11, 19, 23, 31)]:
    print(f"  q={qs}: {mp.nstr(fit0([(q, KP[q]) for q in qs]), 8)}")
print("  1 - kappa(P_q)/e times q:", ", ".join(f"{q}:{mp.nstr((1 - KP[q] / e) * q, 6)}" for q in sorted(KP)))
print("\n(d) ratio H(R_q)/H(P_q):", ", ".join(f"{q}:{mp.nstr(mpf(HR[q]) / HP[q], 8)}" for q in sorted(HP) if q in HR))
for qs in [(19, 23, 31), (19, 23, 27, 31), (23, 27, 31), (11, 19, 23, 31), (11, 19, 23, 27, 31)]:
    print(f"  q={qs}: {mp.nstr(fit0([(q, mpf(HR[q]) / HP[q]) for q in qs]), 8)}")
