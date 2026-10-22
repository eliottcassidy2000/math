# Second order of doubly regular tournaments: P_19 against the non-Paley D_19 (drt_results.out).  Exact forest numbers
# (forest_numbers.c, exact by CRT), the four-loop Delta X_3 (two-matrix reduction), TT5, and the anatomy of
# Delta H by forest length; tests whether Delta N_8 = 7168 Delta TT5 (found at n = 23) holds at n = 19.
# Output: drt19_check.out.  mac-mini-2026-10-22-loop4.
import sys, subprocess, itertools
from fractions import Fraction
from math import factorial, comb
import numpy as np
sys.path.insert(0, '.')
import loop_engine4 as LE4
from drt_switch import code_to_S
from drt23_forest_exact import interpolate_mod, crt, PRIMES

n = 19
codes = {}
for line in open('drt_results.out'):
    t = line.split()
    if len(t) > 6 and t[0] == '19': codes[t[1]] = (int(t[2]), t[6])

def forest(code, H):
    vals = {}
    for p in PRIMES:
        r = subprocess.run(['./forest_numbers', str(p)], input=code + '\n', capture_output=True, text=True, check=True)
        vals[p] = {int(a): int(b) for a, b in (l.split() for l in r.stdout.split('\n') if l.strip())}
    co = [interpolate_mod([vals[p][x] for x in range(n)], p) for p in PRIMES]
    f = [crt([co[j][i] for j in range(3)], PRIMES)[0] for i in range(n)]
    assert f[0] == H and sum(f) == factorial(n) and f == f[::-1]
    N = []
    for L in range(n):
        s = sum(fk * sum(comb(k, a) * comb(n - 1 - k, L - a) * (-1) ** (L - a) for a in range(max(0, L - (n - 1 - k)), min(k, L) + 1)) for k, fk in enumerate(f))
        assert s % factorial(n - L) == 0; N.append(s // factorial(n - L))
    return N

def tt5(S):
    c = 0
    for X in itertools.combinations(range(n), 5):
        sc = sorted(int((S[np.ix_(X, X)][i] == 1).sum()) for i in range(5))
        c += sc == [0, 1, 2, 3, 4]
    return c

def x3_series(S, K=128, R=1.5):
    cores = LE4.load_cores(3)
    G2 = [LE4.drt_gamma2(S, 3, i) for i in range(len(cores))]
    zz = R * np.exp(2j * np.pi * np.arange(K) / K)
    vals = []
    for z in list(zz) + [1.0]:
        C = LE4.drt_coeffs(n, z)
        vals.append(sum(np.sum(LE4.reduce_W(LE4.weight_tensor(3, i), C) * G2[i]) / cores[i][2] for i in range(len(cores))))
    vz = np.array(vals[:K])
    return vals[K].real, [(np.mean(vz * np.exp(-2j * np.pi * np.arange(K) * L / K)) / R ** L).real for L in range(n)]

data = {}
for lab, (H, code) in codes.items():
    S = code_to_S(code); assert LE4.check_drt(S)
    data[lab] = dict(H=H, N=forest(code, H), TT5=tt5(S), X3=x3_series(S))
    print(lab, 'H', H, 'TT5', data[lab]['TT5'], 'N_L (even L <= 10):', data[lab]['N'][0:11:2], 'X_3(1) =', data[lab]['X3'][0], flush=True)
P, D = data['P_19'], data['D_19.1']
dTT5 = D['TT5'] - P['TT5']
dN8 = D['N'][8] - P['N'][8]
print(f"Delta TT5 = {dTT5}; Delta N_L for L = 2, 4, 6: {[D['N'][L] - P['N'][L] for L in (2, 4, 6)]}; Delta N_8 = {dN8} = {Fraction(dN8, dTT5) if dTT5 else 'n/a'} x Delta TT5")
g = [Fraction(n ** L, factorial(n) // factorial(n - L)) for L in range(n)]
phiP = [Fraction(P['N'][L], n ** L) for L in range(n)]; phiD = [Fraction(D['N'][L], n ** L) for L in range(n)]
kapP = sum(a * b for a, b in zip(phiP, g))
assert kapP == Fraction(P['H'] * 2 ** (n - 1), factorial(n))
perL = [float((phiD[L] - phiP[L]) * g[L] / kapP) for L in range(n)]
exact = float(Fraction(D['H'], P['H']) - 1)
delta = np.array(D['X3'][1]) - np.array(P['X3'][1])
d8 = delta[8]
print(f"four-loop delta_8 = {d8:+.6e}, exact Delta phi_8 = {float(phiD[8] - phiP[8]):+.6e} (Delta N_8/n^8)")
strict = D['X3'][0] - P['X3'][0]
D1 = 2 * n / (n + 1)
nxt = sum((1 + (j * D1 + j * (j - 1) / 2) / n) * delta[j] for j in range(n))
print(f"exact H(D_19)/H(P_19) - 1 = {exact:+.4e}; strict four-loop Delta X_3(1) = {strict:+.4e} (ratio exact/strict {exact / strict:.3f}); with the n^-4 finite-size term {nxt:+.4e}")
print("per forest length g(L) Delta phi_L / kappa: " + ' '.join(f"L{L}:{perL[L]:+.3e}" for L in range(8, n, 2)) + f" (sum {sum(perL):+.4e})")
