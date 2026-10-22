# HYP-9364 (b): does kappa(R_n) = H(R_n) 2^(n-1)/n! tend to cosh(1)/cos(1)?  Exact H(R_n) for odd n <= 33
# (ie_ham64 / lattice_bridges_20261021 audit a10_rot.out / R33.out).  Polynomial extrapolation in 1/n on the
# last k points, for kappa and for the ratio kappa/rho with rho(R_n) = ((n+1)^n + (n-1)^n)/(2 Re (n+i)^n)
# (THM-4755 closed form).  mac-mini-2026-10-22-sierp.
from mpmath import mp, mpf, cosh, cos, factorial, matrix, lu_solve, re, mpc
mp.dps = 40
H = {3: 3, 5: 15, 7: 175, 9: 3267, 11: 93027, 13: 3711175, 15: 198464295, 17: 13689269499,
     19: 1184212824763, 21: 125547534942879, 23: 16011537490557279, 25: 2418453569285650675,
     27: 426990256001863046739, 29: 87133267870433699742519, 31: 20351875726012821643802839,
     33: 5394851076280218488551446699, 35: 1610822294345031078388661430315}
import sys
for line in open(sys.argv[1]) if len(sys.argv) > 1 else []:
    if 'H=' in line:
        n = int(line.split('n=')[1].split()[0]); H[n] = int(line.split('H=')[1].split()[0])
L = cosh(1) / cos(1)
def kap(n): return mpf(H[n]) * mpf(2) ** (n - 1) / factorial(n)
def rho(n): return (mpf(n + 1) ** n + mpf(n - 1) ** n) / (2 * re(mpc(n, 1) ** n))
def extrap(ns, f):
    k = len(ns); A = matrix(k, k); b = matrix(k, 1)
    for i, n in enumerate(ns):
        for j in range(k): A[i, j] = mpf(1) / mpf(n) ** j
        b[i] = f(n)
    return lu_solve(A, b)[0]
ns = sorted(H)
print('target cosh(1)/cos(1) =', mp.nstr(L, 15))
print(' n   kappa(R_n)          rho(R_n)            kappa/rho        n(1-kappa/rho)')
for n in ns:
    k, r = kap(n), rho(n)
    print(f"{n:3d}  {mp.nstr(k, 15):18s}  {mp.nstr(r, 15):18s}  {mp.nstr(k / r, 12):15s}  {mp.nstr(n * (1 - k / r), 10)}")
print('\nextrapolation of kappa(R_n) in 1/n, last k points:')
for k in range(2, 9):
    print(f"  k={k}  ns={ns[-k:]}  L_fit = {mp.nstr(extrap(ns[-k:], kap), 12)}   (fit - target = {mp.nstr(extrap(ns[-k:], kap) - L, 4)})")
print('\nextrapolation of kappa/rho in 1/n, last k points (principle predicts 1):')
for k in range(2, 9):
    v = extrap(ns[-k:], lambda n: kap(n) / rho(n))
    print(f"  k={k}  ns={ns[-k:]}  fit = {mp.nstr(v, 12)}")
print('\nextrapolation of y_n = n(1 - kappa/rho), last k points:')
for k in range(2, 8):
    v = extrap(ns[-k:], lambda n: n * (1 - kap(n) / rho(n)))
    print(f"  k={k}  fit y_inf = {mp.nstr(v, 10)}")

# The unstable high-order fits suggest a parity term (-1)^((n-1)/2); fit each residue class mod 4 separately.
for r in (1, 3):
    sub = [n for n in ns if n % 4 == r and n >= 9]
    print(f'\nresidue n = {r} mod 4: {sub}')
    for k in range(2, len(sub) + 1):
        a = extrap(sub[-k:], kap); b = extrap(sub[-k:], lambda n: kap(n) / rho(n))
        print(f"  k={k}  kappa_fit = {mp.nstr(a, 12)} (diff {mp.nstr(a - L, 3)})   kappa/rho_fit = {mp.nstr(b, 10)}")
# second differences of y_n, to see the parity wobble
ys = [(n, n * (1 - kap(n) / rho(n))) for n in ns if n >= 9]
d1 = [(ys[i + 1][0], ys[i + 1][1] - ys[i][1]) for i in range(len(ys) - 1)]
d2 = [(d1[i + 1][0], d1[i + 1][1] - d1[i][1]) for i in range(len(d1) - 1)]
print('\nsecond differences of y_n:', [(n, mp.nstr(v, 4)) for n, v in d2])

# Size of the parity wobble: interpolate the smooth fit of one class (cubic in 1/n through its 4 points nearest
# the target) at the other class's n and compare, relative to kappa.
def interp(pts, f, x):
    from mpmath import mpf
    tot = 0
    for i, a in enumerate(pts):
        w = mpf(1)
        for j, b in enumerate(pts):
            if j != i: w *= (mpf(1) / x - mpf(1) / b) / (mpf(1) / a - mpf(1) / b)
        tot += w * f(a)
    return tot
print('\nparity wobble  w_n = (kappa(n) - smooth fit of the other class at n)/kappa(n):')
for n in ns:
    if n < 13: continue
    other = [m for m in ns if m % 4 != n % 4 and m >= 9]
    near = sorted(other, key=lambda m: abs(m - n))[:4]
    if max(near) < n or min(near) > n: continue      # interpolation only
    w = (kap(n) - interp(sorted(near), kap, n)) / kap(n)
    print(f"  n={n:3d} (d={(n-1)//2} {'even' if (n-1)//2 % 2 == 0 else 'odd '})  w_n = {mp.nstr(w, 4)}   n^3 w_n = {mp.nstr(w * n**3, 4)}   n^4 w_n = {mp.nstr(w * n**4, 4)}")
