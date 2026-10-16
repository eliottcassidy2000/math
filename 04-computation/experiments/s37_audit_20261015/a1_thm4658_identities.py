#!/usr/bin/env python3
"""S37 audit (independent referee), part 1: THM-4658 statements 1-5 checked with independent code.

Parity tau(n) mod 2 for T(n) = n/2, (3n+1)/2 is computed by a sequential descent-below-n scan into a byte array
(a different scheme from the session's octave-parallel bit array), validated against Python big-int trajectories.
Checks (all N <= 2^22 unless stated):
  S1  class-sum recursion C(A;N) = -C(A/2; floor(N/2)) - C((3A+1)/2; 3 ceil(N/2)) for every singleton A mod 3^q,
      q = 0..4 (A = Z for q = 0) and random unions mod 27 / 81;
  S2  H(4N) = H(N) + H_(2,2)(3N) - H_(2,8)(6N) and the three one-octave relations; identity (2) for 120 channels;
  S3  E_(4X) = 4E_X + 8C_X + 4D_X + R_X with |R_X| <= 24 sqrt(X E_(4X)) + 48X, the descent inequality and the
      chain bound, for X = 2^k (k <= 20) and random non-dyadic X; exact R_X also via its epsilon decomposition;
  S4b excess identity E~(A;N) = -E~(A/2;floor(N/2)) - E~(A';3ceil(N/2)) + Haar(A) E~(2+3Z;3ceil(N/2)) (integer-scaled);
  S4a the identity H(3m) = -[S(3m) + 3S(2m) + 3S(m)] (exact) which shows that a main term c N^rho of S forces
      H = -3c(2/3)^rho D(rho) N^rho + o(N^sigma): the S-main-term hypothesis alone does NOT give D(rho) = 0;
  S5  finite-level matrices: f -> f(2.) is a permutation (unitary), V^T V = I for V f(y) = 1_(y=2(3)) f((2y-1)/3),
      phi(sigma) = 4^-sigma (1 + 3^(2 sigma - 1)) convex with phi = 1 exactly at 1/2 and 1.
Runtime target < 1 minute.
"""
import time, math, random
import numpy as np
import numba as nb
import sympy as sp

t00 = time.time()
LOGN = 22
NMAXCHK = 1 << LOGN
NT = 6 * NMAXCHK + 64        # need H up to 6N
LS = 3 * (NMAXCHK // 2) + 8  # length needed for floor(N/2) and 3 ceil(N/2)


@nb.njit(cache=False)
def parity_seq(nt):
    p = np.zeros(nt + 1, np.uint8)
    for n in range(2, nt + 1):
        x = n
        s = 0
        while x >= n:
            if x & 1:
                x = (3 * x + 1) >> 1
            else:
                x >>= 1
            s += 1
        p[n] = (s + p[x]) & 1
    return p


def tau_direct(n):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        s += 1
    return s


p = parity_seq(NT)
print(f'parity bytes to {NT}: {time.time() - t00:.1f} s', flush=True)
rng = random.Random(4658)
bad = 0
for _ in range(3000):
    n = rng.randint(1, NT)
    bad += (tau_direct(n) & 1) != p[n]
for n in range(1, 2000):
    bad += (tau_direct(n) & 1) != p[n]
print(f'validation against big-int trajectories (3000 random + n < 2000): {bad} mismatches', flush=True)
assert bad == 0

t = (1 - 2 * p.astype(np.int64))
t[0] = 0


def cum(m, r, Lc):
    z = np.zeros(Lc + 1, np.int64)
    r0 = r % m
    z[r0::m] = t[r0:Lc + 1:m]
    return np.cumsum(z)


def cumset(m, rs, Lc):
    out = np.zeros(Lc + 1, np.int64)
    for r in rs:
        out += cum(m, r, Lc)
    return out


N = np.arange(0, NMAXCHK + 1, dtype=np.int64)
fl = N // 2
ce3 = 3 * ((N + 1) // 2)

# ---------------- S1 ----------------
t1 = time.time()
nchk = 0
fails = 0
for q in range(0, 5):
    m = 3 ** q
    m1 = 3 ** (q + 1)
    inv2 = pow(2, -1, m) if m > 1 else 0
    inv2b = pow(2, -1, m1)
    cm = {r: cum(m, r, LS) for r in range(m)}
    cm1 = {}
    for r in range(m):
        Ah = (r * inv2) % m if m > 1 else 0
        Ap = ((3 * r + 1) * inv2b) % m1
        assert Ap % 3 == 2
        if Ap not in cm1:
            cm1[Ap] = cum(m1, Ap, LS)
        lhs = cm[r][N]
        rhs = -cm[Ah][fl] - cm1[Ap][ce3]
        fails += int(np.count_nonzero(lhs != rhs))
        nchk += len(N)
    del cm, cm1
for trial in range(6):
    q = 3 if trial < 3 else 4
    m = 3 ** q
    m1 = 3 ** (q + 1)
    rs = sorted(rng.sample(range(m), rng.randint(2, m - 1)))
    inv2 = pow(2, -1, m)
    inv2b = pow(2, -1, m1)
    lhs = cumset(m, rs, LS)[N]
    rhs = -cumset(m, [(r * inv2) % m for r in rs], LS)[fl] - cumset(m1, [((3 * r + 1) * inv2b) % m1 for r in rs], LS)[ce3]
    fails += int(np.count_nonzero(lhs != rhs))
    nchk += len(N)
print(f'S1 class-sum recursion: {nchk} (A, N) checks (all 121 singletons mod 3^q, q <= 4, plus 6 random unions), '
      f'{fails} failures  [{time.time() - t1:.1f} s]', flush=True)


def Hq(q, r, Lc):
    m = 3 ** q
    return 3 * cum(m, r, Lc) - (cum(3 ** (q - 1), r, Lc) if q > 1 else cum(1, 0, Lc))


# ---------------- S2 ----------------
H = Hq(1, 2, NT)
H11 = Hq(1, 1, NT)
H10 = Hq(1, 0, NT)
H22 = Hq(2, 2, NT)
H28 = Hq(2, 8, NT)
H25 = Hq(2, 5, NT)
S = cum(1, 0, NT)
f2 = np.count_nonzero(H[4 * N] != H[N] + H22[3 * N] - H28[6 * N])
f2a = np.count_nonzero(H[2 * N] != -H11[N] - H28[3 * N])
f2b = np.count_nonzero(H11[2 * N] != -H[N] - H22[3 * N])
f2c = np.count_nonzero(H10[2 * N] != -H10[N] - H25[3 * N])
print(f'S2 two-octave identity failures {f2}; one-octave relations failures {f2a}, {f2b}, {f2c} '
      f'(N = 0..2^{LOGN})', flush=True)
del H11, H10, H22, H28, H25
t2 = time.time()
fails2 = 0
for q in range(1, 5):
    m = 3 ** q
    inv2 = pow(2, -1, m)
    inv2b = pow(2, -1, 3 * m)
    for r in range(m):
        e = (r * inv2) % m
        o = ((3 * r + 1) * inv2b) % (3 * m)
        fails2 += int(np.count_nonzero(Hq(q, r, LS)[N] != -Hq(q, e, LS)[fl] - Hq(q + 1, o, LS)[ce3]))
print(f'    wavelet identity (2) for all 120 channels q <= 4: {fails2} failures  [{time.time() - t2:.1f} s]', flush=True)

# ---------------- S3 ----------------


def energies(X):
    n = np.arange(X, 2 * X)
    EX = int(np.sum(H[n] ** 2))
    E4 = int(np.sum(H[4 * X:8 * X] ** 2))
    Dl = H[4 * n] - H[n]
    CX = int(np.sum(H[n] * Dl))
    DX = int(np.sum(Dl * Dl))
    R = E4 - 4 * EX - 8 * CX - 4 * DX
    eps = np.stack([H[4 * n + j] - H[4 * n] for j in range(4)], axis=1)
    sig = eps.sum(axis=1)
    R2 = int(2 * np.sum(H[4 * n] * sig) + np.sum(eps * eps))
    return EX, E4, CX, DX, R, R2, int(np.abs(eps).max()), int(np.sum(eps * eps))


print('S3 energy transport, X = 2^k: k, e_X, c_X, delta_X, e_4X, (e+2c+delta)/4, R/(16X^2), |R|/bound, '
      'descent lhs/rhs, max|eps|, sum eps^2/X', flush=True)
s3fail = 0
for k in range(0, 21):
    X = 1 << k
    EX, E4, CX, DX, R, R2, me, se = energies(X)
    assert R == R2
    bound = 24 * math.sqrt(X * E4) + 48 * X
    lhs = math.sqrt(E4)
    rhs = 2 * math.sqrt(EX) + 2 * math.sqrt(DX) + 4 * math.sqrt(3 * X)
    ok = abs(R) <= bound and lhs <= rhs + 1e-9 and me <= 4 and se <= 48 * X
    s3fail += (not ok)
    if k >= 8:
        print(f'  {k:2d} {EX / X**2:.4f} {CX / X**2:+.4f} {DX / X**2:.4f} {E4 / (16 * X * X):.4f} '
              f'{(EX + 2 * CX + DX) / (4 * X * X):.4f} {R / (16 * X * X):+.6f} {abs(R) / bound:.4f} '
              f'{lhs / rhs:.4f} {me} {se / X:.2f}', flush=True)
print(f'  dyadic k <= 20: {s3fail} violations (bound on R, descent, |eps| <= 4, sum eps^2 <= 48X)', flush=True)
nd_fail = 0
Xs = [3, 5, 7, 100, 999, 12345, 54321, 100003, 262143, 300001, 500000, 777777, 1000003] + \
     [rng.randint(2, 1 << 20) for _ in range(20)]
for X in Xs:
    EX, E4, CX, DX, R, R2, me, se = energies(X)
    ok = R == R2 and abs(R) <= 24 * math.sqrt(X * E4) + 48 * X and \
        math.sqrt(E4) <= 2 * math.sqrt(EX) + 2 * math.sqrt(DX) + 4 * math.sqrt(3 * X) + 1e-9
    nd_fail += (not ok)
print(f'  non-dyadic X ({len(Xs)} values up to 2^20): {nd_fail} violations', flush=True)
chain_fail = 0
nchain = 0
for X0 in [1, 2, 3, 5, 6, 7]:
    Xj = [X0 * 4 ** j for j in range(20) if X0 * 4 ** j <= (1 << 20)]
    vals = []
    for X in Xj:
        EX, E4, CX, DX, R, R2, me, se = energies(X)
        vals.append((X, EX / X ** 2, DX / X ** 2))
    for J in range(1, len(vals)):
        lhs = math.sqrt(vals[J][1])
        rhs = 2 ** (-J) * math.sqrt(vals[0][1]) + sum(2 ** (-(J - j)) * (math.sqrt(vals[j][2]) + 2 * math.sqrt(3 / vals[j][0]))
                                                       for j in range(J))
        chain_fail += lhs > rhs + 1e-12
        nchain += 1
print(f'  chain bound sqrt(e_XJ) <= 2^-J sqrt(e_X0) + sum 2^-(J-j)(sqrt(delta_Xj) + 2 sqrt(3/Xj)): {nchain} chains, '
      f'{chain_fail} violations', flush=True)

# ---------------- S4b ----------------
t4 = time.time()
f4 = 0
S_LS = S[:LS + 1]
C2_LS = cum(3, 2, LS)
for q in range(1, 4):
    m = 3 ** q
    m1 = 3 * m
    inv2 = pow(2, -1, m)
    inv2b = pow(2, -1, m1)
    cm = {r: cum(m, r, LS) for r in range(m)}
    for r in range(m):
        Ah = (r * inv2) % m
        Ap = ((3 * r + 1) * inv2b) % m1
        cAp = cum(m1, Ap, LS)
        lhs = 3 * m * cm[r][N] - 3 * S_LS[N]
        rhs = -(3 * m * cm[Ah][fl] - 3 * S_LS[fl]) - (3 * m * cAp[ce3] - S_LS[ce3]) + (3 * C2_LS[ce3] - S_LS[ce3])
        f4 += int(np.count_nonzero(lhs != rhs))
    del cm
print(f'S4b excess identity (integer-scaled by 3^(q+1)), all 39 singletons mod 3^q, q <= 3: {f4} failures  '
      f'[{time.time() - t4:.1f} s]', flush=True)
print(f'    C(2+3Z;N) - S(N)/3 = H(N)/3: {int(np.count_nonzero(3 * C2_LS[N] - S_LS[N] != H[N]))} failures', flush=True)

# ---------------- S4a / Mercator cancellation identity ----------------
mm = np.arange(1, NMAXCHK + 1)
ident = H[3 * mm] + S[3 * mm] + 3 * S[2 * mm] + 3 * S[mm]
print(f'S4a identity H(3m) = -[S(3m) + 3 S(2m) + 3 S(m)] for m <= 2^{LOGN}: {int(np.count_nonzero(ident))} failures', flush=True)
yy = np.arange(1, 3 * NMAXCHK)
bnd = H[yy] + S[yy] + 3 * S[(2 * yy) // 3] + 3 * S[yy // 3]
print(f'    boundary term at general y < 3*2^{LOGN}: max |H(y) + S(y) + 3S(floor(2y/3)) + 3S(floor(y/3))| = {int(np.abs(bnd).max())}',
      flush=True)
rho = sp.symbols('rho')
Dsym = 1 + 2 ** (-rho) + sp.Rational(3, 2) ** rho / 3
expr = sp.simplify(sp.expand(-(1 + 3 * sp.Rational(2, 3) ** rho + 3 * sp.Rational(1, 3) ** rho) + 3 * sp.Rational(2, 3) ** rho * Dsym, power_base=True))
num = [complex(sp.N((-(1 + 3 * sp.Rational(2, 3) ** rho + 3 * sp.Rational(1, 3) ** rho) + 3 * sp.Rational(2, 3) ** rho * Dsym).subs(rho, v)))
       for v in [0.3, 0.5 + 2j, 1.7 - 0.4j]]
print(f'    check -(1 + 3(2/3)^rho + 3(1/3)^rho) + 3(2/3)^rho D(rho) == 0 at sample points: {[abs(z) < 1e-12 for z in num]}  '
      f'(so a main term c N^rho of S gives H = -3c(2/3)^rho D(rho) N^rho + o(N^sigma))', flush=True)
print(f'    D(sigma) for real sigma > 0 has no zeros: D(0.5) = {float(Dsym.subs(rho, 0.5)):.6f}, '
      f'D(1) = {float(Dsym.subs(rho, 1)):.6f}', flush=True)

# ---------------- S5 ----------------
s5 = []
for Q in range(1, 7):
    M = 3 ** Q
    Mp = 3 ** (Q - 1)
    perm = (2 * np.arange(M)) % M
    U_ok = len(set(perm.tolist())) == M
    V = np.zeros((M, Mp), np.int64)
    for y in range(M):
        if y % 3 == 2:
            V[y, ((2 * y - 1) // 3) % Mp] = 1
    s5.append((Q, U_ok, bool(np.array_equal(V.T @ V, np.eye(Mp, dtype=np.int64)))))
print(f'S5 finite levels (Q, f(2.) permutation, V^T V = I  <=> <Vf,Vg> = <f,g>/3 for normalised Haar): {s5}', flush=True)
sig = sp.symbols('sigma', real=True)
phi = 4 ** (-sig) * (1 + 3 ** (2 * sig - 1))
print(f'    phi(1/2) = {sp.nsimplify(phi.subs(sig, sp.Rational(1, 2)))}, phi(1) = {sp.nsimplify(phi.subs(sig, 1))}, '
      f'phi(3/4) = {float(phi.subs(sig, 0.75)):.6f}, phi(0) = {float(phi.subs(sig, 0)):.6f}, phi(2) = {float(phi.subs(sig, 2)):.6f}',
      flush=True)
fphi = sp.lambdify(sig, phi - 1, 'mpmath')
grid_s = np.linspace(-3, 4, 7001)
vals_s = [float(fphi(float(x))) for x in grid_s]
roots = []
import mpmath
for i in range(len(grid_s) - 1):
    if vals_s[i] == 0 or vals_s[i] * vals_s[i + 1] < 0:
        roots.append(float(mpmath.findroot(fphi, (float(grid_s[i]), float(grid_s[i + 1])), solver='bisect')))
neg_between = all(v < 0 for x, v in zip(grid_s, vals_s) if 0.5 + 1e-9 < x < 1 - 1e-9)
pos_outside = all(v > 0 for x, v in zip(grid_s, vals_s) if x < 0.5 - 1e-9 or x > 1 + 1e-9)
print(f'    roots of phi = 1 on [-3, 4]: {[round(r, 12) for r in roots]}; phi < 1 on (1/2,1): {neg_between}; '
      f'phi > 1 outside: {pos_outside}; phi\'\' = {sp.simplify(sp.diff(phi, sig, 2))} > 0', flush=True)
same = sp.simplify(phi - (2 ** (-2 * sig) + sp.Rational(3, 2) ** (2 * sig) / 3)) == 0
print(f'    phi(sigma) = f(2 sigma) with f(s) = 2^-s + (3/2)^s/3 (total-variation branch norms, = 1 at s = 1, 2): {same}', flush=True)
print(f'total {time.time() - t00:.1f} s')
