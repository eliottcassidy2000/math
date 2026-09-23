#!/usr/bin/env python3
"""
collatz_procgen_20260923_trunk_rh_dirichlet.py

Archimedean side of the trunk/RH lane (session collatz-procgen-20260922, 2026-09-23).

Sections:
  D1  the trunk Dirichlet series D_T(s) = sum_{i>=1} T_i^(-s), T_i = (4^i-1)/3:
      rigorous zero-free half-plane Re s > sigma_0 (sigma_0 < 1/2), meromorphic
      continuation, poles on the lines Re s = -j, and the zeros found in a box
  D2  the complex trunk T(z) = (4^z-1)/3 vanishes at z = i*pi*k/log 2; under
      s = 1-2z these are the zeros of 1 - 2^(1-s) on Re s = 1 (eta's extra zeros)
  D3  Davenport-Heilbronn control: functional equation verified, zeros OFF the
      critical line located and certified by the argument principle
  D4  Epstein control Z(s) = sum' (m^2+5n^2)^(-s): decomposition and functional
      equation verified; zeros off the line searched in a box
  D5  the backward E-game's own dynamical zeta function 1/(1 - W(3^-s; theta))
      of the canonical escape Psi (the only zeta function in which the hostile
      point 1/2 enters, via its prices): dimension check and resonance locations

Numerics: numpy float64 with an Euler-Maclaurin Hurwitz zeta (validated against
mpmath), mpmath for validation and polishing.  Runtime a few minutes; memory
< 300 MB.
"""
import math
import time
import cmath
import numpy as np
import mpmath as mp

T0 = time.time()
LOG2, LOG3, LOG4 = math.log(2), math.log(3), math.log(4)


def hdr(t):
    print()
    print('=' * 78)
    print(t)
    print('=' * 78)


# ----------------------------------------------------------------------------
# numpy Hurwitz zeta by Euler-Maclaurin (vectorized over s)
# ----------------------------------------------------------------------------
_B2K = [1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66, -691 / 2730, 7 / 6, -3617 / 510,
        43867 / 798, -174611 / 330, 854513 / 138, -236364091 / 2730]


def hurwitz(s, a, N=None, K=12):
    """zeta(s, a) for complex array s, 0 < a <= 1, via Euler-Maclaurin."""
    s = np.asarray(s, dtype=complex)
    if N is None:
        N = int(max(60, np.max(np.abs(s.imag)) / 2 + 60))
    n = np.arange(N) + a
    logn = np.log(n)
    tot = np.zeros(s.shape, dtype=complex)
    # sum_{n<N} (n+a)^(-s), chunked
    flat = s.ravel()
    out = np.zeros(flat.shape, dtype=complex)
    CH = int(max(1, min(256, 1.5e6 // max(1, flat.size))))
    for i0 in range(0, N, CH):
        ln = logn[i0:i0 + CH]
        out += np.exp(-np.outer(flat, ln)).sum(axis=1)
    tot = out.reshape(s.shape)
    x = N + a
    lx = math.log(x)
    tot = tot + np.exp((1 - s) * lx) / (s - 1) + 0.5 * np.exp(-s * lx)
    # Bernoulli correction terms
    poch = s.copy()                  # s (s+1) ... (s+2k-2)
    fact = 2.0                       # (2k)!
    for k in range(1, K + 1):
        term = _B2K[k - 1] / fact * poch * np.exp((-s - 2 * k + 1) * lx)
        tot = tot + term
        poch = poch * (s + 2 * k - 1) * (s + 2 * k)
        fact = fact * (2 * k + 1) * (2 * k + 2)
    return tot


def dirichlet_L(s, chi, q):
    """L(s, chi) for a character given as a list chi[0..q-1]"""
    s = np.asarray(s, dtype=complex)
    tot = np.zeros(s.shape, dtype=complex)
    for a in range(1, q + 1):
        c = chi[a % q]
        if c != 0:
            tot = tot + c * hurwitz(s, a / q)
    return tot * np.exp(-s * math.log(q))


def winding(vals):
    """total change of argument / 2pi along a closed polyline of values"""
    ph = np.angle(vals)
    d = np.diff(np.concatenate([ph, ph[:1]]))
    d = (d + np.pi) % (2 * np.pi) - np.pi
    return d.sum() / (2 * np.pi), np.max(np.abs(d))


def rect_count(F, s0, s1, t0, t1, dt=0.02, ds=0.01):
    """argument-principle zero count of F in [s0,s1] x [t0,t1]"""
    bottom = s0 + (np.arange(int((s1 - s0) / ds) + 1) * ds) + 1j * t0
    right = s1 + 1j * (t0 + np.arange(int((t1 - t0) / dt) + 1) * dt)
    top = (s1 - np.arange(int((s1 - s0) / ds) + 1) * ds) + 1j * t1
    left = s0 + 1j * (t1 - np.arange(int((t1 - t0) / dt) + 1) * dt)
    path = np.concatenate([bottom, right, top, left])
    w, mx = winding(F(path))
    return w, mx


# ----------------------------------------------------------------------------
hdr('D0  validation of the Euler-Maclaurin Hurwitz zeta against mpmath')
mp.mp.dps = 25
errs = []
for s in [complex(0.5, 14.1), complex(0.8, 85.7), complex(0.6, 180.3), complex(1.3, 250.0),
          complex(-0.4, 30.0), complex(2.0, 0.5)]:
    for a in (0.2, 0.4, 0.6, 0.8, 1.0, 0.05):
        ref = complex(mp.zeta(mp.mpc(s.real, s.imag), mp.mpf(a)))
        val = complex(hurwitz(np.array([s]), a)[0])
        errs.append(abs(val - ref) / max(1, abs(ref)))
print(f'max relative error over 36 test points (t up to 250): {max(errs):.2e}')

# ----------------------------------------------------------------------------
hdr('D1  the trunk Dirichlet series D_T(s) = sum_{i>=1} T_i^(-s)')
mp.mp.dps = 40
Ti = [mp.mpf(4 ** i - 1) / 3 for i in range(1, 400)]


def g_tail(sig, legal=False):
    """sum_{i>=2} T_i^-sig (legal: only 3 does not divide i), with rigorous tail bound"""
    tot = mp.mpf(0)
    for i in range(2, 300):
        if legal and i % 3 == 0:
            continue
        tot += Ti[i - 1] ** (-sig)
    # tail i >= 300: T_i >= (5/16) 4^i, so sum <= (16/5)^sig 4^(-300 sig)/(1-4^-sig)
    tail = (mp.mpf(16) / 5) ** sig * mp.mpf(4) ** (-300 * sig) / (1 - mp.mpf(4) ** (-sig))
    return tot, tail


def solve_sigma0(legal=False):
    lo, hi = mp.mpf('0.05'), mp.mpf(2)
    for _ in range(80):
        mid = (lo + hi) / 2
        tot, tail = g_tail(mid, legal)
        if tot > 1:
            lo = mid
        else:
            hi = mid
    tot, tail = g_tail(hi, legal)
    return hi, tot, tail


s0, tot, tail = solve_sigma0()
s0L, totL, tailL = solve_sigma0(True)
print(f'sigma_0 (all trunk numbers)       = {mp.nstr(s0, 12)}   (sum_{{i>=2}} T_i^-sigma_0 = 1; tail < {mp.nstr(tail, 3)})')
print(f'sigma_0 (legal exits, 3 does not divide i) = {mp.nstr(s0L, 12)}')
print('PROVED: for Re s > sigma_0,  |D_T(s)| >= 1 - sum_{i>=2} T_i^(-Re s) > 0.  Since sigma_0 < 1/2,')
print('        D_T has NO zeros on the critical line Re s = 1/2 (nor anywhere in Re s > sigma_0).')
for sig in (0.5, 0.6, 1.0):
    t, _ = g_tail(mp.mpf(sig))
    print(f'   margin at Re s = {sig}: 1 - sum_(i>=2) T_i^-sigma = {mp.nstr(1 - t, 8)}')


def DT_cont(s, J=80):
    """D_T(s) = 3^s sum_{j>=0} (s)_j/j! / (4^(s+j) - 1)  (meromorphic continuation)"""
    s = np.asarray(s, dtype=complex)
    tot = np.zeros(s.shape, dtype=complex)
    coef = np.ones(s.shape, dtype=complex)
    for j in range(J):
        tot = tot + coef / (np.exp((s + j) * LOG4) - 1)
        coef = coef * (s + j) / (j + 1)
    return np.exp(s * LOG3) * tot


def DT_direct(s, I=200):
    return sum(complex(Ti[i - 1] ** (-mp.mpc(s.real, s.imag))) for i in range(1, I))


for s in [complex(2, 0), complex(1, 3), complex(0.7, 10), complex(0.3, 25)]:
    a, b = DT_direct(s), complex(DT_cont(np.array([s]))[0])
    print(f'   continuation check at s = {s}: direct {a:.12f}  formula {b:.12f}  |diff| = {abs(a-b):.1e}')
k1 = complex(0, math.pi / LOG2)
eps = 1e-7
res = complex(DT_cont(np.array([k1 + eps]))[0]) * eps
print(f'   pole at s = i*pi/log 2 = {k1:.6f}: residue ~ {res:.8f}; predicted 3^s/log 4 = '
      f'{cmath.exp(k1 * LOG3) / LOG4:.8f}')
# zeros in the box [-3.9, 1.2] x [0.15, 30] via cell windings on a grid
ds, dt = 0.02, 0.02
sig = np.arange(-3.9, 1.2 + 1e-9, ds)
tt = np.arange(0.15, 30 + 1e-9, dt)
S, TT = np.meshgrid(sig, tt)
V = DT_cont(S + 1j * TT)
ph = np.angle(V)


def wrap(d):
    return (d + np.pi) % (2 * np.pi) - np.pi


w = (wrap(ph[:-1, 1:] - ph[:-1, :-1]) + wrap(ph[1:, 1:] - ph[:-1, 1:]) +
     wrap(ph[1:, :-1] - ph[1:, 1:]) + wrap(ph[:-1, :-1] - ph[1:, :-1])) / (2 * np.pi)
wi = np.rint(w).astype(int)
zc = np.argwhere(wi == 1)
pc = np.argwhere(wi == -1)
roots = []
for (r, c) in zc:
    z0 = mp.mpc(sig[c] + ds / 2, tt[r] + dt / 2)
    mp.mp.dps = 30

    def f(z):
        tot = mp.mpf(0)
        coef = mp.mpf(1)
        for j in range(120):
            tot += coef / (mp.power(4, z + j) - 1)
            coef = coef * (z + j) / (j + 1)
        return mp.power(3, z) * tot
    try:
        z = mp.findroot(f, z0, tol=1e-25, maxsteps=50)
        if abs(f(z)) < 1e-20:
            roots.append(complex(z))
    except Exception:
        pass
roots = sorted(set((round(z.real, 10), round(z.imag, 10)) for z in roots), key=lambda x: x[1])
print(f'   box [-3.9,1.2] x [0.15,30]: {len(zc)} zero cells, {len(pc)} pole cells '
      f'(poles expected at Re s = 0,-1,-2,-3, Im s = pi k/log 2: {4 * int(30 / (math.pi / LOG2))})')
print(f'   {len(roots)} zeros polished (|D_T| < 1e-20 at 30 digits); real parts range '
      f'[{min(r[0] for r in roots):.4f}, {max(r[0] for r in roots):.4f}]' if roots else '   no zeros')
for j in range(0, len(roots), 3):
    print('      ' + '   '.join(f'{r[0]:+.8f}{r[1]:+.8f}i' for r in roots[j:j + 3]))
near_half = [r for r in roots if abs(r[0] - 0.5) < 0.05]
print(f'   zeros with |Re s - 1/2| < 0.05: {len(near_half)}  (the proof above forbids any with Re s > {mp.nstr(s0, 4)})')

# ----------------------------------------------------------------------------
hdr('D2  complex zeros of the trunk = zeros of the Euler-type factor 1 - 2^(1-s) on Re s = 1')
mp.mp.dps = 30
for k in range(1, 6):
    z = mp.mpc(0, mp.pi * k / mp.log(2))
    Tz = (mp.power(4, z) - 1) / 3
    s = 1 - 2 * z
    print(f'   k={k}: T(i pi k/log2) = {mp.nstr(abs(Tz), 3)};  s = 1 - 2z = {mp.nstr(s, 12)};  '
          f'|eta(s)| = {mp.nstr(abs(mp.altzeta(s)), 3)},  |zeta(s)| = {mp.nstr(abs(mp.zeta(s)), 8)}')
print('   => the trunk vanishes exactly where eta(s) = (1-2^(1-s)) zeta(s) has its non-zeta zeros, on Re s = 1.')
print('      Any vertical line can be moved to Re s = 1/2 by an affine change of variable; this is not content.')

# ----------------------------------------------------------------------------
hdr('D3  control: the Davenport-Heilbronn function (functional equation, no Euler product)')
q = 5
chi = [0, 1, 1j, -1j, -1]          # chi(2) = i, chi(3) = -i, chi(4) = -1: odd character mod 5
chib = [np.conj(c) for c in chi]
tau = sum(chi[a] * cmath.exp(2j * math.pi * a / q) for a in range(1, q))
W = tau / (1j * math.sqrt(q))
kappa_num = math.tan(cmath.phase(W) / 2)
kappa_closed = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
print(f'   root number W(chi) = tau/(i sqrt5) = {W:.12f} (|W| = {abs(W):.12f})')
print(f'   kappa from (1+i kappa)/(1-i kappa) = W: {kappa_num:.15f};  closed form '
      f'(sqrt(10-2sqrt5)-2)/(sqrt5-1) = {kappa_closed:.15f}')
kap = kappa_num
ca = [((1 - 1j * kap) / 2) * chi[a] + ((1 + 1j * kap) / 2) * chib[a] for a in range(q)]


def DH(s):
    return dirichlet_L(s, ca, q)


def DH_mp(z):
    return mp.power(5, -z) * sum(mp.mpc(ca[a].real, ca[a].imag) * mp.zeta(z, mp.mpf(a) / 5)
                                 for a in range(1, 5))


def LamDH(z):
    return (mp.mpf(5) / mp.pi) ** ((z + 1) / 2) * mp.gamma((z + 1) / 2) * DH_mp(z)


mp.mp.dps = 25
for z in [mp.mpc(0.3, 7), mp.mpc(0.9, 40), mp.mpc(2, 1)]:
    a, b = LamDH(z), LamDH(1 - z)
    print(f'   functional equation Lambda(s) = Lambda(1-s) at s = {mp.nstr(z, 4)}: '
          f'|diff|/|Lambda| = {mp.nstr(abs(a - b) / abs(a), 3)}')
print('   coefficients: a_n = Re-part combination, a_1 = 1, not multiplicative:',
      f'a_2 = {ca[2]:.4f}, a_3 = {ca[3]:.4f}, a_4 = {ca[4]:.4f}, a_2*a_2 = {ca[2]*ca[2]:.4f}')
wcount, mx = rect_count(DH, 0.52, 1.6, 1.0, 200.0, dt=0.01, ds=0.005)
print(f'   zeros of DH in [0.52, 1.6] x [1, 200] (argument principle): {wcount:.4f} '
      f'(max phase step {mx:.3f} rad)')
# localize: sub-rectangles of height 10
found = []
for t0 in range(1, 200, 10):
    wc, _ = rect_count(DH, 0.52, 1.6, float(t0), float(min(t0 + 10, 200)), dt=0.01, ds=0.005)
    n = int(round(wc))
    if n > 0:
        # grid minima of |DH| inside, then polish
        sg = np.arange(0.53, 1.6, 0.01)
        tg = np.arange(t0, min(t0 + 10, 200), 0.01)
        SS, TT2 = np.meshgrid(sg, tg)
        A = np.abs(DH(SS + 1j * TT2))
        idx = np.argsort(A.ravel())[:60]
        cand = []
        for ii in idx:
            z0 = complex(SS.ravel()[ii], TT2.ravel()[ii])
            try:
                z = mp.findroot(DH_mp, mp.mpc(z0.real, z0.imag), tol=1e-18)
            except Exception:
                continue
            z = complex(z)
            if 0.52 < z.real < 1.6 and t0 <= z.imag <= t0 + 10 and \
                    all(abs(z - c) > 1e-6 for c in cand):
                cand.append(z)
            if len(cand) >= n:
                break
        found += cand
        print(f'     t in [{t0},{t0+10}]: {n} zero(s):',
              ', '.join(f'{z.real:.6f}{z.imag:+.6f}i' for z in cand))
print(f'   => {len(found)} zeros with 0.52 < Re s < 1.6 and t <= 200; by the functional equation each has a')
print('      mirror at 1 - conj(rho): the DH function satisfies the zeta-type functional equation and')
print('      VIOLATES the Riemann hypothesis.  (Zeros ON the line are also plentiful; not counted here.)')

# ----------------------------------------------------------------------------
hdr('D4  control: the Epstein zeta function of x^2 + 5y^2 (class number 2)')
chi_m4 = [0, 1, 0, -1]
chi_5 = [0, 1, -1, -1, 1]
chi_m20 = [0] * 20
for n in range(20):
    chi_m20[n] = chi_m4[n % 4] * chi_5[n % 5]
chi_1 = [1]


def Ep(s):
    return (hurwitz(s, 1.0) * dirichlet_L(s, chi_m20, 20) +
            dirichlet_L(s, chi_m4, 4) * dirichlet_L(s, chi_5, 5))


# direct lattice sum at s = 2 (with an integral tail estimate)
Rm = 3000
m = np.arange(-Rm, Rm + 1)
tot = 0.0
for nn in range(-Rm // 2, Rm // 2 + 1):
    v = m.astype(float) ** 2 + 5.0 * nn * nn
    v = v[v > 0]
    tot += np.sum(v ** -2.0)
print(f'   direct sum_(|m|<=3000,|n|<=1500)\' (m^2+5n^2)^-2 = {tot:.10f};  '
      f'zeta(2)L(2,chi_-20)+L(2,chi_-4)L(2,chi_5) = {Ep(np.array([2+0j]))[0].real:.10f}')
mp.mp.dps = 20


def Ep_mp(z):
    def L(chi, qq):
        return mp.power(qq, -z) * sum(chi[a % qq] * mp.zeta(z, mp.mpf(a) / qq)
                                      for a in range(1, qq + 1) if chi[a % qq] != 0)
    return mp.zeta(z) * L(chi_m20, 20) + L(chi_m4, 4) * L(chi_5, 5)


def LamEp(z):
    return (mp.sqrt(20) / (2 * mp.pi)) ** z * mp.gamma(z) * Ep_mp(z)


for z in [mp.mpc(0.3, 7), mp.mpc(0.8, 33)]:
    a, b = LamEp(z), LamEp(1 - z)
    print(f'   functional equation Lambda(s) = Lambda(1-s) at s = {mp.nstr(z, 4)}: '
          f'|diff|/|Lambda| = {mp.nstr(abs(a - b) / abs(a), 3)}')
wcount, mx = rect_count(Ep, 0.52, 1.6, 1.0, 300.0, dt=0.01, ds=0.005)
print(f'   zeros of the Epstein zeta in [0.52, 1.6] x [1, 300] (argument principle): {wcount:.4f} '
      f'(max phase step {mx:.3f})')
if round(wcount) > 0:
    for t0 in range(1, 300, 10):
        wc, _ = rect_count(Ep, 0.52, 1.6, float(t0), float(t0 + 10), dt=0.01, ds=0.005)
        n = int(round(wc))
        if n > 0:
            sg = np.arange(0.53, 1.6, 0.01)
            tg = np.arange(t0, t0 + 10, 0.01)
            SS, TT2 = np.meshgrid(sg, tg)
            A = np.abs(Ep(SS + 1j * TT2))
            idx = np.argsort(A.ravel())[:60]
            cand = []
            for ii in idx:
                z0 = complex(SS.ravel()[ii], TT2.ravel()[ii])
                try:
                    z = complex(mp.findroot(Ep_mp, mp.mpc(z0.real, z0.imag), tol=1e-15))
                except Exception:
                    continue
                if 0.52 < z.real < 1.6 and t0 <= z.imag <= t0 + 10 and all(abs(z - c) > 1e-6 for c in cand):
                    cand.append(z)
                if len(cand) >= n:
                    break
            print(f'     t in [{t0},{t0+10}]: {n} zero(s):',
                  ', '.join(f'{z.real:.6f}{z.imag:+.6f}i' for z in cand))

# ----------------------------------------------------------------------------
hdr('D5  the E-game\'s own zeta: dynamical zeta of the canonical escape Psi')
print('   Psi has full branches (h,k), h in {1,1/2}, k >= 1, each a bijection of its shell onto Z_3^x,')
print('   expanding by 3^k, with price rho_h(k); rho_(1/2)(k) = 2 rho_1(k) (the 1/2-thread pays one extra')
print('   doubling).  For the weighted full shift the Ruelle zeta is 1/(1 - W), W(x;theta) =')
print('   (1 + 2^theta) sum_k rho_1(k)^theta x^k, x = 3^-s  (rho_1(1)=1/3, rho_1(2)=4/9, rho_1(k) =')
print('   2^floor(k log2 3)/3^k = 2^-{k log2 3} for k>=3).  Its poles are the "resonances".')
KMAXW = 4000
kk = np.arange(1, KMAXW + 1)
rho1 = np.empty(KMAXW)
rho1[0] = 1 / 3
rho1[1] = 4 / 9
for k in range(3, KMAXW + 1):
    Kf = (3 ** k).bit_length() - 1          # floor(k log2 3), exact
    rho1[k - 1] = math.exp(Kf * LOG2 - k * LOG3)


def Wser(x, th, K=KMAXW):
    return (1 + 2 ** th) * np.sum(rho1[:K] ** th * x ** kk[:K])


def s_of_theta(th):
    lo, hi = 1e-9, 0.999
    for _ in range(200):
        mid = (lo + hi) / 2
        if Wser(mid, th) < 1:
            lo = mid
        else:
            hi = mid
    x0 = (lo + hi) / 2
    return -math.log(x0) / LOG3, x0


ths = np.linspace(0, 3, 301)
svals = [s_of_theta(t)[0] for t in ths]
i_min = int(np.argmin(svals))
lo, hi = ths[max(0, i_min - 2)], ths[min(300, i_min + 2)]
for _ in range(60):
    m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
    if s_of_theta(m1)[0] < s_of_theta(m2)[0]:
        hi = m2
    else:
        lo = m1
thstar = (lo + hi) / 2
sstar, x0 = s_of_theta(thstar)
print(f'   min_theta s(theta) = {sstar:.5f} at theta* = {thstar:.4f}  (endgame note: dim_H Bad_Psi <= 0.64150)')
# Lundberg exponent: Haar mass s = 1, solve W(1/3; theta) = 1
lo, hi = 1.0, 20.0
for _ in range(100):
    mid = (lo + hi) / 2
    if Wser(1 / 3, mid) < 1:
        lo = mid
    else:
        hi = mid
thL = (lo + hi) / 2
print(f'   Lundberg exponent: W(1/3; theta_L) = 1 at theta_L = {thL:.4f}  (endgame note: 8.6434)')
for th, name in [(thstar, 'theta* (dimension)'), (thL, 'theta_L (Lundberg)'), (1.0, 'theta = 1 (prices)')]:
    s0v, x0v = s_of_theta(th)
    out = []
    for K in (150, 300):
        c = (1 + 2 ** th) * rho1[:K] ** th
        poly = np.concatenate([[-1.0], c])[::-1]      # highest degree first: c_K x^K + ... + c_1 x - 1
        rts = np.roots(poly)
        out.append(rts[np.abs(rts) < 0.8])
    stable = []
    for r in out[1]:
        if np.min(np.abs(out[0] - r)) < 1e-7:
            # polish with the long series
            z = complex(r)
            for _ in range(30):
                fz = (1 + 2 ** th) * np.sum(rho1 ** th * z ** kk) - 1
                dfz = (1 + 2 ** th) * np.sum(rho1 ** th * kk * z ** (kk - 1))
                z = z - fz / dfz
            stable.append(z)
    stable = sorted(stable, key=abs)
    print(f'   {name}: leading pole x_0 = {x0v:.6f} (s = {s0v:.5f}); stable poles with |x| < 0.8 '
          f'({len(stable)}), as |x| and Re s = -log|x|/log 3:')
    for z in stable[:10]:
        print(f'      x = {z.real:+.6f}{z.imag:+.6f}i   |x| = {abs(z):.6f}   Re s = {-math.log(abs(z))/LOG3:+.5f}'
              f'   (Ihara-type critical value s_0/2 = {s0v/2:.5f})')

print()
print(f'[dirichlet] done in {time.time() - T0:.1f} s')
