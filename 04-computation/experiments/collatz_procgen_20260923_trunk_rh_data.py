#!/usr/bin/env python3
"""
collatz_procgen_20260923_trunk_rh_data.py

Candidate "bridges" built from Collatz / E-SCC data, tested for RH-type
structure against controls (session collatz-procgen-20260922, trunk/RH lane,
2026-09-23).

Data (m <= N = 10^6 for T2/T3; the T1 exponents use sigma and the Psi statistics to 10^7):
  d(m)   number of canonical-escape (Psi) steps until the value drops below m
         (m prime to 3; the backward E-game of E-SCC Q2)
  L1(m)  number of expensive links (1/2-transfers with precision k >= 3) used
  H(m)   number of hostile landings (steps taken from a point with k >= 3,
         i.e. m = 1 or 14 mod 27 at that step)
  sig(m) Collatz shortcut stopping time (first j with T^j(m) < m), m >= 2
  tau(m) Collatz shortcut total stopping time (steps to reach 1)
  trunk  indicator of the trunk numbers T_i = (4^i-1)/3, i >= 1

Tests:
  T1  growth exponent of the summatory fluctuation S(x) = sum_{m<=x} (a_m - mean)
      = abscissa of convergence of the fluctuation Dirichlet series (path fit on
      [1e4,1e7]) and the variance-time (Hurst) exponent of block sums; null =
      20 shuffles; controls Moebius, Liouville, random signs, periodic 3-adic
      sequences.  Heuristic prediction for digit-determined statistics:
      alpha = H = (dimension of the exceptional set)/2
  T2  functional-equation detector: does Lambda(s) = gamma-factor * sum a_n n^-s
      satisfy Lambda(s) = eps Lambda(1-s) for some degree-1 (q <= 300, kappa in
      {0,1}) or degree-2 (weight w <= 12, level q <= 150) gamma factor?
      Equivalent theta-modularity relations are fitted by least squares on
      t in [1/12, 12] with polar terms allowed; the residual is normalized by
      the part of theta(1/t) orthogonal to the polar terms (fits where that part
      vanishes are degenerate and skipped); controls zeta, L(chi_-4), L(chi_5),
      Delta, E_4, 1e-6 perturbations (sensitivity) and random signs
  T3  zeros of the Dirichlet polynomials sum_{n<=1000} a_n n^-s in the box
      [-1.5, 2] x [2, 40]: distribution of real parts; controls zeta_N, mu, random

Runtime ~5 minutes, memory < 500 MB (10^7 arrays processed in chunks).
"""
import math
import time
import numpy as np

T0 = time.time()
N = 10 ** 6
rng = np.random.default_rng(20260923)


def hdr(t):
    print()
    print('=' * 78)
    print(t)
    print('=' * 78)


# ----------------------------------------------------------------------------
hdr('C0  data: Psi-descent statistics (E-SCC backward game) and Collatz stopping times, m <= 10^6')
KSTAR = [0, 2] + [(3 ** (s + 1)).bit_length() - 1 for s in range(2, 400)]   # K*(s)
d = np.zeros(N + 1, dtype=np.int32)
L1 = np.zeros(N + 1, dtype=np.int32)
H = np.zeros(N + 1, dtype=np.int32)
maxd = (0, 0)
for m in range(2, N + 1):
    if m % 3 == 0:
        continue
    x = m
    steps = l1 = hl = 0
    while True:
        if x % 3 == 1:
            u = x - 1
            k = 0
            while u % 3 == 0:
                u //= 3
                k += 1
            if k >= 3:
                hl += 1
            x = u << KSTAR[k - 1]
        else:
            w = 2 * x - 1
            k = 0
            while w % 3 == 0:
                w //= 3
                k += 1
            if k >= 3:
                hl += 1
                l1 += 1
            x = w << KSTAR[k - 1]
        steps += 1
        if x < m:
            break
        if steps > 2000:
            raise RuntimeError(m)
    d[m], L1[m], H[m] = steps, l1, hl
    if steps > maxd[0]:
        maxd = (steps, m)
units = np.array([(m % 3 != 0) and m >= 2 for m in range(N + 1)])
print(f'Psi descent: all {units.sum()} m <= 10^6 prime to 3 descend; mean d = {d[units].mean():.4f}, '
      f'max d = {maxd[0]} at m = {maxd[1]}; mean L1 = {L1[units].mean():.4f}; mean H = {H[units].mean():.4f}')
print('   (cross-check with the endgame lane: mean Psi steps on m = 14 mod 27 up to 1e11 was 2.2848)')
m14 = np.arange(14, N + 1, 27)
print(f'   mean d on m = 14 (mod 27), m <= 10^6: {d[m14].mean():.4f}')
sig = np.zeros(N + 1, dtype=np.int32)
tau = np.zeros(N + 1, dtype=np.int32)
for m in range(2, N + 1):
    x = m
    j = 0
    while x >= m:
        x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
        j += 1
    sig[m] = j
    tau[m] = j + tau[x]
print(f'Collatz shortcut: mean stopping time sigma = {sig[2:].mean():.4f} (max {sig.max()} at m = {int(sig.argmax())});'
      f' mean total stopping time tau = {tau[2:].mean():.3f} (max {tau.max()} at m = {int(tau.argmax())})')
trunk = np.zeros(N + 1)
for i in range(1, 12):
    t = (4 ** i - 1) // 3
    if t <= N:
        trunk[t] = 1
# arithmetic controls
mu = np.ones(N + 1, dtype=np.int8)
mu[0] = 0
isprime = np.ones(N + 1, dtype=bool)
isprime[:2] = False
for p in range(2, int(N ** 0.5) + 1):
    if isprime[p]:
        isprime[p * p::p] = False
for p in np.nonzero(isprime)[0]:
    mu[p::p] *= -1
    if p * p <= N:
        mu[p * p::p * p] = 0
# Liouville via Omega (count prime factors with multiplicity)
Om = np.zeros(N + 1, dtype=np.int8)
for p in np.nonzero(isprime)[0]:
    pk = p
    while pk <= N:
        Om[pk::pk] += 1
        pk *= p
lam = np.where(Om % 2 == 0, 1, -1).astype(np.int8)
lam[0] = 0
del Om, isprime
allm = np.ones(N + 1, dtype=bool)
allm[:2] = False
print(f'controls built: Moebius (M(10^6) = {int(mu[1:].sum())}), Liouville (L(10^6) = {int(lam[1:].sum())})')

# ----------------------------------------------------------------------------
hdr('T1  fluctuation exponents of the summatory functions S(x) = sum_{m<=x} (a_m - mean)')
print('    alpha = growth exponent of S = abscissa of convergence of sum (a_m - mean) m^-s.')
print('    RH <=> alpha(Moebius) = 1/2 (Littlewood); every weakly dependent sequence also gives 1/2.')
N7 = 10 ** 7


def psi_descent_vec(lo, hi):
    """vectorized Psi descent for lo <= m < hi (m prime to 3, m >= 2): steps, L1, H"""
    m = np.arange(lo, hi, dtype=np.int64)
    keep = (m % 3 != 0) & (m >= 2)
    m = m[keep]
    x = m.copy()
    st = np.zeros(len(m), np.int16)
    l1 = np.zeros(len(m), np.int16)
    hl = np.zeros(len(m), np.int16)
    act = np.arange(len(m))
    KS = np.array(KSTAR, dtype=np.int64)
    while len(act):
        xa = x[act]
        c1 = (xa % 3 == 1)
        num = np.where(c1, xa - 1, 2 * xa - 1)
        k = np.zeros(len(act), np.int64)
        z = (num % 3 == 0)
        while z.any():
            num = np.where(z, num // 3, num)
            k += z
            z = (num % 3 == 0)
        y = num << KS[k - 1]
        big = k >= 3
        hl[act] += big
        l1[act] += big & ~c1
        st[act] += 1
        x[act] = y
        act = act[y >= m[act]]
    return m, st, l1, hl


def sigma_vec(lo, hi):
    m = np.arange(max(lo, 2), hi, dtype=np.int64)
    x = m.copy()
    st = np.zeros(len(m), np.int16)
    act = np.arange(len(m))
    while len(act):
        xa = x[act]
        xa = np.where(xa % 2 == 0, xa // 2, (3 * xa + 1) // 2)
        x[act] = xa
        st[act] += 1
        act = act[xa >= m[act]]
    return m, st


d7 = np.zeros(N7 + 1, np.int16)
L17 = np.zeros(N7 + 1, np.int16)
H7 = np.zeros(N7 + 1, np.int16)
s7 = np.zeros(N7 + 1, np.int16)
CH7 = 10 ** 6
for lo in range(0, N7 + 1, CH7):
    hi = min(N7 + 1, lo + CH7)
    mm_, a_, b_, c_ = psi_descent_vec(lo, hi)
    d7[mm_], L17[mm_], H7[mm_] = a_, b_, c_
    mm_, a_ = sigma_vec(lo, hi)
    s7[mm_] = a_
assert np.array_equal(d7[:N + 1], d.astype(np.int16)) and np.array_equal(s7[:N + 1], sig.astype(np.int16))
units7 = np.zeros(N7 + 1, bool)
units7[2:] = (np.arange(2, N7 + 1) % 3 != 0)
all7 = np.zeros(N7 + 1, bool)
all7[2:] = True
print(f'    vectorized recomputation to 10^7 agrees with the loop code on m <= 10^6; max Psi steps to 10^7:'
      f' {int(d7.max())} (m = {int(d7.argmax())}); max sigma to 10^7: {int(s7.max())} (m = {int(s7.argmax())})')
# Moebius / Liouville to 10^7
mu7 = np.ones(N7 + 1, dtype=np.int8)
mu7[0] = 0
ip = np.ones(N7 + 1, dtype=bool)
ip[:2] = False
for p in range(2, int(N7 ** 0.5) + 1):
    if ip[p]:
        ip[p * p::p] = False
Om7 = np.zeros(N7 + 1, dtype=np.int8)
for p in np.nonzero(ip)[0]:
    mu7[p::p] *= -1
    if p * p <= N7:
        mu7[p * p::p * p] = 0
    pk = p
    while pk <= N7:
        Om7[pk::pk] += 1
        pk *= p
lam7 = np.where(Om7 % 2 == 0, 1, -1).astype(np.int8)
lam7[0] = 0
del Om7, ip
XS = np.unique(np.round(np.logspace(4, 7, 31)).astype(np.int64))
LB = [2 ** j for j in range(4, 19)]          # block lengths (powers of 2: neutral for the 3-adic data)


def path_exponents(f):
    """memory-lean: fits log max|S| and log rms(S) against log x on XS"""
    S = np.cumsum(f[1:], dtype=np.float64)
    A = np.abs(S)
    np.maximum.accumulate(A, out=A)
    Mv = A[XS - 1].copy()
    np.multiply(S, S, out=A)
    np.cumsum(A, out=A)
    Rv = np.sqrt(A[XS - 1] / XS)
    del A, S
    lx = np.log(XS)
    return (np.polyfit(lx, np.log(Mv + 1e-300), 1)[0],
            np.polyfit(lx, np.log(Rv + 1e-300), 1)[0])


def block_vars(f):
    g = f[1:]
    v = []
    for L in LB:
        nb = len(g) // L
        bs = g[:nb * L].reshape(nb, L).sum(axis=1)
        v.append(bs.var())
    return np.array(v)


def hurst(f, lo=0, hi=None):
    v = block_vars(f)
    lb = np.log(LB)[lo:hi]
    return 0.5 * np.polyfit(lb, np.log(v[lo:hi] + 1e-300), 1)[0]


def hurst_split(f):
    v = np.log(block_vars(f) + 1e-300)
    lb = np.log(LB)
    return (0.5 * np.polyfit(lb[:8], v[:8], 1)[0], 0.5 * np.polyfit(lb[7:], v[7:], 1)[0])


def centered7(a, mask):
    f = a.astype(np.float64)
    mbar = f[mask].mean()
    f -= mbar
    f[~mask] = 0.0
    f[0] = 0.0
    return f, mbar


print(f'{"sequence (m <= 10^7)":44s} {"mean":>7s} {"a(max)":>7s} {"a(rms)":>7s} {"H":>6s} |'
      f' {"shuffle null: a(max)":>21s} {"a(rms)":>13s} {"H":>13s}')
NSH = 20
SPLIT = []
h95 = -(0.63092975357 * math.log2(0.63092975357) + 0.36907024643 * math.log2(0.36907024643))
for name, a, mask in [('Psi steps d(m)', d7, units7), ('L1 links', L17, units7),
                      ('hostile landings H', H7, units7), ('Collatz stopping time sigma', s7, all7)]:
    f, mbar = centered7(a, mask)
    am, ar = path_exponents(f)
    hh = hurst(f)
    SPLIT.append((name, hurst_split(f)))
    nulls = []
    vals = f[mask]
    for _ in range(NSH):
        rng.shuffle(vals)
        f[mask] = vals            # reuse the buffer: zeros off the support stay zero
        nulls.append(path_exponents(f) + (hurst(f),))
    del f, vals
    nulls = np.array(nulls)
    mu_, sd_ = nulls.mean(axis=0), nulls.std(axis=0)
    print(f'{name:44s} {mbar:7.4f} {am:7.3f} {ar:7.3f} {hh:6.3f} |'
          f'  {mu_[0]:6.3f} +- {sd_[0]:5.3f}   {mu_[1]:5.3f} +- {sd_[1]:4.3f}  {mu_[2]:5.3f} +- {sd_[2]:4.3f}')
for name, a in [('Moebius mu (RH <=> alpha = 1/2)', mu7), ('Liouville lambda', lam7)]:
    f = a.astype(np.float64)
    am, ar = path_exponents(f)
    print(f'{name:44s} {"":7s} {am:7.3f} {ar:7.3f} {hurst(f):6.3f}')
r7 = rng.choice([-1.0, 1.0], size=N7 + 1)
r7[0] = 0
am, ar = path_exponents(r7)
print(f'{"random signs":44s} {"":7s} {am:7.3f} {ar:7.3f} {hurst(r7):6.3f}')
per7 = np.zeros(N7 + 1)
per7[1:] = (np.arange(1, N7 + 1) % 27 == 1) - 1 / 27
am, ar = path_exponents(per7)
print(f'{"1_(m = 1 mod 27) - 1/27 (periodic 3-adic)":44s} {"":7s} {am:7.3f} {ar:7.3f} {hurst(per7):6.3f}')
vv = np.zeros(N7 + 1)
t3 = np.arange(N7 + 1, dtype=np.int64)
t3[0] = 1
z = (t3 % 3 == 0)
while z.any():
    vv += z
    t3 = np.where(z, t3 // 3, t3)
    z = (t3 % 3 == 0)
vv[1:] -= 0.5
vv[0] = 0
am, ar = path_exponents(vv)
print(f'{"v_3(m) - 1/2 (3-adic digit function)":44s} {"":7s} {am:7.3f} {ar:7.3f} {hurst(vv):6.3f}')
del t3, vv, per7, r7
print('  scale dependence of H (block lengths 2^4..2^11 vs 2^11..2^18):')
for name, (hlo, hhi) in SPLIT:
    print(f'     {name:40s} H_low = {hlo:.3f}   H_high = {hhi:.3f}')
# exact counts N_k of parity classes mod 2^k with no multiplicative descent within k shortcut steps
cnt = {0: 1}
NK = [1]
for j in range(1, 41):
    new = {}
    for a_, c_ in cnt.items():
        for b_ in (0, 1):
            if 3 ** (a_ + b_) > 2 ** j:
                new[a_ + b_] = new.get(a_ + b_, 0) + c_
    cnt = new
    NK.append(sum(cnt.values()))
rate = [math.log2(NK[k + 1] / NK[k]) for k in range(40)]
print(f'  Collatz: exact N_k (classes mod 2^k without descent) N_22 = {NK[22]}, N_26 = {NK[26]} (= the choice-ladder counts);')
print(f'  local growth rate log2(N_(k+1)/N_k) averaged over k = 4..10: {np.mean(rate[4:11]):.3f}, k = 11..18: '
      f'{np.mean(rate[11:19]):.3f}, k = 30..39: {np.mean(rate[30:40]):.3f}  (-> h(log_3 2) = {h95:.4f})')
print(f'  predicted H = rate/2: {np.mean(rate[4:11])/2:.3f} (low), {np.mean(rate[11:19])/2:.3f} (high), limit {h95/2:.4f}.')
print('  Psi: the endgame lane counts 24, 274, 3050, 31438 alive classes at n = 7, 11, 15, 19 give local rates')
print(f'  {math.log(274/24,3)/4:.3f}, {math.log(3050/274,3)/4:.3f}, {math.log(31438/3050,3)/4:.3f} -> predicted H about '
      f'{math.log(31438/24,3)/12/2:.3f} (limit <= 0.32075).')
print(f'  Heuristic (random placement of the exceptional residue classes; see note): for a statistic that is')
print(f'  determined by the digits consumed, alpha = H = (exceptional dimension)/2.  Collatz: h(log_3 2)/2 = {h95/2:.5f};')
print(f'  Psi descent: dim Bad_Psi / 2 <= 0.64150/2 = 0.32075.  RH plays no role in either number.')
del d7, L17, H7, s7, mu7, lam7, units7, all7

# ----------------------------------------------------------------------------
hdr('T2  functional-equation detector (theta-modularity), degree 1 and degree 2')
NA = 3000
tgrid = np.exp(np.linspace(math.log(1 / 12), math.log(12), 81))


def _fit(y, first, polar):
    """least squares y ~ b0*first + polar @ c; residual normalized by the part of y orthogonal to polar"""
    Pn = polar / np.linalg.norm(polar, axis=0)
    cp, *_ = np.linalg.lstsq(Pn, y, rcond=None)
    yperp = y - Pn @ cp
    ny, np_ = np.linalg.norm(y), np.linalg.norm(yperp)
    if np_ < 1e-9 * ny or np.linalg.norm(first) == 0:
        return None
    A = np.hstack([first[:, None] / np.linalg.norm(first), Pn])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = np.linalg.norm(A @ beta - y) / np_
    eps = beta[0] / np.linalg.norm(first)
    return res, eps


def fe_deg1(a, qmax=300):
    best = (1e9, None)
    ndeg = 0
    a = np.asarray(a[1:NA + 1], dtype=float)
    n = np.arange(1, NA + 1, dtype=float)
    tmin = tgrid[0]
    for kap in (0, 1):
        for q in range(1, qmax + 1):
            nmax = min(NA, int(math.sqrt(45 * q / (math.pi * tmin))) + 2)
            nn = n[:nmax]
            c = a[:nmax] * nn ** kap
            if not np.any(c):
                continue
            th = np.exp(-math.pi * np.outer(tgrid, nn ** 2) / q) @ c
            thi = np.exp(-math.pi * np.outer(1 / tgrid, nn ** 2) / q) @ c
            polar = np.vstack([np.ones_like(tgrid), tgrid ** 0.5, tgrid, tgrid ** 1.5]).T
            r = _fit(thi, tgrid ** (kap + 0.5) * th, polar)
            if r is None:
                ndeg += 1
                continue
            if r[0] < best[0]:
                best = (r[0], (q, kap, r[1]))
    return best, ndeg


def fe_deg2(a, qmax=150, wmax=12):
    best = (1e9, None)
    ndeg = 0
    a = np.asarray(a[1:NA + 1], dtype=float)
    n = np.arange(1, NA + 1, dtype=float)
    ymin = tgrid[0]
    for w in range(1, wmax + 1):
        for q in range(1, qmax + 1):
            rq = math.sqrt(q)
            nmax = min(NA, int(rq / (2 * math.pi * ymin) * (44 + 0.5 * (w - 1) * 8.0)) + 20)
            nn = n[:nmax]
            c = a[:nmax] * nn ** ((w - 1) / 2)
            if not np.any(c):
                continue
            g = np.exp(-2 * math.pi * np.outer(tgrid, nn) / rq) @ c
            gi = np.exp(-2 * math.pi * np.outer(1 / tgrid, nn) / rq) @ c
            polar = np.vstack([np.ones_like(tgrid), tgrid ** w]).T
            r = _fit(gi, tgrid ** w * g, polar)
            if r is None:
                ndeg += 1
                continue
            if r[0] < best[0]:
                best = (r[0], (q, w, r[1]))
    return best, ndeg


def chi_seq(chi, q):
    s = np.zeros(NA + 1)
    for nn in range(1, NA + 1):
        s[nn] = chi[nn % q]
    return s


ones = np.ones(NA + 1)
ones[0] = 0
zeta_pert = ones.copy()
zeta_pert[2] += 1e-6
# Ramanujan tau and sigma_3 via q-expansions
M2 = NA + 1
qs = np.zeros(M2, dtype=object)
# Delta = q prod (1-q^n)^24, integer arithmetic (object) for exactness up to NA terms
prod = [0] * M2
prod[0] = 1
for nn in range(1, M2):
    # multiply by (1 - q^nn)^24 via 24 successive multiplications by (1 - q^nn)
    for _ in range(24):
        for j in range(M2 - 1, nn - 1, -1):
            prod[j] -= prod[j - nn]
taus = np.zeros(NA + 1)
for nn in range(1, NA + 1):
    taus[nn] = float(prod[nn - 1])
print(f'   Ramanujan tau(1..6) = {[int(taus[k]) for k in range(1, 7)]}')
taus_n = taus / np.arange(0, NA + 1, dtype=float).clip(1) ** 5.5
sig3 = np.zeros(NA + 1)
for dd in range(1, NA + 1):
    sig3[dd::dd] += dd ** 3
del_pert = taus.copy()
del_pert[2] *= (1 + 1e-6)

tests = [('CONTROL zeta (a_n = 1)', ones, 1), ('CONTROL zeta with a_2 perturbed by 1e-6', zeta_pert, 1),
         ('CONTROL L(s, chi_-4)', chi_seq([0, 1, 0, -1], 4), 1),
         ('CONTROL L(s, chi_5)', chi_seq([0, 1, -1, -1, 1], 5), 1),
         ('CONTROL random signs', rng.choice([-1.0, 1.0], size=NA + 1), 1),
         ('CONTROL Delta (weight 12), a_n = tau(n)/n^5.5', taus_n, 2),
         ('CONTROL Delta with tau(2) perturbed by 1e-6', del_pert / np.arange(0, NA + 1).clip(1) ** 5.5, 2),
         ('CONTROL E_4, a_n = sigma_3(n)/n^1.5', sig3 / np.arange(0, NA + 1).clip(1) ** 1.5, 2),
         ('trunk indicator', trunk[:NA + 1], 0), ('Psi steps d(m)', d[:NA + 1].astype(float), 0),
         ('L1 links', L1[:NA + 1].astype(float), 0), ('hostile landings H', H[:NA + 1].astype(float), 0),
         ('Collatz stopping time sigma', sig[:NA + 1].astype(float), 0),
         ('Collatz total stopping time tau', tau[:NA + 1].astype(float), 0)]
print(f'{"sequence":48s} {"best deg-1 residual (q,kappa,eps)":>44s} {"best deg-2 residual (q,w,eps)":>44s}')
for name, a, which in tests:
    (r1, nd1) = fe_deg1(a) if which in (0, 1) else ((float('nan'), None), 0)
    (r2, nd2) = fe_deg2(a) if which in (0, 2) else ((float('nan'), None), 0)
    s1 = f'{r1[0]:.1e} {r1[1][0]},{r1[1][1]},{r1[1][2]:+.4f} [{nd1} deg]' if r1[1] else '-'
    s2 = f'{r2[0]:.1e} {r2[1][0]},{r2[1][1]},{r2[1][2]:+.4f} [{nd2} deg]' if r2[1] else '-'
    print(f'{name:48s} {s1:>44s} {s2:>44s}')
print('  Residual = ||fit - theta(1/t)|| / ||part of theta(1/t) orthogonal to the polar terms|| over 81 t in [1/12, 12];')
print('  [n deg] = number of (q, .) pairs skipped because that orthogonal part vanishes (no arithmetic signal).')
print('  A genuine functional equation gives ~1e-12 or less with eps = +-1; one coefficient changed by 1e-6 gives')
print('  a residual orders of magnitude larger; an unrelated sequence gives O(1e-2..1).')

# ----------------------------------------------------------------------------
hdr('T3  zeros of Dirichlet polynomials sum_{n<=1000} a_n n^-s in [-1.5, 2] x [2, 40]')
NP = 1000
ds_, dt_ = 0.025, 0.025
sg = np.arange(-1.5, 2.0 + 1e-9, ds_)
tg = np.arange(2.0, 40.0 + 1e-9, dt_)
SS, TT = np.meshgrid(sg, tg)
Z = SS + 1j * TT
logn = np.log(np.arange(1, NP + 1, dtype=float))
zeta_zeros = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178]


def wrapd(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def dp_zeros(a):
    a = np.asarray(a[1:NP + 1], dtype=float)
    V = np.zeros(Z.shape, dtype=complex)
    for i in np.nonzero(a)[0]:
        V += a[i] * np.exp(-Z * logn[i])
    ph = np.angle(V)
    w = (wrapd(ph[:-1, 1:] - ph[:-1, :-1]) + wrapd(ph[1:, 1:] - ph[:-1, 1:]) +
         wrapd(ph[1:, :-1] - ph[1:, 1:]) + wrapd(ph[:-1, :-1] - ph[1:, :-1])) / (2 * np.pi)
    cells = np.argwhere(np.rint(w).astype(int) != 0)
    nz = np.nonzero(a)[0]
    roots = []
    for (r, c) in cells:
        z = complex(sg[c] + ds_ / 2, tg[r] + dt_ / 2)
        for _ in range(40):
            e = np.exp(-z * logn[nz])
            fz = np.sum(a[nz] * e)
            dfz = -np.sum(a[nz] * logn[nz] * e)
            if dfz == 0:
                break
            z = z - fz / dfz
        roots.append(z)
    return np.array(roots), int(np.rint(w).sum())


print(f'{"sequence":34s} {"#cells":>6s} {"Re s range":>18s} {"|Re-1/2|<0.05":>14s} {"median|Re-1/2|":>15s}'
      f' {"dist to zeta zeros (t<40)":>28s}')
seqs = [('zeta partial sum (a_n = 1)', ones), ('Moebius mu', mu[:NP + 1].astype(float)),
        ('random signs', rng.choice([-1.0, 1.0], size=NP + 1)),
        ('trunk indicator (5 terms)', trunk[:NP + 1]),
        ('Psi steps d(m)', d[:NP + 1].astype(float)), ('L1 links', L1[:NP + 1].astype(float)),
        ('hostile landings H', H[:NP + 1].astype(float)),
        ('Collatz stopping time sigma', sig[:NP + 1].astype(float))]
for name, a, msk in [('Psi steps d(m)', d, units), ('Collatz stopping time sigma', sig, allm)]:
    for rep in range(10):
        b = a[:NP + 1].astype(float).copy()
        idx = np.nonzero(msk[:NP + 1])[0]
        b[idx] = rng.permutation(b[idx])
        seqs.append((f'{name} shuffled #{rep}', b))
SHUF = {}
for name, a in seqs:
    roots, cnt = dp_zeros(a)
    if len(roots) == 0:
        print(f'{name:34s} {0:6d}   no zeros in the box')
        continue
    re = roots.real
    dist = [np.min(np.abs(roots - complex(0.5, g))) for g in zeta_zeros]
    if 'shuffled' in name:
        key = name.split(' shuffled')[0]
        SHUF.setdefault(key, []).append((np.mean(np.abs(re - 0.5) < 0.05), np.median(np.abs(re - 0.5)), dist))
        continue
    print(f'{name:34s} {len(roots):6d} [{re.min():+7.3f},{re.max():+7.3f}] '
          f'{np.mean(np.abs(re - 0.5) < 0.05):14.3f} {np.median(np.abs(re - 0.5)):15.3f}'
          f'   {" ".join(f"{x:.2f}" for x in dist)}')
for key, rows in SHUF.items():
    fr = np.array([r[0] for r in rows])
    md = np.array([r[1] for r in rows])
    ds = np.array([r[2] for r in rows])
    print(f'{key + " (10 mean-matched shuffles)":34s}        {"":18s} {fr.mean():8.3f}+-{fr.std():.3f} {md.mean():9.3f}+-{md.std():.3f}'
          f'   {" ".join(f"{x:.2f}" for x in ds.mean(axis=0))}  (mean)')
    print(f'{"":34s}        {"":18s} {"":14s} {"":15s}   {" ".join(f"{x:.2f}" for x in ds.min(axis=0))}  (min over shuffles)')
print('  The zeros of zeta itself in this window are the 6 points 1/2 + i*gamma listed above (all on the line);')
print('  finite Dirichlet polynomials, including the partial sums of zeta, show no alignment with Re s = 1/2.')
print('  A sequence with mean c contains c*zeta_N(s) (times (1-3^-s) on units), so its zeros are perturbed zeta_N')
print('  zeros; the shuffles show what that alone produces.')

print()
print(f'[data] done in {time.time() - T0:.1f} s')
