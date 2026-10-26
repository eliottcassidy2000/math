# Audit of HYP-9391 (c) / check E.  All inputs are the auditor's own: Y sums (ysums.c) at primes q = 3 mod 4 up to
# 9999 and PARI coefficients (forms_dump.gp).
#  1. Reproduce the lane's q^-5-weighted projection on [307, 4999] (coefficients, Gram, moments, max|m| and share
#     |m| > 2 at q >= 2000, small-value shares, twist and mod-8 tables).
#  2. Projection error: covariance of the fitted known part (white-noise model, residual variance from the fit) at
#     every prime; split into the low-order block {1, q, q^2, a84, |b|^2} and the rest; max / mean over q >= 2000.
#  3. Robustness: (i) snapped fit (q^3, q a84, a86, a46 coefficients rounded to integers, refit the rest);
#     (ii) fit on [307, 1999], evaluate on [2000, 9999]; (iii) fit on [307, 4999], evaluate on [5000, 9999];
#     (iv) fit on [307, 9999].  Statistics recomputed each time.
#  4. Reference laws by quadrature: SU(2), USp(4), SU(2) x SU(2) (tensor), Sym^3 SU(2); moments 2, 4, 6, 8, 12,
#     P(|x| < 0.05), P(|x| < 0.2), P(|x| > 2), P(|x| > 3); one-sample KS distances of each m_i to each law.
# Usage: part_c_analysis.py ysums_lo ysums_hi forms ; writes part_c_analysis.out
import sys, os
import numpy as np
from sympy import jacobi_symbol
HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)

rows = {}
for fn in sys.argv[1:3]:
    for line in open(fn):
        r = list(map(int, line.split()))
        if r[0] % 4 == 3: rows[r[0]] = r[1:7]
forms = {}
for line in open(sys.argv[3]):
    p, x84, x86, x46, nb = map(int, line.split()); forms[p] = (x84, x86, x46, nb)
allq = np.array(sorted(rows), dtype=np.int64)
NAMES = ['1', 'q', 'q^2', 'q^3', 'a84', 'q a84', 'a86', 'a46', '|b|^2']
LOW = [0, 1, 2, 4, 8]

def design(qs):
    qf = qs.astype(float)
    a84 = np.array([forms[q][0] for q in qs], float); a86 = np.array([forms[q][1] for q in qs], float)
    a46 = np.array([forms[q][2] for q in qs], float); b2 = np.array([forms[q][3] for q in qs], float)
    B = np.stack([np.ones_like(qf), qf, qf ** 2, qf ** 3, a84, qf * a84, a86, a46, b2], axis=1)
    Y = np.array([[rows[q][0], rows[q][2], rows[q][4], rows[q][5]] for q in qs], float)   # Y02 Y05 Y09 Y18
    return B, Y

def fit(qs, fixed=None):
    """weighted LS, weight q^-5 (i.e. OLS on values / q^(5/2)); fixed: dict col -> (4,) coefficients held fixed"""
    B, Y = design(qs)
    w = qs.astype(float) ** -2.5
    Bw = B * w[:, None]; Yw = Y * w[:, None]
    free = [j for j in range(B.shape[1]) if not fixed or j not in fixed]
    if fixed:
        for j, cvec in fixed.items(): Yw = Yw - np.outer(Bw[:, j], cvec)
    sc = np.linalg.norm(Bw[:, free], axis=0)
    cf, *_ = np.linalg.lstsq(Bw[:, free] / sc, Yw, rcond=None)
    cf = cf / sc[:, None]
    coef = np.zeros((B.shape[1], 4))
    coef[free] = cf
    if fixed:
        for j, cvec in fixed.items(): coef[j] = cvec
    return coef

def newparts(qs, coef):
    B, Y = design(qs)
    w = qs.astype(float) ** -2.5
    R = (Y - B @ coef) * w[:, None]
    return np.stack([R[:, 0], R[:, 2], R[:, 3], R[:, 1] + R[:, 3]], axis=1), R

def stats(qs, m, label):
    P(f"  [{label}] {len(qs)} primes in [{qs.min()}, {qs.max()}]")
    for j, nm in enumerate(['m1=N(Y02)', 'm2=N(Y09)', 'm3=N(Y18)', 'm4=N(Y05)+N(Y18)']):
        x = m[:, j]
        P(f"    {nm}: E2={np.mean(x**2):.3f} E4={np.mean(x**4):.3f} E6={np.mean(x**6):.2f} max|x|={np.max(np.abs(x)):.3f}"
          f" P(|x|>2)={np.mean(np.abs(x) > 2):.3f} P(|x|<0.05)={np.mean(np.abs(x) < 0.05):.3f} P(|x|<0.2)={np.mean(np.abs(x) < 0.2):.3f}")
    G = m.T @ m / len(qs)
    P(f"    Gram(m1..m4): max |G - I| = {np.max(np.abs(G - np.eye(4))):.3f}; diag {np.round(np.diag(G), 3).tolist()}")

# ---------- 1. reproduce the lane on [307, 4999] ----------
P("1. Lane's projection reproduced (primes q = 3 mod 4 in [307, 4999])")
qA = allq[(allq >= 307) & (allq <= 4999)]
cA = fit(qA)
for i, nm in enumerate(['Y02', 'Y05', 'Y09', 'Y18']):
    P(f"  {nm}: " + ', '.join(f"{n}: {c:+.3f}" for n, c in zip(NAMES, cA[:, i])))
mA, RA = newparts(qA, cA)
GR = RA.T @ RA / len(qA)
P("  Gram of raw new parts (Y02, Y05, Y09, Y18):", np.round(GR, 3).tolist())
stats(qA, mA, 'all 307')
big = qA >= 2000
P(f"  at the {big.sum()} primes >= 2000: max|m_i| = {np.round(np.max(np.abs(mA[big]), axis=0), 3).tolist()}, "
  f"share |m_i| > 2 = {np.round(np.mean(np.abs(mA[big]) > 2, axis=0), 3).tolist()}")
for cls in (3, 7):
    sel = (qA % 8 == cls)
    P(f"  q = {cls} mod 8 ({sel.sum()}): E2 = {np.round(np.mean(mA[sel]**2, axis=0), 3).tolist()}, "
      f"E4 = {np.round(np.mean(mA[sel]**4, axis=0), 3).tolist()}, max|x| = {np.round(np.max(np.abs(mA[sel]), axis=0), 3).tolist()}")
chars = {'(2/p)': 2, '(3/p)': 3, '(5/p)': 5, '(7/p)': 7, '(-7/p)': -7, '(6/p)': 6, '(-6/p)': -6, '(10/p)': 10, '(11/p)': 11, '(13/p)': 13}
for nm, D in chars.items():
    e = np.array([jacobi_symbol(D % int(p), int(p)) for p in qA], float)
    M = np.array([[np.mean(e * mA[:, i] * mA[:, j]) for j in range(4)] for i in range(4)])
    P(f"  twist {nm:7s}: max |E[e m_i m_j]| (all i, j) = {np.max(np.abs(M)):.3f}")

# ---------- 2. projection error ----------
P("2. Projection error of the fitted known part (white-noise model)")
B, Y = design(qA)
w = qA.astype(float) ** -2.5
Bw = B * w[:, None]
sc = np.linalg.norm(Bw, axis=0)
Bs = Bw / sc
Cinv = np.linalg.inv(Bs.T @ Bs)
for i, nm in enumerate(['Y02', 'Y05', 'Y09', 'Y18']):
    s2 = np.sum(RA[:, i] ** 2) / (len(qA) - B.shape[1])
    se = np.sqrt(np.diag(Cinv) * s2) / sc
    P(f"  {nm}: residual variance {s2:.3f}; coefficient SEs: " + ', '.join(f"{n}: {s:.3g}" for n, s in zip(NAMES, se)))
s2 = 1.0
lev_all = np.einsum('ij,jk,ik->i', Bs, Cinv, Bs)
Bl = np.zeros_like(Bs); Bl[:, LOW] = Bs[:, LOW]
lev_low = np.einsum('ij,jk,ik->i', Bl, Cinv, Bl)
Bh = Bs - Bl
lev_high = np.einsum('ij,jk,ik->i', Bh, Cinv, Bh)
for lab, sel in (('q >= 2000', qA >= 2000), ('q < 1000', qA < 1000)):
    P(f"  {lab}: SD of fitted known part (sigma = 1): total mean {np.mean(np.sqrt(lev_all[sel])):.3f} max {np.max(np.sqrt(lev_all[sel])):.3f}; "
      f"low-order block {{1,q,q^2,a84,|b|^2}} mean {np.mean(np.sqrt(lev_low[sel])):.3f} max {np.max(np.sqrt(lev_low[sel])):.3f}; "
      f"other block mean {np.mean(np.sqrt(lev_high[sel])):.3f} max {np.max(np.sqrt(lev_high[sel])):.3f}")

# ---------- 3. robustness ----------
P("3. Robustness of the separation")
snap = {3: np.round(cA[3]), 5: np.round(cA[5]), 6: np.round(cA[6]), 7: np.round(cA[7])}
P(f"  (i) snapped columns q^3, q a84, a86, a46 -> {[(NAMES[j], v.tolist()) for j, v in snap.items()]}")
cS = fit(qA, fixed=snap)
for i, nm in enumerate(['Y02', 'Y05', 'Y09', 'Y18']):
    P(f"    {nm}: " + ', '.join(f"{n}: {c:+.3f}" for n, c in zip(NAMES, cS[:, i])))
mS, _ = newparts(qA, cS)
stats(qA, mS, 'snapped fit, [307,4999]')
P(f"    |m(snapped) - m(lane fit)| at q >= 2000: max {np.round(np.max(np.abs(mS[big] - mA[big]), axis=0), 3).tolist()}, "
  f"rms {np.round(np.sqrt(np.mean((mS[big] - mA[big])**2, axis=0)), 3).tolist()}")
qB = allq[(allq >= 307) & (allq < 2000)]; qC = allq[(allq >= 2000)]
cB = fit(qB)
mC_B, _ = newparts(qC, cB)
cC = fit(qC)
mC_C, _ = newparts(qC, cC)
P("  (ii) fit on [307, 1999] only, new parts evaluated on [2000, 9999] (out of sample):")
stats(qC, mC_B, 'fit<2000 -> eval>=2000')
P(f"    vs. in-sample fit on [2000, 9999]: max |difference| {np.round(np.max(np.abs(mC_B - mC_C), axis=0), 3).tolist()}")
qD = allq[(allq >= 5000)]
mD_A, _ = newparts(qD, cA)
P("  (iii) lane's fit ([307, 4999]) extrapolated to [5000, 9999]:")
stats(qD, mD_A, 'fit[307,4999] -> eval[5000,9999]')
mD_S, _ = newparts(qD, cS)
stats(qD, mD_S, 'snapped fit -> eval[5000,9999]')
qE = allq[allq >= 307]
cE = fit(qE)
mE, RE = newparts(qE, cE)
P("  (iv) fit on all primes [307, 9999]:")
for i, nm in enumerate(['Y02', 'Y05', 'Y09', 'Y18']):
    P(f"    {nm}: " + ', '.join(f"{n}: {c:+.3f}" for n, c in zip(NAMES, cE[:, i])))
stats(qE, mE, 'all [307,9999]')
bigE = qE >= 2000
P(f"    at the {bigE.sum()} primes >= 2000: max|m_i| = {np.round(np.max(np.abs(mE[bigE]), axis=0), 3).tolist()}, "
  f"share |m_i| > 2 = {np.round(np.mean(np.abs(mE[bigE]) > 2, axis=0), 3).tolist()}")
for cls in (3, 7):
    sel = (qE % 8 == cls)
    P(f"    q = {cls} mod 8 ({sel.sum()}): E2 = {np.round(np.mean(mE[sel]**2, axis=0), 3).tolist()}, "
      f"E4 = {np.round(np.mean(mE[sel]**4, axis=0), 3).tolist()}, max|x| = {np.round(np.max(np.abs(mE[sel]), axis=0), 3).tolist()}")
for nm, D in chars.items():
    e = np.array([jacobi_symbol(D % int(p), int(p)) for p in qE], float)
    M = np.array([[np.mean(e * mE[:, i] * mE[:, j]) for j in range(4)] for i in range(4)])
    P(f"    twist {nm:7s}: max |E[e m_i m_j]| (all i, j) = {np.max(np.abs(M)):.3f}")
np.save(os.path.join(os.environ.get('AUDIT_SCRATCH', '/tmp'), 'mE.npy'), np.column_stack([qE, mE]))

# ---------- 4. reference laws ----------
P("4. Reference trace laws (quadrature)")
n = 1200
th = (np.arange(n) + 0.5) * np.pi / n
wsu2 = (2 / np.pi) * np.sin(th) ** 2 * (np.pi / n)
xsu2 = 2 * np.cos(th)
T1, T2 = np.meshgrid(th, th, indexing='ij')
wus = (np.cos(T1) - np.cos(T2)) ** 2 * np.sin(T1) ** 2 * np.sin(T2) ** 2
wus = (wus / wus.sum()).ravel()
xus = (2 * np.cos(T1) + 2 * np.cos(T2)).ravel()
X1, X2 = np.meshgrid(xsu2, xsu2, indexing='ij')
W1, W2 = np.meshgrid(wsu2, wsu2, indexing='ij')
xss = (X1 * X2).ravel(); wss = (W1 * W2).ravel()
xs3 = np.sin(4 * th) / np.sin(th); ws3 = wsu2
laws = {'SU(2)': (xsu2, wsu2), 'USp(4)': (xus, wus), 'SU2xSU2': (xss, wss), 'Sym3 SU(2)': (xs3, ws3)}
for nm, (x, wt) in laws.items():
    wt = wt / wt.sum()
    mom = [np.sum(wt * x ** k) for k in (2, 4, 6, 8, 12)]
    P(f"  {nm:10s}: E2={mom[0]:.3f} E4={mom[1]:.3f} E6={mom[2]:.3f} E8={mom[3]:.2f} E12={mom[4]:.1f} "
      f"P(|x|<0.05)={np.sum(wt * (np.abs(x) < 0.05)):.4f} P(|x|<0.2)={np.sum(wt * (np.abs(x) < 0.2)):.4f} "
      f"P(|x|>2)={np.sum(wt * (np.abs(x) > 2)):.4f} P(|x|>3)={np.sum(wt * (np.abs(x) > 3)):.4f}")
    laws[nm] = (x, wt)
def ks(sample, x, wt):
    o = np.argsort(x); xs, cw = x[o], np.cumsum(wt[o])
    s = np.sort(sample); ecdf_hi = np.arange(1, len(s) + 1) / len(s); ecdf_lo = np.arange(len(s)) / len(s)
    F = np.interp(s, xs, cw)
    return max(np.max(np.abs(ecdf_hi - F)), np.max(np.abs(ecdf_lo - F)))
P("  one-sample KS distance D of m_i (fit on [307, 9999], q >= 307) to each law; 5% critical value ~ 1.36/sqrt(n) = "
  f"{1.36 / np.sqrt(len(qE)):.3f}")
for j, nm in enumerate(['m1', 'm2', 'm3', 'm4']):
    P(f"    {nm}: " + ', '.join(f"{ln}: {ks(mE[:, j], x, wt):.3f}" for ln, (x, wt) in laws.items()))
open(os.path.join(HERE, 'part_c_analysis.out'), 'w').write('\n'.join(out) + '\n')
