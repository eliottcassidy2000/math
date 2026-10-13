"""
collatz_parity_codex_census_crosscheck_20261013.py -- out-of-sample test of the (-1)^tau line model across agents.
The Codex census 05-knowledge/results/collatz_parity_bins_20261013b.json (renewal_shapes_oct11, C++ stream, jump depth 12)
holds 1024 bins per octave over [2^14, 2^34): [lo, hi, count, signed sum, even count, even signed sum].  Our model is
fitted on [2^29, 2^30) from our own tau array (numba, exact to 2^30): nine lines at the predicted effective frequencies
with the predicted decays (THM-4629).  Tests on the Codex octaves 2^30..2^33 (never seen by the fit):
  (1) per line: predicted amplitude a_29 decay^(j-29) and constant phase, against the least-squares amplitude and phase
      of the Codex bin means in octave j;
  (2) the prefix mean A(N) = S(N)/N at the Codex bin edges over 2^30..2^34 (S(2^14) taken from our exact array), against
      the Cesaro shadow of the fitted lines (exact integral formula), including the value at 2^34;
  (3) the Mercator-line amplitude per octave 2^16..2^29 from both censuses (consistency of the two computations).
"""
import json, numpy as np, time
from numba import njit
from mpmath import mp, mpc, exp, findroot, log, pi as mpi

mp.dps = 30
L = log(3) / log(2) - 1
LN2 = np.log(2.0)
LINES = [2.5, 14.5, 26.5, 79.5, 120.5, 173.5, 226.5, 279.5, 332.5]

@njit(cache=True)
def tau16(N):
    tau = np.zeros(N + 8, dtype=np.int16)
    for n in range(2, N + 8):
        m = n; steps = 0
        while m >= n:
            if m & 1:
                m = (3 * m + 1) >> 1
            else:
                m >>= 1
            steps += 1
        tau[n] = steps + tau[m]
    return tau

@njit(cache=True)
def bin_top(tau, J, R):
    s = np.zeros(R); c = np.zeros(R)
    for n in range(1 << J, 1 << (J + 1)):
        k = int((np.log2(n) - J) * R)
        if k >= R: k = R - 1
        s[k] += -1.0 if (tau[n] & 1) else 1.0
        c[k] += 1
    return s / np.maximum(c, 1)

@njit(cache=True)
def prefix_sum_to(tau, N):
    s = 0
    for n in range(1, N + 1):
        s += -1 if (tau[n] & 1) else 1
    return s

t0 = time.time()
tau = tau16(1 << 30)
print(f"tau to 2^30 in {time.time()-t0:.1f}s")
d = json.load(open("../../05-knowledge/results/collatz_parity_bins_20261013b.json"))
b = np.array(d["bins"], dtype=float)
lo, hi, cnt, ssum = b[:, 0], b[:, 1], b[:, 2], b[:, 3]
xc = np.log2((lo + hi - 1) / 2.0)                   # bin centre in log2
yc = ssum / cnt
print(f"Codex bins: {len(b)} over [2^{np.log2(lo[0]):.0f}, 2^{np.log2(hi[-1]):.0f}); prefix_sum at 2^34 = {d['prefix_sum']}")

roots = []
for f in LINES:
    s0 = 2j * mpi * f
    w = findroot(lambda w: exp(-(s0 + w)) + exp((s0 + w) * L) + 2, mpc(-0.01, 0))
    roots.append((f, f + float(w.imag) / (2 * np.pi), float(w.real)))

# fit on our own top octave
J, R = 29, 768
y = bin_top(tau, J, R); x = J + (np.arange(R) + 0.5) / R
cols = []
for f, fe, g in roots:
    env = np.exp(g * (x - 29)); cols.append(env * np.cos(2 * np.pi * fe * x)); cols.append(env * np.sin(2 * np.pi * fe * x))
cols.append(np.ones_like(x))
coef, *_ = np.linalg.lstsq(np.column_stack(cols), y, rcond=None)

def ls(xs, ys, f):
    A = np.column_stack([np.cos(2 * np.pi * f * xs), np.sin(2 * np.pi * f * xs), np.ones_like(xs)])
    c, *_ = np.linalg.lstsq(A, ys, rcond=None)
    return float(np.hypot(c[0], c[1])), float(np.arctan2(-c[1], c[0]))

print("\n# (1) lines on the Codex octaves 2^30..2^33 (never seen by the fit): predicted amplitude / phase vs Codex least squares")
print(" line    f_eff     j   pred amp   Codex amp   pred phase  Codex phase")
for i, (f, fe, g) in enumerate(roots):
    c, s = coef[2 * i], coef[2 * i + 1]; a29 = np.hypot(c, s); ph = np.arctan2(-s, c)
    if a29 < 0.03: continue
    for j in (30, 31, 32, 33):
        sel = (xc >= j) & (xc < j + 1)
        am, pm = ls(xc[sel], yc[sel], fe)
        print(f" {f:6.1f} {fe:9.4f}  {j}   {a29 * np.exp(g * (j - 29)):.4f}     {am:.4f}      {ph:+.3f}      {pm:+.3f}")

print("\n# (2) prefix means over 2^30..2^34 from the Codex census vs the Cesaro shadow of our nine lines")
S14 = prefix_sum_to(tau, 1 << 14)
S_edges = S14 + np.concatenate([[0.0], np.cumsum(ssum)])   # prefix sums at lo[0], lo[1], ..., hi[-1]
N_edges = np.concatenate([lo, [hi[-1]]])
def model_A(N):
    tot = coef[-1]
    for i, (f, fe, g) in enumerate(roots):
        c, s = coef[2 * i], coef[2 * i + 1]
        p = 1 + g / LN2 + 1j * 2 * np.pi * fe / LN2
        tot += np.exp(-29 * g) * ((c - 1j * s) * (N ** p - 1) / (p * N)).real
    return tot
sel = N_edges >= 2 ** 30
Ns = N_edges[sel][::4]; Aex = S_edges[sel][::4] / Ns
Am = np.array([model_A(float(N)) for N in Ns])
print(f" points {len(Ns)} (256 per octave): exact rms {np.sqrt((Aex**2).mean()):.6f}, max |A| {np.abs(Aex).max():.6f}; model rms {np.sqrt((Am**2).mean()):.6f}; rms(exact - model) {np.sqrt(((Aex-Am)**2).mean()):.6f}")
xs = np.log2(Ns); top = xs >= 33
print(f" exact prefix mean on 2^33..2^34: amplitude at 26.5036 = {ls(xs[top], Aex[top], 26.5036)[0]:.6f} (model {ls(xs[top], Am[top], 26.5036)[0]:.6f})")
for k in (30, 31, 32, 33, 34):
    i = int(np.argmin(np.abs(Ns - 2.0 ** k)))
    print(f" N = 2^{k}: exact {Aex[i]:+.6f}   model {Am[i]:+.6f}")
print(f" check: exact A(2^34) from bins = {S_edges[-1] / N_edges[-1]:+.6f}, Codex prefix_sum/2^34 = {d['prefix_sum'] / 2.0**34:+.6f}")

print("\n# (3) Mercator-line amplitude per octave, Codex census vs ours (ours: 768 bins, sinc-corrected)")
for j in range(16, 30):
    sel = (xc >= j) & (xc < j + 1); am, pm = ls(xc[sel], yc[sel], 26.5036)
    yo = bin_top(tau, j, 768); xo = j + (np.arange(768) + 0.5) / 768; ao, po = ls(xo, yo, 26.5036)
    print(f" 2^{j}: Codex {am / np.sinc(26.5036 / 1024):.4f} ({pm:+.3f})   ours {ao / np.sinc(26.5036 / 768):.4f} ({po:+.3f})")
print("done")
