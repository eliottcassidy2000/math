"""
collatz_parity_fast_modes_20261013.py -- the untested half of THM-4629 statement 2: the class-specific modes of the
3-adic closure decay by the factor 1/2 per octave (|e^(-s)| = 2) at the frequencies alpha/m - l/(2 3^(i-1)).  For the
parity of tau (alpha = 1, m = 2) at depth i = 1 the predicted fast modes sit at f in Z/2, at depth 2 at f in Z/6.
Measured: D_r(u) = m_r(u) - m(u), the class mean of (-1)^tau over n = r (mod 3) (and mod 9) minus the global bin mean,
per octave 2^8..2^29, with the bin count R_j chosen so that every class bin has >= 48 integers; the rms of D over the
bins, the sampling-noise floor sqrt((1 - m^2)/count), the noise-corrected signal rms, and its octave-to-octave ratio
(model: 1/2).  Also the prefix-form discrepancy of the Codex note, eps_M = [3S(2M) + 3S(M) + S(3M)]/(6M), at 64 points
per octave over 2^8..2^28, and its fitted power M^(-1+beta).
"""
import numpy as np, time
from numba import njit

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
def class_bins(tau, J, R, M):
    s = np.zeros((M, R)); c = np.zeros((M, R))
    for n in range(1 << J, 1 << (J + 1)):
        k = int((np.log2(n) - J) * R)
        if k >= R: k = R - 1
        r = n % M
        s[r, k] += -1.0 if (tau[n] & 1) else 1.0
        c[r, k] += 1
    return s, c

@njit(cache=True)
def prefix_at(tau, Ns):
    out = np.zeros(len(Ns)); s = 0; k = 0
    for n in range(1, Ns[-1] + 1):
        s += -1 if (tau[n] & 1) else 1
        while k < len(Ns) and n == Ns[k]:
            out[k] = s; k += 1
    return out

t0 = time.time()
tau = tau16(1 << 30)
print(f"tau to 2^30 in {time.time()-t0:.1f}s")

def ls_amp(x, y, f):
    if f == 0:
        return abs(y.mean())
    A = np.column_stack([np.cos(2 * np.pi * f * x), np.sin(2 * np.pi * f * x), np.ones_like(x)])
    cc, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(np.hypot(cc[0], cc[1]))

for M in (3, 9):
    print(f"\n## classes mod {M}: discrepancy D_r = class mean - global mean of (-1)^tau, per octave")
    print(" octave  R   count/bin   rms(D) all classes   noise floor   signal rms   ratio to previous octave | LS amplitude of D_2 (mod 3) or D_1 (mod 9) at f = 0, 1/2, 1, 3/2 (and 1/6, 1/3 for mod 9)")
    prev = None
    for J in range(8, 30):
        R = max(1, min(64, (1 << J) // (M * 48)))
        R = int(2 ** int(np.log2(R)))
        s, c = class_bins(tau, J, R, M)
        m_r = s / np.maximum(c, 1); glob = s.sum(axis=0) / c.sum(axis=0)
        D = m_r - glob[None, :]
        cnt = c.mean()
        noise2 = ((1 - glob ** 2)[None, :] / np.maximum(c, 1)).mean()
        rms2 = (D ** 2).mean()
        sig = np.sqrt(max(rms2 - noise2, 0.0))
        ratio = sig / prev if (prev and sig > 0) else float("nan")
        x = J + (np.arange(R) + 0.5) / R
        row = D[2] if M == 3 else D[1]
        fs = (0, 0.5, 1, 1.5) if M == 3 else (0, 1/6, 1/3, 0.5, 1, 1.5)
        amps = " ".join(f"{ls_amp(x, row, f):.4f}" for f in fs) if R >= 4 else "(R < 4)"
        print(f" 2^{J:2d}  {R:3d}  {cnt:9.0f}   {np.sqrt(rms2):.5f}            {np.sqrt(noise2):.5f}     {sig:.5f}      {ratio:.3f} | {amps}", flush=True)
        prev = sig if sig > 0 else prev

print("\n## prefix discrepancy eps_M = [3S(2M) + 3S(M) + S(3M)]/(6M) (Codex note, identity (7)), 64 points per octave, 2^8..2^28")
Ms = np.array(sorted(set(int(round(2.0 ** (8 + k / 64))) for k in range(20 * 64 + 1))), dtype=np.int64)
allN = np.array(sorted(set(np.concatenate([Ms, 2 * Ms, 3 * Ms]).tolist())), dtype=np.int64)
S = dict(zip(allN.tolist(), prefix_at(tau, allN).tolist()))
eps = np.array([(3 * S[2 * M] + 3 * S[M] + S[3 * M]) / (6.0 * M) for M in Ms.tolist()])
xM = np.log2(Ms.astype(float))
print(" log2 M : rms of eps over the octave, max |eps|, sign changes")
for J in range(8, 28):
    sel = (xM >= J) & (xM < J + 1)
    e = eps[sel]
    print(f" {J:2d}: rms {np.sqrt((e**2).mean()):.3e}  max {np.abs(e).max():.3e}  sign changes {int((np.diff(np.sign(e)) != 0).sum())}")
# power-law fit of the octave rms: rms ~ M^(-(1-beta))
rms = np.array([np.sqrt((eps[(xM >= J) & (xM < J + 1)] ** 2).mean()) for J in range(12, 28)])
slope = np.polyfit(np.arange(12, 28), np.log2(rms), 1)[0]
print(f" fitted log2 rms(eps) slope over octaves 12..27: {slope:.3f} per octave (model 1/M: -1.000; effective beta = {1 + slope:.3f})")
print("done")
