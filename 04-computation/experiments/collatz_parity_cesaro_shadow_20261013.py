"""
collatz_parity_cesaro_shadow_20261013.py -- the ordinary (Cesaro) mean of the parity of tau as the shadow of the
log-periodic lines.  A(N) = (1/N) sum_(n<=N) (-1)^tau(n).  For a line e^(g(x-29)) [c cos(2 pi f x) + s sin(2 pi f x)]
in x = log2 n, the prefix mean is Re[(c - i s) (N^(1+g/ln2+i om) - 1) / ((1 + g/ln2 + i om) N)] e^(-29 g), om = 2 pi f/ln 2:
the Cesaro mean of a log-periodic oscillation is the oscillation divided by sqrt(1 + om^2) (phase lag arctan om).
Fit (c, s) of every predicted tau-parity line jointly on the top octave [2^29, 2^30) at 1/768 octave, with the
predicted effective frequency and decay (THM-4629 roots); compare the exact A(N) at 256 points per octave over
2^24..2^30 with the model; predict A(2^34) (Codex census, collatz_signed_resonance_pair_mass_20261013b.md: 0.001452).
"""
import numpy as np, time
from numba import njit
from mpmath import mp, mpc, exp, findroot, log, pi as mpi

mp.dps = 30
L = log(3) / log(2) - 1
LINES = [2.5, 14.5, 26.5, 79.5, 120.5, 173.5, 226.5, 279.5, 332.5]
LN2 = np.log(2.0)

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
def prefix_means(tau, Ns):
    A = np.zeros(len(Ns)); s = 0; k = 0
    for n in range(1, Ns[-1] + 1):
        s += -1 if (tau[n] & 1) else 1
        while k < len(Ns) and n == Ns[k]:
            A[k] = s / n; k += 1
    return A

@njit(cache=True)
def bin_top(tau, J, R):
    s = np.zeros(R); c = np.zeros(R)
    for n in range(1 << J, 1 << (J + 1)):
        k = int((np.log2(n) - J) * R)
        if k >= R: k = R - 1
        s[k] += -1.0 if (tau[n] & 1) else 1.0
        c[k] += 1
    return s / np.maximum(c, 1)

t0 = time.time()
tau = tau16(1 << 30)
print(f"tau to 2^30 in {time.time()-t0:.1f}s")

# roots
roots = []
for f in LINES:
    s0 = 2j * mpi * f
    w = findroot(lambda w: exp(-(s0 + w)) + exp((s0 + w) * L) + 2, mpc(-0.01, 0))
    roots.append((f, f + float(w.imag) / (2 * np.pi), float(w.real)))

# joint fit on the top octave
J, R = 29, 768
y = bin_top(tau, J, R)
x = J + (np.arange(R) + 0.5) / R
cols = []
for f, fe, g in roots:
    env = np.exp(g * (x - 29))
    cols.append(env * np.cos(2 * np.pi * fe * x)); cols.append(env * np.sin(2 * np.pi * fe * x))
cols.append(np.ones_like(x))
M = np.column_stack(cols)
coef, *_ = np.linalg.lstsq(M, y, rcond=None)
resid = y - M @ coef
print(f"joint fit of {len(roots)} lines + constant on [2^29, 2^30): rms(y) = {np.sqrt((y**2).mean()):.4f}, rms(residual) = {np.sqrt((resid**2).mean()):.4f}, constant = {coef[-1]:+.5f}")
print("line   eff f      decay/oct   amplitude   phase     Cesaro factor 1/sqrt(1+om^2)   Cesaro amplitude at 2^29")
for i, (f, fe, g) in enumerate(roots):
    c, s = coef[2 * i], coef[2 * i + 1]
    a = np.hypot(c, s); ph = np.arctan2(-s, c)
    om = 2 * np.pi * fe / LN2
    print(f"{f:6.1f} {fe:9.4f}   {np.exp(g):.4f}     {a:.4f}   {ph:+.3f}     {1/np.sqrt(1+om*om):.6f}                  {a/np.sqrt(1+om*om):.6f}")

def model_A(N):
    # Cesaro mean of the fitted modal sum, exact integral from 1 to N
    tot = coef[-1]
    for i, (f, fe, g) in enumerate(roots):
        c, s = coef[2 * i], coef[2 * i + 1]
        om = 2 * np.pi * fe / LN2
        p = 1 + g / LN2 + 1j * om
        tot += np.exp(-29 * g) * ((c - 1j * s) * (N ** p - 1) / (p * N)).real
    return tot

# exact prefix means at 256 points per octave over 2^24..2^30
J0, J1, RR = 24, 30, 256
Ns = np.array([int(round(2.0 ** (J0 + k / RR))) for k in range((J1 - J0) * RR + 1)], dtype=np.int64)
A = prefix_means(tau, Ns)
Am = np.array([model_A(float(N)) for N in Ns])
xs = np.log2(Ns.astype(np.float64))
print(f"\nexact prefix means over 2^24..2^30 ({len(Ns)} points): rms {np.sqrt((A**2).mean()):.6f}, max |A| {np.abs(A).max():.6f}; model rms {np.sqrt((Am**2).mean()):.6f}; rms(exact - model) {np.sqrt(((A-Am)**2).mean()):.6f}")
# oscillation frequency of the exact prefix means on the top octave
def ls_amp(x, y, f):
    B = np.column_stack([np.cos(2 * np.pi * f * x), np.sin(2 * np.pi * f * x), np.ones_like(x)])
    cc, *_ = np.linalg.lstsq(B, y, rcond=None)
    return float(np.hypot(cc[0], cc[1]))
sel = xs >= 29
fs = np.arange(0.5, 128, 0.25)
amps = np.array([ls_amp(xs[sel], A[sel], f) for f in fs])
top = np.argsort(-amps)[:6]
print("top-octave spectrum of the exact prefix mean A(N): " + " ".join(f"({fs[i]:.2f},{amps[i]:.6f})" for i in top))
print("\n N = 2^k:   exact A(N)      model A(N)")
for k in range(J0, J1 + 1):
    i = (k - J0) * RR
    print(f" 2^{k}:   {A[i]:+.6f}     {Am[i]:+.6f}")
print("\nprediction beyond the data (extrapolated decays): N = 2^31..2^34 model A(N) and envelope")
for k in (31, 32, 33, 34):
    N = 2.0 ** k
    env = sum(np.hypot(coef[2*i], coef[2*i+1]) * np.exp(g * (k - 29)) / np.sqrt(1 + (2*np.pi*fe/LN2)**2) for i, (f, fe, g) in enumerate(roots))
    print(f" 2^{k}: model {model_A(N):+.6f}, envelope {env:.6f}")
# fine scan of the model near 2^34 and the Codex value
Nf = 2.0 ** np.linspace(33.9, 34.1, 2001)
vals = np.array([model_A(N) for N in Nf])
print(f"model over 2^33.9..2^34.1: min {vals.min():+.6f} max {vals.max():+.6f}; Codex exact census at 2^34: +0.001452")
np.save("collatz_parity_cesaro_shadow_20261013_prefix.npy", np.vstack([Ns.astype(np.float64), A, Am]))
print("done")
