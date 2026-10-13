"""
collatz_parity_three_characters_20261013.py -- the three step-count parities of a Collatz orbit and their
log-periodic spectra, against the twisted transfer law.

Functionals of n (T(n) = n/2 or (3n+1)/2; tau = T-steps to 1; omega = odd T-steps to 1; tau + omega = steps of the
Collatz map C(n) = 3n+1 / n/2 to 1):
  t_tau   = (-1)^tau       transfer  psi(u) = -1/2 psi(u-1) - 1/2 psi(u+theta-1)   resonances: f half-integer, ||f theta|| small
  t_omega = (-1)^omega     transfer  psi(u) = +1/2 psi(u-1) - 1/2 psi(u+theta-1)   resonances: f integer,       ||f theta - 1/2|| small
  t_C     = (-1)^(tau+om)  transfer  psi(u) = -1/2 psi(u-1) + 1/2 psi(u+theta-1)   resonances: f half-integer, ||f theta - 1/2|| small
(theta = log2 3; u = log2 n).  Mode e^{(2 pi i f + w)u}: the characteristic equations
  tau:   e^{-s} + e^{s(theta-1)} = -2,   omega: e^{-s} - e^{s(theta-1)} = 2,   C: -e^{-s} + e^{s(theta-1)} = 2,  s = 2 pi i f + w.
Part 1: exact roots (mpmath continuation from the resonance point) at the predicted frequencies: drift Im w / 2 pi, decay Re w.
Part 2: per-octave least-squares spectra (1/R octave bins, 2^16..2^NL) of the three parities, peaks and candidate amplitudes.
Part 3: top-octave scan at 1/768 octave for the tau-parity up to f = 340 (the 665-mode at 332.5) and of omega up to 160.
Part 4: direct look at the 26.5 mode: bin means over two periods near 2^29.
Run: python3 collatz_parity_three_characters_20261013.py <NL> <R>
"""
import sys, time
import numpy as np
from numba import njit
from mpmath import mp, mpc, exp, findroot, log, pi as mpi

NL = int(sys.argv[1]) if len(sys.argv) > 1 else 30
R = int(sys.argv[2]) if len(sys.argv) > 2 else 64
THETA = np.log2(3.0)
mp.dps = 30
L = log(3) / log(2) - 1

def root_tau(f):
    # continuation in the phase mismatch from the exact resonance: e^{-s} + e^{s L} = -2, s = 2 pi i f + w
    q = round(2 * f)                      # odd
    def F(w, fr):
        s = 2j * mpi * fr + w
        return exp(-s) + exp(s * L) + 2
    # start from the resonant frequency where the phases align: solve for f0 with q theta nearest even integer
    w = mpc(0)
    f0 = q / 2
    for fr in np.linspace(f0, f0, 1):
        pass
    # continuation: deform from the exactly resonant fictitious theta to the true one is implicit; just solve directly from w=0
    return findroot(lambda w: F(w, f0), mpc(-0.01, 0))

def generic_root(eq):
    try:
        return findroot(eq, mpc(-0.01, 0.0))
    except Exception:
        return None

print('# Part 1: predicted modes and exact roots (Re w = decay per octave, Im w/2pi = frequency drift)')
rows = []
preds = {
 'tau':   [0.5, 2.5, 9.5, 14.5, 26.5, 332.5],
 'omega': [1, 6, 18, 47, 153],
 'C':     [0.5, 1.5, 3.5, 8.5, 20.5, 73.5, 179.5],
}
eqs = {
 'tau':   lambda s: exp(-s) + exp(s * L) + 2,
 'omega': lambda s: exp(-s) - exp(s * L) - 2,
 'C':     lambda s: -exp(-s) + exp(s * L) - 2,
}
pred_table = {}
for name, fs in preds.items():
    for f in fs:
        s0 = 2j * mpi * f
        w = generic_root(lambda w: eqs[name](s0 + w))
        if w is None:
            print(f'  {name:5s} f={f:7.2f}: no root found'); continue
        drift = float(w.imag) / (2 * np.pi); decay = float(w.real)
        pred_table[(name, f)] = (f + drift, decay)
        print(f'  {name:5s} f={f:7.2f}: effective f = {f + drift:9.4f}, decay {decay:+.5f}/octave (factor {np.exp(decay):.4f})')

@njit(cache=True)
def tau_omega(N):
    tau = np.zeros(N + 8, dtype=np.int16); om = np.zeros(N + 8, dtype=np.int16)
    for n in range(2, N + 8):
        m = n; steps = 0; odd = 0
        while m >= n:
            if m & 1:
                m = (3 * m + 1) >> 1; odd += 1
            else:
                m >>= 1
            steps += 1
        tau[n] = steps + tau[m]; om[n] = odd + om[m]
    return tau, om

@njit(cache=True)
def bin_means3(tau, om, J0, J1, R):
    nb = (J1 - J0) * R
    s = np.zeros((3, nb)); c = np.zeros(nb)
    for n in range(1 << J0, 1 << J1):
        k = int((np.log2(n) - J0) * R)
        if k >= nb: k = nb - 1
        a = tau[n] & 1; b = om[n] & 1
        s[0, k] += -1.0 if a else 1.0
        s[1, k] += -1.0 if b else 1.0
        s[2, k] += -1.0 if (a ^ b) else 1.0
        c[k] += 1
    return s / np.maximum(c, 1), c

t0 = time.time()
N = 1 << NL
tau, om = tau_omega(N)
print(f'\ntau, omega to 2^{NL} in {time.time()-t0:.1f}s', flush=True)
J0, J1 = 16, NL
ms, cnts = bin_means3(tau, om, J0, J1, R)
xs = J0 + (np.arange(ms.shape[1]) + 0.5) / R

def ls_amp(x, y, f):
    A = np.column_stack([np.cos(2 * np.pi * f * x), np.sin(2 * np.pi * f * x), np.ones_like(x)])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return np.hypot(c[0], c[1])

freqs = np.arange(0.5, R / 2, 0.05)
names = ['tau', 'omega', 'C']
print(f'\n# Part 2: per-octave spectra (1/{R} octave bins), top 5 peaks and amplitudes at the predicted effective frequencies')
for idx, name in enumerate(names):
    cand = [(f, pred_table[(name, f)][0]) for f in preds[name] if (name, f) in pred_table and pred_table[(name, f)][0] < R / 2]
    print(f'\n## {name}: candidates ' + ' '.join(f'{f}->{fe:.3f}' for f, fe in cand))
    for j in (16, 18, 20, 22, 24, 26, 28, J1 - 1):
        sel = (xs >= j) & (xs < j + 1); x = xs[sel]; y = ms[idx][sel]
        amps = np.array([ls_amp(x, y, f) for f in freqs])
        order = np.argsort(-amps); peaks = []
        for i in order:
            if all(abs(freqs[i] - g) > 0.4 for g, _ in peaks):
                peaks.append((freqs[i], amps[i]))
            if len(peaks) >= 5: break
        ca = [ls_amp(x, y, fe) for _, fe in cand]
        print(f' 2^{j:2d}: ' + ' '.join(f'({f:.2f},{a:.4f})' for f, a in peaks) + ' | pred ' + ' '.join(f'{a:.4f}' for a in ca), flush=True)

print('\n# Part 3: top-octave scan at 1/768 octave (frequencies to 340) for tau and omega parities')
@njit(cache=True)
def bin_means_top(tau, om, J, R):
    s = np.zeros((2, R)); c = np.zeros(R)
    for n in range(1 << J, 1 << (J + 1)):
        k = int((np.log2(n) - J) * R)
        if k >= R: k = R - 1
        s[0, k] += -1.0 if (tau[n] & 1) else 1.0
        s[1, k] += -1.0 if (om[n] & 1) else 1.0
        c[k] += 1
    return s / np.maximum(c, 1), c
RT = 768
mt, ct = bin_means_top(tau, om, J1 - 1, RT)
xt = (J1 - 1) + (np.arange(RT) + 0.5) / RT
hf = np.arange(0.5, 340, 0.25)
for idx, name in enumerate(['tau', 'omega']):
    amps = np.array([ls_amp(xt, mt[idx], f) for f in hf])
    order = np.argsort(-amps); peaks = []
    for i in order:
        if all(abs(hf[i] - g) > 1.0 for g, _ in peaks):
            peaks.append((hf[i], amps[i]))
        if len(peaks) >= 10: break
    print(f' {name} top octave 2^{J1-1}: ' + ' '.join(f'({f:.2f},{a:.4f})' for f, a in peaks))
    for f in ([26.5036, 332.5, 332.5 + 0.0001] if name == 'tau' else [6.0, 47.0, 153.0]):
        print(f'    amplitude at f={f:.4f}: {ls_amp(xt, mt[idx], f):.4f}')
    # fine scan around 332.5 / 153
    centre = 332.5 if name == 'tau' else 153.0
    fine = np.arange(centre - 1.0, centre + 1.0, 0.01)
    fa = np.array([ls_amp(xt, mt[idx], f) for f in fine])
    print(f'    fine scan near {centre}: best f = {fine[int(np.argmax(fa))]:.3f} amplitude {fa.max():.4f}')

print('\n# Part 4: the 26.5-mode seen directly: tau-parity means in 1/212-octave bins over two periods near 2^29 (period 2/53 octave)')
@njit(cache=True)
def bin_means_range(tau, lo, hi, nb):
    s = np.zeros(nb); c = np.zeros(nb)
    xl = np.log2(lo); xh = np.log2(hi)
    for n in range(lo, hi):
        k = int((np.log2(n) - xl) / (xh - xl) * nb)
        if k >= nb: k = nb - 1
        s[k] += -1.0 if (tau[n] & 1) else 1.0
        c[k] += 1
    return s / np.maximum(c, 1)
lo = 1 << 29; hi = int(lo * 2 ** (4 / 53))
mr = bin_means_range(tau, lo, hi, 16)
print('   ' + ' '.join(f'{v:+.3f}' for v in mr))
np.save('collatz_parity_three_characters_20261013_bins.npy', np.vstack([xs, ms, cnts]))
print('done')
