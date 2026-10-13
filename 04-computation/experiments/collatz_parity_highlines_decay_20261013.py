"""
collatz_parity_highlines_decay_20261013.py -- decay per octave of the high resonance lines of the three parities of
(tau, omega) against the exact roots of the twisted transfer law (THM-4629).  Bins of 1/768 octave over octaves 20..29
(tau, omega exact to 2^30); least-squares amplitude at each predicted effective frequency per octave; fitted decay =
exp(slope of log amplitude vs octave) over the octaves with amplitude > 0.003; predicted decay = exp(Re w).
Lines: a SELECTION of near-returns (those visible in the top-octave scans), not the complete census; the census with
|Q log2 3 - P| < 0.07, Q <= 768 is tested line by line in collatz_parity_census_m2_m3_20261013.py (audit finding 3).
"""
import numpy as np, time
from numba import njit
from mpmath import mp, mpc, exp, findroot, log, pi as mpi

mp.dps = 30
L = log(3) / log(2) - 1
EQS = {
 "tau":   lambda s: exp(-s) + exp(s * L) + 2,
 "omega": lambda s: exp(-s) - exp(s * L) - 2,
 "C":     lambda s: -exp(-s) + exp(s * L) - 2,
}
LINES = {
 "tau":   [2.5, 14.5, 26.5, 79.5, 120.5, 173.5, 226.5, 279.5, 332.5],
 "omega": [6, 18, 35, 47, 59, 100, 112, 153, 194, 206, 247, 259, 365],
 "C":     [3.5, 8.5, 20.5, 32.5, 73.5, 126.5, 179.5],
}

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
def bin3(tau, om, J0, J1, R):
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
    return s / np.maximum(c, 1)

t0 = time.time()
tau, om = tau_omega(1 << 30)
J0, J1, R = 20, 30, 768
ms = bin3(tau, om, J0, J1, R)
xs = J0 + (np.arange(ms.shape[1]) + 0.5) / R
print(f"tau, omega to 2^30 and 1/{R}-octave bins over 2^{J0}..2^{J1} in {time.time()-t0:.1f}s")

def ls_amp(x, y, f):
    A = np.column_stack([np.cos(2 * np.pi * f * x), np.sin(2 * np.pi * f * x), np.ones_like(x)])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(np.hypot(c[0], c[1]))

TH = np.log2(3.0)
print("parity  f      (Q,P)      defect   eff f      pred decay  fitted  | amplitudes per octave 2^20 .. 2^29")
for idx, name in enumerate(["tau", "omega", "C"]):
    for f in LINES[name]:
        Q = int(round(2 * f)); P = int(round(Q * TH)); de = Q * TH - P
        s0 = 2j * mpi * f
        try:
            w = findroot(lambda w: EQS[name](s0 + w), mpc(-0.01, 0))
        except Exception:
            print(f"{name:5s} f={f:6.1f}: no root"); continue
        fe = f + float(w.imag) / (2 * np.pi); dec = float(np.exp(float(w.real)))
        amps = np.array([ls_amp(xs[(xs >= j) & (xs < j + 1)], ms[idx][(xs >= j) & (xs < j + 1)], fe) for j in range(J0, J1)])
        ok = amps > 0.003
        if ok.sum() >= 3:
            cfit = np.polyfit(np.arange(J0, J1)[ok], np.log(amps[ok]), 1); fit = float(np.exp(cfit[0]))
        else:
            fit = float("nan")
        print(f"{name:5s} {f:6.1f} ({Q:4d},{P:5d}) {de:+.4f} {fe:9.4f}   {dec:.4f}     {fit:.4f}  | " + " ".join(f"{a:.4f}" for a in amps), flush=True)
print("done")
