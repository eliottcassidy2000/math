"""
collatz_parity_classes_mod9_20261013.py -- data test of the 3-adic silence (THM-4629 statement 3): the parity means
restricted to the classes n = r mod 9 (and their aggregates mod 3) carry the same lines, with the same amplitudes, as the
global mean; the class-specific modes decay by 1/2 per octave and are invisible at 2^29.
Top octave [2^29, 2^30), 768 bins per class; least-squares amplitude at the predicted effective frequencies of the three
parities; report the global amplitude, the nine class amplitudes and the max relative spread.
"""
import numpy as np, time
from numba import njit

TH = np.log2(3.0)
LINES = {
 "tau":   [2.4549, 14.4689, 26.5036, 79.5105, 120.4762, 173.4809, 226.4865, 279.4929, 332.5001],
 "omega": [6.0203, 46.9850, 99.9912, 152.9982, 193.9697, 259.0121],
 "C":     [8.4606, 20.4822, 73.4880, 126.4946, 179.5019],
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
def bin_classes(tau, om, J, R, M):
    s = np.zeros((3, M, R)); c = np.zeros((M, R))
    for n in range(1 << J, 1 << (J + 1)):
        k = int((np.log2(n) - J) * R)
        if k >= R: k = R - 1
        r = n % M
        a = tau[n] & 1; b = om[n] & 1
        s[0, r, k] += -1.0 if a else 1.0
        s[1, r, k] += -1.0 if b else 1.0
        s[2, r, k] += -1.0 if (a ^ b) else 1.0
        c[r, k] += 1
    return s, c

t0 = time.time()
tau, om = tau_omega(1 << 30)
J, R, M = 29, 768, 9
s, c = bin_classes(tau, om, J, R, M)
print(f"tau, omega to 2^30 and class bins in {time.time()-t0:.1f}s; min count per class bin {int(c.min())}")
x = J + (np.arange(R) + 0.5) / R

def ls_amp(y, f):
    A = np.column_stack([np.cos(2 * np.pi * f * x), np.sin(2 * np.pi * f * x), np.ones_like(x)])
    cc, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(np.hypot(cc[0], cc[1])), float(np.arctan2(-cc[1], cc[0]))

for idx, name in enumerate(["tau", "omega", "C"]):
    glob = s[idx].sum(axis=0) / c.sum(axis=0)
    print(f"\n## {name}: line | global amplitude | class amplitudes r = 0..8 mod 9 | max |class - global| / global | max phase spread (rad)")
    for f in LINES[name]:
        ag, pg = ls_amp(glob, f)
        cls = [ls_amp(s[idx, r] / np.maximum(c[r], 1), f) for r in range(M)]
        amps = np.array([a for a, _ in cls]); phs = np.array([p for _, p in cls])
        dph = np.abs(np.angle(np.exp(1j * (phs - pg)))).max()
        print(f" {f:9.4f} | {ag:.4f} | " + " ".join(f"{a:.4f}" for a in amps) + f" | {np.abs(amps - ag).max() / ag:.3f} | {dph:.3f}")
    # mod-3 aggregates
    print("   mod 3 aggregates (r = 0, 1, 2): " + " ; ".join(
        f"{f:.1f}: " + " ".join(f"{ls_amp(s[idx, r::3].sum(axis=0) / np.maximum(c[r::3].sum(axis=0), 1), f)[0]:.4f}" for r in range(3))
        for f in LINES[name][:4]))
print("done")
