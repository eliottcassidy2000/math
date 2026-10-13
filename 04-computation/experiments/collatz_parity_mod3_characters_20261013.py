"""
collatz_parity_mod3_characters_20261013.py -- the general character law tested at m = 3.
Character chi_{alpha,beta}(n) = exp(2 pi i (alpha tau(n) + beta omega(n))/3). Predicted modes: f = Q/3 for the near-returns
2^P ~ 3^Q with Q = alpha, P = -beta (mod 3). From the list (Q,P): (53,84): Q=2,P=0 -> chi_{2,0} at f=53/3=17.667;
(147,233): Q=0,P=2 -> chi_{0,1} (beta=1: P=-1=2 ok) at f=49; (306,485): Q=0,P=2 -> chi_{0,1} at f=102;
(559,886): Q=1,P=1 -> chi_{1,2} (beta=2: P=-2=1 ok) at f=186.33; (241,382): Q=1,P=1 -> chi_{1,2} at f=80.33;
(665,1054): Q=2,P=1 -> chi_{2,2} at f=221.67; (359,569): Q=2,P=2 -> chi_{2,1} at f=119.67; (200,317): Q=2,P=2 -> chi_{2,1} at f=66.67;
(41,65): Q=2,P=2 -> chi_{2,1} at f=13.667; (5,8): Q=2,P=2 -> chi_{2,1} at f=1.667; (12,19): Q=0,P=1 -> chi_{0,2} at f=4;
(94,149): Q=1,P=2 -> chi_{1,1} at f=31.33; (29,46): Q=2,P=1 -> chi_{2,2} at f=9.667; (7,11): Q=1,P=2 -> chi_{1,1} at f=2.333;
(17,27): Q=2,P=0 -> chi_{2,0} at f=5.667.
Top octave [2^29, 2^30), 768 bins; least-squares amplitude of the complex bin means at the predicted f and the scan peaks.
"""
import numpy as np, time
from numba import njit

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
def bin_chars(tau, om, J, R):
    # complex means of chi_{a,b} for all 9 characters of (Z/3)^2
    re = np.zeros((9, R)); im = np.zeros((9, R)); c = np.zeros(R)
    for n in range(1 << J, 1 << (J + 1)):
        k = int((np.log2(n) - J) * R)
        if k >= R: k = R - 1
        t = tau[n] % 3; o = om[n] % 3
        for a in range(3):
            for b in range(3):
                ph = 2 * np.pi * ((a * t + b * o) % 3) / 3
                re[3 * a + b, k] += np.cos(ph); im[3 * a + b, k] += np.sin(ph)
        c[k] += 1
    return re / np.maximum(c, 1), im / np.maximum(c, 1)

t0 = time.time()
tau, om = tau_omega(1 << 30)
print(f'tau, omega to 2^30 in {time.time()-t0:.1f}s')
J, R = 29, 768
re, im = bin_chars(tau, om, J, R)
x = J + (np.arange(R) + 0.5) / R

def amp(z, f):
    # complex least squares: z(x) ~ A e^{2 pi i f x} + B e^{-2 pi i f x} + C ; report |A| and |B|
    E = np.exp(2j * np.pi * f * x)
    M = np.column_stack([E, np.conj(E), np.ones_like(E)])
    coef, *_ = np.linalg.lstsq(M, z, rcond=None)
    return abs(coef[0]), abs(coef[1])

pred = {(2, 0): [17.667, 5.667], (0, 1): [49.0, 102.0], (1, 2): [186.333, 80.333], (2, 2): [221.667, 9.667],
        (2, 1): [119.667, 66.667, 13.667, 1.667], (0, 2): [4.0], (1, 1): [31.333, 2.333]}
fs = np.arange(0.5, 380, 0.25)
for a in range(3):
    for b in range(3):
        if (a, b) == (0, 0):
            continue
        z = re[3 * a + b] + 1j * im[3 * a + b]
        amps = np.array([max(amp(z, f)) for f in fs])
        order = np.argsort(-amps); peaks = []
        for i in order:
            if all(abs(fs[i] - g) > 1.0 for g, _ in peaks):
                peaks.append((fs[i], amps[i]))
            if len(peaks) >= 8: break
        pr = pred.get((a, b), [])
        print(f'chi_({a},{b}): peaks ' + ' '.join(f'({f:.2f},{v:.4f})' for f, v in peaks))
        if pr:
            print(f'          predicted ' + ' '.join(f'{f:.3f}:{max(amp(z, f)):.4f}' for f in pr))
print('done')
