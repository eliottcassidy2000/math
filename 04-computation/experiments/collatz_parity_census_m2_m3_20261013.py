"""
collatz_parity_census_m2_m3_20261013.py -- audit repair (THM-4629 audit, findings 3, 6, 9): the COMPLETE near-return census
for the three parities (m = 2) and for the (Z/3)^2 characters (m = 3), tested line by line over ten octaves with the sign
of the frequency retained (positive-frequency coefficient only), and the C-parity scan near the Nyquist frequency 32 of the
old 64-bin analysis at 768 bins per octave.
Census: all (Q, P) with |Q theta - P| < 0.07 and Q <= 768 (m = 2, f = Q/2 <= 384) or Q <= 1152 (m = 3, f = Q/3 <= 384).
For each pair the law assigns the character (Q mod m, -P mod m); the complex bin means of all m^2 characters over
[2^20, 2^30) at 1/768 octave are fitted with e^{g (x-29)} e^{+-2 pi i f_eff x} and a constant (g, drift from the exact root
with theta_f = (Q theta - P)/m); report the amplitude in the predicted character and in the others.
"""
import numpy as np, time
from numba import njit
from mpmath import mp, mpc, exp, findroot, log, pi as mpi

mp.dps = 25
TH = np.log2(3.0); D = 2.0 - TH; L = log(3) / log(2) - 1
MOD = 6

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
def hist_octave(tau, om, J, R, MOD):
    H = np.zeros((R, MOD, MOD))
    for n in range(1 << J, 1 << (J + 1)):
        k = int((np.log2(n) - J) * R)
        if k >= R: k = R - 1
        H[k, tau[n] % MOD, om[n] % MOD] += 1.0
    return H

t0 = time.time()
tau, om = tau_omega(1 << 30)
R = 768
Hoct = {J: hist_octave(tau, om, J, R, MOD) for J in range(16, 30)}
print(f"tau, omega to 2^30 and (tau, omega) mod 6 histograms for octaves 16..29 in {time.time()-t0:.1f}s", flush=True)

def root_w(theta_f):
    w = findroot(lambda w: exp(-w) + exp(2j * mpi * theta_f) * exp(w * L) - 2, mpc(-0.001, 0))
    return float(w.real), float(w.imag) / (2 * np.pi)

def char_means(m, J0, J1):
    nb = (J1 - J0) * R
    Hm = np.zeros((nb, m, m)); cnt = np.zeros(nb)
    for J in range(J0, J1):
        H = Hoct[J]; sl = slice((J - J0) * R, (J - J0 + 1) * R)
        cnt[sl] = H.sum(axis=(1, 2))
        for t in range(MOD):
            for o in range(MOD):
                Hm[sl, t % m, o % m] += H[:, t, o]
    tt = np.arange(m); W = np.exp(2j * np.pi * np.outer(tt, tt) / m)
    z = np.einsum("at,bo,rto->abr", W, W, Hm) / cnt[None, None, :]
    x = J0 + (np.arange(nb) + 0.5) / R
    return z, x

names2 = {(1, 0): "tau", (0, 1): "omega", (1, 1): "C"}
for m, Qmax in ((2, 768), (3, 1152)):
    J0, J1 = 20, 30
    z, x = char_means(m, J0, J1); nb = len(x)
    zz = z.reshape(m * m, nb)
    print(f"\n## m = {m}: complete census |Q theta - P| < 0.07, Q <= {Qmax}; ten-octave fit 2^20..2^30; amplitude at 2^29")
    print("   (Q,P)   defect   char    f_eff    decay  | predicted  others (max, which)  verdict")
    carried = weak = fail = 0; weak_slow = []
    for Q in range(1, Qmax + 1):
        P = int(round(Q * TH)); d = Q * TH - P
        if abs(d) >= 0.07: continue
        a, b = Q % m, (-P) % m
        if (a, b) == (0, 0):
            continue        # the trivial character (class density, HYP-9230) is not a function of (tau, omega)
        g, drift = root_w(d / m)
        f = Q / m + drift
        env = np.exp(g * (x - 29)); E = env * np.exp(2j * np.pi * f * x)
        M = np.column_stack([E, np.conj(E), np.ones(nb)])
        coef, *_ = np.linalg.lstsq(M, zz.T, rcond=None)
        A = np.abs(coef[0]); A2 = A.reshape(m, m)
        pred = A2[a, b]
        others = [(A2[i, j], (i, j)) for i in range(m) for j in range(m) if (i, j) not in ((a, b), (0, 0))]
        mx, who = max(others)
        if pred < 0.01:
            tag = "WEAK"; weak += 1
            if np.exp(g) > 0.95: weak_slow.append((Q, P, f, np.exp(g), pred, mx))
        elif mx > 1.5 * pred:
            tag = "FAIL"; fail += 1
        else:
            tag = "CARRIED"; carried += 1
        cname = names2.get((a, b), "") if m == 2 else f"({a},{b})"
        print(f" ({Q:4d},{P:5d}) {d:+.4f} {cname:6s} {f:8.3f}  {np.exp(g):.4f} |  {pred:.4f}     {mx:.4f} {who}   {tag}", flush=True)
    print(f" m = {m}: CARRIED {carried}, WEAK {weak}, FAIL {fail}")
    if weak_slow:
        print(f" predicted slow lines (decay > 0.95) with amplitude below 0.01 at m = {m}:")
        for Q, P, f, dec, pred, mx in weak_slow:
            print(f"   ({Q},{P}) f = {f:.3f} decay {dec:.4f}: predicted-character amplitude {pred:.4f}, max other {mx:.4f}")

print("\n## C parity near the old Nyquist frequency 32, 768 bins per octave, least squares on single octaves")
for J in (16, 20, 29):
    z, x = char_means(2, J, J + 1)
    y = z[1, 1].real      # (-1)^(tau+omega)
    best = []
    for f in np.arange(28.0, 36.0, 0.05):
        A = np.column_stack([np.cos(2 * np.pi * f * x), np.sin(2 * np.pi * f * x), np.ones(len(x))])
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        best.append((float(np.hypot(c[0], c[1])), f))
    best.sort(reverse=True)
    print(f" 2^{J}: top amplitudes in [28, 36]: " + " ".join(f"({f:.2f},{a:.4f})" for a, f in best[:5]) + f"; at 32.00: {[a for a, f in best if abs(f-32.0)<0.026][0]:.4f}; at 32.523 (line (65,103)): {[a for a, f in best if abs(f-32.5)<0.026][0]:.4f}")
print("done")
