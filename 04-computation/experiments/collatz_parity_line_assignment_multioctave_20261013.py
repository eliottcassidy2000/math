"""
collatz_parity_line_assignment_multioctave_20261013.py -- the line-by-line character test over ten octaves.
Same rule as collatz_parity_line_assignment_20261013.py (THM-4629, statement 3): the near-return (Q, P) is carried by the
character (Q mod m, -P mod m) at f_eff = Q/m + theta_f/(2 - theta).  Here the complex bin means of every character are
fitted over [2^20, 2^30) at 1/768 octave (7680 points, resolution 0.1 cycle per octave) with the basis
e^{g (x - 29)} e^{+-2 pi i f_eff x} (g = Re w of the exact root, continued from 0) and a constant; the amplitude is reported at
2^29.  Verdicts as before: CARRIED (predicted amplitude > 0.01 and >= 3 x the median over the other characters),
WEAK (< 0.01), FAIL (another character exceeds the predicted one by more than 1.5 x).  Lines of different characters within
0.15 cycle per octave of each other are still noted.
"""
import numpy as np, time
from numba import njit
from mpmath import mp, mpc, exp, findroot, log, pi as mpi

mp.dps = 25
TH = np.log2(3.0); D = 2.0 - TH; L = log(3) / log(2) - 1
MOD = 120
MS = (4, 5, 6, 8, 10, 12)

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
J0, J1, R = 20, 30, 768
nb = (J1 - J0) * R
x = J0 + (np.arange(nb) + 0.5) / R
# reduced histograms per m: Hm[m] of shape (nb, m, m)
Hm = {m: np.zeros((nb, m, m)) for m in MS}
cnt = np.zeros(nb)
for J in range(J0, J1):
    H = hist_octave(tau, om, J, R, MOD)
    sl = slice((J - J0) * R, (J - J0 + 1) * R)
    cnt[sl] = H.sum(axis=(1, 2))
    for m in MS:
        for t in range(MOD):
            for o in range(MOD):
                Hm[m][sl, t % m, o % m] += H[:, t, o]
    del H
print(f"tau, omega to 2^30 and reduced histograms over 2^{J0}..2^{J1} in {time.time()-t0:.1f}s", flush=True)

NR = []
for Q in range(1, 2301):
    P = int(round(Q * TH)); d = Q * TH - P
    if abs(d) < 0.07:
        NR.append((Q, P, d))

# exact roots for the tau-parity-type equation are not what we need: for a character (a, b) mod m the root equation is
# e^{-w} + e^{2 pi i theta_f} e^{w(theta-1)} = 2 with theta_f = d/m (statement 3): compute Re w once per (Q, m)
def root_w(theta_f):
    try:
        w = findroot(lambda w: exp(-w) + exp(2j * mpi * theta_f) * exp(w * L) - 2, mpc(-0.001, 0))
        return float(w.real), float(w.imag) / (2 * np.pi)
    except Exception:
        return None

summary = {}
for m in MS:
    tt = np.arange(m); W = np.exp(2j * np.pi * np.outer(tt, tt) / m)
    z = np.einsum("at,bo,rto->abr", W, W, Hm[m]) / cnt[None, None, :]
    zz = z.reshape(m * m, nb)
    Qmax = {4: 1536, 5: 1920}.get(m, 2300)
    carried = weak = fail = 0; rows = []
    for Q, P, d in NR:
        if Q > Qmax: continue
        a, b = Q % m, (-P) % m
        if (a, b) == (0, 0): continue
        rw = root_w(d / m)
        if rw is None: continue
        g, drift = rw
        f = Q / m + drift
        env = np.exp(g * (x - 29))
        E = env * np.exp(2j * np.pi * f * x)
        M = np.column_stack([E, np.conj(E), np.ones(nb)])
        coef, *_ = np.linalg.lstsq(M, zz.T, rcond=None)
        A = np.abs(coef[0])
        pred = A[a * m + b]
        others = np.delete(A, a * m + b)
        med = float(np.median(others)); mx = float(others.max())
        ia, ib = divmod(int(np.delete(np.arange(m * m), a * m + b)[int(np.argmax(others))]), m)
        if pred < 0.01:
            tag = "WEAK"; weak += 1
        elif mx > 1.5 * pred:
            tag = "FAIL"; fail += 1
        elif pred >= 3 * med:
            tag = "CARRIED"; carried += 1
        else:
            tag = "FAIL"; fail += 1
        nbrs = [(Q2, P2) for Q2, P2, d2 in NR if Q2 != Q and abs(Q2 / m - Q / m) < 0.15]
        rows.append((Q, P, d, a, b, f, np.exp(g), pred, med, mx, ia, ib, tag, nbrs))
    summary[m] = (carried, weak, fail)
    print(f"\n## m = {m}: CARRIED {carried}, WEAK {weak}, FAIL {fail}")
    print("   (Q,P)   defect  char   f_eff    decay   |A| pred   median others  max others (char)  verdict  [lines of other characters within 0.15]")
    for Q, P, d, a, b, f, dec, pred, med, mx, ia, ib, tag, nbrs in rows:
        print(f" ({Q:4d},{P:4d}) {d:+.4f} ({a:2d},{b:2d}) {f:8.3f}  {dec:.4f}   {pred:.4f}    {med:.4f}        {mx:.4f} ({ia},{ib})   {tag:7s} {nbrs if nbrs else ''}", flush=True)
print("\nsummary m: (carried, weak, fail):", summary)
print("done")
