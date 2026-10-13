"""
collatz_parity_line_assignment_20261013.py -- direct test of the character law (THM-4629, statement 3) for m = 4, 5, 6, 8, 10, 12:
for every near-return 2^P ~ 3^Q (|Q theta - P| < 0.07, Q <= 2300 so that f = Q/m stays below the Nyquist 384 at m >= 6;
for m = 4, 5 only Q <= 1536, 1920) the law predicts that exactly the character (alpha, beta) = (Q mod m, -P mod m) carries a
line at f_eff = Q/m + theta_f/(2 - theta), theta_f = (Q theta - P)/m, and that no other character does.
Top octave [2^29, 2^30) at 1/768 octave; complex bin means of all m^2 characters from the histogram of (tau mod 120, omega mod 120);
complex least squares at f_eff: A e^{2 pi i f x} + B e^{-2 pi i f x} + C.  Report per line: |A| in the predicted character,
the median and the maximum of |A| over the other m^2 - 1 characters (and which one attains the maximum), and the same for the
conjugate test (|B| in the character (-alpha, -beta)).  Verdict per line: CARRIED if the predicted amplitude exceeds 0.01 and
is at least 3 x the median of the others; WEAK if below 0.01 (initial amplitude too small to test); FAIL if another character
exceeds the predicted one by more than 1.5 x (a neighbouring line of another character within 1 cycle is noted).
"""
import numpy as np, time
from numba import njit

TH = np.log2(3.0); D = 2.0 - TH
MOD = 120

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
def hist_top(tau, om, J, R, MOD):
    H = np.zeros((R, MOD, MOD))
    for n in range(1 << J, 1 << (J + 1)):
        k = int((np.log2(n) - J) * R)
        if k >= R: k = R - 1
        H[k, tau[n] % MOD, om[n] % MOD] += 1.0
    return H

t0 = time.time()
tau, om = tau_omega(1 << 30)
J, R = 29, 768
H = hist_top(tau, om, J, R, MOD)
cnt = H.sum(axis=(1, 2))
x = J + (np.arange(R) + 0.5) / R
print(f"tau, omega to 2^30 and the (tau, omega) mod {MOD} histogram of the top octave in {time.time()-t0:.1f}s")

# near-returns
NR = []
for Q in range(1, 2301):
    P = int(round(Q * TH)); d = Q * TH - P
    if abs(d) < 0.07:
        NR.append((Q, P, d))
print(f"{len(NR)} near-returns with |Q theta - P| < 0.07 and Q <= 2300")

def char_means_all(m):
    # reduce the histogram to (tau mod m, omega mod m), then 2D DFT over the two residues: z[alpha, beta, bin]
    Hm = np.zeros((R, m, m))
    for t in range(MOD):
        for o in range(MOD):
            Hm[:, t % m, o % m] += H[:, t, o]
    tt = np.arange(m)
    W = np.exp(2j * np.pi * np.outer(tt, tt) / m)          # W[alpha, t]
    z = np.einsum("at,bo,rto->abr", W, W, Hm) / cnt[None, None, :]
    return z

def fitA(z, f):
    E = np.exp(2j * np.pi * f * x)
    M = np.column_stack([E, np.conj(E), np.ones(R)])
    coef, *_ = np.linalg.lstsq(M, z.T, rcond=None)   # z: (K, R) -> coef: (3, K)
    return np.abs(coef[0]), np.abs(coef[1])

results = {}
for m in (4, 5, 6, 8, 10, 12):
    z = char_means_all(m)                            # (m, m, R)
    zz = z.reshape(m * m, R)
    Qmax = {4: 1536, 5: 1920}.get(m, 2300)
    carried = weak = fail = 0; rows = []
    for Q, P, d in NR:
        if Q > Qmax: continue
        a, b = Q % m, (-P) % m
        if (a, b) == (0, 0): continue
        f = Q / m + (d / m) / D
        A, B = fitA(zz, f)
        A2 = A.reshape(m, m); B2 = B.reshape(m, m)
        pred = A2[a, b]
        others = np.delete(A, a * m + b)
        med = float(np.median(others)); mx = float(others.max()); imx = int(np.argmax(others))
        ia, ib = divmod(int(np.delete(np.arange(m * m), a * m + b)[imx]), m)
        predB = B2[(-a) % m, (-b) % m]
        if pred < 0.01:
            tag = "WEAK"; weak += 1
        elif mx > 1.5 * pred:
            tag = "FAIL"; fail += 1
        elif pred >= 3 * med:
            tag = "CARRIED"; carried += 1
        else:
            tag = "FAIL"; fail += 1
        # neighbouring line of another character within 1 cycle?
        nb = [(Q2, P2) for Q2, P2, d2 in NR if Q2 != Q and abs(Q2 / m - Q / m) < 1.0]
        rows.append((Q, P, d, a, b, f, pred, predB, med, mx, ia, ib, tag, nb))
    results[m] = (carried, weak, fail)
    print(f"\n## m = {m}: lines CARRIED by the predicted character {carried}, WEAK (predicted amplitude < 0.01) {weak}, FAIL {fail}")
    print("   (Q,P)   defect  char   f_eff    |A| predicted  |B| conj   median others  max others (char)  verdict  [neighbouring lines of other characters within 1 cycle]")
    for Q, P, d, a, b, f, pred, predB, med, mx, ia, ib, tag, nb in rows:
        print(f" ({Q:4d},{P:4d}) {d:+.4f} ({a},{b}) {f:8.3f}   {pred:.4f}       {predB:.4f}    {med:.4f}        {mx:.4f} ({ia},{ib})   {tag:7s} {nb if nb else ''}")
print("\nsummary m: (carried, weak, fail):", results)
print("done")
