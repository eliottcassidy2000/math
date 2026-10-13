"""
collatz_parity_characters_modM_20261013.py -- the general character law tested where no temperament reading exists:
all characters chi_(alpha,beta)(n) = exp(2 pi i (alpha tau(n) + beta omega(n))/m) for m = 4, 5, 6 (and m = 2, 3 as controls),
top octave [2^29, 2^30) at 1/768 octave.  The complex bin means z(x) are fitted at every f on a 0.25-grid as
A e^{2 pi i f x} + B e^{-2 pi i f x} + C; the POSITIVE-frequency coefficient A tests the rule for (alpha, beta):
f = Q/m with Q = alpha, P = round(Q theta) = -beta (mod m) and small defect; the NEGATIVE-frequency coefficient B tests
the conjugate character (-alpha, -beta).  Peaks (local maxima, exclusion radius 1.0, amplitude > 0.015) are classified:
HIT (predicted near-return of the right class, |defect| < 0.07), SIDELOBE (within 1.3..1.7 of a larger peak of the same
list), or MISS.  Also: the amplitude and phase of the Mercator line 26.5036 of (-1)^tau octave by octave from 2^7.
"""
import numpy as np, time
from numba import njit

TH = np.log2(3.0)
MOD = 60                      # lcm of 3, 4, 5 (and 2, 6): histogram of (tau mod 60, omega mod 60) per bin

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

@njit(cache=True)
def tau_bins(tau, J, R):
    s = np.zeros(R); c = np.zeros(R)
    for n in range(1 << J, 1 << (J + 1)):
        k = int((np.log2(n) - J) * R)
        if k >= R: k = R - 1
        s[k] += -1.0 if (tau[n] & 1) else 1.0
        c[k] += 1
    return s / np.maximum(c, 1)

t0 = time.time()
tau, om = tau_omega(1 << 30)
J, R = 29, 768
H = hist_top(tau, om, J, R, MOD)
cnt = H.sum(axis=(1, 2))
x = J + (np.arange(R) + 0.5) / R
print(f"tau, omega to 2^30 and the (tau, omega) mod {MOD} histogram of the top octave in {time.time()-t0:.1f}s")

tt = np.arange(MOD)[:, None]; oo = np.arange(MOD)[None, :]
def char_means(m, a, b):
    ph = np.exp(2j * np.pi * ((a * tt + b * oo) % m) / m)        # shape MOD x MOD
    return (H * ph[None, :, :]).sum(axis=(1, 2)) / cnt

fs = np.arange(0.5, 380.0, 0.25)
E = np.exp(2j * np.pi * np.outer(fs, x))                           # len(fs) x R
ones = np.ones(R)
def fit_all(z):
    # for every f: least squares z ~ A E_f + B conj(E_f) + C ; solve 3x3 normal equations per f (vectorised)
    Ef = E; Eb = np.conj(E)
    # Gram matrices
    g11 = (np.conj(Ef) * Ef).sum(axis=1); g12 = (np.conj(Ef) * Eb).sum(axis=1); g13 = np.conj(Ef).sum(axis=1)
    g22 = (np.conj(Eb) * Eb).sum(axis=1); g23 = np.conj(Eb).sum(axis=1); g33 = np.full(len(fs), float(R))
    r1 = (np.conj(Ef) * z[None, :]).sum(axis=1); r2 = (np.conj(Eb) * z[None, :]).sum(axis=1); r3 = np.full(len(fs), z.sum())
    G = np.stack([np.stack([g11, g12, g13], -1), np.stack([np.conj(g12), g22, g23], -1), np.stack([np.conj(g13), np.conj(g23), g33], -1)], -2)
    rhs = np.stack([r1, r2, r3], -1)
    sol = np.linalg.solve(G, rhs[..., None])[..., 0]
    return np.abs(sol[:, 0]), np.abs(sol[:, 1])

def peaks_of(amp, thr=0.015, rad=1.0, nmax=12):
    order = np.argsort(-amp); out = []
    for i in order:
        if amp[i] < thr: break
        if all(abs(fs[i] - g) > rad for g, _ in out):
            out.append((fs[i], amp[i]))
        if len(out) >= nmax: break
    return out

def classify(peaks, m, a, b):
    res = []
    for f, v in peaks:
        Q = int(round(m * f)); P = int(round(Q * TH)); d = Q * TH - P
        hit = (abs(d) < 0.07) and (Q % m == a % m) and (P % m == (-b) % m) and abs(m * f - Q) < 0.4
        side = any((1.3 <= abs(f - g) <= 1.7) and w > v for g, w in peaks)
        res.append((f, v, Q, P, d, "HIT" if hit else ("SIDE" if side else "MISS")))
    return res

summary = {}
for m in (2, 3, 4, 5, 6):
    hits = sides = misses = 0; miss_list = []
    print(f"\n## m = {m}: characters (alpha, beta), positive-frequency peaks [A] test (alpha, beta); negative-frequency peaks [B] test (-alpha, -beta)")
    for a in range(m):
        for b in range(m):
            if (a, b) == (0, 0):
                continue
            z = char_means(m, a, b)
            A, B = fit_all(z)
            for lab, amp, (aa, bb) in (("A", A, (a, b)), ("B", B, ((-a) % m, (-b) % m))):
                if m == 2 and lab == "B":
                    continue
                res = classify(peaks_of(amp), m, aa, bb)
                for f, v, Q, P, d, tag in res:
                    if tag == "HIT": hits += 1
                    elif tag == "SIDE": sides += 1
                    else:
                        misses += 1; miss_list.append((m, a, b, lab, f, v, Q, P, d))
                txt = " ".join(f"{f:.2f}:{v:.3f}{'*' if tag == 'HIT' else ('~' if tag == 'SIDE' else '?')}" for f, v, Q, P, d, tag in res)
                print(f" chi({a},{b}) [{lab} -> ({aa},{bb})]: {txt}")
    summary[m] = (hits, sides, misses)
    print(f" m = {m}: HIT {hits}, SIDELOBE {sides}, MISS {misses}")
    for mm, a, b, lab, f, v, Q, P, d in miss_list:
        print(f"   MISS m={mm} chi({a},{b})[{lab}] f={f:.2f} amp={v:.3f} nearest (Q,P)=({Q},{P}) defect {d:+.4f} Q mod m={Q % mm} P mod m={P % mm}")
print("\nsummary (m: hits, sidelobes, misses):", summary)

print("\n## Mercator line of (-1)^tau octave by octave (least squares at f = 26.5036; R bins per octave)")
for j in range(7, 30):
    Rj = 64 if j < 14 else 768
    y = tau_bins(tau, j, Rj); xj = j + (np.arange(Rj) + 0.5) / Rj
    Aj = np.column_stack([np.cos(2 * np.pi * 26.5036 * xj), np.sin(2 * np.pi * 26.5036 * xj), np.ones(Rj)])
    c, *_ = np.linalg.lstsq(Aj, y, rcond=None)
    amp = np.hypot(c[0], c[1]) / (np.sinc(26.5036 / Rj))        # sinc-corrected for the bin width
    print(f" 2^{j:2d} (R={Rj:3d}): amplitude {amp:.4f}, phase {np.arctan2(-c[1], c[0]):+.3f}, rms(bin means) {np.sqrt((y**2).mean()):.3f}")
print("done")
