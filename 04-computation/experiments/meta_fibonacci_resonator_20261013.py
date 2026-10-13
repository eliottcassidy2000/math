"""
meta_fibonacci_resonator_20261013.py -- follow-up to meta_fibonacci_one_sided_walks_20261013.py.
(1) FINITE-EXACT: the deviation identities, checked to machine precision for n <= 2^22:
      Conway:  d(n) = (p/n) d(p) + (q/n) d(q),  p = a(n-1), q = n - a(n-1)
      Q:       d(n) = (p/n) d(p) + (q/n) d(q) - (1/2)[(1-1/n) d(n-1) + (1-2/n) d(n-2)] + 3/(4n),  p = n - Q(n-1), q = n - Q(n-2)
(2) The feedback resonator of Q: the homogeneous part d(n) + d(n-1)/2 + d(n-2)/2 = 0 has roots z = (-1 +- i sqrt 7)/4,
    |z| = 1/sqrt 2, period 2 pi/arg z = 3.253; the Yule-Walker autocorrelation of an AR(2) with phi1 = phi2 = -1/2 is
    r = (-1/3, -1/3, +1/3, 0, -1/6, +1/12, +1/24, -1/16, ...).  Measured: the autocorrelation of d_Q(n) at lags 1..12 in
    octaves 12, 16, 20, 24, 28 (after removing the octave mean), and the same for Conway; the dominant period from the
    power spectrum of a 2^16 window in each octave.
(3) Conway hump convergence: A_k = max d over [2^k, 2^(k+1)), A_k sqrt(k), the position of the maximum, and the profile at
    16 points for k = 8, 12, 16, 20, 24, 28; Q: the same statistics plus the rms.
(4) Where p, q fall: the fraction of n in generation k whose backward indices p, q lie in generation k-1 (both), and the
    distribution of log2(n/p), log2(n/q) around 1 (mean, std) -- the one-sided walk.
"""
import numpy as np, time
from numba import njit

N = 1 << 30

@njit(cache=True)
def conway(N):
    a = np.zeros(N + 1, dtype=np.int32); a[1] = 1; a[2] = 1
    for n in range(3, N + 1):
        a[n] = a[a[n - 1]] + a[n - a[n - 1]]
    return a

@njit(cache=True)
def hofq(N):
    q = np.zeros(N + 1, dtype=np.int32); q[1] = 1; q[2] = 1
    for n in range(3, N + 1):
        i = n - q[n - 1]; j = n - q[n - 2]
        if i < 1 or j < 1:
            q[0] = n; break
        q[n] = q[i] + q[j]
    return q

@njit(cache=True)
def check_identities(a, q, N):
    err_a = 0.0; err_q = 0.0
    for n in range(3, N + 1):
        p = a[n - 1]; qq = n - a[n - 1]
        d = a[n] / n - 0.5
        rhs = (p / n) * (a[p] / p - 0.5) + (qq / n) * (a[qq] / qq - 0.5)
        e = abs(d - rhs)
        if e > err_a: err_a = e
        p2 = n - q[n - 1]; q2 = n - q[n - 2]
        dq = q[n] / n - 0.5
        rhs2 = (p2 / n) * (q[p2] / p2 - 0.5) + (q2 / n) * (q[q2] / q2 - 0.5) - 0.5 * ((1 - 1.0 / n) * (q[n - 1] / (n - 1) - 0.5) + (1 - 2.0 / n) * (q[n - 2] / (n - 2) - 0.5)) + 3.0 / (4 * n)
        e2 = abs(dq - rhs2)
        if e2 > err_q: err_q = e2
    return err_a, err_q

@njit(cache=True)
def gen_stats(x, J, kind):
    lo = 1 << J; hi = 1 << (J + 1)
    both_prev = 0; cnt = 0
    for n in range(lo, hi):
        if kind == 0:
            p = x[n - 1]; qq = n - x[n - 1]
        else:
            p = n - x[n - 1]; qq = n - x[n - 2]
        if (lo // 2) <= p < lo and (lo // 2) <= qq < lo:
            both_prev += 1
        cnt += 1
    return both_prev / cnt

t0 = time.time()
a = conway(N); q = hofq(N)
print(f"sequences to 2^30 in {time.time()-t0:.1f}s")
ea, eq = check_identities(a, q, 1 << 22)
print(f"(1) identities checked for n <= 2^22: max |error| Conway {ea:.2e}, Q {eq:.2e}  (float64 round-off level)")

print("\n(2) autocorrelation of the deviation at lags 1..12 (octave mean removed); AR(2) prediction for Q: -0.333 -0.333 +0.333 0 -0.167 +0.083 +0.042 -0.062 ...")
for kind, name, x in ((1, "Q", q), (0, "Conway", a)):
    for J in (12, 16, 20, 24, 28):
        lo, hi = 1 << J, 1 << (J + 1)
        d = x[lo:hi].astype(np.float64) / np.arange(lo, hi, dtype=np.float64) - 0.5
        d = d - d.mean()
        v = (d * d).mean()
        r = [float((d[:-l] * d[l:]).mean() / v) for l in range(1, 13)]
        # dominant period from the power spectrum of the first 2^16 samples (detrended by a running mean of 64)
        w = d[: 1 << 16]
        w = w - np.convolve(w, np.ones(64) / 64, mode="same")
        P = np.abs(np.fft.rfft(w)) ** 2
        freqs = np.fft.rfftfreq(len(w))
        sel = freqs > 1.0 / 32
        fpk = freqs[sel][int(np.argmax(P[sel]))]
        print(f" {name:6s} 2^{J}: r = " + " ".join(f"{v:+.3f}" for v in r) + f" | dominant short period {1/fpk:.3f} (n units)")

print("\n(3) hump statistics per generation: Conway A_k = max d, A_k sqrt(k), x of max; Q max d, min d, rms")
for J in range(4, 30):
    lo, hi = 1 << J, 1 << (J + 1)
    nn = np.arange(lo, hi, dtype=np.float64)
    da = a[lo:hi] / nn - 0.5; dq = q[lo:hi] / nn - 0.5
    ia = int(np.argmax(da))
    print(f" k={J:2d}: Conway A={da[ia]:.5f} A sqrt(k)={da[ia]*np.sqrt(J):.4f} x={1+ia/lo:.4f} rms={np.sqrt((da**2).mean()):.5f} | Q max={dq.max():+.5f} min={dq.min():+.5f} rms={np.sqrt((dq**2).mean()):.5f}")
print("\n Conway profile d(x) at 16 points x = 1 + (i+0.5)/16, for k = 8, 12, 16, 20, 24, 28:")
for J in (8, 12, 16, 20, 24, 28):
    lo, hi = 1 << J, 1 << (J + 1)
    da = a[lo:hi] / np.arange(lo, hi, dtype=np.float64) - 0.5
    prof = da.reshape(16, -1).mean(axis=1)
    print(f"  k={J:2d}: " + " ".join(f"{v:+.4f}" for v in prof))
print(" Q profile (16-point means) for the same k:")
for J in (8, 12, 16, 20, 24, 28):
    lo, hi = 1 << J, 1 << (J + 1)
    dq = q[lo:hi] / np.arange(lo, hi, dtype=np.float64) - 0.5
    prof = dq.reshape(16, -1).mean(axis=1)
    print(f"  k={J:2d}: " + " ".join(f"{v:+.4f}" for v in prof))

print("\n(4) fraction of n in generation k with both backward indices in generation k-1")
for J in (8, 12, 16, 20, 24, 28):
    print(f" k={J}: Conway {gen_stats(a, J, 0):.4f}, Q {gen_stats(q, J, 1):.4f}")
print("done")
