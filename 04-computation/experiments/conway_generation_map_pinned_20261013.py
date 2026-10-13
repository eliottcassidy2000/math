"""
conway_generation_map_pinned_20261013.py -- the Conway continuum generation map with the exact boundary conditions.
Exact facts (from the recursion): a(2^k) = 2^(k-1) and, by induction on j, a(2^k + j) = 2^(k-1) + j for all small j
(a(2^k + j) = a(a(2^k + j - 1)) + a(2^k + j - a(2^k + j - 1)) = a(2^(k-1) + j - 1) + a(2^(k-1) + 1) while the pattern holds
in generation k-1), so d(2^k + j) = (j/2)/(2^k + j): the hump starts with slope 1/2 in x and the boundary value is 0;
at the right end a is non-decreasing with increments 0 or 1 and a(2^(k+1)) = 2^k, so d(2^(k+1) - 1) ~ 0 too.
Part 1 (FINITE-EXACT): j*(k) = the largest j with a(2^k + i) = 2^(k-1) + i for all i <= j, for k = 4..29.
Part 2: the map of conway_generation_map_20261013.py with both ends pinned to 0 after every step, iterated from the
measured generation 14 to k = 400: A_k, A_k sqrt(k), x_max, and the 1/A^2 fits; the limiting profile.
"""
import numpy as np, time
from numba import njit

N = 1 << 30
G = 4096

@njit(cache=True)
def conway(N):
    a = np.zeros(N + 1, dtype=np.int32); a[1] = 1; a[2] = 1
    for n in range(3, N + 1):
        a[n] = a[a[n - 1]] + a[n - a[n - 1]]
    return a

@njit(cache=True)
def profile(a, J, R):
    s = np.zeros(R); c = np.zeros(R)
    lo = 1 << J
    for n in range(lo, lo << 1):
        k = int((n - lo) * R / lo)
        if k >= R: k = R - 1
        s[k] += a[n] / n - 0.5; c[k] += 1
    return s / np.maximum(c, 1)

@njit(cache=True)
def jstar(a, k):
    base = 1 << (k - 1); lo = 1 << k; j = 0
    while lo + j + 1 < (lo << 1) and a[lo + j + 1] == base + j + 1:
        j += 1
    return j

def interp(h, y):
    pos = (y - 1.0) * G - 0.5
    i = int(np.floor(pos))
    if i < 0:
        return h[0] * (pos + 1.0)
    if i >= G - 1:
        return h[G - 1] * (G - 0.5 - pos)
    f = pos - i
    return (1 - f) * h[i] + f * h[i + 1]

def cgm_step(h, pin):
    hn = np.zeros(G)
    for i in range(G):
        x = 1.0 + (i + 0.5) / G
        d = h[i]
        yp = x * (1 + 2 * d); ym = x * (1 - 2 * d)
        vp = interp(hn, yp / 2) if yp >= 2.0 else interp(h, yp)
        vm = interp(h, ym) if ym >= 1.0 else interp(h, ym * 2)
        hn[i] = (0.5 + d) * vp + (0.5 - d) * vm
    if pin:
        hn[0] = 0.0; hn[G - 1] = 0.0
    return hn

t0 = time.time()
a = conway(N)
print(f"a(n) to 2^30 in {time.time()-t0:.1f}s")
print("\n# Part 1: j*(k) = length of the exact initial run a(2^k + j) = 2^(k-1) + j")
print(" " + " ".join(f"k={k}:{jstar(a, k)}" for k in range(4, 30)))

xs = 1 + (np.arange(G) + 0.5) / G
h14 = profile(a, 14, G)
measured = {J: profile(a, J, G) for J in (20, 25, 29)}
def stats(h):
    i = int(np.argmax(h)); return h[i], xs[i], float(np.sqrt((h ** 2).mean()))
print("\n# Part 2: pinned map (both ends 0) from generation 14")
h = h14.copy(); rows = []
for J in range(15, 401):
    h = cgm_step(h, True)
    A, x, r = stats(h); rows.append((J, A, x, r))
    if J in (20, 25, 29, 40, 60, 80, 100, 150, 200, 300, 400):
        extra = ""
        if J in measured:
            mA, mx, mr = stats(measured[J]); extra = f"   measured A={mA:.5f} ratio {A/mA:.4f}"
        print(f" k={J:3d}: A={A:.6f} A sqrt(k)={A*np.sqrt(J):.5f} x_max={x:.4f} rms={r:.6f} rms sqrt(k)={r*np.sqrt(J):.5f}{extra}", flush=True)
J = np.array([r[0] for r in rows]); A = np.array([r[1] for r in rows])
for lo_, hi_ in ((20, 60), (60, 150), (150, 400)):
    sel = (J >= lo_) & (J <= hi_)
    sl, ic = np.polyfit(J[sel], 1 / A[sel] ** 2, 1)
    print(f" fit 1/A^2 = {sl:.3f} k + {ic:.1f} over [{lo_}, {hi_}]: C = {1/np.sqrt(sl):.5f}, k0 = {-ic/sl:.2f}")
print(" pinned profile at k = 400 scaled by sqrt(400), 16 points:")
print("  " + " ".join(f"{v*20:+.4f}" for v in h.reshape(16, -1).mean(axis=1)))
print(" measured generation 29 scaled by sqrt(29):")
print("  " + " ".join(f"{v*np.sqrt(29):+.4f}" for v in measured[29].reshape(16, -1).mean(axis=1)))
print("done")
