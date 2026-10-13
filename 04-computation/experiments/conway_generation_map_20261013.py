"""
Lagged numerical model for Conway's generation profile; no discrete-model
error or asymptotic rate theorem is implemented here.

The exact deviation identity conserves child mass. Its selector is
s_n=a(n-1)/n-1/2=h_k(x)-(a(n)-a(n-1))/n, with x=n/2^k.
Both child arguments are in the CLOSED previous arc. This code instead
uses the lagged h_(k-1)(x) to predict the next profile and includes an
interpolation rule. Replacing the current selector by that lagged value
requires a separate, unproved error estimate; it is not justified by the
one-step O(1/n) slow-growth bound. The Taylor calculation applies formally
to smooth profiles of this lagged model. Fitted inverse-square-root decay
and comparison with generations15..29 remain numerical evidence.

Focused correction and exact closed-arc proof:
05-knowledge/results/conway_tail_barrier_20261013c.md, section6.
The algorithm and historical output are unchanged by this clarification.
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

def interp(h, y):
    # h on the grid x_i = 1 + (i + 0.5)/G ; linear interpolation, h(1) = h(2) = 0
    pos = (y - 1.0) * G - 0.5
    i = int(np.floor(pos))
    if i < 0:
        return h[0] * (pos + 1.0)          # from h(1) = 0
    if i >= G - 1:
        return h[G - 1] * (G - 0.5 - pos)  # to h(2) = 0
    f = pos - i
    return (1 - f) * h[i] + f * h[i + 1]

def cgm_step(h):
    hn = np.zeros(G)
    for i in range(G):
        x = 1.0 + (i + 0.5) / G
        d = h[i]
        yp = x * (1 + 2 * d); ym = x * (1 - 2 * d)
        vp = interp(hn, yp / 2) if yp >= 2.0 else interp(h, yp)     # same generation (already computed, since yp/2 < x)
        vm = interp(h, ym) if ym >= 1.0 else interp(h, ym * 2)        # (ym < 1 would be the generation after; does not occur for h > 0)
        hn[i] = (0.5 + d) * vp + (0.5 - d) * vm
    return hn

t0 = time.time()
a = conway(N)
print(f"a(n) to 2^30 in {time.time()-t0:.1f}s")
measured = {J: profile(a, J, G) for J in range(10, 30)}
xs = 1 + (np.arange(G) + 0.5) / G

def stats(h):
    i = int(np.argmax(h)); p16 = h.reshape(16, -1).mean(axis=1)
    return h[i], xs[i], float(np.sqrt((h ** 2).mean())), p16

print("\n## CGM iterated from the measured generation 14 (4096 points) against the measured generations")
print(" k  | measured: max   x_max   rms    | CGM: max    x_max   rms    | corr(16-pt profiles)  max ratio CGM/measured")
h = measured[14].copy()
for J in range(15, 30):
    h = cgm_step(h)
    mA, mx, mr, mp = stats(measured[J]); cA, cx, cr, cp = stats(h)
    cc = np.corrcoef(mp, cp)[0, 1]
    print(f" {J:2d} | {mA:.5f} {mx:.4f} {mr:.5f} | {cA:.5f} {cx:.4f} {cr:.5f} | {cc:+.5f}  {cA / mA:.4f}", flush=True)

print("\n## CGM iterated from generation 14 to k = 400: A_k, A_k sqrt(k), x_max, and the fitted 1/A^2 slope")
h = measured[14].copy()
rows = []
for J in range(15, 401):
    h = cgm_step(h)
    A, x, r, p = stats(h)
    rows.append((J, A, x, r))
    if J in (20, 30, 40, 60, 80, 100, 150, 200, 300, 400):
        print(f" k={J:3d}: A={A:.6f}  A sqrt(k)={A*np.sqrt(J):.5f}  x_max={x:.4f}  rms={r:.6f}  rms sqrt(k)={r*np.sqrt(J):.5f}")
J = np.array([r[0] for r in rows]); A = np.array([r[1] for r in rows])
for lo_, hi_ in ((20, 60), (100, 200), (200, 400)):
    sel = (J >= lo_) & (J <= hi_)
    sl, ic = np.polyfit(J[sel], 1 / A[sel] ** 2, 1)
    print(f" fit 1/A^2 = {sl:.3f} k + {ic:.1f} over k in [{lo_}, {hi_}]: C = 1/sqrt(slope) = {1/np.sqrt(sl):.5f}, k0 = {-ic/sl:.2f}")
print(" measured 1/A^2 over k = 10..29 (from the a array):")
Jm = np.arange(10, 30); Am = np.array([measured[J].max() for J in Jm])
sl, ic = np.polyfit(Jm, 1 / Am ** 2, 1)
print(f"  fit 1/A^2 = {sl:.3f} k + {ic:.1f}: C = {1/np.sqrt(sl):.5f}, k0 = {-ic/sl:.2f}")
print(" limiting shape (CGM at k = 400, scaled by sqrt(k), 16-point profile):")
A, x, r, p = stats(h)
print("  " + " ".join(f"{v*np.sqrt(400):+.4f}" for v in p))
print(" measured generation 29 scaled by sqrt(29):")
A, x, r, p = stats(measured[29])
print("  " + " ".join(f"{v*np.sqrt(29):+.4f}" for v in p))
print("done")
