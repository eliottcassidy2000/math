"""
hofstadter_q_envelope_node_20261013.py -- where is the calm node of the Q jitter envelope (HYP-9259, statement 4)?
The rms of d(n) = Q(n)/n - 1/2 over 128 sub-intervals of each generation [2^k, 2^(k+1)), k = 22..29; the sub-interval of
the minimum, a parabolic interpolation of its position in x = n/2^k, the value there, and the position of the maximum;
the candidates sqrt 2 = 1.41421 (the midpoint of the generation in log scale) and 1.5 (the old spike position).
"""
import numpy as np, time
from numba import njit

N = 1 << 30
B = 128

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
def env(q, J, B):
    s2 = np.zeros(B); c = np.zeros(B); lo = 1 << J
    for n in range(lo, lo << 1):
        k = int((n - lo) * B / lo)
        if k >= B: k = B - 1
        d = q[n] / n - 0.5
        s2[k] += d * d; c[k] += 1
    return np.sqrt(s2 / np.maximum(c, 1))

t0 = time.time(); q = hofq(N); print(f"Q to 2^30 in {time.time()-t0:.1f}s")
print(" k : node sub-interval [x_lo, x_hi]   parabolic x_min   rms at node   rms max (at x)   node/max   | rms at the bins containing sqrt2 and 1.5")
for J in range(22, 30):
    e = env(q, J, B); xs = 1 + (np.arange(B) + 0.5) / B
    i = int(np.argmin(e[8:-8])) + 8
    # parabolic interpolation on log rms
    y0, y1, y2 = np.log(e[i - 1]), np.log(e[i]), np.log(e[i + 1])
    denom = (y0 - 2 * y1 + y2)
    off = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
    xmin = xs[i] + off / B
    j = int(np.argmax(e))
    b_s2 = int((np.sqrt(2) - 1) * B); b_15 = int(0.5 * B)
    print(f" {J:2d}: [{1 + i / B:.4f}, {1 + (i + 1) / B:.4f}]   {xmin:.4f}          {e[i]:.5f}      {e[j]:.5f} ({xs[j]:.3f})   {e[i] / e[j]:.3f}   | {e[b_s2]:.5f}  {e[b_15]:.5f}")
print(" sqrt 2 = 1.41421; the 128-bin profile at 2^29 around the node (x from 1.33 to 1.50):")
e = env(q, 29, B); xs = 1 + (np.arange(B) + 0.5) / B
sel = (xs > 1.33) & (xs < 1.50)
print("  " + " ".join(f"{x:.3f}:{v:.4f}" for x, v in zip(xs[sel], e[sel])))
print("done")
