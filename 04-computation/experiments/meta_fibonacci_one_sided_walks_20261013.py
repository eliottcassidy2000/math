"""
meta_fibonacci_one_sided_walks_20261013.py -- Hofstadter's Q sequence (Q(1) = Q(2) = 1, Q(n) = Q(n - Q(n-1)) + Q(n - Q(n-2)))
and the Hofstadter-Conway sequence (a(1) = a(2) = 1, a(n) = a(a(n-1)) + a(n - a(n-1))) read with the Collatz instrument of
THM-4629: both are two-branch backward recursions whose branches both contract by about one octave (p, q ~ n/2), the
one-sided case of the transfer law, so the log-periodic spectrum of the deviation d(n) = x(n)/n - 1/2 is predicted to be a
comb at INTEGER frequencies (generations = octaves), with no irrational near-return lines.
Exact identities (from the definitions, with p, q the two backward indices):
  Conway:  d(n) = (p/n) d(p) + (q/n) d(q),            p = a(n-1), q = n - a(n-1), p + q = n           (conservative averaging)
  Q:       d(n) = (p/n) d(p) + (q/n) d(q) - (1/2)[(1-1/n) d(n-1) + (1-2/n) d(n-2)] + 3/(4n),  p = n - Q(n-1), q = n - Q(n-2)
           (the backward indices overshoot by the previous deviations: a damped feedback z^2 + z/2 + 1/2 = 0, |z| = 1/sqrt 2, period 3.25).
Computed to 2^30: per octave the bin means of d_a, d_Q (R = min(1024, 2^(J-6)) bins), the comb amplitudes at f = 1..6 and a
scan for non-integer peaks, the hump amplitude (rms and max) and its position, the self-similarity between consecutive
octaves, the step statistics log2(n/p), log2(n/q) (mean, std), the spikes Q(3 2^k) and a(3 2^k), the values at 20 and 1024,
the parities, and the Conway record (last n with a(n)/n >= 0.55).
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
            q[0] = n
            break
        q[n] = q[i] + q[j]
    return q

@njit(cache=True)
def octave_stats(x, J, R, kind):
    # kind 0: Conway (p = x[n-1], q = n - x[n-1]); kind 1: Q (p = n - x[n-1], q = n - x[n-2])
    s = np.zeros(R); sp = np.zeros(R); c = np.zeros(R)
    m1 = 0.0; v1 = 0.0; m2 = 0.0; v2 = 0.0; cnt = 0
    dmax = -1.0; nmax = 0; dmin = 1.0; nmin = 0
    for n in range(1 << J, 1 << (J + 1)):
        k = int((np.log2(n) - J) * R)
        if k >= R: k = R - 1
        d = x[n] / n - 0.5
        s[k] += d; c[k] += 1
        sp[k] += -1.0 if (x[n] & 1) else 1.0
        if d > dmax: dmax = d; nmax = n
        if d < dmin: dmin = d; nmin = n
        if kind == 0:
            p = x[n - 1]; qq = n - x[n - 1]
        else:
            p = n - x[n - 1]; qq = n - x[n - 2]
        l1 = np.log2(n / p); l2 = np.log2(n / qq)
        m1 += l1; v1 += l1 * l1; m2 += l2; v2 += l2 * l2; cnt += 1
    return s / np.maximum(c, 1), sp / np.maximum(c, 1), c, m1 / cnt, np.sqrt(max(v1 / cnt - (m1 / cnt) ** 2, 0.0)), m2 / cnt, np.sqrt(max(v2 / cnt - (m2 / cnt) ** 2, 0.0)), dmax, nmax, dmin, nmin

def ls_amp(x, y, f):
    A = np.column_stack([np.cos(2 * np.pi * f * x), np.sin(2 * np.pi * f * x), np.ones_like(x)])
    cc, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(np.hypot(cc[0], cc[1]))

t0 = time.time()
a = conway(N)
print(f"Conway a(n) to 2^30 in {time.time()-t0:.1f}s; a(20) = {a[20]}, a(1024) = {a[1024]}, a(1489) = {a[1489]}, a(2^30) = {a[N]}", flush=True)
t0 = time.time()
q = hofq(N)
print(f"Hofstadter Q(n) to 2^30 in {time.time()-t0:.1f}s; defined throughout: {q[0] == 0}; Q(20) = {q[20]}, Q(1024) = {q[1024]}, Q(1536) = {q[1536]}, Q(2^30) = {q[N]}", flush=True)

print("\n## values at 20 and 1024, and the spikes at 3 2^k")
print(" k : Q(3 2^k)  2^(k+1)  | a(3 2^k)  | a(2^k)  Q(2^k)")
for k in range(1, 29):
    n = 3 << k
    print(f" {k:2d}: {q[n]:10d} {1 << (k+1):9d}  | {a[n]:10d} | {a[1 << k]:9d} {q[1 << k]:9d}")
print(f" a(20)/20 = {a[20]/20:.4f}, Q(20)/20 = {q[20]/20:.4f}, a(1024)/1024 = {a[1024]/1024:.4f}, Q(1024)/1024 = {q[1024]/1024:.4f}")
# Conway record: last n with a(n)/n >= 0.55, and the maximum of a(n)/n per generation
nn = np.arange(1, 1 << 24, dtype=np.float64)
ratio = a[1:1 << 24] / nn
big = np.nonzero(ratio >= 0.55)[0]
print(f" Conway: last n < 2^24 with a(n)/n >= 0.55: n = {big[-1] + 1} (Mallows: 1488 is the last, 1489 the prize number); number of such n: {len(big)}")
print(" generation k: argmax of a(n)/n in [2^k, 2^(k+1)), its value, position x = n/2^k; and the same for Q")
for k in range(3, 24):
    lo, hi = 1 << k, 1 << (k + 1)
    ra = a[lo:hi] / np.arange(lo, hi, dtype=np.float64); ia = int(np.argmax(ra))
    rq = q[lo:hi] / np.arange(lo, hi, dtype=np.float64); iq = int(np.argmax(rq)); jq = int(np.argmin(rq))
    print(f"  k={k:2d}: Conway n={lo+ia:9d} ratio {ra[ia]:.5f} x={1+ia/lo:.4f} | Q max n={lo+iq:9d} ratio {rq[iq]:.5f} x={1+iq/lo:.4f}, Q min ratio {rq[jq]:.5f} at x={1+jq/lo:.4f}")

print("\n## per octave: comb amplitudes of d = x(n)/n - 1/2 at integer f (least squares), rms and max of d, step statistics, parity mean")
hdr = " J    R   | rms(d)   max d (at x)      | f=1      f=2      f=3      f=4      f=5      f=6   | non-integer peaks (f, amp) | steps: mean1 std1 mean2 std2 | parity mean"
prof = {0: {}, 1: {}}
for kind, name, x in ((0, "Conway", a), (1, "Q", q)):
    print(f"\n### {name}")
    print(hdr)
    for J in range(8, 30):
        R = int(min(1024, 1 << (J - 6)))
        md, mp, c, m1, s1, m2, s2, dmax, nmax, dmin, nmin = octave_stats(x, J, R, kind)
        xs = J + (np.arange(R) + 0.5) / R
        prof[kind][J] = md
        comb = [ls_amp(xs, md, f) for f in range(1, 7)]
        fs = np.arange(0.25, 8.0, 0.05)
        amps = np.array([ls_amp(xs, md, f) for f in fs])
        peaks = []
        for i in np.argsort(-amps):
            if abs(fs[i] - round(fs[i])) < 0.3: continue
            if all(abs(fs[i] - g) > 0.5 for g, _ in peaks):
                peaks.append((fs[i], amps[i]))
            if len(peaks) >= 2: break
        print(f" {J:2d} {R:5d} | {np.sqrt((md**2).mean()):.5f}  {dmax:+.5f} ({nmax / (1 << J):.3f}) | " + " ".join(f"{v:.5f}" for v in comb) + " | " + " ".join(f"({f:.2f},{v:.5f})" for f, v in peaks) + f" | {m1:.4f} {s1:.4f} {m2:.4f} {s2:.4f} | {mp.mean():+.4f}", flush=True)

print("\n## self-similarity: correlation between the hump of octave J and of octave J+1 (profiles at min(R, 16) points); rms of the profile; Conway rms x sqrt(J)")
for kind, name in ((0, "Conway"), (1, "Q")):
    for J in range(8, 29):
        h1 = prof[kind][J]; h2 = prof[kind][J + 1]
        P = min(len(h1), len(h2), 16)
        r1 = h1.reshape(P, -1).mean(axis=1); r2 = h2.reshape(P, -1).mean(axis=1)
        cc = np.corrcoef(r1, r2)[0, 1]
        print(f" {name:6s} J={J:2d}: corr(J, J+1) {cc:+.3f}, profile rms {np.sqrt((r1**2).mean()):.5f}, ratio of profile rms J+1/J {np.sqrt((r2**2).mean()) / np.sqrt((r1**2).mean()):.3f}, rms*sqrt(J) {np.sqrt((r1**2).mean()) * np.sqrt(J):.4f}")

print("\n## envelope of the Q jitter: rms of the pointwise deviation over 16 sub-intervals of each octave (x = 1 + (i+0.5)/16), and the sub-interval of the maximum")
@njit(cache=True)
def env16(x, J):
    s2 = np.zeros(16); c = np.zeros(16)
    lo = 1 << J
    for n in range(lo, lo << 1):
        k = int((n - lo) * 16 / lo)
        if k >= 16: k = 15
        d = x[n] / n - 0.5
        s2[k] += d * d; c[k] += 1
    return np.sqrt(s2 / np.maximum(c, 1))
for J in range(8, 30):
    e = env16(q, J)
    print(f" Q 2^{J:2d}: " + " ".join(f"{v:.4f}" for v in e) + f" | max at x={1 + (int(np.argmax(e)) + 0.5) / 16:.3f}, min/max {e.min() / e.max():.3f}, octave rms {np.sqrt((e**2).mean()):.5f}")
for J in (8, 12, 16, 20, 24, 28):
    e = env16(a, J)
    print(f" Conway 2^{J:2d}: " + " ".join(f"{v:.4f}" for v in e))
print("done")
