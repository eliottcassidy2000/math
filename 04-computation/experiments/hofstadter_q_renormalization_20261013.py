"""
hofstadter_q_renormalization_20261013.py -- why the Q jitter is scale-invariant (HYP-9259, statement 4).
Exact identity: d(n) + (1/2)[(1-1/n) d(n-1) + (1-2/n) d(n-2)] = s(n) + 3/(4n),  s(n) = (p/n) d(p) + (q/n) d(q),
p = n - Q(n-1), q = n - Q(n-2).  So d = H * s (up to O(1/n)) with H(omega) = 1/(1 + e^{-i omega}/2 + e^{-2 i omega}/2),
and s(n) is the previous generation read at p, q ~ n/2: a copy of d stretched by the factor two in n.  Hence the
renormalization fixed point S_d(omega) = |H(omega)|^2 S_s(omega), S_s(omega) ~ S_d(2 omega) for omega < pi/2, which has
no scale in it.  Tests on the octave [2^22, 2^23) (and [2^18, 2^19) for scale invariance):
  (1) acf of s at lags 1..12 against acf of d at the half lags (stretching: acf_s(2 l) ~ acf_d(l));
  (2) the averaged periodogram S_d(omega) (64 segments of 2^16) against |H(omega)|^2 S_d(2 omega) for omega < pi/2,
      and the shape of S_d on (pi/2, pi] against |H(omega)|^2 times a constant floor;
  (3) the infinite-product prediction S(omega) ~ prod_(m>=0, 2^m omega <= pi) |H(2^m omega)|^2 against the measured
      periodogram shape, with the resonance peak at omega_0 = 1.93 (period 3.25) and its subharmonic at omega_0/2.
"""
import numpy as np, time
from numba import njit

N = 1 << 24

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
def source(q, lo, hi):
    s = np.zeros(hi - lo); d = np.zeros(hi - lo); e = np.zeros(hi - lo)
    for n in range(lo, hi):
        p = n - q[n - 1]; qq = n - q[n - 2]
        dn = q[n] / n - 0.5
        sn = (p / n) * (q[p] / p - 0.5) + (qq / n) * (q[qq] / qq - 0.5)
        d[n - lo] = dn; s[n - lo] = sn
        e[n - lo] = dn + 0.5 * ((1 - 1.0 / n) * (q[n - 1] / (n - 1) - 0.5) + (1 - 2.0 / n) * (q[n - 2] / (n - 2) - 0.5)) - sn - 3.0 / (4 * n)
    return d, s, e

t0 = time.time()
q = hofq(N)
print(f"Q to 2^24 in {time.time()-t0:.1f}s")

def acf(x, L=12):
    x = x - x.mean(); v = (x * x).mean()
    return np.array([(x[:-l] * x[l:]).mean() / v for l in range(1, L + 1)])

def H2(w):
    return 1.0 / np.abs(1 + 0.5 * np.exp(-1j * w) + 0.5 * np.exp(-2j * w)) ** 2

for J in (18, 22):
    lo, hi = 1 << J, 1 << (J + 1)
    d, s, e = source(q, lo, hi)
    print(f"\n## octave 2^{J}: identity residual max {np.abs(e).max():.2e}; var(d) {d.var():.3e}, var(s) {s.var():.3e}, ratio var(d)/var(s) {d.var()/s.var():.3f}")
    ad, As = acf(d), acf(s)
    print(" acf_d(l), l=1..12: " + " ".join(f"{v:+.3f}" for v in ad))
    print(" acf_s(l), l=1..12: " + " ".join(f"{v:+.3f}" for v in As))
    print(" stretching test acf_s(2l) vs acf_d(l): " + " ".join(f"[{As[2*l-1]:+.3f}|{ad[l-1]:+.3f}]" for l in range(1, 7)))
    # averaged periodogram
    seg = 1 << 16; nseg = len(d) // seg
    dd = (d - d.mean())[: nseg * seg].reshape(nseg, seg)
    Sd = (np.abs(np.fft.rfft(dd, axis=1)) ** 2).mean(axis=0); Sd /= Sd.mean()
    w = np.fft.rfftfreq(seg) * 2 * np.pi
    def S_at(om):
        return float(np.interp(om, w, Sd))
    print(" periodogram S_d(omega) normalised to mean 1; the renormalization test S_d(omega) / [|H(omega)|^2 S_d(2 omega)] for omega < pi/2:")
    for om in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.966, 1.05, 1.2, 1.35, 1.5):
        print(f"   omega={om:.3f}: S={S_at(om):.3f}  |H|^2={H2(om):.3f}  S(2 omega)={S_at(2*om):.3f}  ratio S/(|H|^2 S(2 omega)) = {S_at(om) / (H2(om) * S_at(2*om)):.3f}")
    print(" upper half-band: S_d(omega) / |H(omega)|^2 (should be a flat floor if the stretched source is empty there):")
    for om in (1.6, 1.7, 1.8, 1.9, 1.93, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.14):
        print(f"   omega={om:.2f}: S={S_at(om):.3f}  |H|^2={H2(om):.3f}  S/|H|^2={S_at(om)/H2(om):.3f}")
    # infinite-product prediction (floor = 1 above pi/2), compared with the measured S on a grid, both normalised at omega_0 = 1.93
    def prod_pred(om):
        p = 1.0; x = om
        while x <= np.pi:
            p *= H2(x); x *= 2
        return p
    grid = np.array([0.3, 0.5, 0.7, 0.966, 1.2, 1.5, 1.93, 2.2, 2.6, 3.0])
    pm = np.array([prod_pred(o) for o in grid]); sm = np.array([S_at(o) for o in grid])
    print(" infinite-product prediction vs measured (both normalised at omega = 1.93):")
    print("   omega: " + " ".join(f"{o:6.3f}" for o in grid))
    print("   pred : " + " ".join(f"{v:6.3f}" for v in pm / pm[6]))
    print("   meas : " + " ".join(f"{v:6.3f}" for v in sm / sm[6]))
    i0 = int(np.argmax(Sd[w > 1.0])); print(f" measured peak in (1, pi]: omega = {w[w > 1.0][i0]:.3f} (period {2*np.pi/w[w>1.0][i0]:.3f}); predicted resonance omega_0 = {np.arctan2(np.sqrt(7), -1):.3f} (period {2*np.pi/np.arctan2(np.sqrt(7), -1):.3f})")
print("done")
