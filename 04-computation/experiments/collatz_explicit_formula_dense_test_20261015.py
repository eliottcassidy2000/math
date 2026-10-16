#!/usr/bin/env python3
"""Dense test of the explicit formula for the Collatz parity sum (HYP-9285) on exact data through 2^LOG (default 27).

S(N) = sum_(n<=N) t(n), t(n) = (-1)^tau(n), is sampled exactly at 2048 log-uniform points per octave on [2^16, 2^LOG)
(dense enough that no zero below Im s = 8000, i.e. about 880 cycles per octave, aliases: the Nyquist limit is 1024
cycles per octave).  Least squares for S(N) ~ Re sum_rho c_rho N^rho over the zeros rho of D(s) = 1 + 2^(-s) + (3/2)^s/3
with Re rho > RMIN and 0 < Im rho < T is fitted on [2^18, 2^(LOG-2)) and tested on [2^(LOG-2), 2^LOG); reported:
relative rms errors and the rms of the holdout residual divided by sqrt(N), against the rms of H(N)/sqrt(N),
H(N) = sum_(n<=N) (3[n = 2 mod 3] - 1) t(n).  If S = (modal sum) + R exactly, then H = -[R(N) + 3R(2N/3) + 3R(N/3)];
for a remainder uncorrelated across the three scales this predicts rms(R/sqrt N) ~ rms(H/sqrt N)/sqrt(10).
NUMERICAL.  Session opus-2026-10-15-S37.  Usage: python collatz_explicit_formula_dense_test_20261015.py LOG out.json
"""
import sys, os, json, time
import importlib.util
import numpy as np
import mpmath as mp
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('cen', os.path.join(HERE, 'collatz_octave_memory_census_20261015.py'))
cen = importlib.util.module_from_spec(spec)
_argv = sys.argv
sys.argv = ['x', '24', '20', 'unused.json', 'unused.npz']
spec.loader.exec_module(cen)
sys.argv = _argv
LOG = int(sys.argv[1]) if len(sys.argv) > 1 else 27
OUT = sys.argv[2] if len(sys.argv) > 2 else 'explicit_dense.json'
PER = 2048
TMAX = 8000.0


@nb.njit
def prefix_at(P, pts):
    out = np.zeros(pts.shape[0], dtype=np.int64)
    out3 = np.zeros((pts.shape[0], 3), dtype=np.int64)
    S = 0
    C3 = np.zeros(3, dtype=np.int64)
    gi = 0
    for n in range(1, pts[-1] + 1):
        t = 1 - 2 * int((P[n >> 6] >> np.uint64(n & 63)) & np.uint64(1))
        S += t
        C3[n % 3] += t
        while gi < pts.shape[0] and pts[gi] == n:
            out[gi] = S
            for a in range(3):
                out3[gi, a] = C3[a]
            gi += 1
    return out, out3


def D(s):
    return 1 + mp.power(2, -s) + mp.power(mp.mpf(3) / 2, s) / 3


def zeros_upto(T, rmin=0.5):
    zs = []
    for t in np.arange(0.5, T, 0.2):
        for sg in (0.55, 0.7, 0.85, 0.95, 0.995):
            try:
                z = mp.findroot(D, mp.mpc(sg, t))
            except Exception:
                continue
            if rmin < z.real < 1 and 0 < z.imag < T and all(abs(z - w) > 1e-7 for w in zs):
                zs.append(z)
    return sorted(zs, key=lambda z: float(z.imag))


if __name__ == '__main__':
    t0 = time.time()
    P = cen.parity_bits(LOG)
    pts = np.unique(np.round(2.0 ** (16 + np.arange(0, (LOG - 16) * PER) / PER)).astype(np.int64))
    pts = pts[pts < (1 << LOG)]
    S, C3 = prefix_at(P, pts)
    S = S.astype(np.float64)
    H = (3 * C3[:, 2] - C3.sum(axis=1)).astype(np.float64)
    N = pts.astype(np.float64)
    ZS = zeros_upto(TMAX)
    print(f'{len(pts)} samples; {len(ZS)} zeros of D with Re > 1/2, Im < {TMAX:.0f} ({time.time() - t0:.0f} s)', flush=True)
    lnN = np.log(N)
    fit = (N >= 2 ** 18) & (N < 2 ** (LOG - 2))
    te = N >= 2 ** (LOG - 2)
    hsq = float(np.sqrt(((H[te] / np.sqrt(N[te])) ** 2).mean()))
    rows = []
    for T in (300, 1000, 3000, 5000, 8000):
        for rmin in (0.65, 0.5):
            zs = [z for z in ZS if float(z.imag) < T and float(z.real) > rmin]
            cols = []
            for z in zs:
                w = np.exp(float(z.real) * lnN) * np.exp(1j * float(z.imag) * lnN)
                cols += [w.real, -w.imag]
            A = np.stack(cols, axis=1)
            coef, *_ = np.linalg.lstsq(A[fit], S[fit], rcond=None)
            r = S - A @ coef
            row = dict(T=T, rmin=rmin, n_zeros=len(zs),
                       rel_rms_fit=float(np.sqrt((r[fit] ** 2).mean() / (S[fit] ** 2).mean())),
                       rel_rms_holdout=float(np.sqrt((r[te] ** 2).mean() / (S[te] ** 2).mean())),
                       holdout_residual_over_sqrtN=float(np.sqrt(((r[te] / np.sqrt(N[te])) ** 2).mean())))
            rows.append(row)
            print(f"Im < {T}, Re > {rmin}: {len(zs)} zeros; relative rms fit {row['rel_rms_fit']:.4f}, holdout "
                  f"{row['rel_rms_holdout']:.4f}; holdout residual / sqrt(N) rms {row['holdout_residual_over_sqrtN']:.3f} "
                  f"(H / sqrt(N) rms {hsq:.3f}; H / sqrt(10 N) {hsq / np.sqrt(10):.3f})", flush=True)
    json.dump(dict(logmax=LOG, samples_per_octave=PER, fit=[2 ** 18, 2 ** (LOG - 2)], holdout=[2 ** (LOG - 2), 2 ** LOG],
                   H_over_sqrtN_rms_holdout=hsq, zeros=[[str(z.real), str(z.imag)] for z in ZS], rows=rows),
              open(OUT, 'w'), indent=1)
