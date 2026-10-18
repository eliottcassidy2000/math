#!/usr/bin/env python3
"""S42 audit a5: re-derivation of the P13-P14 statistics (own code) and the checks of their reading.

Sampling as in the session: N = 2^24 + j*STEP < 2^28 (STEP = 2^12 for Y = 100, 2^14 for Y = 1000); z_u = c_Y(u)/sqrt(N)
over the seeds u <= Y.  Reproduces s_u, F_eff, model sd sqrt(2 sum s_u^2), observed sd, raw pooled kurtosis.  Adds:
  - centred pooled kurtosis (each normalised coefficient centred at its own sample mean, scaled by its own sd);
  - the raw 4th moment a centred-or-not Gaussian with the observed per-seed mean and sd would give (mu^4+6mu^2sig^2+3sig^4);
  - Var(E) = diag + off with diag = sum_u Var(z_u^2), off = sum_{u != v} Cov(z_u^2, z_v^2) (sign of the cross part);
  - the model sd of independent Gaussians with the OBSERVED means and sds: sqrt(sum (2 sig^4 + 4 mu^2 sig^2));
  - per-octave means of E/N and lag-1 autocorrelation of the sampled E/N series.
"""
import os, sys, json, time
import numpy as np
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from a2_independent_incremental import fmap


@njit(cache=False)
def sample(Nmax, Y, dest, sgn, Nstart, step):
    c = np.zeros(Y + 1, dtype=np.int64)
    ns = (Nmax - 1 - Nstart) // step + 1
    Z = np.zeros((ns, Y + 1))
    Ns = np.zeros(ns, dtype=np.int64)
    j = 0
    for N in range(1, Nmax):
        for t in range(3):
            m = 4 * N - 3 if t == 0 else (4 * N - 1 if t == 1 else 2 * N - 1)
            r = m % 3
            w = 2 if r == 2 else -1
            amt = (-3 * (0 if r == 0 else (1 if r == 1 else -1)) - w) if t == 2 else w
            i = (m - 1) // 2
            c[dest[i]] += sgn[i] * amt
        if N >= Nstart and (N - Nstart) % step == 0:
            q = np.sqrt(N)
            for u in range(Y + 1):
                Z[j, u] = c[u] / q
            Ns[j] = N
            j += 1
    return Z[:j], Ns[:j]


def main():
    Nmax = 2 ** 28; M = 4 * Nmax
    out = {}
    for Y, step in ((100, 2 ** 12), (1000, 2 ** 14)):
        t0 = time.time()
        dest, sgn = fmap(M, Y)
        Zall, Ns = sample(Nmax, Y, dest, sgn, 2 ** 24, step)
        del dest, sgn
        seeds = np.array([u for u in range(1, Y + 1) if u % 2 == 1 and u % 8 != 5])
        Z = Zall[:, seeds]
        su = (Z ** 2).mean(axis=0)
        m = su.sum(); sd_model = np.sqrt(2 * (su ** 2).sum()); Feff = m ** 2 / (su ** 2).sum()
        E = (Z ** 2).sum(axis=1)
        order = np.argsort(-su); cum = np.cumsum(su[order]) / m
        top = order[:int(np.searchsorted(cum, 0.9)) + 1]
        W = Z[:, top] / np.sqrt(su[top])
        kurt_raw = float((W ** 4).mean() / (W ** 2).mean() ** 2)
        mu = W.mean(axis=0); sig = W.std(axis=0)
        Wc = (W - mu) / sig
        kurt_centred = float((Wc ** 4).mean())
        gauss_raw = float(np.mean(mu ** 4 + 6 * mu ** 2 * sig ** 2 + 3 * sig ** 4))   # pooled raw 4th moment of N(mu, sig^2)
        # variance decomposition of E over all seeds
        Q = Z ** 2
        varE = float(E.var())
        diag = float(Q.var(axis=0).sum())
        off = varE - diag
        muZ = Z.mean(axis=0); sZ = Z.std(axis=0)
        sd_model_obs = float(np.sqrt(np.sum(2 * sZ ** 4 + 4 * muZ ** 2 * sZ ** 2)))
        # correlations as in the session (Pearson of normalised coefficients over samples)
        Cm = np.corrcoef(W.T); offc = Cm[~np.eye(len(top), dtype=bool)]
        # octave means and autocorrelation
        octs = {}
        for k in range(24, 28):
            sel = (Ns >= 2 ** k) & (Ns < 2 ** (k + 1))
            octs[k] = float(E[sel].mean())
        e0 = E - E.mean()
        ac1 = float((e0[1:] * e0[:-1]).mean() / e0.var())
        res = dict(Y=Y, samples=int(len(Ns)), F_Y=int(len(seeds)), F_eff=float(Feff), model_mean=float(m), model_sd=float(sd_model),
                   observed_mean=float(E.mean()), observed_sd=float(E.std()), observed_max=float(E.max()), top90=int(len(top)),
                   kurtosis_raw_pooled=kurt_raw, kurtosis_centred_pooled=kurt_centred,
                   raw_4th_moment_of_gaussians_with_observed_means=gauss_raw,
                   var_E=varE, diag_sum_var_zsq=diag, offdiag_sum_cov_zsq=off,
                   model_sd_gaussian_with_observed_means=sd_model_obs,
                   mean_abs_corr=float(np.abs(offc).mean()), max_abs_corr=float(np.abs(offc).max()),
                   frac_positive_corr=float((offc > 0).mean()),
                   max_abs_normalised_mean=float(np.abs(mu).max()), octave_means=octs, lag1_autocorr_E=ac1,
                   max_minus_mean_over_sd=float((E.max() - E.mean()) / E.std()))
        out[Y] = res
        print('Y = %d [%.0f s]' % (Y, time.time() - t0))
        for k, v in res.items():
            print('   %s: %s' % (k, v))
    json.dump({str(k): v for k, v in out.items()}, open(os.path.join(HERE, 'a5_statistics_check.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
