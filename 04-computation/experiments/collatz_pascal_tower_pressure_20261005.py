#!/usr/bin/env python3
"""
Pressure function and tail structure of the Pascal towers at moderate level (opus, 2026-10-05).

For the uniform-start q-tower (and the i.i.d.-digit model) the second moment is EXACTLY E_R |f_N(0)|^2 = 3^-N
(Fourier note 4k), while the typical (Lyapunov) rate is below 3^(-1/2): the difference is a Jensen gap.  This script
measures, from S independent samples at a moderate level N, the pressure
    P_N(s) = (1/N) log E_R |f_N(0)|^(2s),   s in [0, 1.5],
whose exact value at s = 1 is -log 3 (a sampling control), whose derivative at 0 is 2 lambda, and whose shape
between says where the q = 3 tower and the other towers differ; plus the share of the second moment carried by the
top quantiles of |f_N|^2 (how heavy the tail is), and quantiles of (1/N) log|f_N|.

Usage: python collatz_pascal_tower_pressure_20261005.py [N=100] [samples=200000] [workers=10] [models=tower3,tower5,tower7,tower9,iid]
"""
import sys, os, math, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import collatz_pascal_tower_lyapunov_20261005 as T


def chunk(args):
    model, q, N, seeds, Mfac = args
    out = np.zeros((len(seeds), 2))
    Mmax = Mfac * N
    for i, s in enumerate(seeds):
        bits, corr = T.make_bits(model, N, Mmax, q=q, u=1, seed=s)
        L = np.zeros(N)
        T.columns_lyap(bits, corr, Mmax, N, L)
        out[i, 0] = L[N - 1]
        out[i, 1] = L[N // 2 - 1]
    return out


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    S = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    models = sys.argv[4].split(",") if len(sys.argv) > 4 else ["tower3", "tower5", "tower7", "tower9", "iid"]
    import multiprocessing as mp
    lines = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); lines.append(s)
    P(f"pressure experiment: N={N}, samples={S}, Mmax=5N; exact control: P(1) = -log 3 = {-math.log(3):.5f}, i.e. (1/N) log E|f|^2 = -1.09861")
    svals = [0.1, 0.25, 0.5, 0.75, 0.9, 1.0, 1.25, 1.5]
    t0 = time.time()
    for mdl in models:
        if mdl.startswith("tower"):
            model, q = "tower", int(mdl[5:])
        else:
            model, q = mdl, 3
        jobs = []
        cs = 2000
        for a in range(0, S, cs):
            jobs.append((model, q, N, list(range(10_000_000 + a, 10_000_000 + min(a + cs, S))), 5))
        with mp.Pool(workers) as pool:
            res = pool.map(chunk, jobs)
        res = np.concatenate(res)
        LN = res[:, 0]; LH = res[:, 1]
        # pressure: (1/N) log mean exp(2 s L)
        row = []
        for s in svals:
            x = 2 * s * LN
            mx = x.max()
            Ps = (mx + math.log(np.mean(np.exp(x - mx)))) / N
            row.append(Ps)
        lam = LN.mean() / N
        lam_half = (LN - LH).mean() / (N - N // 2)
        se = (LN - LH).std(ddof=1) / math.sqrt(len(LN)) / (N - N // 2)
        # share of the second moment from the top quantiles
        w = np.exp(2 * (LN - LN.max())); w /= w.sum()
        ws = np.sort(w)[::-1]
        share = {p: ws[: max(1, int(p * len(ws)))].sum() for p in (0.001, 0.01, 0.1)}
        qs = np.quantile(LN / N, [0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999])
        P(f"{mdl}: [{time.time()-t0:.0f}s] lambda_N (mean L_N/N) = {lam:.5f}; (L_N - L_N/2)/(N/2) = {lam_half:.5f} +- {se:.5f} (rate {math.exp(lam_half):.5f})")
        P("   P_N(s) at s=" + ", ".join(f"{s}" for s in svals) + ":  " + ", ".join(f"{v:.4f}" for v in row) + f"   [P(1) exact -1.0986]")
        P(f"   second-moment share of top 0.1% / 1% / 10% of samples: {share[0.001]:.3f} / {share[0.01]:.3f} / {share[0.1]:.3f}")
        P("   quantiles of L_N/N (0.1%,1%,10%,50%,90%,99%,99.9%): " + ", ".join(f"{v:.4f}" for v in qs))
        P(f"   Var(L_N)/N = {LN.var()/N:.4f}; Var(L_N - L_N/2)/(N/2) = {(LN-LH).var()/(N-N//2):.4f}; skew of L_N: {((LN-LN.mean())**3).mean()/LN.std()**3:.3f}")
    with open(__file__.replace(".py", f"_N{N}.out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
