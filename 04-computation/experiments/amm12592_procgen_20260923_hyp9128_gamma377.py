#!/usr/bin/env python3
"""THM-4467 extension: subordination certificate for Lambda(gamma) > 0 at gamma = 377/1000 (C* >= 1.377).

Same method as the gamma = 3/8 certificate (amm12592_procgen_20260923_polya_capacity.py [6]):
 * Riemann map of W_gamma from a lightning Laplace solve (numerical, used only to BUILD the certificate);
 * certificate = explicit polynomial P(z) = sum_k b_k z^k with real float64 coefficients, P(0) = 0, b_1 > 1;
 * check: P(e^{it}) lies in W_gamma for all t (FFT grid + explicit Lipschitz and rounding bounds; the
   grid is the union of two interleaved 2^L grids), hence R(W_gamma,0) >= P'(0) > 1 by Schwarz's lemma,
   hence Lambda(gamma) > 0 and, by THM-4467, no fair extractor has T(n) <= (1+gamma) n + D.
Usage: python3 amm12592_procgen_20260923_hyp9128_gamma377.py [--gamma-num 377 --gamma-den 1000] [--rr 0.99975]
       python3 amm12592_procgen_20260923_hyp9128_gamma377.py --verify-only CERT.json
"""
from __future__ import annotations
import argparse, hashlib, json, math, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from amm12592_procgen_20260923_polya_capacity import solve_W, T_member_w   # lightning solve on W_gamma

REPO = Path(__file__).resolve().parents[2]
DEFAULT_CERT = REPO / "05-knowledge" / "results" / "amm12592_procgen_20260923_hyp9128_gamma377_certificate.json"


def riemann_inverse_coeffs(g, rr, M, npol=60, nlight=36, steps=160, iters=6):
    lam0, lt, wc = solve_W(g, npol, nlight)
    R = math.exp(lam0)
    H0 = lt.H(0.0)[0]
    c = lt.c
    npol_ = lt.npol
    poles = lt.poles

    def H_and_dH(w):
        zz = (w - lt.zc) / lt.sc
        # polynomial part by Horner, with derivative
        Hp = np.zeros_like(w)
        dHp = np.zeros_like(w)
        for k in range(npol_, -1, -1):
            dHp = dHp * zz + Hp
            Hp = Hp * zz + c[k]
        dHp = dHp / lt.sc
        Hr = np.zeros_like(w)
        dHr = np.zeros_like(w)
        for j, (pj, d) in enumerate(poles):
            cj = c[npol_ + 1 + j]
            inv = 1.0 / (w - pj)
            Hr += cj * d * inv
            dHr -= cj * d * inv * inv
        return Hp + Hr - 1j * H0.imag, dHp + dHr

    t = 2 * math.pi * np.arange(M) / M

    def damped(w, z, n_it):
        # damped Newton for Phi(w) = z that never leaves W_gamma (spurious roots beyond the reentrant corner are
        # possible because the lightning H continues analytically into the notch)
        for _ in range(n_it):
            Hw, dH = H_and_dH(w)
            e = np.exp(-Hw)
            F = w * e - z
            step = F / (e * (1 - w * dH))
            lam = np.ones(len(w))
            for _k in range(30):
                wn = w - lam * step
                Hn, _ = H_and_dH(wn)
                bad = (np.abs(wn * np.exp(-Hn) - z) > np.abs(F)) | (T_member_w(wn, g) >= 1)
                if not bad.any():
                    break
                lam = np.where(bad, lam / 2, lam)
            w = np.where(bad, w, wn)
        return w

    w = None
    for s in np.linspace(0, rr, steps + 1)[1:]:
        z = s * np.exp(1j * t)
        if w is None:
            w = z * R
        w = damped(w, z, iters)
    # polishing at the final radius
    z = rr * np.exp(1j * t)
    w = damped(w, z, 60)
    Hw, _ = H_and_dH(w)
    resid = float(np.max(np.abs(w * np.exp(-Hw) - rr * np.exp(1j * t))))
    coeff = np.fft.fft(w) / M
    b = coeff[: M // 2].real.copy()
    b[0] = 0.0
    return b, R, lt.resid, resid


def verify(b, g, logM=22, eval_allow=1e-9):
    """P(e^{it}) in W_gamma for all t in [0, pi] (real coefficients): two interleaved grids of 2^logM points."""
    M = 1 << logM
    n = len(b)
    assert b[0] == 0.0 and b[1] > 1.0
    S1 = float(np.sum(np.arange(n) * np.abs(b))) * (1 + 1e-12)
    worst = (-1.0, None)
    maxT = -1.0
    minsep = np.inf
    for shift in (0.0, 0.5):
        pad = np.zeros(M, dtype=complex)
        pad[:n] = b * np.exp(1j * 2 * math.pi * shift * np.arange(n) / M)
        vals = np.fft.ifft(pad) * M          # P(e^{2 pi i (j+shift)/M})
        half = M // 2 + 1
        wv = vals[:half]
        del vals, pad
        p = (1 - np.sqrt(1 - 4 * wv)) / 2
        one_m_2p = np.abs(1 - 2 * p)
        ap, aq = np.abs(p), np.abs(1 - p)
        m_, S_ = np.minimum(ap, aq), ap + aq
        Tj = m_ * S_ ** g
        delta = eval_allow + (math.pi / (2 * M)) * S1      # combined grid spacing 2 pi / (2M): half-spacing pi/(2M)
        assert np.all(4 * delta <= one_m_2p ** 2 / 2)
        dp = 2 * delta / one_m_2p
        lip = (S_ + 2 * dp) ** g + 2 * g * (m_ + dp)
        bound = Tj + dp * lip + 1e-12
        k = int(np.argmax(bound))
        if bound[k] > worst[0]:
            worst = (float(bound[k]), (shift, k))
        maxT = max(maxT, float(Tj.max()))
        minsep = min(minsep, float(one_m_2p.min()))
        del p, one_m_2p, ap, aq, m_, S_, Tj, dp, lip, bound, wv
    return {"b1": float(b[1]), "S1": S1, "grid": 2 * M, "delta": eval_allow + (math.pi / (2 * M)) * S1,
            "max_T_grid": maxT, "max_T_bound": worst[0], "min|1-2p|": minsep, "certified": worst[0] < 1.0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gamma-num", type=int, default=377)
    ap.add_argument("--gamma-den", type=int, default=1000)
    ap.add_argument("--rr", type=float, default=0.99975)
    ap.add_argument("--M", type=int, default=1 << 15)
    ap.add_argument("--logM", type=int, default=22)
    ap.add_argument("--cert", type=Path, default=DEFAULT_CERT)
    ap.add_argument("--verify-only", type=Path, default=None)
    args = ap.parse_args()
    if args.verify_only:
        data = json.loads(args.verify_only.read_text())
        b = np.array([float.fromhex(h) for h in data["coefficients_hex"]])
        res = verify(b, data["gamma_num"] / data["gamma_den"], args.logM)
        print("verify-only:", res)
        return 0 if res["certified"] else 1
    g = args.gamma_num / args.gamma_den
    t0 = time.time()
    print(f"=== THM-4467 extension: subordination certificate at gamma = {args.gamma_num}/{args.gamma_den} ===")
    b, R, lres, nres = riemann_inverse_coeffs(g, args.rr, args.M)
    print(f"  lightning R(W,0) = {R:.9f} (boundary residual {lres:.1e}); inverse-map Newton residual {nres:.1e}; "
          f"rr = {args.rr}; degree {len(b) - 1}; b1 = {b[1]:.9f}  [{time.time() - t0:.0f}s]", flush=True)
    res = verify(b, g, args.logM)
    print(f"  grid 2x2^{args.logM}, delta = {res['delta']:.3e}, max T on grid = {res['max_T_grid']:.9f}, "
          f"max certified bound = {res['max_T_bound']:.9f}, min|1-2p| = {res['min|1-2p|']:.4f}, S1 = {res['S1']:.4f}")
    print(f"  certified: {res['certified']}  [{time.time() - t0:.0f}s]")
    if res["certified"]:
        data = {"gamma_num": args.gamma_num, "gamma_den": args.gamma_den,
                "description": "P(z)=sum_k b_k z^k; P(0)=0, b_1>1, P(|z|=1) in W_gamma (THM-4467 subordination)",
                "coefficients_hex": [float(x).hex() for x in b]}
        txt = json.dumps(data, separators=(",", ":")) + "\n"
        args.cert.parent.mkdir(parents=True, exist_ok=True)
        args.cert.write_text(txt)
        print(f"  certificate written: {args.cert.name} sha256={hashlib.sha256(txt.encode()).hexdigest()}")
        print(f"  => R(W_gamma,0) >= b1 = {b[1]:.9f} > 1  =>  C* >= 1 + {args.gamma_num}/{args.gamma_den} = {1 + g}")
    return 0 if res["certified"] else 1


if __name__ == "__main__":
    sys.exit(main())
