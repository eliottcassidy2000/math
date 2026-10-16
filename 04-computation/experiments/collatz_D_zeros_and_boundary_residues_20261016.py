#!/usr/bin/env python3
"""Zeros of D(s) = 1 + 2^(-s) + (3/2)^s/3 in the strip 1/2 < Re s < 1 up to height T, and floating-point estimates of
the residues of F(s) = sum t(n) n^(-s) at all of them from the finite boundary numerator of Theorem B3.

(1) Zeros: vectorised Newton iteration from a lattice of starting points, deduplicated, then counted against the
    argument principle on rectangles [0.5, 1] x [kH, (k+1)H] (exact integer winding numbers; NUMERICAL).
(2) Residues: for each zero rho, r_hat(rho) = -(3/2)^rho numer(rho)/D'(rho), where numer is the boundary numerator of
    collatz_pole_certificate_general_N_20261015.py evaluated in complex128 from the exact signed block moments through
    N = 2^LOGN (first-order Taylor expansion per 2^17 block, as in the certificate), WITHOUT the future tail and the
    quadrature error, which the Arb certificate bounds.  These are point estimates; the conditional intervals are in
    the S37 certificate outputs.  Also reported: |numer(sigma + i t)| on vertical lines (growth in t), and the
    amplitudes |2 r_hat / rho| that the modes contribute to S(N)/N^(Re rho).
Session opus-2026-10-16-S38.  Usage: python collatz_D_zeros_and_boundary_residues_20261016.py moments.npz census.json
       T out.json
"""
import sys, json, math, time
import numpy as np

MOM, CEN, T, OUT = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
L2, L15 = math.log(2.0), math.log(1.5)


def D(s):
    return 1 + np.exp(-s * L2) + np.exp(s * L15) / 3


def Dp(s):
    return -L2 * np.exp(-s * L2) + L15 * np.exp(s * L15) / 3


def newton(s, it=60):
    for _ in range(it):
        s = s - D(s) / Dp(s)
    return s


def find_zeros(T):
    re0 = np.array([0.5005, 0.502, 0.505, 0.51, 0.52, 0.58, 0.64, 0.7, 0.76, 0.82, 0.88, 0.94, 0.98, 0.995, 0.9995])
    im0 = np.arange(0.05, T + 1, 0.08)
    S0 = (re0[:, None] + 1j * im0[None, :]).ravel()
    Z = newton(S0)
    ok = np.isfinite(Z) & (np.abs(D(Z)) < 1e-12) & (Z.real > 0.5) & (Z.real < 1.0) & (Z.imag > 0) & (Z.imag < T)
    Z = Z[ok]
    Z = Z[np.lexsort((Z.real, np.round(Z.imag, 6)))]
    out = []
    for z in Z:
        if not out or abs(z - out[-1]) > 1e-7:
            if all(abs(z - w) > 1e-7 for w in out[-20:]):
                out.append(z)
    out = np.array(sorted(out, key=lambda z: z.imag))
    return out


def winding(s0, s1, t0, t1, n=4000):
    a = np.linspace(s0, s1, n) + 1j * t0
    b = s1 + 1j * np.linspace(t0, t1, n)
    c = np.linspace(s1, s0, n) + 1j * t1
    d = s0 + 1j * np.linspace(t1, t0, n)
    path = np.concatenate([a, b, c, d])
    ph = np.unwrap(np.angle(D(path)))
    return (ph[-1] - ph[0]) / (2 * np.pi)


if __name__ == '__main__':
    t0 = time.time()
    Z = find_zeros(T)
    # argument-principle count on [0.5, 1] x [k*500, (k+1)*500]
    counts = []
    mism = 0
    for k in range(int(math.ceil(T / 500))):
        a, b = k * 500 + 1e-3, min((k + 1) * 500, T) - 1e-3
        w = winding(0.5 + 1e-9, 1.0 - 1e-9, a, b, n=int(20 * (b - a)) + 2000)
        c = int(((Z.imag > a) & (Z.imag < b)).sum())
        counts.append((a, b, round(w, 3), c))
        if abs(w - c) > 0.01:
            mism += 1
    # dense local search in every rectangle whose count disagrees with its winding number
    if mism:
        extra = []
        for (a, b, w, c) in counts:
            if abs(w - c) > 0.01:
                S0 = (np.linspace(0.5005, 0.9995, 200)[:, None] + 1j * np.arange(a, b, 0.02)[None, :]).ravel()
                Zl = newton(S0, it=100)
                okl = np.isfinite(Zl) & (np.abs(D(Zl)) < 1e-11) & (Zl.real > 0.5) & (Zl.real < 1) & (Zl.imag > a) & (Zl.imag < b)
                for z in Zl[okl]:
                    if all(abs(z - w_) > 1e-7 for w_ in list(Z) + extra):
                        extra.append(z)
        Z = np.array(sorted(list(Z) + extra, key=lambda z: z.imag))
        mism = 0
        counts2 = []
        for (a, b, w, c) in counts:
            c2 = int(((Z.imag > a) & (Z.imag < b)).sum())
            counts2.append((a, b, w, c2))
            if abs(w - c2) > 0.01:
                mism += 1
        counts = counts2
    print(f'{len(Z)} zeros with 1/2 < Re < 1, 0 < Im < {T:.0f} ({time.time() - t0:.1f} s); argument-principle rectangles '
          f'with a count mismatch: {mism} of {len(counts)}', flush=True)
    d = np.load(MOM)
    v = d['moments'].astype(np.float64)
    N = int(d['N'])
    B = int(d['B'])
    cen = json.load(open(CEN))
    tN = int(cen['t_N'])
    L = (2 * N - 1) // 3
    J = L // 2
    start_long = (J + 1 + B - 1) // B
    start_short = (L + 1 + B - 1) // B
    jj = np.arange(start_long, len(v))
    lo = jj.astype(np.float64) * B
    center = (2 * lo + B - 1) / 2
    sm = v[jj, 0]
    mu1 = (2 * v[jj, 1] - (B - 1) * sm) / 2
    lc = np.log(center)
    short_mask = jj >= start_short

    def numer(s):
        pw = np.exp(-s * lc)
        val = pw * (sm - s * mu1 / center)
        last = tN * np.exp(-s * math.log(N))
        longsum = val.sum() + last
        shortsum = val[short_mask].sum() + last
        return np.exp(s * math.log(2 / 3)) * (shortsum + np.exp(-s * L2) * longsum)

    rows = []
    for z in Z:
        nm = numer(z)
        r = -np.exp(z * L15) * nm / Dp(z)
        rows.append(dict(re=float(z.real), im=float(z.imag), r_re=float(r.real), r_im=float(r.imag), abs_r=float(abs(r)),
                         amp_in_S=float(2 * abs(r) / abs(z)), abs_numer=float(abs(nm)), abs_Dp=float(abs(Dp(z)))))
    print(f'residue estimates at {len(rows)} zeros from block moments through 2^{int(math.log2(N))} ({time.time() - t0:.1f} s)',
          flush=True)
    lines = {}
    for sig in (0.55, 0.7, 0.85, 1.0):
        ts = np.arange(10.0, T, 37.0)
        vals = np.array([abs(numer(sig + 1j * t)) for t in ts])
        q = [float(np.quantile(vals[(ts >= a) & (ts < b)], 0.5)) for a, b in ((10, 1000), (1000, 4000), (4000, 10000), (10000, T))
             if ((ts >= a) & (ts < b)).any()]
        lines[str(sig)] = dict(median_abs_numer_by_height_band=q, max=float(vals.max()))
        print(f'|numer({sig} + it)| median over t in [10,1e3), [1e3,4e3), [4e3,1e4), [1e4,T): {[round(x, 4) for x in q]}', flush=True)
    top = sorted(rows, key=lambda r: -r['amp_in_S'] * (2.0 ** (40 * (r['re'] - 1))))[:12]
    print('largest contributions to S(2^40)/2^40 (|2r/rho| 2^(40(Re rho - 1))):')
    for r in top:
        print(f"   {r['re']:.5f} + {r['im']:9.3f}i  |r| = {r['abs_r']:.4g}  |2r/rho| = {r['amp_in_S']:.4g}  "
              f"weight at 2^40: {r['amp_in_S'] * 2.0 ** (40 * (r['re'] - 1)):.4g}", flush=True)
    json.dump(dict(T=T, n_zeros=len(rows), argument_principle=counts, N=N, zeros_and_residues=rows, vertical_lines=lines),
              open(OUT, 'w'), indent=1)
