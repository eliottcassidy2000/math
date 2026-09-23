#!/usr/bin/env python3
"""AMM 12592 (Glazer's time-limited fair coin): a UNIFORM lower bound C* >= 1 + gamma_1 for every exactly
fair deterministic extractor, via Polya's capacity theorem applied to the p <-> 1-p quotient of the
spine functions.  procgen lane 2026-09-23.

Mathematical chain (proved in the companion note, sections 3-4):
  * spine normal form (THM-2966): F(p)=sum_m p^m q W_m(p) in Z[[p]], G(u)=sum_m u^m (1-u) V_m(1-u) in Z[[u]],
    F(p)+G(1-p)=1/2 on (0,1); Bernstein boxes give |W_m(p)| <= (|p|+|1-p|)^{d_m};
  * if T(n) <= (1+g) n + D then F converges on Omega0(g)={|p| (|p|+|1-p|)^g < 1} and continues to
    U_g = {min(|p|,|1-p|) (|p|+|1-p|)^g < 1}, simply connected and symmetric under s: p -> 1-p;
  * Delta=F-G is s-even, Sigma=F+G-1/2 is s-odd; phi(w)=Delta(p), psi(w)=Sigma(p)/(1-2p) (w=p(1-p))
    lie in Z[[w]] (resp. -1/2+Z[[w]]) and are analytic on W_g = w(U_g);
  * Polya (1928): an integer power series analytic on a domain of conformal radius > 1 is rational;
  * rational F contradicts F(p)+G(1-p)=1/2 by Fatou's lemma + Gauss's lemma (denominator parity).
  * log R(W_g,0) = V(U_g,0) + g_{U_g}(0,1) =: Lambda(g)  (2:1 proper map w: U_g -> W_g).
Hence: Lambda(g) > 0  ==>  no fair extractor has T(n) <= (1+g) n + D.

This script:
  [1] validates the Laplace solver against closed forms at g=0
      (R(U_0,0)=4 sqrt6/9, g_{U_0}(0,1)=(1/2)ln2, R(W_0,0)=8 sqrt3/9);
  [2] tabulates V(g), G(g), Lambda(g) (two independent solvers: on U, and directly on W);
  [3] bisects gamma_0 (one-point Polya, V=0) and gamma_1 (quotient, Lambda=0);
  [4] natural-boundary lemma: max_{W_g} |w| = r_pi (1+r_pi); W_g inside the unit disk iff g >= gamma*;
  [5] RIGOROUS two-disk certificate at g = 7/20 (mpmath interval arithmetic, closed form);
  [6] RIGOROUS subordination certificate at g = 3/8: explicit polynomial P, P(0)=0, P'(0)>1,
      P(unit circle) inside W_{3/8} (FFT grid + explicit Lipschitz/rounding bounds) => R(W_{3/8},0) > 1;
  [7] consistency: Lambda(gamma*) < 0 and Lambda(1) < 0 (constructions exist there);
  [8] parity sanity: p(w) = (1-sqrt(1-4w))/2 = sum_j w^(2^j) (mod 2).
Numerical parts [2]-[4] are floating point with a-posteriori residuals (NUMERICAL); [5],[6] are
certificates (interval arithmetic / explicit error bounds).
Usage: python3 amm12592_procgen_20260923_polya_capacity.py [--quick] [--cert PATH] [--verify-only]
"""
from __future__ import annotations
import argparse, hashlib, json, math, sys, time
from pathlib import Path
import numpy as np
import mpmath as mp

REPO = Path(__file__).resolve().parents[2]
DEFAULT_CERT = REPO / "05-knowledge" / "results" / "amm12592_procgen_20260923_subordination_certificate.json"
PHI = (1 + 5 ** 0.5) / 2
GSTAR = math.log(PHI) / math.log(5 ** 0.5)          # log_5(phi^2) = 0.59798743566544...


# ----------------------------------------------------------------------------- geometry
def r_boundary(theta: float, g: float) -> float:
    """t > 0 with t (t + |1 - t e^{i theta}|)^g = 1 (strictly increasing in t)."""
    lo, hi = 0.0, 1.0
    f = lambda t: t * (t + abs(1 - t * complex(math.cos(theta), math.sin(theta)))) ** g - 1
    while f(hi) < 0:
        hi *= 2
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def corner(g: float) -> complex:
    rc = 2.0 ** (-g / (1 + g))           # |p|=|1-p|=rc,  rc (2 rc)^g = 1
    return complex(0.5, math.sqrt(rc * rc - 0.25))


def left_arc(g: float, n_arc: int = 900) -> tuple[np.ndarray, complex]:
    """boundary of U_g with Re p <= 1/2: from corner a through -r_pi to conj(a) (clustered at corners)."""
    a = corner(g)
    thc = math.atan2(a.imag, a.real)
    u = np.linspace(0, 1, n_arc)
    u = np.unique(np.clip(np.concatenate([np.exp(-np.linspace(14, 0.5, 90)), u]), 0, 1))
    th = thc + (math.pi - thc) * u
    up = np.array([r_boundary(t, g) * np.exp(1j * t) for t in th])
    return np.concatenate([up, np.conj(up[::-1])]), a


def T_member(p: np.ndarray, g: float) -> np.ndarray:
    ap, aq = np.abs(p), np.abs(1 - p)
    return np.minimum(ap, aq) * (ap + aq) ** g


def T_member_w(w: np.ndarray, g: float) -> np.ndarray:
    p = (1 - np.sqrt(1 - 4 * np.asarray(w, dtype=complex))) / 2
    return T_member(p, g)


# ----------------------------------------------------------------------------- lightning Laplace solver
class Lightning:
    """harmonic h = Re H, H = polynomial + simple poles clustered at corners (Gopal-Trefethen style)."""

    def __init__(self, Z, data, corners, zc, sc, npol=50, nlight=30, sigma=4.0, Lfar=0.7):
        self.zc, self.sc, self.npol = zc, sc, npol
        self.poles = []
        for cp, direc in corners:
            for j in range(1, nlight + 1):
                d = Lfar * math.exp(-sigma * (math.sqrt(nlight) - math.sqrt(j)))
                self.poles.append((cp + direc * d, d))
        B = self.cbasis(Z)
        A = np.hstack([B.real, -B.imag])
        coef, *_ = np.linalg.lstsq(A, data, rcond=None)
        n = B.shape[1]
        self.c = coef[:n] + 1j * coef[n:]
        self.resid = float(np.max(np.abs(A @ coef - data)))

    def cbasis(self, z):
        z = np.atleast_1d(np.asarray(z, dtype=complex))
        cols = [((z - self.zc) / self.sc) ** k for k in range(self.npol + 1)]
        cols += [d / (z - pj) for pj, d in self.poles]
        return np.array(cols).T

    def H(self, z):
        return self.cbasis(z) @ self.c


def solve_U(g, npol=50, nlight=30):
    Zl, a = left_arc(g)
    Z = np.concatenate([Zl, 1 - Zl])
    lt = Lightning(Z, np.log(np.abs(Z)), [(a, 1j), (a.conjugate(), -1j)], 0.5, 1.3, npol, nlight, Lfar=0.6)
    V = lt.H(0.0)[0].real
    G = lt.H(1.0)[0].real
    return V, G, lt.resid


def solve_W(g, npol=50, nlight=30):
    Zl, a = left_arc(g)
    Wb = Zl * (1 - Zl)
    wc = (a * (1 - a)).real
    wmin = Wb.real.min()
    lt = Lightning(Wb, np.log(np.abs(Wb)), [(complex(wc, 0), 1.0)], (wc + wmin) / 2, (wc - wmin) / 2 * 1.1,
                   npol, nlight, Lfar=0.8)
    H0 = lt.H(0.0)[0]
    return H0.real, lt, wc


def bisect(fun, lo, hi, it=46):
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if fun(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ----------------------------------------------------------------------------- two-disk closed form
def two_disk_Lambda_float(c, rho):
    """V = D(c,rho) U D(1-c,rho):  log R(w(V),0) = log[ tan(A) (1/4+s^2) / (kappa s) ]."""
    s = math.sqrt(rho * rho - (0.5 - c) ** 2)
    b0 = math.atan(s / (0.5 + rho - c))
    kap = math.pi / (2 * math.pi - 4 * b0)
    A = 2 * kap * (math.atan(2 * s) - b0)
    return math.log(math.tan(A) * (0.25 + s * s) / (kap * s)), math.log(math.sin(A) * (0.25 + s * s) / (kap * s))


def two_disk_certificate(gnum, gden, cnum, cden, rnum, rden, dps=40):
    iv = mp.iv
    iv.dps = dps
    g = iv.mpf(gnum) / gden
    c = iv.mpf(cnum) / cden
    rho = iv.mpf(rnum) / rden
    one = iv.mpf(1)
    # (i) containment of the circle c + rho e^{it}, t in [0,pi] (conjugation symmetric), in Omega0(g);
    #     Omega0 is star-shaped w.r.t. 0 and 0 lies in the disk, so the circle suffices.
    todo = [(mp.mpf(0), mp.mpf(mp.pi) + mp.mpf(10) ** -30)]
    npieces, maxup = 0, mp.mpf(0)
    while todo:
        lo, hi = todo.pop()
        t = iv.mpf([lo, hi])
        x = c + rho * iv.cos(t)
        y = rho * iv.sin(t)
        az = iv.sqrt(x ** 2 + y ** 2)
        a1 = iv.sqrt((1 - x) ** 2 + y ** 2)
        val = az * iv.exp(g * iv.log(az + a1))
        if val.b < 1:
            npieces += 1
            maxup = max(maxup, mp.mpf(val.b))
        else:
            if hi - lo < mp.mpf(10) ** -12:
                return {"containment": False, "where": (float(lo), float(hi))}
            mid = (lo + hi) / 2
            todo += [(lo, mid), (mid, hi)]
    # 0 inside the disk, overlap of the two disks
    assert (c * c).b < (rho * rho).a
    assert (c + rho).a > (one - c - rho).b
    # (ii) Lambda for the two-disk union, interval evaluation of the closed form
    s = iv.sqrt(rho * rho - (iv.mpf(1) / 2 - c) ** 2)
    b0 = iv.atan2(s, iv.mpf(1) / 2 + rho - c)
    kap = iv.pi / (2 * iv.pi - 4 * b0)
    A = 2 * kap * (iv.atan2(2 * s, one) - b0)
    lam = iv.log(iv.tan(A) * (iv.mpf(1) / 4 + s * s) / (kap * s))
    return {"containment": True, "pieces": npieces, "max_upper_T_on_circle": float(maxup),
            "Lambda_interval": (float(lam.a), float(lam.b)), "Lambda_positive": bool(lam.a > 0)}


# ----------------------------------------------------------------------------- subordination certificate
def build_subordination_polynomial(g, rr=0.9985, M=4096, ncut=2048, npol=60, nlight=36):
    lam0, lt, wc = solve_W(g, npol, nlight)
    R = math.exp(lam0)
    Hn = lambda z: lt.H(z) - 1j * lt.H(0.0)[0].imag
    t = 2 * math.pi * np.arange(M) / M
    w = None
    for s in np.linspace(0, rr, 400)[1:]:
        z = s * np.exp(1j * t)
        if w is None:
            w = z * R
        for _ in range(8):
            Hw = Hn(w)
            h = 1e-7
            dH = (Hn(w + h) - Hn(w - h)) / (2 * h)
            w = w - (w * np.exp(-Hw) - z) / (np.exp(-Hw) * (1 - w * dH))
    coeff = np.fft.fft(w) / M
    b = coeff[:ncut].real.copy()
    b[0] = 0.0
    return b, R, lt.resid


def verify_subordination(b: np.ndarray, g: float, logM=21, eval_allow=1e-9):
    """Rigorous-by-bounds check that P(e^{it}) lies in W_g for all t, where P = sum_k b_k z^k (real b).
    Grid of M=2^logM points; every t is within pi/M of a grid point; |P(e^{it})-P(e^{it_j})| <= (pi/M) S1,
    S1 = sum k|b_k|; FFT evaluation error <= eval_allow (the standard FFT bound gives < 1e-10 here).
    For w within delta of w_j = P(e^{it_j}):  root displacement |dp| <= 2 delta/|1-2p_j|  (if 4delta<=|1-2p_j|^2/2),
    and |dT| <= |dp| (S_max^g + 2 g m_max) with S=|p|+|1-p|>=1, m=min(|p|,|1-p|)."""
    M = 1 << logM
    n = len(b)
    assert b[0] == 0.0 and b[1] > 1.0
    S1 = float(np.sum(np.arange(n) * np.abs(b))) * (1 + 1e-12)
    pad = np.zeros(M, dtype=complex)
    pad[:n] = b
    vals = np.fft.ifft(pad) * M                     # P(e^{2 pi i j/M}), j=0..M-1
    half = M // 2 + 1                                 # t in [0, pi] suffices (real coefficients, W conj-symmetric)
    wv = vals[:half]
    p = (1 - np.sqrt(1 - 4 * wv)) / 2
    one_m_2p = np.abs(1 - 2 * p)
    ap, aq = np.abs(p), np.abs(1 - p)
    m = np.minimum(ap, aq)
    S = ap + aq
    Tj = m * S ** g
    delta = eval_allow + (math.pi / M) * S1
    assert np.all(4 * delta <= one_m_2p ** 2 / 2)
    dp = 2 * delta / one_m_2p
    lip = (S + 2 * dp) ** g + 2 * g * (m + dp)
    bound = Tj + dp * lip + 1e-12
    return {"b1": float(b[1]), "S1": S1, "grid": M, "delta": delta, "max_T_grid": float(Tj.max()),
            "max_T_bound": float(bound.max()), "min|1-2p|": float(one_m_2p.min()),
            "certified": bool(bound.max() < 1.0)}


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--cert", type=Path, default=DEFAULT_CERT)
    ap.add_argument("--verify-only", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    print("=== AMM12592 procgen 2026-09-23: Polya-capacity lower bound for the uniform constant C* ===")
    print(f"gamma* = log_5(phi^2) = {GSTAR:.15f}   (C_* = {1 + GSTAR:.15f})")

    if args.verify_only:
        data = json.loads(args.cert.read_text())
        b = np.array([float.fromhex(h) for h in data["coefficients_hex"]])
        res = verify_subordination(b, data["gamma_num"] / data["gamma_den"])
        print("[6] verify-only:", res)
        return 0 if res["certified"] else 1

    # [1] validation at gamma = 0
    V0, G0, res0 = solve_U(0.0)
    L0 = solve_W(0.0)[0]
    exV, exG, exL = math.log(4 * math.sqrt(6) / 9), 0.5 * math.log(2), math.log(8 * math.sqrt(3) / 9)
    tdL, tdV = two_disk_Lambda_float(0.0, 1.0)
    print("\n[1] validation at gamma=0 (U_0 = D(0,1) u D(1,1)):")
    print(f"    V=log R(U,0): lightning {V0:.10f}  closed form log(4sqrt6/9) {exV:.10f}  two-disk formula {tdV:.10f}")
    print(f"    g_U(0,1)    : lightning {G0:.10f}  closed form (1/2)ln2      {exG:.10f}")
    print(f"    Lambda=log R(W,0): on-U {V0 + G0:.10f}  direct-on-W {L0:.10f}  closed form log(8sqrt3/9) {exL:.10f}  two-disk {tdL:.10f}")
    print(f"    boundary residual of the U-solve: {res0:.2e}")
    assert abs(V0 - exV) < 1e-6 and abs(G0 - exG) < 1e-6 and abs(L0 - exL) < 1e-6

    # [2] table
    print("\n[2] table: gamma, V=log R(U,0), G=g_U(0,1), Lambda=V+G (on U), Lambda (direct on W), r_pi")
    grid = [0.0, 0.05, 0.0955, 0.1, 0.2, 0.3, 0.35, 0.37, 0.375, 0.3775, 0.38, 0.4, 0.5, GSTAR, 0.7, 1.0]
    for g in grid:
        V, G, res = solve_U(g)
        Lw = solve_W(g)[0]
        print(f"    {g:8.5f}  V={V:+.8f}  G={G:.8f}  Lambda={V + G:+.8f}  direct={Lw:+.8f}  "
              f"r_pi={r_boundary(math.pi, g):.6f}  resid={res:.1e}")

    # [3] thresholds
    print("\n[3] thresholds (bisection):")
    for (npol, nl) in ([(50, 30)] if args.quick else [(40, 24), (50, 30), (70, 40)]):
        g0 = bisect(lambda g: solve_U(g, npol, nl)[0], 0.05, 0.15)
        g1 = bisect(lambda g: solve_W(g, npol, nl)[0], 0.30, 0.45)
        print(f"    npol={npol} nlight={nl}: gamma_0 (one-point Polya) = {g0:.9f}   gamma_1 (quotient) = {g1:.9f}"
              f"   => C* >= {1 + g1:.7f}")

    # [4] natural-boundary lemma
    print("\n[4] max_{closure W_g}|w| = r_pi(1+r_pi) (attained at p=-r_pi);  <= 1  iff  g >= gamma*:")
    for g in [0.3, 0.5, 0.59, GSTAR, 0.61, 0.8]:
        Zl, a = left_arc(g, 1500)
        wmax = float(np.max(np.abs(Zl * (1 - Zl))))
        rp = r_boundary(math.pi, g)
        print(f"    g={g:.6f}: sampled max|w|={wmax:.10f}  r_pi(1+r_pi)={rp * (1 + rp):.10f}")
    mp.mp.dps = 50
    ph = (1 + mp.sqrt(5)) / 2
    gs = mp.log(ph) / mp.log(mp.sqrt(5))
    print(f"    exact: (1/phi)(1+2/phi)^gamma* - 1 = {mp.nstr((1 / ph) * (1 + 2 / ph) ** gs - 1, 5)}  "
          f"(1+2/phi = sqrt5, sqrt5^gamma* = phi), and r_pi(1+r_pi)=1 at r_pi=1/phi")

    # [5] two-disk interval certificate
    print("\n[5] rigorous two-disk certificate (mpmath interval arithmetic):")
    best = None
    for c in np.linspace(0.0, 0.12, 25):
        # largest rho with the circle inside Omega0 (float scan) -- for reporting only
        lo, hi = 0.0, 1.2
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            tt = np.linspace(0, math.pi, 3000)
            z = c + mid * np.exp(1j * tt)
            if np.max(np.abs(z) * (np.abs(z) + np.abs(1 - z)) ** 0.35) < 1:
                lo = mid
            else:
                hi = mid
        L2, _ = two_disk_Lambda_float(c, lo)
        if best is None or L2 > best[0]:
            best = (L2, c, lo)
    print(f"    float scan at gamma=0.35: best two-disk Lambda = {best[0]:+.6f} at c={best[1]:.4f}, rho={best[2]:.6f}")
    cert = two_disk_certificate(7, 20, 13, 200, 39, 50)
    print(f"    gamma=7/20, c=13/200, rho=39/50: {cert}")
    assert cert["containment"] and cert["Lambda_positive"]
    print("    => R(W_{7/20},0) >= R(w(D(c,rho) u D(1-c,rho)),0) > 1  =>  C* >= 27/20 = 1.35 (PROVED, given Polya/Fatou)")

    # [6] subordination certificate at gamma = 3/8
    print("\n[6] subordination certificate at gamma = 3/8:")
    g = 3 / 8
    b, R, resid = build_subordination_polynomial(g)
    res = verify_subordination(b, g)
    print(f"    Riemann-map estimate R(W_3/8,0) = {R:.8f} (lightning residual {resid:.1e})")
    print(f"    certificate polynomial: degree {len(b) - 1}, P(0)=0, P'(0)=b1={res['b1']:.8f}, sum k|b_k|={res['S1']:.6f}")
    print(f"    grid 2^{int(math.log2(res['grid']))}, delta={res['delta']:.3e}, max T on grid={res['max_T_grid']:.8f}, "
          f"max certified bound={res['max_T_bound']:.8f}, min|1-2p|={res['min|1-2p|']:.4f}")
    print(f"    certified P(unit circle) in W_3/8: {res['certified']}")
    assert res["certified"]
    data = {"gamma_num": 3, "gamma_den": 8, "description": "P(z)=sum_k b_k z^k; P(0)=0, b_1>1, P(|z|=1) in W_gamma",
            "coefficients_hex": [float(x).hex() for x in b]}
    txt = json.dumps(data, separators=(",", ":")) + "\n"
    args.cert.parent.mkdir(parents=True, exist_ok=True)
    args.cert.write_text(txt)
    print(f"    certificate written: {args.cert.name}  sha256={hashlib.sha256(txt.encode()).hexdigest()}")
    res2 = verify_subordination(np.array([float.fromhex(h) for h in data['coefficients_hex']]), 3 / 8)
    assert res2["certified"]
    print("    => by Schwarz's lemma R(W_{3/8},0) >= P'(0) > 1  =>  C* >= 11/8 = 1.375 (PROVED, given Polya/Fatou)")

    # [7] consistency with constructions
    print("\n[7] consistency (the obstruction must NOT fire where constructions exist):")
    for gg, name in [(GSTAR, "gamma* (Long / golden block construction)"), (1.0, "gamma=1 (classical 2n)")]:
        V, G, _ = solve_U(gg)
        print(f"    {name}: V={V:+.6f}, Lambda={V + G:+.6f}  -> both negative: {V < 0 and V + G < 0}")
        assert V < 0 and V + G < 0

    # [8] parity of p(w)
    print("\n[8] p(w) = sum_n Cat(n-1) w^n ; Cat(n-1) odd  iff  n is a power of two (checked n <= 4096):")
    ok = True
    c = 1  # Cat(0)
    for n in range(1, 4097):
        # Cat(n-1) parity via Kummer: odd iff n is a power of 2
        odd = (math.comb(2 * n - 2, n - 1) // n) % 2 == 1
        ok &= (odd == (n & (n - 1) == 0))
    print(f"    verified: {ok}   => every complement-symmetric fair extractor has phi(w) = sum_j w^(2^j) (mod 2)")
    print(f"\nAll checks passed. [{time.time() - t0:.1f}s]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
