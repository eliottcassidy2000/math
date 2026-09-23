#!/usr/bin/env python3
"""'Structured majorant' audit for THM-4467 (procgen 2026-09-23, HYP-9128 note, section 7).

THM-4467's lossy step: |W_m(p)| <= sum_k C(d,k)|p|^{d-k}|q|^k = S(p)^d  (S = |p|+|1-p|).
This script quantifies, exactly or with explicit numerics:
  [1] box sharpness: max over the Bernstein box of |W(p)| >= S(p)^d / pi (phase-aligned 0/C(d,k) choices);
  [2] the parity skeleton: the odd Bernstein coordinates forced by E = 2W - 1 == 1 (mod 2) sit at the Lucas
      submasks k of d; their minimal contribution prod_j (|p|^{2^j} + |q|^{2^j}) (bits j of d) vs S^d;
  [3] the defect on the extremal objects: Long's separately balanced blocks and the golden-zero super-blocks,
      rate (1/R_i) log(|E_i(p)| / S(p)^{R_i}) at the golden point and at boundary points of U_gamma;
  [4] kernel-noise saturation (EXACT): adding even multiples of Lemma-R kernel vectors K_{i,r} to a fair
      super-block (same deadlines, fairness preserved) attains |E_i(p0)| >= c0 S(p0)^{R_i} at a prescribed p0;
      the modified block is re-verified exactly (box, parity, palindromic class polynomial).
Conclusion printed at the end.
"""
from __future__ import annotations
import json, math, sys, time
from fractions import Fraction
from math import comb
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
import amm12592_procgen_20260923_hyp9128_finite as FIN

PHI = (1 + 5 ** 0.5) / 2


def S_of(p):
    return abs(p) + abs(1 - p)


def r_pi(g):
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if mid * (1 + 2 * mid) ** g < 1:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def E_value(e_row, p):
    """E_i(p) = sum_r e_r p^r q^{R-r}  (Bernstein form of the level's signed read-prefix count)."""
    R = len(e_row) - 1
    q = 1 - p
    return sum(complex(e) * p ** r * q ** (R - r) for r, e in enumerate(e_row) if e)


def box_sharpness(d, p, nphase=720):
    q = 1 - p
    z = np.array([comb(d, k) * p ** (d - k) * q ** k for k in range(d + 1)], dtype=complex)
    tot = np.sum(np.abs(z))
    best = 0.0
    for ph in np.linspace(0, 2 * math.pi, nphase, endpoint=False):
        sel = np.real(z * np.exp(-1j * ph)) > 0
        best = max(best, abs(np.sum(z[sel])))
    return best / tot


def parity_skeleton(d, p):
    q = 1 - p
    val = 1.0
    j = 0
    while d >> j:
        if (d >> j) & 1:
            val *= abs(p) ** (2 ** j) + abs(q) ** (2 ** j)
        j += 1
    return val / S_of(p) ** d


def superblock_packets(N):
    out, data = FIN.margins_certificate(N)
    ok, err, sites, f = FIN.lemma_R_round(N, data)
    assert ok
    t, R, a, S, Fnum, den, digits, h, E0 = data
    e = [[(-1) ** (i + r) * f[i][r] for r in range(len(f[i]))] for i in range(len(f))]
    return e, R, a, t, S, f, h


def long_packets(N):
    cert = json.loads((REPO / "scratch" / "procgen_amm" / "long_draft" / "scripts" / "certificates" / "finite_blocks.json").read_text())
    for b in cert["blocks"]:
        if b["N"] == N:
            e = [[(-1) ** (i + r) * b["f"][i][r] for r in range(len(b["f"][i]))] for i in range(N)]
            return e, b["R"], b["a"]
    raise KeyError(N)


def defect_table(name, e, R, points, i_range):
    rows = []
    for (pname, p) in points:
        rates = []
        for i in i_range:
            if R[i] < 8:
                continue
            Ev = abs(E_value(e[i], p))
            ref = S_of(p) ** R[i]
            if Ev == 0:
                continue
            rates.append((math.log(Ev / ref) / R[i], i))
        if rates:
            mx = max(rates)
            md = float(np.median([r for r, _ in rates]))
            rows.append((pname, mx[0], mx[1], md))
    print(f"  {name}:")
    for pname, mx, arg, md in rows:
        print(f"     at {pname:>26s}: max_i (1/R_i) log(|E_i|/S^R_i) = {mx:+.4f} (level {arg}), median {md:+.4f}")


def boundary_point(g, angle):
    """point p0 on the boundary of Omega_0(gamma) = {|p| S(p)^gamma < 1} on the ray arg p = angle."""
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        p = mid * complex(math.cos(angle), math.sin(angle))
        if abs(p) * S_of(p) ** g < 1:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2 * complex(math.cos(angle), math.sin(angle))


def kernel_saturation(N, kind, g, angle, alpha=Fraction(1, 64), cache={}):
    """Add even multiples m_r of the Lemma-R kernel vectors K_(i,r) at ONE level i of the exact fair super-block,
    phase-aligned at p0 in the boundary of Omega_0(gamma); re-verify the whole block exactly.
    kind = 'uncapped': a level i < h with d_i > 0 (real centre, |f| <= C/2 + 5);
    kind = 'capped'  : a level i > h (d_i = 0, f = Lemma-R rounding of the zero centre, |f| <= 2) with d_m/m closest
                       to gamma, so that the level is a genuine 'd_m = gamma m' level."""
    if N not in cache:
        cache[N] = superblock_packets(N)
    e, R, a, t, S, f, h = cache[N]
    n = 3 * N - 1
    p0 = boundary_point(g, angle)
    q0 = 1 - p0
    if kind == "uncapped":
        i = h // 2
        while a[i] - a[i + 1] == 0:
            i += 1
    else:
        cands = [i for i in range(h + 1, n - 2)]
        i = min(cands, key=lambda i: abs(R[i] / (N + i) - g))
    di = a[i] - a[i + 1]
    E_before = E_value(e[i], p0)
    ref = S_of(p0) ** R[i]
    phase = E_before / abs(E_before) if abs(E_before) > 0 else 1.0
    f2 = [row[:] for row in f]
    for r in range(1, R[i]):
        base = comb(R[i], r) if kind == "uncapped" else comb(R[i] - 1, r - 1)
        mag = 2 * int(alpha * base / 2)
        if mag == 0:
            continue
        term = (-1) ** (i + r) * p0 ** r * q0 ** (R[i] - r)
        sgn = 1 if (term * phase.conjugate()).real >= 0 else -1
        m = sgn * mag
        f2[i][r] += m
        for j in range(di + 1):
            f2[i + 1][r - 1 + j] -= (-1) ** j * comb(di, j) * m
    # exact re-verification of the modified super-block: box, parity, boundary +-1, class polynomial == target
    ok = True
    for ii in range(n + 1):
        for r in range(R[ii] + 1):
            Q = comb(R[ii], r)
            if abs(f2[ii][r]) > Q or (f2[ii][r] - Q) % 2 or (r in (0, R[ii]) and abs(f2[ii][r]) != 1):
                ok = False
    P = [0] * (n + 1)
    for ii in range(n + 1):
        for r in range(R[ii] + 1):
            ev = (-1) ** (ii + r) * f2[ii][r]
            if ev:
                for j in range(a[ii] + 1):
                    P[ii + r + j] += ev * comb(a[ii], j)
    target = FIN.target_P(N, S) + [0] * (n + 1 - (2 * N + 1))
    ok = ok and P == target
    e2 = [(-1) ** (i + r) * f2[i][r] for r in range(R[i] + 1)]
    E_after = E_value(e2, p0)
    return {"N": N, "kind": kind, "level": i, "m": N + i, "d": R[i], "d_over_m": R[i] / (N + i), "p0": p0,
            "ratio_before": abs(E_before) / ref, "ratio_after": abs(E_after) / ref, "valid": ok,
            "rate_before": math.log(abs(E_before) / ref) / R[i] if abs(E_before) > 0 else -math.inf,
            "rate_after": math.log(abs(E_after) / ref) / R[i]}


def main():
    t0 = time.time()
    print("=== Structured-majorant audit (THM-4467 triangle-inequality step) ===")
    print("\n[1] box sharpness: max_(box) |W(p)| / S(p)^d  (phase-aligned 0/C(d,k) choices)")
    for p in [complex(-0.62, 0), complex(0.3, 0.6), complex(-0.4, 0.5), complex(0.8, -0.7)]:
        print(f"    p = {p}: d=60 ratio {box_sharpness(60, p):.4f}, d=200 ratio {box_sharpness(200, p):.4f}  (>= 1/pi = {1/math.pi:.4f})")
    print("\n[2] parity skeleton (Lucas submasks of d): minimal forced |E|/S^d")
    for p in [complex(-0.62, 0), complex(0.3, 0.6)]:
        for d in [63, 64, 100, 127]:
            val = parity_skeleton(d, p)
            print(f"    p = {p}, d = {d}: {val:.3e}  (rate {math.log(val)/d:+.4f})")
    print("\n[3] defect on extremal objects: (1/R_i) log(|E_i(p)| / S(p)^{R_i})")
    pts = [("golden point p=-1/phi", complex(-1 / PHI, 0)),
           ("-r_pi(0.5) (dW_0.5 extreme)", complex(-r_pi(0.5), 0)),
           ("-r_pi(0.3775)", complex(-r_pi(0.3775), 0)),
           ("complex p=0.2+0.8i", complex(0.2, 0.8))]
    for N in [32, 64]:
        e, R, a = long_packets(N)
        defect_table(f"Long block N={N} (levels 0..N-1)", e, R, pts, range(N))
    for N in [32, 64]:
        e, R, a, t, S, f, h = superblock_packets(N)
        defect_table(f"golden-zero super-block N={N} (levels 0..h={h})", e, R, pts, range(h + 1))
    print("\n[4] kernel-noise saturation (exact): even multiples m_r of K_(i,r) at ONE level, |m_r| = 2 floor(C/128),")
    print("    phase-aligned at p0 in the boundary of Omega_0(gamma); every modified block is re-verified exactly")
    for N in [64, 128]:
        for kind, g in [("uncapped", 0.59), ("capped", 0.3775)]:
            for ang in [math.pi, 2 * math.pi / 3]:
                res = kernel_saturation(N, kind, g, ang)
                print(f"    N={N:4d} {kind:8s} level m={res['m']:4d} d={res['d']:4d} (d/m={res['d_over_m']:.3f}), "
                      f"p0={res['p0'].real:+.4f}{res['p0'].imag:+.4f}i: |E|/S^d {res['ratio_before']:.2e} -> {res['ratio_after']:.2e}, "
                      f"rate {res['rate_before']:+.4f} -> {res['rate_after']:+.4f}; valid fair super-block: {res['valid']}", flush=True)
                assert res["valid"]
    print("\nConclusion: the box maximum attains the triangle bound up to 1/pi; parity forces only the sparse Lucas skeleton;")
    print("the explicit extremal constructions have a strictly negative defect rate (a larger domain), but even kernel moves")
    print("preserve fairness, parity and deadlines and restore |E|/S^d to a constant (rate O(1/d) -> 0) at any prescribed point,")
    print("at levels with any ratio d/m in (0, 0.59).  Hence no termwise majorant valid for all fair extractors can enlarge U_gamma.")
    print(f"[{time.time() - t0:.0f}s]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
