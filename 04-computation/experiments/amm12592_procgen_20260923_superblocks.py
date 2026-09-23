#!/usr/bin/env python3
"""AMM 12592 construction side: cross-annulus cancellation inside SUPER-BLOCKS.  procgen lane 2026-09-23.

A ratio-B super-block is the set of critical values [N, BN) decided by the common horizon BN and balanced
(heads = tails in every Hamming layer) as a WHOLE, not annulus by annulus.  B = 2 is the separately
balanced class (THM-3007/3009, Long's draft).  With Long's complement pairing 0^L 1 w <-> 1^L 0 wbar,
balance is equivalent to: the 0-spine class polynomial
      P_+(x) = sum_i x^i (1+x)^{a_i} E_i(x),  |[x^r]E_i| <= C(R_i,r),  [x^r]E_i = C(R_i,r) (mod 2),
is PALINDROMIC of degree (B-2)N with P_+(0) = +-1.  Equivalently the super-block's 0-spine deviation is
w^N S_N(w), w = pq, S_N in Z[w], deg S_N <= (B-2)N/2, S_N(0) = +-1: the signed HANDOFF STATE.
(B = 2 forces S_N = +-1; that is exactly why separately balanced schemes are lacunary and golden.)

Parts
 [A] exact CP-SAT: rho_B(N) = min over super-blocks of max_{N<=L<BN} T(L)/L (bisection over ratios t/L;
     monotone feasibility; UNKNOWN = solver timeout, reported, never counted as feasible);
 [B] the handoff states S_N of the optima (all vanish at the golden point w=-1);
 [C] asymptotic NECESSARY thresholds for the kappa-families S_N ~ (1+w)^(kappa N): natural boundary (NB)
     and single-block evaluation (EV) inequalities (grid computation, NUMERICAL);
 [D] exact fold thresholds: Long's real-center recursion with target F_0(u) = ((1+u)/2)^(N-1) S((1-u^2)/4),
     sufficient l1 criterion (levels >= 1) + exact Krawtchouk box at level 0, exact integers;
 [E] INTEGER super-blocks: fold center + Long's lattice rounding (his function, imported unmodified from the
     fetched draft), re-verified by independent exact code; plus an exact rational certificate (Long's
     evaluation inequality) that NO separately balanced block [N,2N) reaches the same ratio;
 [F] the chained ILP family sum_j w^(4^j) S_(4^j)(w), S_N = -(1+w)(1-w^N)/(1-w), equals
     -((1+w)/(1-w)) sum_k (-1)^k w^(2^k): a Hadamard gap series, so that family cannot beat gamma*.
Usage: python3 amm12592_procgen_20260923_superblocks.py [--quick] [--tlimit 120] [--draft DIR]
"""
from __future__ import annotations
import argparse, importlib.util, math, sys, time
from fractions import Fraction
from math import comb
from pathlib import Path

import numpy as np
import sympy as sp
from ortools.sat.python import cp_model

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DRAFT = REPO / "scratch" / "procgen_amm" / "long_draft"
GSTAR = math.log((1 + 5 ** 0.5) / 2) / math.log(5 ** 0.5)


# ============================================================================= [A] exact CP-SAT
def profile_floor(N, B, C):
    return [min(B * N, max(N + i + 1, math.floor(C * (N + i)))) for i in range((B - 1) * N)]


def profile_ceil(N, B, C):
    return [min(B * N, max(N + i + 1, math.ceil(C * (N + i)))) for i in range((B - 1) * N)]


def solve_block(N, B, t, tlimit=120.0, workers=2):
    H, D = B * N, (B - 2) * N
    ntail = H - N - 1
    m = cp_model.CpModel()
    coef = [dict() for _ in range(ntail + 1)]
    const = [0] * (ntail + 1)
    ev = {}
    for i in range((B - 1) * N):
        L = N + i
        Ri, ai = t[i] - L - 1, H - t[i]
        for r in range(Ri + 1):
            Q = comb(Ri, r)
            par = Q % 2
            y = m.NewIntVar((-Q - par) // 2, (Q - par) // 2, f"y_{i}_{r}")
            ev[(i, r)] = (y, par)
            for j in range(ai + 1):
                l = i + r + j
                cf = comb(ai, j)
                coef[l][y] = coef[l].get(y, 0) + 2 * cf
                const[l] += par * cf
    ex = lambda l: sum(c * v for v, c in coef[l].items()) + const[l]
    for l in range(ntail + 1):
        if l > D:
            m.Add(ex(l) == 0)
        elif l < D - l:
            m.Add(ex(l) - ex(D - l) == 0)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = tlimit
    s.parameters.num_workers = workers
    st = s.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return "FEASIBLE", {k: 2 * s.Value(v) + par for k, (v, par) in ev.items()}
    return ("INFEASIBLE" if st == cp_model.INFEASIBLE else "UNKNOWN"), None


def class_poly(N, B, t, e_of):
    """e_of[(i,r)] = e_{i,r} (0-spine packet coefficients, x-convention)."""
    H = B * N
    c = [0] * (H - N)
    for (i, r), val in e_of.items():
        L = N + i
        Ri, ai = t[i] - L - 1, H - t[i]
        Q = comb(Ri, r)
        assert abs(val) <= Q and (val - Q) % 2 == 0, ("box/parity", i, r)
        for j in range(ai + 1):
            c[i + r + j] += val * comb(ai, j)
    return c


def is_balanced(N, B, c):
    D = (B - 2) * N
    return all(c[l] == 0 for l in range(D + 1, len(c))) and all(c[l] == c[D - l] for l in range(D + 1)) and abs(c[0]) == 1


def min_ratio(N, B, tlimit):
    H = B * N
    cands = sorted({Fraction(tt, L) for L in range(N, H) for tt in range(L + 1, H + 1)})
    lo, hi = 0, len(cands) - 1
    st, sol = solve_block(N, B, profile_floor(N, B, cands[hi]), tlimit)
    best, log, infeas = (cands[hi], sol), [], None
    while lo < hi:
        mid = (lo + hi) // 2
        t = profile_floor(N, B, cands[mid])
        st, sol = solve_block(N, B, t, tlimit)
        log.append((str(cands[mid]), st))
        if sol is not None:
            assert is_balanced(N, B, class_poly(N, B, t, sol))
            hi, best = mid, (cands[mid], sol)
        else:
            if st == "INFEASIBLE":
                infeas = cands[mid] if infeas is None else max(infeas, cands[mid])
            lo = mid + 1
    C, sol = best
    t = profile_floor(N, B, C)
    return max(Fraction(t[i], N + i) for i in range(len(t))), t, sol, log, infeas


# ============================================================================= [B] handoff state
def state_polynomial(N, B, c):
    D = (B - 2) * N
    p, w = sp.symbols("p w")
    q = 1 - p
    assert all(ci == 0 for ci in c[D + 1:])
    Sp = sp.expand(sum(ci * p ** l * q ** (D - l) for l, ci in enumerate(c[:D + 1])))
    deg = D // 2
    s = sp.symbols(f"s0:{deg + 1}")
    eqs = sp.Poly(sp.expand(sum(s[k] * (p * q) ** k for k in range(deg + 1)) - Sp), p).all_coeffs()
    sol_s = sp.solve(eqs, s, dict=True)[0]
    S = sp.expand(sum(sol_s[s[k]] * w ** k for k in range(deg + 1)))
    return sp.factor(S), S.subs(w, -1)


# ============================================================================= [C] asymptotic necessary thresholds
def kappa_thresholds(kappas, grid=401):
    import cmath
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from amm12592_procgen_20260923_polya_capacity import left_arc
    X, Y = np.meshgrid(np.linspace(-3, 3, grid), np.linspace(-3, 3, grid))
    Wg = (X + 1j * Y).ravel()
    Wg = Wg[np.abs(Wg) > 1e-6]
    p1 = (1 - np.sqrt(1 - 4 * Wg + 0j)) / 2
    p2 = 1 - p1
    vs = np.linspace(1, 4, 241)

    def rhs(p, g):
        lp, lS = np.log(np.abs(p)), np.log(np.abs(p) + np.abs(1 - p))
        d = np.minimum(g * vs, 4 - vs)
        return np.max(vs[None, :] * lp[:, None] + d[None, :] * lS[:, None], axis=1)

    def lhs(w, k):
        return np.log(np.abs(w)) + k * np.log(np.abs(1 + w) + 1e-300) + (1 - k) * np.maximum(0, np.log(np.abs(w)))

    def ev_ok(g, k):
        return np.all(lhs(Wg, k) <= np.minimum(rhs(p1, g), rhs(p2, g)) + 1e-12)

    def nb_ok(g, k):
        Zl, _ = left_arc(g, 700)
        return np.all(lhs(Zl * (1 - Zl), k) < 0)

    def thr(test, k):
        lo, hi = 0.3, 1.2
        for _ in range(24):
            mid = 0.5 * (lo + hi)
            lo, hi = (lo, mid) if test(mid, k) else (mid, hi)
        return hi

    return [(k, thr(nb_ok, k), thr(ev_ok, k)) for k in kappas]


# ============================================================================= [D] exact fold with target
def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return out


def kraw(R, t):
    c, prev = [1], 0
    for r in range(R):
        nxt, rem = divmod((R - 2 * t) * c[-1] - (R - r + 1) * prev, r + 1)
        assert rem == 0
        prev = c[-1]
        c.append(nxt)
    return c


def family(N, j):
    """S = (1+w)^(2^j) (1-w^N)/(1-w^(2^j)); j<0 gives S=1.  S(0)=1, deg<=N, S = 1+w^N (mod 2) for j>=0."""
    if j < 0:
        return [1]
    m = 2 ** j
    return poly_mul([comb(m, k) for k in range(m + 1)], [1 if k % m == 0 else 0 for k in range(N - m + 1)])


def F0_numer(N, S):
    d = len(S) - 1
    tot, pw = [0], [1]
    for k in range(d + 1):
        if S[k]:
            term = [S[k] * 4 ** (d - k) * v for v in pw]
            tot += [0] * max(0, len(term) - len(tot))
            for jj, v in enumerate(term):
                tot[jj] += v
        pw = poly_mul(pw, [1, 0, -1])
    return poly_mul([comb(N - 1, k) for k in range(N)], tot), 2 ** (N - 1) * 4 ** d


def fold(N, Fnum, den, t, want_center=False):
    """Long's fold with general target.  Returns (max_{i>=1} l1 mass, level-0 box ratio, centers or None)."""
    F = np.array(Fnum, dtype=object)
    R = [t[i] - (N + i) - 1 for i in range(len(t))]
    worstM, box0, centers = Fraction(0), Fraction(0), []
    for i in range(len(R)):
        nz = np.nonzero(F != 0)[0]
        Ri = R[i]
        if len(nz) == 0:
            if want_center:
                centers.append([Fraction(0)] * (Ri + 1))
            continue
        deg = int(nz[-1])
        F = F[:deg + 1]
        if deg <= Ri:
            W, F = F.copy(), np.array([0], dtype=object)
        else:
            W = np.concatenate([F[:Ri], np.array([F[Ri:].sum()], dtype=object)])
            suf = np.cumsum(F[:Ri:-1])[::-1]
            F = np.concatenate([np.zeros(Ri, dtype=object), -2 * suf])
        if i == 0 or want_center:
            e = [0] * (Ri + 1)
            for tt, v in enumerate(W):
                if v:
                    kc = kraw(Ri, tt)
                    for r in range(Ri + 1):
                        e[r] += v * kc[r]
            if i == 0 and Ri >= 1:
                box0 = max(Fraction(abs(e[r]), den * comb(Ri, r)) for r in range(1, Ri + 1))
            if want_center:   # Long's convention u_{i,r} = (-1)^{i+r} e_{i,r}
                centers.append([Fraction((-1) ** (i + r) * e[r], den) for r in range(Ri + 1)])
        if i >= 1:
            worstM = max(worstM, Fraction(int(np.abs(W).sum()), den))
    assert not np.any(F != 0), "fold did not terminate"
    return worstM, box0, (centers, R) if want_center else None


def fold_threshold(N, S, lo=1.50, hi=1.64, it=12):
    Fnum, den = F0_numer(N, S)
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        M, b0, _ = fold(N, Fnum, den, profile_ceil(N, 4, mid))
        lo, hi = (lo, mid) if (M <= 1 and b0 <= 1) else (mid, hi)
    return hi


# ============================================================================= [E] integer super-blocks
def load_long(draft: Path):
    p = draft / "scripts" / "verify_glazer_critical.py"
    if not p.exists():
        return None
    spec = importlib.util.spec_from_file_location("long_verifier", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def integer_superblock(long, N, j, C):
    B = 4
    S = family(N, j)
    t = profile_floor(N, B, C)
    Fnum, den = F0_numer(N, S)
    _, _, (u, R) = fold(N, Fnum, den, t, want_center=True)
    a = [B * N - ti for ti in t]
    f, err = long.integral_rounding(len(t), a, R, u)
    e_of = {(i, r): (-1) ** (i + r) * f[i][r] for i in range(len(t)) for r in range(len(f[i]))}
    c = class_poly(N, B, t, e_of)
    tgt = [0] * len(c)
    for k, sk in enumerate(S):
        for jj in range(2 * N - 2 * k + 1):
            tgt[k + jj] += sk * comb(2 * N - 2 * k, jj)
    eff = max(Fraction(t[i], N + i) for i in range(len(t)))
    return is_balanced(N, B, c), c == tgt, eff, float(err), sum(len(x) for x in f)


def ratio2_evaluation_certificate(N, Ceff):
    """exact: min over rational th in [0.30,0.46] of sum_i th^i (1-th)^{a_i} (1+th)^{R_i} for the ratio-2 profile
    min(2N, floor(Ceff L)); a value < 1 certifies that no separately balanced block [N,2N) has max T/L <= Ceff."""
    best = None
    for k in range(300, 462, 2):
        th = Fraction(k, 1000)
        tot = Fraction(0)
        for i in range(N):
            L = N + i
            tt = min(2 * N, max(L + 1, math.floor(Ceff * L)))
            tot += th ** i * (1 - th) ** (2 * N - tt) * (1 + th) ** (tt - L - 1)
        if best is None or tot < best[0]:
            best = (tot, th)
    return best


# ============================================================================= [F] chained family
def chained_family_check(K=8):
    w = sp.symbols("w")
    order = 4 ** (K // 2)
    lhs, Nn = 0, 1
    while Nn < order:
        lhs += -w ** Nn * (1 + w) * sum(w ** k for k in range(Nn))
        Nn *= 4
    lac = sum((-1) ** k * w ** (2 ** k) for k in range(0, 2 * K) if 2 ** k < order)
    rhs = sp.series(-(1 + w) / (1 - w) * lac, w, 0, order).removeO()
    diff = sp.expand(lhs - rhs)
    return all(diff.coeff(w, k) == 0 for k in range(order)), order


# ============================================================================= [G] steepest-descent asymptotics of the fold
def steepest_threshold(kappa, s1=False, L=9.0, n=721, nv=36, it=14):
    """Asymptotic SUFFICIENT threshold for the fold with target rate
       log|(1+u)/2| + kappa log|(5-u^2)/4| + (1-kappa) max(0, log|(1-u^2)/4|)   (S=1: first term only).
    c^{(i)}_s = (-2)^i (1/2 pi i) oint F_0(u) u^{-s-1} (u-1)^{-i} du over any contour around {0,1};
    exponent <= bottleneck level of Psi = rate - v log|(1-u)/2| - sigma log|u| between {0,1} and infinity
    (grid percolation, NUMERICAL).  Condition: level <= 0 for all v in [0, v_h], sigma = min(g(1+v), 3-v)."""
    from scipy import ndimage
    xs = np.linspace(-L, L, n)
    X, Y = np.meshgrid(xs, xs)
    U = X + 1j * Y
    eps = 1e-300
    LY, LU = np.log(np.abs((1 - U) / 2) + eps), np.log(np.abs(U) + eps)
    src = (np.abs(U) < 0.03) | (np.abs(U - 1) < 0.03)
    border = np.zeros(U.shape, bool)
    border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    if s1:
        RT, degfac = np.log(np.abs((1 + U) / 2) + eps), 1.0
    else:
        W = (1 - U ** 2) / 4
        RT = (np.log(np.abs((1 + U) / 2) + eps) + kappa * np.log(np.abs((5 - U ** 2) / 4) + eps)
              + (1 - kappa) * np.maximum(0, np.log(np.abs(W) + eps)))
        degfac = 3.0

    def cstar(v, sig):
        Psi = RT - v * LY - sig * LU
        Psi[src] = np.inf
        lo, hi = -5.0, 5.0
        for _ in range(22):
            c = 0.5 * (lo + hi)
            lab, _ = ndimage.label(Psi > c)
            sl = np.unique(lab[src]); sl = sl[sl > 0]
            if np.intersect1d(sl, np.unique(lab[border])).size:
                lo = c
            else:
                hi = c
        return hi

    def ok(g):
        vh = (degfac - g) / (1 + g)
        return all(cstar(v, min(g * (1 + v), 3 - v)) <= 1e-3 for v in np.linspace(0, vh, nv))

    lo, hi = 0.45, 0.75
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if ok(mid) else (mid, hi)
    return hi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--tlimit", type=float, default=120.0)
    ap.add_argument("--draft", type=Path, default=DEFAULT_DRAFT)
    ap.add_argument("--fold-nmax", type=int, default=1024)
    args = ap.parse_args()
    t0 = time.time()
    print("=== AMM12592 procgen 2026-09-23: super-blocks (cross-annulus cancellation) ===")

    print("\n[A] exact CP-SAT  rho_B(N) := min over complement-symmetric super-blocks [N,BN) of max T(L)/L")
    specs = [(2, 4), (2, 8), (4, 4), (8, 2)] if args.quick else [(2, 2), (2, 4), (2, 8), (2, 16), (4, 2), (4, 4), (4, 8), (8, 2), (8, 4)]
    found = {}
    for B, N in specs:
        ts = time.time()
        eff, t, sol, log, infeas = min_ratio(N, B, args.tlimit)
        unk = [c for c, s in log if s == "UNKNOWN"]
        found[(B, N)] = (t, sol)
        print(f"  B={B} N={N:3d}: best {eff} = {float(eff):.5f}; largest ratio proved INFEASIBLE {infeas}; "
              f"timeouts at {unk}  [{time.time() - ts:.0f}s]")
        print(f"        T(L), L={N}..{B * N - 1}: {t}")

    print("\n[B] handoff states S_N (super-block 0-spine deviation = w^N S_N(w)):")
    for (B, N), (t, sol) in found.items():
        if B > 2:
            Sf, Sm1 = state_polynomial(N, B, class_poly(N, B, t, sol))
            print(f"  B={B} N={N}: S_N(w) = {Sf}    S_N(-1) = {Sm1}")

    print("\n[C] asymptotic necessary thresholds for S_N ~ (1+w)^(kappa N) (1-w^N)/(1-w^(kappa N)), B=4 (NUMERICAL grid):")
    for k, nb, evv in kappa_thresholds([0.0, 1 / 16, 1 / 8, 0.15, 0.25] if not args.quick else [0.0, 1 / 8]):
        print(f"  kappa={k:.4f}: natural boundary C >= {1 + nb:.4f};  single-block evaluation C >= {1 + evv:.4f}   (C_* = {1 + GSTAR:.4f})")

    print("\n[D] exact fold thresholds C (profile ceil(C L), l1 masses <= 1, exact level-0 box), B=4:")
    Ns = [64, 128, 256] if args.quick else [n for n in [64, 128, 256, 512, 1024] if n <= args.fold_nmax]
    for N in Ns:
        row = []
        for name, j in [("S=1", -1), ("k=1/32", N.bit_length() - 6), ("k=1/16", N.bit_length() - 5), ("k=1/8", N.bit_length() - 4)]:
            if j >= -1 and (j < 0 or 2 ** j <= N):
                row.append(f"{name}:{fold_threshold(N, family(N, j)):.4f}")
        print(f"  N={N:5d}: " + "  ".join(row), flush=True)

    print("\n[E] integer super-blocks (fold center + Long's lattice rounding, independent re-verification):")
    long = load_long(args.draft)
    if long is None:
        print("  Long draft not available; skipped.")
    else:
        for N, j, C in [(64, 2, Fraction(53, 34)), (128, 4, Fraction(83, 53)), (256, 4, Fraction(157, 100))]:
            ok, eqt, eff, err, sites = integer_superblock(long, N, j, C)
            cert = ratio2_evaluation_certificate(N, eff)
            print(f"  [{N},{4 * N}) kappa={2 ** j / N}: balanced={ok} (target reproduced {eqt}), sites={sites}, "
                  f"max|f-u|={err:.3f}, max T(L)/L = {eff} = {float(eff):.5f}")
            print(f"     separately balanced [{N},{2 * N}) at the same ratio: Long's evaluation sum = {float(cert[0]):.5f} < 1 "
                  f"at theta={cert[1]}  => impossible: {cert[0] < 1}")
            assert ok and eqt and cert[0] < 1

    print("\n[G] asymptotic SUFFICIENT fold thresholds by optimal-contour (steepest-descent) bounds, B=4 (NUMERICAL grid):")
    print(f"  S=1 (golden family): C <= {1 + steepest_threshold(0, s1=True):.4f}   (exact value C_* = {1 + GSTAR:.4f}; grid error ~1e-3)")
    for k in ([1 / 16] if args.quick else [1 / 32, 1 / 16, 0.09, 1 / 8]):
        print(f"  kappa={k:.4f}: C <= {1 + steepest_threshold(k):.4f}", flush=True)

    print("\n[F] chained ILP family sum_j w^(4^j) S_(4^j)(w), S_N = -(1+w)(1-w^N)/(1-w):")
    ok, order = chained_family_check()
    print(f"  = -((1+w)/(1-w)) sum_k (-1)^k w^(2^k)  modulo w^{order}: {ok}  (Hadamard gaps: natural boundary |w|=1)")
    print(f"\nDone. [{time.time() - t0:.0f}s]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
