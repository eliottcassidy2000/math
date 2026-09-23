#!/usr/bin/env python3
"""AMM 12592 (Glazer time-limited fair coin) -- independent audit of C. D. Long's draft
"A Critical Linear Bound with Logarithmic Savings for Glazer's Time-Limited Fair-Coin Problem"
(GitHub long-mathematics/deterministic-von-neumann-fair-extractor, draft dated 2026-09-20).

What this script does (all exact integer / rational arithmetic unless stated):
  A. records the SHA-256 of the fetched draft files (expected verifier digest
     07801d62aedf91b87e122bd96f53a9bd450d64dd4e265caaff364f3da0dced61);
  B. imports Long's verifier ONLY to obtain his packet arrays f_{i,r} for dyadic N
     (his fold + lattice-rounding construction), then re-checks them with code written
     from the statements alone:
       - capacity |f| <= C(R,r) and parity f = C(R,r) (mod 2), boundary packets = +-1;
       - the packet identity  sum_i z^i (1-z)^{a_i} f_i(z) = 1  by an O(N^2) Horner
         recursion in exact integers (different algorithm from his radix evaluation);
       - the deadline profile T(L) = 2N - a_i against min(2N, ceil(C_* L) - s_N), with
         floor(beta m) recomputed by 80-digit mpmath and a separation margin check;
       - the two headline bounds T(L) <= ceil(C_* L) and T(L) <= C_* L - log_5 L + 6;
       - consistency with his own lower bound (Cor. 11.2): sum_i sqrt5^(t_i - C_*(N+i))
         >= 1 + theta_* must hold for every constructed block;
  C. the lacunary identity (new observation, this session): for every separately balanced
     block the 0-spine deviation is exactly +-(pq)^N, so the symmetrized spine function is
     phi(w) = sum_j eps_j w^(2^j); checked here as an exact polynomial identity in p for
     the blocks N <= NMAX_LAC;
  D. corpus cross-walk numbers: Long's per-block slack against the repository's
     "gamma* floor profile" T(L) = L + 1 + floor(gamma* L) (THM-3029/3302/3330 D0 language).
Usage:  python3 amm12592_procgen_20260923_long_audit.py [--draft DIR] [--nmax 512]
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, math, sys
from math import comb
from pathlib import Path

import mpmath as mp

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DRAFT = REPO / "scratch" / "procgen_amm" / "long_draft"
EXPECTED_VERIFIER_SHA = "07801d62aedf91b87e122bd96f53a9bd450d64dd4e265caaff364f3da0dced61"

mp.mp.dps = 80
PHI = (1 + mp.sqrt(5)) / 2
BETA = 2 * mp.log(PHI) / mp.log(5)          # gamma* = C_* - 1
CSTAR = 1 + BETA
THETA = PHI ** -2


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def floor_beta(m: int) -> int:
    x = BETA * m
    f = int(mp.floor(x))
    assert x - f > mp.mpf(10) ** -60 and (f + 1) - x > mp.mpf(10) ** -60, "floor too close"
    return f


def ceil_cstar(m: int) -> int:
    # C_* m = m + beta m irrational
    return m + floor_beta(m) + 1


def s_shift(N: int) -> int:
    j = 0
    while 5 ** (j + 1) <= N:
        j += 1
    return max(0, j - 3)


def load_long(draft: Path):
    spec = importlib.util.spec_from_file_location("long_verifier", draft / "scripts" / "verify_glazer_critical.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["long_verifier"] = mod
    spec.loader.exec_module(mod)
    return mod


def horner_identity(N, a, f):
    """Independent exact check of sum_i z^i (1-z)^{a_i} f_i(z) == 1 via
    V_i := sum_{j<=i} z^j f_j(z) (1-z)^{a_j - a_i},  V_{i+1} = V_i (1-z)^{a_i - a_{i+1}} + z^{i+1} f_{i+1}."""
    V = [0] * (N + 1)
    for r, c in enumerate(f[0]):
        V[r] += c
    for i in range(1, N):
        drop = a[i - 1] - a[i]
        assert drop >= 0
        for _ in range(drop):                  # multiply by (1 - z)
            for k in range(N, 0, -1):
                V[k] -= V[k - 1]
        for r, c in enumerate(f[i]):
            V[i + r] += c
    assert a[N - 1] == 0
    return V[0] == 1 and all(v == 0 for v in V[1:])


def audit_block(N, a, R, f):
    out = {}
    n = N - 1
    # shape
    assert len(a) == N and len(R) == N and len(f) == N
    sites = 0
    min_margin = None
    for i in range(N):
        assert a[i] + R[i] == n - i and a[i] >= 0 and R[i] >= 0
        assert len(f[i]) == R[i] + 1
        for r, c in enumerate(f[i]):
            Q = comb(R[i], r)
            assert abs(c) <= Q and (c - Q) % 2 == 0, (N, i, r)
            sites += 1
            if 0 < r < R[i]:
                mg = Q - abs(c)
                min_margin = mg if min_margin is None else min(min_margin, mg)
        assert abs(f[i][0]) == 1 and abs(f[i][R[i]]) == 1
    out["sites"] = sites
    out["min_interior_margin"] = min_margin
    out["identity"] = horner_identity(N, a, f)
    # deadlines
    sN = s_shift(N)
    worst_ceil = -10 ** 9
    worst_log = mp.mpf(-10) ** 9
    lower_sum = mp.mpf(0)
    slack_vs_floor_profile = []
    for i in range(N):
        L = N + i
        T = 2 * N - a[i]
        want = min(2 * N, ceil_cstar(L) - sN)
        assert T == want, (N, i, T, want)
        assert T <= ceil_cstar(L)
        worst_ceil = max(worst_ceil, T - ceil_cstar(L))
        worst_log = max(worst_log, T - (CSTAR * L - mp.log(L) / mp.log(5)))
        lower_sum += mp.sqrt(5) ** (T - CSTAR * L)
        # repository floor profile L+1+floor(gamma* L) (= ceil(C_* L)); slack D0 = T - that
        if ceil_cstar(L) <= 2 * N:
            slack_vs_floor_profile.append(T - (L + 1 + floor_beta(L)))
    out["s_N"] = sN
    out["max(T-ceil(C*L))"] = worst_ceil
    out["max(T-(C*L-log5 L))"] = float(worst_log)
    out["log_bound_ok(<6)"] = bool(worst_log < 6)
    out["sum sqrt5^(t-C*L) (must be >= 1+theta*)"] = float(lower_sum)
    out["lower_bound_consistent"] = bool(lower_sum >= 1 + THETA)
    out["D0 vs floor profile (max over unsaturated levels)"] = max(slack_vs_floor_profile) if slack_vs_floor_profile else None
    return out


# ---------- C. lacunary identity of the 0-spine deviation ----------
def poly_mul(A, B):
    out = [0] * (len(A) + len(B) - 1)
    for i, x in enumerate(A):
        if x:
            for j, y in enumerate(B):
                if y:
                    out[i + j] += x * y
    return out


def binom_poly(e, sgn):  # (1 + sgn*p)^e  as list in p
    return [comb(e, k) * (sgn ** k) for k in range(e + 1)]


def spine_deviation_block(N, a, R, f):
    """0-spine deviation of block N as a polynomial in p:
       sum_i p^{N+i} q * sum_r e_{i,r} p^r q^{R_i - r}  (unread bits contribute (p+q)^{a_i}=1),
       with e_{i,r} = (-1)^{i+r} f_{i,r} (Long's sign convention)."""
    tot = [0]
    for i in range(N):
        L = N + i
        Ei = [0] * (R[i] + 1)
        for r in range(R[i] + 1):
            e = (-1) ** (i + r) * f[i][r]
            if e:
                term = poly_mul([0] * r + [1], binom_poly(R[i] - r, -1))
                for k, v in enumerate(term):
                    Ei[k] += e * v
        term = poly_mul([0] * L + [1, -1], Ei)  # p^L q E_i(p)
        if len(term) > len(tot):
            tot += [0] * (len(term) - len(tot))
        for k, v in enumerate(term):
            tot[k] += v
    while len(tot) > 1 and tot[-1] == 0:
        tot.pop()
    return tot


def check_lacunary(N, a, R, f):
    dev = spine_deviation_block(N, a, R, f)
    # expected +-(p q)^N = +-p^N (1-p)^N
    target = poly_mul([0] * N + [1], binom_poly(N, -1))
    if dev == target:
        return +1
    if dev == [-x for x in target]:
        return -1
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draft", type=Path, default=DEFAULT_DRAFT)
    ap.add_argument("--nmax", type=int, default=512)
    ap.add_argument("--nmax-lac", type=int, default=128)
    args = ap.parse_args()
    draft = args.draft
    print("=== AMM12592 procgen 2026-09-23: audit of Long's draft ===")
    files = ["deterministic_von_neumann_fair_extractor.tex", "scripts/verify_glazer_critical.py",
             "scripts/supplementary_checks.py", "scripts/certificates/finite_blocks.json",
             "scripts/verification_output.txt"]
    if not (draft / files[1]).exists():
        print(f"Long draft not found under {draft}; fetch it first (see run script). Skipping audit.")
        return 0
    for fn in files:
        print(f"  sha256 {sha(draft / fn)}  {fn}")
    vs = sha(draft / files[1])
    print(f"  verifier digest matches manuscript/README: {vs == EXPECTED_VERIFIER_SHA}")
    assert vs == EXPECTED_VERIFIER_SHA

    long = load_long(draft)
    # certificate blocks (N<128) -- checked from the JSON, not from his construction
    cert = json.loads((draft / files[3]).read_text())
    print("\n[B1] distributed finite certificate (N = 2..64), independent checker:")
    for blk in cert["blocks"]:
        N, a, R, f = blk["N"], blk["a"], blk["R"], blk["f"]
        res = audit_block(N, a, R, f)
        print(f"  N={N:4d}: sites={res['sites']:6d} identity={res['identity']} "
              f"margin={res['min_interior_margin']} s_N={res['s_N']} "
              f"sum_lb={res['sum sqrt5^(t-C*L) (must be >= 1+theta*)']:.4f} ok={res['lower_bound_consistent']}")
        assert res["identity"] and res["lower_bound_consistent"]
    print("\n[B2] Long's construction (fold + lattice rounding) regenerated for N = 128 .. nmax,"
          " re-checked by independent code:")
    N = 128
    while N <= args.nmax:
        a, R = long.profile(N)
        u, _ = long.real_center(N, a, R)
        f, err = long.integral_rounding(N, a, R, u)
        res = audit_block(N, a, R, f)
        print(f"  N={N:5d}: sites={res['sites']:7d} identity={res['identity']} margin={res['min_interior_margin']} "
              f"s_N={res['s_N']} max(T-ceil)={res['max(T-ceil(C*L))']} "
              f"max(T-(C*L-log5L))={res['max(T-(C*L-log5 L))']:.3f} "
              f"lb_sum={res['sum sqrt5^(t-C*L) (must be >= 1+theta*)']:.4f} "
              f"D0_vs_floor_profile={res['D0 vs floor profile (max over unsaturated levels)']} max|f-u|={float(err):.4f}")
        assert res["identity"] and res["log_bound_ok(<6)"] and res["lower_bound_consistent"]
        N *= 2
    print("\n[C] lacunary 0-spine deviation: block N contributes exactly eps_N (pq)^N")
    signs = []
    for blk in cert["blocks"]:
        s = check_lacunary(blk["N"], blk["a"], blk["R"], blk["f"])
        signs.append((blk["N"], s))
        assert s != 0
    N = 128
    while N <= args.nmax_lac:
        a, R = long.profile(N)
        u, _ = long.real_center(N, a, R)
        f, _ = long.integral_rounding(N, a, R, u)
        s = check_lacunary(N, a, R, f)
        signs.append((N, s))
        assert s != 0
        N *= 2
    print("  (N, eps_N):", signs, "  plus the L=1 rule 01->H,10->T giving +pq.")
    print("  => 2F(p) - p = sum_j eps_j (pq)^(2^j): phi(w) is lacunary (Hadamard gaps), natural boundary |w|=1.")
    print("\nAll audit checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
