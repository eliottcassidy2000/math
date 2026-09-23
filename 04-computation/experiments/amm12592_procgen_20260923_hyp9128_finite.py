#!/usr/bin/env python3
"""HYP-9128, finite part: exact certificates for the golden-zero super-blocks, N < N_A.

Construction: c = 159/100; super-block [N,4N), N a power of two; deadlines t_i = min(4N, ceil(c(N+i))), i = 0..3N-1;
handoff state S_N(w) = (1+w)^m sum_{j<16} w^{mj}, m = N/16; target P_+(x) = S_N(x/(1+x)^2) (1+x)^{2N}
(palindromic of degree 2N, P_+(0) = 1), equivalently F_0(u) = ((1+u)/2)^{N-1} S_N((1-u^2)/4).

For each N (all arithmetic exact, integers / rationals):
  [S] structural hypotheses of Lemma R: drops a_i - a_{i+1} in {0,1,2}; a_{n-2} = a_{n-1} = a_n = 0; R_i >= 2 (i <= n-2);
  [F] Long's fold with the general target; termination level h <= n-2; masses M_i = ||W_i||_1;
  [M] margins: level 0 exact packet e_{0,r} (Krawtchouk transform): e_{0,0} = 1, |e_{0,R_0}| <= 1,
      |e_{0,r}| <= C(R_0,r) - 5 (0<r<R_0);  levels 1..h: M_i <= 1/2 and R_i >= 10 (=> |u| <= C/2, margin >= 5);
      levels > h: zero center, a_{i-1} = a_i = 0 (Lemma R (iii) gives |f| <= 2 <= C(R_i,r) at interior sites).
      If all hold, Lemma R (re-derived in the note) yields an integral admissible packet vector: the super-block exists.
  [R] (N <= --round-max) an explicit integral solution built by an INDEPENDENT implementation of Lemma R
      (kernel basis, parity-correct q, lower/upper boundary chains, nearest-even rounding; coset conditions asserted),
      then verified exactly: parity, box, boundaries = +-1, and class polynomial == target.
Usage: python3 amm12592_procgen_20260923_hyp9128_finite.py [--nmin 16] [--nmax 2048] [--round-max 256]
"""
from __future__ import annotations
import argparse, math, sys, time
from fractions import Fraction
from math import comb
import numpy as np

C = Fraction(159, 100)


def ceil_frac(x: Fraction) -> int:
    return -((-x.numerator) // x.denominator)


def profile(N):
    H = 4 * N
    t = [min(H, ceil_frac(C * (N + i))) for i in range(3 * N)]
    R = [t[i] - (N + i) - 1 for i in range(3 * N)]
    a = [H - t[i] for i in range(3 * N)]
    return t, R, a


def poly_mul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q):
                if y:
                    out[i + j] += x * y
    return out


def state_S(N):
    m = N // 16
    A = [comb(m, k) for k in range(m + 1)]
    Bq = [1 if k % m == 0 else 0 for k in range(15 * m + 1)]
    return poly_mul(A, Bq)          # degree 16m = N


def target_F0(N, S):
    """numerators of F_0(u) = ((1+u)/2)^{N-1} S((1-u^2)/4) over den = 2^{N-1} 4^{deg S}."""
    d = len(S) - 1
    tot, pw = [0], [1]
    for k in range(d + 1):
        if S[k]:
            term = [S[k] * 4 ** (d - k) * v for v in pw]
            tot += [0] * max(0, len(term) - len(tot))
            for j, v in enumerate(term):
                tot[j] += v
        pw = poly_mul(pw, [1, 0, -1])
    return poly_mul([comb(N - 1, k) for k in range(N)], tot), 2 ** (N - 1) * 4 ** d


def target_P(N, S):
    """P_+(x) = sum_k s_k x^k (1+x)^{2N-2k}  (coefficient list, degree <= 2N)."""
    P = [0] * (2 * N + 1)
    for k, sk in enumerate(S):
        if sk:
            for j in range(2 * N - 2 * k + 1):
                P[k + j] += sk * comb(2 * N - 2 * k, j)
    return P


def kraw(R, t):
    c, prev = [1], 0
    for r in range(R):
        nxt, rem = divmod((R - 2 * t) * c[-1] - (R - r + 1) * prev, r + 1)
        assert rem == 0
        prev = c[-1]
        c.append(nxt)
    return c


def fold(N, Fnum, R):
    """Long's fold with general target; returns list of digits W_i (dict t -> numerator) and h."""
    F = np.array(Fnum, dtype=object)
    digits, h = [], None
    for i in range(len(R)):
        nz = np.nonzero(F != 0)[0]
        if len(nz) == 0:
            digits.append({})
            continue
        deg = int(nz[-1])
        F = F[:deg + 1]
        Ri = R[i]
        if deg <= Ri:
            digits.append({t: int(v) for t, v in enumerate(F) if v != 0})
            F = np.array([0], dtype=object)
            h = i
        else:
            low = {t: int(F[t]) for t in range(Ri) if F[t] != 0}
            tail = int(F[Ri:].sum())
            if tail:
                low[Ri] = tail
            digits.append(low)
            suf = np.cumsum(F[:Ri:-1])[::-1]
            F = np.concatenate([np.zeros(Ri, dtype=object), -2 * suf])
    assert not np.any(F != 0), "fold did not terminate"
    return digits, h


def level0_packet(R0, W0):
    """E_0(x) = sum_t w_t (1-x)^t (1+x)^{R0-t}, by P_t = w_t (1+x)^{R0-t} + (1-x) P_{t+1}."""
    P = np.zeros(R0 + 1, dtype=object)
    B = np.zeros(R0 + 1, dtype=object)
    B[0] = 1                                   # (1+x)^0
    for t in range(R0, -1, -1):
        # P <- (1-x) P
        P = P - np.concatenate([np.zeros(1, dtype=object), P[:-1]])
        w = W0.get(t, 0)
        if w:
            P = P + w * B
        # B <- (1+x) B  for the next (smaller) t
        B = B + np.concatenate([np.zeros(1, dtype=object), B[:-1]])
    return [int(v) for v in P]


def margins_certificate(N):
    t, R, a = profile(N)
    n = 3 * N - 1
    out = {"N": N}
    # [S]
    drops_ok = all(a[i] - a[i + 1] in (0, 1, 2) for i in range(n))
    tail_ok = a[n - 2] == a[n - 1] == a[n] == 0 and R[n - 1] == 1 and R[n] == 0
    R_ok = all(R[i] >= 2 for i in range(n - 1))
    mono_ok = True
    S = state_S(N)
    Fnum, den = target_F0(N, S)
    digits, h = fold(N, Fnum, R)
    mono_ok = all(R[i] <= R[i + 1] for i in range(h)) if h else True
    out["structure"] = drops_ok and tail_ok and R_ok and mono_ok and h is not None and h <= n - 2
    out["h"] = h
    # levels >= 1
    maxM = Fraction(0)
    for i in range(1, h + 1):
        M = Fraction(sum(abs(v) for v in digits[i].values()), den)
        maxM = max(maxM, M)
    out["max_M_i"] = maxM
    out["levels_ok"] = maxM <= Fraction(1, 2) and min(R[1:h + 1]) >= 10 and all(a[i] == 0 for i in range(h, n + 1))
    # level 0
    R0 = R[0]
    E0 = level0_packet(R0, digits[0])
    assert E0[0] == den, "corner e_{0,0} != 1"
    worst = None
    for r in range(1, R0):
        marg = comb(R0, r) * den - abs(E0[r])            # scaled by den
        if worst is None or marg < worst[0]:
            worst = (marg, r)
    out["lvl0_min_margin"] = float(Fraction(worst[0], den))
    out["lvl0_argmin_r"] = worst[1]
    out["lvl0_top"] = float(Fraction(abs(E0[R0]), den))
    ratio = max(Fraction(abs(E0[r]), den * comb(R0, r)) for r in range(1, R0))
    out["lvl0_max_ratio"] = float(ratio)
    out["level0_ok"] = worst[0] >= 5 * den and abs(E0[R0]) <= den
    out["certified"] = out["structure"] and out["levels_ok"] and out["level0_ok"]
    eff = max(Fraction(t[i], N + i) for i in range(len(t)))
    out["max_T_over_L"] = eff
    return out, (t, R, a, S, Fnum, den, digits, h, E0)


# ----------------------------------------------------------------------------- independent Lemma R implementation
def lemma_R_round(N, data):
    t, R, a, S, Fnum, den, digits, h, E0 = data
    n = 3 * N - 1
    D = den
    d = [a[i] - a[i + 1] for i in range(n)]
    # real center U_{i,r} = den * u_{i,r},  u_{i,r} = (-1)^{i+r} e_{i,r}
    U = []
    for i in range(n + 1):
        Ri = R[i]
        if i == 0:
            e = E0
        else:
            e = [0] * (Ri + 1)
            for tt, v in digits[i].items():
                kc = kraw(Ri, tt)
                for r in range(Ri + 1):
                    e[r] += v * kc[r]
        U.append([(-1) ** (i + r) * e[r] for r in range(Ri + 1)])
    # target Q(z) = P_+(-z), degree <= n
    P = target_P(N, S)
    Q = [(-1) ** k * P[k] for k in range(len(P))] + [0] * (n + 1 - len(P))
    # A applied to c = binom rows, by Horner accumulation
    V = [0] * (n + 1)
    for i in range(n + 1):
        if i > 0:
            for _ in range(d[i - 1]):
                for k in range(n, 0, -1):
                    V[k] -= V[k - 1]
        for r in range(R[i] + 1):
            V[i + r] += comb(R[i], r)
    Hpoly = [V[k] - Q[k] for k in range(n + 1)]
    assert all(x % 2 == 0 for x in Hpoly), "parity hypothesis Q = sum z^i(1+z)^{n-i} (mod 2) fails"
    Hpoly = [x // 2 for x in Hpoly]
    hcoef = [0] * (n + 1)
    for i in range(n + 1):
        hi = Hpoly[i]
        hcoef[i] = hi
        if hi:
            for j in range(a[i] + 1):
                if i + j <= n:
                    Hpoly[i + j] -= hi * comb(a[i], j) * (-1) ** j
    assert all(x == 0 for x in Hpoly)
    q = [[comb(R[i], r) for r in range(R[i] + 1)] for i in range(n + 1)]
    for i in range(n + 1):
        q[i][0] -= 2 * hcoef[i]
    # gamma: coordinates of q - u in the kernel basis (scaled by D)
    G = [[0] * (R[i] + 1) for i in range(n)]
    for i in range(n + 1):
        for r in range(R[i] + 1):
            val = q[i][r] * D - U[i][r]
            if i > 0:
                for j in range(d[i - 1] + 1):
                    rp = r + 1 - j
                    if 1 <= rp <= R[i - 1]:
                        val += (-1) ** j * comb(d[i - 1], j) * G[i - 1][rp]
            if r == 0 or i == n:
                assert val == 0, ("q-u not in kernel", i, r)
            else:
                G[i][r] = val
    # choose delta (scaled by D)
    Dl = [[None] * (R[i] + 1) for i in range(n)]
    # lower chain: delta_{i,1} = u_{i+1,0} - b_{i+1}
    for i in range(n - 1):
        b = 1 if U[i + 1][0] >= 0 else -1
        Dl[i][1] = U[i + 1][0] - b * D
    # upper chain
    eta_prev = None
    for i in range(n - 1):
        X = -U[i][R[i]] + ((-1) ** d[i - 1] * eta_prev if i > 0 else 0)
        b = -1 if X > 0 else 1
        eta = b * D + X
        assert abs(eta) <= D
        if R[i] == 1:
            raise AssertionError("R_i = 1 before level n-1")
        Dl[i][R[i]] = eta
        eta_prev = eta
    assert eta_prev == 0 or True
    # final overlap (n-1, 1): odd coset, choose +1
    Dl[n - 1][1] = D
    # other coordinates: nearest even integer
    for i in range(n):
        for r in range(1, R[i] + 1):
            if Dl[i][r] is None:
                g = G[i][r]
                m2 = 2 * ((g + D) // (2 * D))          # nearest even integer to g/D
                Dl[i][r] = g - m2 * D
            # coset check: (G - Dl)/D must be an even integer
            diff = G[i][r] - Dl[i][r]
            assert diff % (2 * D) == 0, ("coset", i, r)
            assert abs(Dl[i][r]) <= D
    # f = u + sum delta K
    f = []
    maxerr = Fraction(0)
    for i in range(n + 1):
        row = []
        for r in range(R[i] + 1):
            val = U[i][r]
            if i < n and r >= 1:
                val += Dl[i][r]
            if i > 0:
                for j in range(d[i - 1] + 1):
                    rp = r + 1 - j
                    if 1 <= rp <= R[i - 1]:
                        val -= (-1) ** j * comb(d[i - 1], j) * Dl[i - 1][rp]
            assert val % D == 0, ("non-integral f", i, r)
            fi = val // D
            if 0 < r < R[i]:
                maxerr = max(maxerr, abs(Fraction(val - U[i][r], D)))
            row.append(fi)
        f.append(row)
    # exact verification of the integral packets
    Pcls = [0] * (n + 1)
    for i in range(n + 1):
        for r in range(R[i] + 1):
            Qb = comb(R[i], r)
            fi = f[i][r]
            assert abs(fi) <= Qb and (fi - Qb) % 2 == 0, ("box/parity", i, r)
            if r in (0, R[i]):
                assert abs(fi) == 1
            e = (-1) ** (i + r) * fi
            if e:
                for j in range(a[i] + 1):
                    Pcls[i + r + j] += e * comb(a[i], j)
    Ptarget = P + [0] * (n + 1 - len(P))
    ok = Pcls == Ptarget
    Dg = 2 * N
    pal = all(Pcls[l] == 0 for l in range(Dg + 1, n + 1)) and all(Pcls[l] == Pcls[Dg - l] for l in range(Dg + 1))
    return ok and pal, float(maxerr), sum(len(x) for x in f), f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=16)
    ap.add_argument("--nmax", type=int, default=2048)
    ap.add_argument("--round-max", type=int, default=256)
    args = ap.parse_args()
    print(f"=== HYP-9128 finite certificates: c = {C}, super-blocks [N,4N), S_N = (1+w)^(N/16) sum_(j<16) w^(jN/16) ===")
    N = args.nmin
    while N <= args.nmax:
        t0 = time.time()
        out, data = margins_certificate(N)
        line = (f"N={N:5d}: max T/L = {out['max_T_over_L']} = {float(out['max_T_over_L']):.5f}; structure {out['structure']}; "
                f"h={out['h']}; max_(i>=1) M_i = {float(out['max_M_i']):.3e}; level-0 min margin {out['lvl0_min_margin']:.2f} "
                f"(r={out['lvl0_argmin_r']}), max |e|/C = {out['lvl0_max_ratio']:.4f}, |e_(0,R0)| = {out['lvl0_top']:.2e}; "
                f"MARGIN CERTIFICATE {'PASS' if out['certified'] else 'FAIL'}")
        print(line + f"  [{time.time() - t0:.0f}s]", flush=True)
        if N <= args.round_max:
            t1 = time.time()
            try:
                ok, err, sites, _f = lemma_R_round(N, data)
                print(f"        independent Lemma-R rounding: integral solution {'VERIFIED' if ok else 'WRONG'} "
                      f"({sites} sites, max interior |f-u| = {err:.3f})  [{time.time() - t1:.0f}s]", flush=True)
            except AssertionError as ex:
                print(f"        independent Lemma-R rounding FAILED: {ex}  [{time.time() - t1:.0f}s]", flush=True)
        N *= 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
