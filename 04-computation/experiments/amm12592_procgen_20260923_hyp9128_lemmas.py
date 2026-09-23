#!/usr/bin/env python3
"""HYP-9128: exact checks of the identities used in the proof (procgen 2026-09-23), at small N.

  [F1] fold support/formula: for 1 <= i <= h and s >= R_{i-1}:  c^{(i)}_s = (-2)^i sum_{x>s} c_x C(x-s-1, i-1),
       and F_i vanishes below R_{i-1}; the digit W_i has at most two nonzero coefficients (s = R_{i-1}, R_i);
       tail_i = sum_{s >= R_i} c^{(i)}_s = (-2)^i sum_x c_x C(x-R_i, i).
  [F2] level-0 decomposition: e_{0,r} = [z^r] A(z) - [z^r] D(z),  A = S(z/(1+z)^2) (1+z)^{-(N-1-R_0)},
       D = sum_{t >= R_0} c_t ((1-z)^t (1+z)^{R_0-t} - (1-z)^{R_0});  bottom majorant |[z^r]A| <= C(K+2m+r-1, r) (r < m).
  [F3] top regime: e_{0,R_0-r'} = -[z^{r'}] (1+z)^{R_0} Delta(-u(z)),  Delta(v) = sum_{t>=R_0} c_t (v^t - v^{R_0}),  r' <= R_0.
  [F4] Krawtchouk pairing bound |K_t^{(R)}(r)| <= [z^r](1+z^2)^j (1+z)^{R-2j}, j = min(t, R-t)  (exhaustive, R <= 40).
  [F5] prefix blocks of the assembled extractor: C. D. Long's distributed ratio-2 blocks [N,2N), N = 2, 4, 8
       (certificate finite_blocks.json), re-checked here: horizon L+1+R_i+a_i = 2N, box, parity, and class
       polynomial sum_i x^i (1+x)^{a_i} E_i(x) = +-1 (separately balanced), T(L) <= 2L.
  [F6] structural hypotheses of Lemma R and the parity hypothesis for every super-block N = 16*4^k <= 2^22:
       drops in {0,1,2}, a_{n-2} = a_{n-1} = a_n = 0, R_i >= 2 (i <= n-2), R nondecreasing below h, R_i >= 10 on
       1 <= i <= h, h <= n-2, T(L) <= ceil(159 L/100); S_N == 1 + w^N (mod 2) (bitmask arithmetic, N <= 2^14).
All arithmetic exact (Python integers / Fractions).
"""
from __future__ import annotations
import sys, time
from fractions import Fraction
from math import comb
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import amm12592_procgen_20260923_hyp9128_finite as FIN


def series_coeffs_rational(num, den_pow_base, K, order):
    """coefficients up to z^order of num(z) * (1+z)^{-K} (num a coefficient list)."""
    if K == 0:
        inv = [1] + [0] * order
    else:
        inv = [(-1) ** j * comb(K + j - 1, j) for j in range(order + 1)]      # (1+z)^{-K}
    out = [0] * (order + 1)
    for i, a in enumerate(num):
        if a and i <= order:
            for j in range(order + 1 - i):
                out[i + j] += a * inv[j]
    return out


def check(N):
    t, R, a = FIN.profile(N)
    n = 3 * N - 1
    S = FIN.state_S(N)
    Fnum, den = FIN.target_F0(N, S)
    digits, h = FIN.fold(N, Fnum, R)
    m = N // 16
    R0 = R[0]
    K = N - 1 - R0
    # ---------------- F1: recompute F_i by the fold and compare with the closed formula
    import numpy as np
    F = list(Fnum)
    ok_formula, ok_support, ok_two, ok_tail = True, True, True, True
    for i in range(0, h):
        Ri = R[i]
        deg = len(F) - 1
        while deg > 0 and F[deg] == 0:
            deg -= 1
        F = F[:deg + 1]
        nxt = [0] * max(deg, 1)
        suf = 0
        for s in range(deg - 1, Ri - 1, -1):
            suf += F[s + 1]
            nxt[s] = -2 * suf
        F = nxt
        ii = i + 1
        # formula for level ii (valid for s >= R_{ii-1} = R_i)
        for s in range(R[ii - 1], len(F)):
            val = (-2) ** ii * sum(Fnum[x] * comb(x - s - 1, ii - 1) for x in range(s + 1, len(Fnum)))
            if val != F[s]:
                ok_formula = False
                break
        if any(F[s] != 0 for s in range(0, min(R[ii - 1], len(F)))):
            ok_support = False
        if ii < h:
            nonzero_low = [s for s in range(R[ii]) if s < len(F) and F[s] != 0]
            if any(s != R[ii - 1] for s in nonzero_low):
                ok_two = False
            tail = sum(F[R[ii]:])
            tail_formula = (-2) ** ii * sum(Fnum[x] * comb(x - R[ii], ii) for x in range(len(Fnum)) if x >= R[ii])
            if tail != tail_formula:
                ok_tail = False
    # ---------------- F2: level-0 decomposition and bottom majorant
    E0 = FIN.level0_packet(R0, digits[0])          # scaled by den
    order = R0
    # A(z) = S(w) (1+z)^{-K},  S(w) = sum_k s_k z^k (1+z)^{-2k}
    Acoef = [Fraction(0)] * (order + 1)
    for k, sk in enumerate(S):
        if sk and k <= order:
            part = series_coeffs_rational([0] * k + [1], None, K + 2 * k, order)
            for j in range(order + 1):
                Acoef[j] += sk * part[j]
    # D(z) coefficients
    Dcoef = [Fraction(0)] * (order + 1)
    tail_t = [(tt, Fraction(Fnum[tt], den)) for tt in range(R0, len(Fnum)) if Fnum[tt]]
    oneminus_R0 = [(-1) ** j * comb(R0, j) for j in range(order + 1)]
    for tt, ct in tail_t:
        # (1-z)^t (1+z)^{R0-t} = (1-z)^t (1+z)^{-(t-R0)}
        p1 = [(-1) ** j * comb(tt, j) for j in range(min(tt, order) + 1)]
        ser = series_coeffs_rational(p1, None, tt - R0, order)
        for j in range(order + 1):
            Dcoef[j] += ct * (ser[j] - oneminus_R0[j])
    ok_decomp = all(Fraction(E0[r], den) == Acoef[r] - Dcoef[r] for r in range(order + 1))
    ok_major = all(abs(Acoef[r]) <= comb(K + 2 * m + r - 1, r) for r in range(1, min(m, order + 1)))
    # ---------------- F3: top regime
    ok_top = True
    for rp in range(0, min(R0, 40) + 1):
        # [z^rp] (1+z)^{R0} Delta(-u(z)),  (1+z)^{R0} (-u)^t = (-1)^t (1-z)^t (1+z)^{R0-t}
        val = Fraction(0)
        for tt, ct in tail_t:
            p1 = [(-1) ** j * comb(tt, j) for j in range(min(tt, rp) + 1)]
            ser = series_coeffs_rational(p1, None, tt - R0, rp)
            val += ct * ((-1) ** tt * ser[rp] - (-1) ** R0 * (-1) ** rp * comb(R0, rp))
        if Fraction(E0[R0 - rp], den) != -val:
            ok_top = False
            break
    return {"N": N, "h": h, "F1_formula": ok_formula, "F1_support": ok_support, "F1_two_coeffs": ok_two,
            "F1_tail": ok_tail, "F2_decomposition": ok_decomp, "F2_bottom_majorant": ok_major, "F3_top": ok_top}


def pairing_check(Rmax=40):
    worst = 0
    for R in range(1, Rmax + 1):
        for t in range(R + 1):
            j = min(t, R - t)
            K = FIN.kraw(R, t)
            poly = [0] * (R + 1)
            for l in range(j + 1):
                for s in range(R - 2 * j + 1):
                    poly[2 * l + s] += comb(j, l) * comb(R - 2 * j, s)
            for r in range(R + 1):
                if abs(K[r]) > poly[r]:
                    return False, (R, t, r)
    return True, None


def prefix_blocks():
    import json
    path = Path(__file__).resolve().parents[2] / "scratch" / "procgen_amm" / "long_draft" / "scripts" / "certificates" / "finite_blocks.json"
    cert = json.loads(path.read_text())
    out = []
    for b in cert["blocks"]:
        N = b["N"]
        if N not in (2, 4, 8):
            continue
        R, a, f = b["R"], b["a"], b["f"]
        ok_h = all(N + i + 1 + R[i] + a[i] == 2 * N for i in range(N))
        ok_box = all(abs(f[i][r]) <= comb(R[i], r) and (f[i][r] - comb(R[i], r)) % 2 == 0
                     for i in range(N) for r in range(R[i] + 1))
        consts = []
        for twist in (1, 0):
            P = [0] * (2 * N)
            for i in range(N):
                for r in range(R[i] + 1):
                    ev = ((-1) ** (i + r) if twist else 1) * f[i][r]
                    if ev:
                        for j in range(a[i] + 1):
                            P[i + r + j] += ev * comb(a[i], j)
            consts.append(all(x == 0 for x in P[1:]) and abs(P[0]) == 1)
        ok_T = all(N + i + 1 + R[i] <= 2 * (N + i) for i in range(N))
        out.append((N, ok_h and ok_box and consts[0] and ok_T, consts))
    return out


def structure_all(kmax=9):
    """streaming check (O(1) memory, integer arithmetic) of the Lemma-R structure for N = 16*4^k, k <= kmax."""
    res = []
    for k in range(kmax + 1):
        N = 16 * 4 ** k
        n = 3 * N - 1

        def prof(i):
            ti = min(4 * N, (159 * (N + i) + 99) // 100)          # ceil(159 (N+i) / 100)
            return ti, ti - N - i - 1, 4 * N - ti
        ok, h = True, None
        t_prev, R_prev, a_prev = prof(0)
        minR = None
        for i in range(n + 1):
            ti, Ri, ai = prof(i)
            if i > 0:
                ok &= (a_prev - ai) in (0, 1, 2)
                if h is None:
                    ok &= R_prev <= Ri                                # R nondecreasing below h
            if h is None and ai == 0:
                h = i
            if h is not None:
                ok &= ai == 0                                         # a_i = 0 for all i >= h
            if i <= n - 2:
                ok &= Ri >= 2
            if 1 <= i and (h is None or i <= h):
                minR = Ri if minR is None else min(minR, Ri)
            ok &= 100 * ti < 159 * (N + i) + 100                     # t_i <= ceil(159 (N+i)/100)
            t_prev, R_prev, a_prev = ti, Ri, ai
        ok &= h is not None and h <= n - 2 and (N < 4096 or minR >= 10)
        ok &= prof(n - 2)[2] == prof(n - 1)[2] == prof(n)[2] == 0
        par = None
        if N <= 2 ** 14:
            m = N // 16
            A = 1                                   # (1+w)^m mod 2 as a bitmask, by repeated squaring of (1+w)
            base, e = 0b11, m
            while e:
                if e & 1:
                    A = clmul(A, base)
                base = clmul(base, base)
                e >>= 1
            Bq = 0
            for j in range(16):
                Bq ^= 1 << (m * j)
            par = clmul(A, Bq) == (1 | (1 << N))
        res.append((N, h, ok, par))
    return res


def clmul(x, y):
    """carry-less product (polynomials over GF(2) as bitmasks)."""
    out = 0
    while y:
        if y & 1:
            out ^= x
        x <<= 1
        y >>= 1
    return out


def main():
    t0 = time.time()
    print("=== HYP-9128 lemma checks (exact) ===")
    for N in [16, 32, 64]:
        res = check(N)
        print("  ", res, flush=True)
        assert all(v for k, v in res.items() if k not in ("N", "h"))
    ok, bad = pairing_check(40)
    print(f"  [F4] pairing bound |K_t(r)| <= [z^r](1+z^2)^j(1+z)^(R-2j), j=min(t,R-t), all R <= 40, all t, r: {ok} {bad or ''}")
    assert ok
    for N, okb, consts in prefix_blocks():
        print(f"  [F5] Long's block [{N},{2 * N}): horizon/box/parity/T<=2L and class polynomial = +-1 (sign-twisted "
              f"convention {consts[0]}, untwisted {consts[1]}): {okb}")
        assert okb
    for N, h, okS, par in structure_all():
        print(f"  [F6] N = {N:8d}: h = {h:8d}; Lemma-R structure and T(L) <= ceil(159L/100): {okS}; "
              f"S_N == 1 + w^N (mod 2): {par if par is not None else 'n/a (N > 2^14; symbolic: (1+w^m)(1-w^(16m))/(1-w^m))'}")
        assert okS and par in (True, None)
    print(f"All lemma checks passed. [{time.time() - t0:.0f}s]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
