#!/usr/bin/env python3
"""collatz_dyadic_energy_breaks_20261013.py -- the dyadic energy of the signed mod-3 parity discrepancy H, reduced to the
break process of the stopping-time parity (session collatz-posets-zeta5-20260927, opus S28, 2026-10-13).

 Objects (Codex amplitude note, collatz_amplitude_boundary_20261013c.md): t(n) = (-1)^tau(n), tau = strict first-hit length of
 1 under T(n) = n/2, (3n+1)/2; w_n = 3 [n = 2 mod 3] - 1; H(N) = sum_(n<=N) w_n t(n); W(n) = 0, -1, +1 for n = 0, 1, 2 mod 3;
 identity (6): H(N) = W(N) t(N) - sum_(n<N) W(n) (t(n+1) - t(n)).
 (1) Parity of tau for 1 <= n <= 2^27 (numba; descend below n and look up).
 (2) Breaks beta(n) = [t(n+1) != t(n)]: per octave the count B_k, the counts by residue of n mod 3, the signed current
     J_k = sum_(breaks in octave) W(n) t(n+1) (so Delta H = -2 J_k), max |H|, the normalized energy int_X^(2X) H^2 dx / X^2
     (Codex: 1.77 at 2^24 ... 0.095 at 2^33), the break density B_k/2^k and its octave ratio.
 (3) Whiteness of the signed break current: autocorrelation of the sequence W(n_i) t(n_i + 1) over consecutive breaks n_i
     (lags 1..8) per octave; the ratio J_k^2 / (number of breaks with W != 0) (a white current has ratio ~ 1 on average).
 (4) The exact first-step identities for breaks: beta(2m) = [t(m) != t(3m+2)], beta(2m-1) = [t(3m-1) != t(m)] (verified),
     the even/odd split of breaks per octave, and the octave increment of H as the exact transfer
     Delta H_(k+1) = sum_m W(2m)(t(3m+2) - t(m)) + sum_m W(2m-1)(t(m) - t(3m-1)) over m in octave k (verified).
Usage: python3 collatz_dyadic_energy_breaks_20261013.py [log2 N]
"""
import sys, time, math
import numpy as np
from numba import njit

T0 = time.time()
LOGN = int(sys.argv[1]) if len(sys.argv) > 1 else 27
N = 1 << LOGN


@njit(cache=True)
def parity_table(N):
    """p[n] = tau(n) mod 2 for 1 <= n <= N, tau = number of T-steps to first reach 1 (tau(1) = 0)."""
    p = np.zeros(N + 2, dtype=np.int8)
    p[1] = 0
    for n in range(2, N + 1):
        x = n; k = 0
        while x >= n:
            if x & 1:
                x = (3 * x + 1) >> 1
            else:
                x >>= 1
            k += 1
            if x == 1:
                break
        # x < n now (or x == 1): tau(n) = k + tau(x)
        p[n] = (k + p[x]) & 1
    return p


@njit(cache=True)
def octave_stats(p, N, LOGN):
    """returns per-octave arrays: B (breaks), B1, B2 (breaks at n = 1, 2 mod 3), J (signed current), maxH, energy (sum H^2 over
    the octave), Hend (H at 2^(k+1)), Beven (breaks at even n), and the whiteness lag products."""
    nk = LOGN
    B = np.zeros(nk + 1, dtype=np.int64); B1 = np.zeros(nk + 1, dtype=np.int64); B2 = np.zeros(nk + 1, dtype=np.int64)
    J = np.zeros(nk + 1, dtype=np.int64); maxH = np.zeros(nk + 1, dtype=np.int64); energy = np.zeros(nk + 1, dtype=np.float64)
    Hend = np.zeros(nk + 1, dtype=np.int64); Beven = np.zeros(nk + 1, dtype=np.int64)
    lagsum = np.zeros((nk + 1, 9), dtype=np.float64); lagcnt = np.zeros((nk + 1, 9), dtype=np.int64)
    hist = np.zeros(9, dtype=np.int64)   # last signed break values (ring)
    hpos = 0; nhist = 0
    H = 0; k = 0
    for n in range(1, N + 1):
        r = n % 3
        w = 2 if r == 2 else -1
        t = 1 - 2 * p[n]
        H += w * t
        if n >= (1 << (k + 1)):
            k += 1
        if H > maxH[k]:
            maxH[k] = H
        if -H > maxH[k]:
            maxH[k] = -H
        energy[k] += float(H) * float(H)
        if n < N:
            tn1 = 1 - 2 * p[n + 1]
            if tn1 != t:
                B[k] += 1
                if (n & 1) == 0:
                    Beven[k] += 1
                Wn = 0
                if r == 1:
                    Wn = -1; B1[k] += 1
                elif r == 2:
                    Wn = 1; B2[k] += 1
                if Wn != 0:
                    v = Wn * tn1
                    J[k] += v
                    # whiteness: products with the previous signed breaks (lags 1..8)
                    for l in range(1, 9):
                        if nhist >= l:
                            idx = (hpos - l) % 9
                            lagsum[k, l] += v * hist[idx]; lagcnt[k, l] += 1
                    hist[hpos] = v; hpos = (hpos + 1) % 9
                    if nhist < 9:
                        nhist += 1
        if n == (1 << (k + 1)) - 1 or n == N:
            Hend[k] = H
    return B, B1, B2, J, maxH, energy, Hend, Beven, lagsum, lagcnt


@njit(cache=True)
def check_identities(p, N):
    """beta(2m) = [t(m) != t(3m+2)], beta(2m-1) = [t(3m-1) != t(m)] for 2m+1 <= N and 3m+2 <= N; and the exact octave transfer."""
    bad = 0; checked = 0
    for m in range(2, N // 3 - 1):
        tm = p[m]; t3p = p[3 * m + 2]; t3m = p[3 * m - 1]
        # break at 2m: t(2m+1) != t(2m)
        b_even = (p[2 * m + 1] != p[2 * m]); pred_even = (tm != t3p)
        b_odd = (p[2 * m] != p[2 * m - 1]); pred_odd = (t3m != tm)
        if b_even != pred_even or b_odd != pred_odd:
            bad += 1
        checked += 1
    return bad, checked


@njit(cache=True)
def transfer_check(p, N, k):
    """Delta H over octave k+1 = sum_(m in octave k) W(2m)(t(3m+2) - t(m)) + W(2m-1)(t(m) - t(3m-1)), with the boundary terms."""
    lo = 1 << k; hi = 1 << (k + 1)
    direct = 0
    for n in range(2 * lo, 2 * hi):
        r = n % 3; w = 2 if r == 2 else -1
        direct += w * (1 - 2 * p[n])
    # Delta H = H(2hi - 1) - H(2lo - 1) = W(2hi-1) t(2hi-1) - W(2lo-1) t(2lo-1) - sum_(n = 2lo-1)^(2hi-2) W(n)(t(n+1) - t(n)),
    # and the n-sum splits into n = 2m (m in [lo, hi)) with t(2m+1) - t(2m) = -(t(3m+2) - t(m)) and n = 2m-1 (m in [lo, hi))
    # with t(2m) - t(2m-1) = -(t(m) - t(3m-1)).
    via = Wof(2 * hi - 1) * (1 - 2 * p[2 * hi - 1]) - Wof(2 * lo - 1) * (1 - 2 * p[2 * lo - 1])
    for m in range(lo, hi):
        tm = 1 - 2 * p[m]
        t3p = 1 - 2 * p[3 * m + 2]
        t3m = 1 - 2 * p[3 * m - 1]
        via += Wof(2 * m) * (t3p - tm) + Wof(2 * m - 1) * (tm - t3m)
    return direct, via


@njit(cache=True)
def Wof(n):
    r = n % 3
    if r == 0:
        return 0
    return -1 if r == 1 else 1


@njit(cache=True)
def identity6_check(p, N):
    """H(N) = W(N) t(N) - sum_(n<N) W(n)(t(n+1) - t(n)) for all N < Nmax."""
    H = 0; S = 0; bad = 0
    for n in range(1, N):
        r = n % 3; w = 2 if r == 2 else -1
        t = 1 - 2 * p[n]
        H += w * t
        Wn = 0 if r == 0 else (-1 if r == 1 else 1)
        if Wn * t - S != H:
            bad += 1
        S += Wn * ((1 - 2 * p[n + 1]) - t)
    return bad


def main():
    print("== parity table to 2^%d ==" % LOGN, flush=True)
    p = parity_table(N)
    print(" done in %.0fs; t(1..12) = %s" % (time.time() - T0, [1 - 2 * int(p[n]) for n in range(1, 13)]), flush=True)
    B, B1, B2, J, maxH, energy, Hend, Beven, lagsum, lagcnt = octave_stats(p, N, LOGN)
    print("== (2) per-octave break census, current, energy ==")
    print(" k | B_k (breaks) | B_k/2^k | ratio B_k/B_(k-1) | even-n share | B1 | B2 | J_k | J_k^2/(B1+B2) | H(2^(k+1)) | max|H| | max|H|/sqrt(2^(k+1)) | energy/X^2")
    for k in range(8, LOGN):
        X = 1 << k
        nb = B1[k] + B2[k]
        print(" %2d | %9d | %.5f | %.3f | %.3f | %8d | %8d | %7d | %.3f | %8d | %7d | %.3f | %.4f" % (
            k, B[k], B[k] / X, (B[k] / B[k - 1]) if B[k - 1] else 0, Beven[k] / B[k] if B[k] else 0, B1[k], B2[k], J[k],
            (J[k] ** 2 / nb) if nb else 0, Hend[k], maxH[k], maxH[k] / math.sqrt(2 * X), energy[k] / (X * X)))
    # break density exponent
    ks = np.arange(14, LOGN); dens = np.array([B[k] / (1 << k) for k in ks])
    slope = np.polyfit(ks, np.log2(dens), 1)[0]
    print(" break density B_k/2^k over octaves %d..%d: log2-slope %.4f per octave (B(X) ~ X^(%.3f))" % (ks[0], ks[-1], slope, 1 + slope))
    print("== (3) whiteness of the signed break current (autocorrelation over consecutive breaks with W != 0) ==")
    for k in (16, 20, 24, LOGN - 1):
        rho = [(lagsum[k, l] / lagcnt[k, l]) if lagcnt[k, l] else 0.0 for l in range(1, 9)]
        print(" octave %2d: rho(1..8) = %s" % (k, [round(x, 4) for x in rho]))
    print("== (4) exact identities ==")
    print(" identity (6) H(N) = W(N)t(N) - sum_(n<N) W(n)(t(n+1) - t(n)): violations for N < 2^%d: %d" % (LOGN, identity6_check(p, N)))
    bad, checked = check_identities(p, N)
    print(" beta(2m) = [t(m) != t(3m+2)], beta(2m-1) = [t(3m-1) != t(m)]: %d violations in %d checks" % (bad, checked))
    for k in (12, 16, 20):
        d, v = transfer_check(p, N, k)
        print(" octave transfer Delta H_(%d) direct = %d, via the previous octave = %d: %s" % (k + 1, d, v, d == v))
    print("== (5) the twisted-residue formula rho(l) = (-1)^l E[W(n_i) W(n_(i+l))], the gap law of breaks, and a renewal closure ==")
    gaps, cls, cls_next, twist, direct = gap_census(p, N, 1 << 22, N)
    tot = gaps.sum()
    print(" breaks in [2^22, 2^27): %d; gap law g = 1..12: %s; P(g = 0 mod 3) = %.4f, 1 mod 3: %.4f, 2 mod 3: %.4f" % (
        tot, [round(gaps[g] / tot, 4) for g in range(1, 13)], gaps[3::3].sum() / tot, gaps[1::3].sum() / tot, gaps[2::3].sum() / tot))
    print(" class of break position n mod 3: %s; transition matrix of consecutive break classes (rows from, cols to 0,1,2):" % (cls / cls.sum()).round(4).tolist())
    M = cls_next / cls_next.sum(axis=1, keepdims=True)
    print("  ", M.round(4).tolist())
    print(" twisted-residue formula gamma(l) = (-1)^l E[W(c_i) W(c_(i+l))] (lag in breaks) at lags 1..12: %s" % [round(x, 4) for x in twist[1:13]])
    print(" direct autocovariance of the signed current v_i = W t(n_i + 1) at the same lags:            %s (identical: the formula is exact)" % [round(x, 4) for x in direct[1:13]])
    g0 = 2.0 / 3.0
    emp_factor = 1 + 2 * sum(direct[1:13]) / g0
    print(" empirical long-run variance factor (gamma(0) + 2 sum_(l<=12) gamma(l))/gamma(0), gamma(0) = 2/3: %.4f" % emp_factor)
    a = gaps[3::3].sum() / tot; b = gaps[1::3].sum() / tot; c = gaps[2::3].sum() / tot
    om = complex(-0.5, 3 ** 0.5 / 2); what = a + b * om + c * om * om
    gm = [((-1) ** l) * (2.0 / 3.0) * (what ** l).real for l in range(1, 13)]
    print(" renewal model with i.i.d. gaps, P(g mod 3 = 0,1,2) = (%.4f, %.4f, %.4f): gamma(l) = (-1)^l (2/3) Re(w^l), w = a + b omega + c omega^2 = %.4f%+.4fi:" % (a, b, c, what.real, what.imag))
    print("   model gamma(1..12) = %s; closed-form variance factor Re[(1 - w)/(1 + w)] = %.4f" % ([round(x, 4) for x in gm], ((1 - what) / (1 + what)).real))
    rho_model = renewal_rho(M, 12)
    print(" Markov-in-residue model (measured 3x3 transitions): gamma(l)/gamma(0) = %s" % [round(x, 4) for x in rho_model])
    fac = ((1 - what) / (1 + what)).real
    print(" conditional energy (renewal model): E[H(N)^2] = gamma(0) x factor x B(N) = (2/3)(%.3f) B(N) = %.3f B(N) ~ %.3f N at break density 0.24; dyadic energy int_X^(2X) H^2 ~ %.3f X^2" % (
        fac, g0 * fac, g0 * fac * 0.24, g0 * fac * 0.24 * 1.5))
    print("== (6) eps_M = A(2M) + A(M)/2 + A(3M)/2 (Codex identity (7)) is -H(3M)/(6M) exactly; its decay under the break model ==")
    S, Hp = prefix_sums(p, N)
    bad = 0; rows = []
    for j in range(8, LOGN - 1):
        M = 1 << j
        eps = S[2 * M] / (2 * M) + S[M] / (2 * M) + S[3 * M] / (6 * M)
        if abs(eps + Hp[3 * M] / (6.0 * M)) > 1e-12:
            bad += 1
        # rms over 64 prefixes per octave [2^j, 2^(j+1))
        acc = 0.0; cnt = 0
        for i in range(64):
            Mi = M + (i * M) // 64
            if 3 * Mi <= N:
                e = S[2 * Mi] / (2 * Mi) + S[Mi] / (2 * Mi) + S[3 * Mi] / (6 * Mi)
                acc += e * e; cnt += 1
        rows.append((j, eps, math.sqrt(acc / cnt) if cnt else float("nan")))
    print(" identity eps_M = -H(3M)/(6M) checked at the dyadic M: %d violations" % bad)
    print(" j | eps_(2^j) | rms over 64 prefixes in the octave | rms x sqrt(M) | model 0.105/sqrt(M) x sqrt(M) = 0.105")
    for j, e, r in rows:
        print(" %2d | %+.3e | %.3e | %.4f" % (j, e, r, r * math.sqrt(1 << j)))
    def slope(j0, j1):
        xs = [j for j, e, r in rows if j0 <= j <= j1]; ys = [math.log2(r) for j, e, r in rows if j0 <= j <= j1]
        return np.polyfit(xs, ys, 1)[0]
    print(" log2-slopes of the octave rms over 2^12..2^22: %.3f; 2^18..2^25: %.3f; 2^22..2^25: %.3f (S27 reported -0.45, -0.7, -1.2 on other ranges: the range dependence is the walk's own fluctuation; the break model predicts -1/2 with a slowly falling prefactor)" % (
        slope(12, 22), slope(18, 25), slope(22, 25)))
    print("total %.0fs" % (time.time() - T0))


@njit(cache=True)
def prefix_sums(p, N):
    S = np.zeros(N + 1, dtype=np.float64); H = np.zeros(N + 1, dtype=np.float64)
    s = 0.0; h = 0.0
    for n in range(1, N + 1):
        t = 1.0 - 2.0 * p[n]
        s += t; S[n] = s
        r = n % 3
        h += (2.0 if r == 2 else -1.0) * t; H[n] = h
    return S, H


@njit(cache=True)
def gap_census(p, N, lo, hi):
    gaps = np.zeros(64, dtype=np.int64); cls = np.zeros(3, dtype=np.int64); cls_next = np.zeros((3, 3), dtype=np.int64)
    # twisted residue products over consecutive breaks with W != 0, lags 1..8
    # gamma(l) = E[v_i v_(i+l)] over ALL breaks (v = 0 on class 0), lag counted in breaks; the parity flips at every break,
    # so v_i v_(i+l) = (-1)^l W(c_i) W(c_(i+l)); the direct product of the actual signed values is kept as a check
    tw_sum = np.zeros(13, dtype=np.float64); dir_sum = np.zeros(13, dtype=np.float64); tw_cnt = np.zeros(13, dtype=np.int64)
    histW = np.zeros(13, dtype=np.int64); histV = np.zeros(13, dtype=np.int64); hpos = 0; nhist = 0
    last = -1; lastc = -1
    for n in range(lo, hi):
        if p[n + 1] != p[n]:
            c = n % 3
            cls[c] += 1
            if last > 0:
                g = n - last
                if g < 64:
                    gaps[g] += 1
                cls_next[lastc, c] += 1
            last = n; lastc = c
            W = 0 if c == 0 else (-1 if c == 1 else 1)
            v = W * (1 - 2 * p[n + 1])
            for l in range(1, 13):
                if nhist >= l:
                    idx = (hpos - l) % 13
                    tw_sum[l] += W * histW[idx] * (1 if l % 2 == 0 else -1); dir_sum[l] += v * histV[idx]; tw_cnt[l] += 1
            histW[hpos] = W; histV[hpos] = v; hpos = (hpos + 1) % 13
            if nhist < 13:
                nhist += 1
    twist = np.zeros(13); direct = np.zeros(13)
    for l in range(1, 13):
        twist[l] = tw_sum[l] / tw_cnt[l] if tw_cnt[l] else 0.0
        direct[l] = dir_sum[l] / tw_cnt[l] if tw_cnt[l] else 0.0
    return gaps, cls, cls_next, twist, direct


def renewal_rho(M, L):
    """Markov chain on break classes {0,1,2} with transition M; the current takes W = 0, -1, +1 on classes 0, 1, 2 and the
    parity twist (-1)^l between the l-th next break; rho(l) = (-1)^l E[W(c_i) W(c_(i+l))] / E[W^2] over pairs with both W != 0,
    which equals the stationary twisted correlation of the W-chain restricted to nonzero classes."""
    import numpy.linalg as la
    w, v = la.eig(M.T)
    pi = np.real(v[:, np.argmin(abs(w - 1))]); pi = pi / pi.sum()
    Wv = np.array([0.0, -1.0, 1.0])
    out = []
    Ml = np.eye(3)
    for l in range(1, L + 1):
        Ml = Ml @ M
        num = sum(pi[a] * Wv[a] * Ml[a, b] * Wv[b] for a in range(3) for b in range(3))
        den = sum(pi[a] * Wv[a] ** 2 * sum(Ml[a, b] * (Wv[b] != 0) for b in range(3)) for a in range(3))
        out.append(((-1) ** l) * num / den)
    return out


if __name__ == "__main__":
    main()
