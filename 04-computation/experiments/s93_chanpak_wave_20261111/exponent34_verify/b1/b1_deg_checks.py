#!/usr/bin/env python3
"""Independent adversarial check of THM-5141 (draft): section "The half split", inequalities (1), (2),
(D1), (D2), and the arithmetic of (B1).

Everything here is written from the definitions in the theorem file, not from the coordinator's scripts.

Objects (w = s1 + s2, s1 = ceil(w/2), s2 = floor(w/2)):
  x = v 2^s1 + y, u = rev_s1(y);   a = u-u', b = y-y', a' = v-v', b' = rev_s2 v - rev_s2 v'
  X = x - x' ,  Xi = rev_w x - rev_w x'
  claimed: X = a' 2^s1 + b,  Xi = a 2^s2 + b',  #pairs with (a,b,a',b') = r_s1(a,b) r_s2(a',b')
  M = max(|X|,|Xi|)/gcd(X,Xi)
  Deg_K(a,b)   = #{(a',b') : |a'|,|b'| < 2^s2, b' odd, M <= K}
  Deg'_K(a',b') = #{(a,b)  : |a|,|b|  < 2^s1, b odd,  M <= K}

Usage: python b1_deg_checks.py WMIN WMAX [wide]
"""
import os
import sys
import math
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
import numpy as np
from numba import njit, prange


def rev(x, w):
    r = 0
    for i in range(w):
        if (x >> i) & 1:
            r |= 1 << (w - 1 - i)
    return r


def r_table(s):
    """R[a+off, b+off] = #{(u,u') in [0,2^s)^2 : u-u' = a, rev_s u - rev_s u' = b}."""
    n = 1 << s
    off = n - 1
    R = np.zeros((2 * n - 1, 2 * n - 1), dtype=np.int64)
    rv = [rev(u, s) for u in range(n)]
    for u in range(n):
        for u2 in range(n):
            R[u - u2 + off, rv[u] - rv[u2] + off] += 1
    return R


@njit(cache=True)
def gcd_(a, b):
    while b != 0:
        a, b = b, a % b
    return a


@njit(parallel=True, cache=True)
def deg_scan(s_out, s_in, R_out, R_in, pmax, qmax, Mcap):
    """Outer points (p,q): |p| <= pmax, q odd, |q| <= qmax.  Inner points (p',q'): |p'|,|q'| < 2^s_in, q' odd.
    The two numbers are  A = p' 2^s_out + q  and  B = p 2^s_in + q'  (for Deg: s_out=s1, s_in=s2, A = X, B = Xi;
    for Deg': s_out=s2, s_in=s1, A = Xi, B = X).   M = max(|A|,|B|)/gcd(A,B).
    Returns per value of p (row index ip = p + pmax):
      maxdeg[ip, K]      = max over q of #{inner : M <= K}
      maxdeg_supp[ip, K] = the same max restricted to outer points in the support of r_{s_out}
      histN[ip, M]       = sum over q and inner points of r_out r_in at height exactly M
      histQ[ip, M]       = sum over q and ALL inner box points of r_out^2 at height exactly M
    """
    n_out = 1 << s_out
    n_in = 1 << s_in
    off_out = n_out - 1
    off_in = n_in - 1
    npv = 2 * pmax + 1
    maxdeg = np.zeros((npv, Mcap), dtype=np.int64)
    maxdeg_supp = np.zeros((npv, Mcap), dtype=np.int64)
    histN = np.zeros((npv, Mcap), dtype=np.int64)
    histQ = np.zeros((npv, Mcap), dtype=np.int64)
    for ip in prange(npv):
        p = ip - pmax
        hist = np.zeros(Mcap, dtype=np.int64)
        for q in range(-qmax, qmax + 1):
            if (q & 1) == 0:
                continue
            ro = 0
            if abs(p) < n_out and abs(q) < n_out:
                ro = R_out[p + off_out, q + off_out]
            for t in range(Mcap):
                hist[t] = 0
            for pp in range(-(n_in - 1), n_in):
                A = pp * n_out + q
                aA = abs(A)
                for qq in range(-(n_in - 1), n_in):
                    if (qq & 1) == 0:
                        continue
                    B = p * n_in + qq
                    aB = abs(B)
                    g = gcd_(aA, aB)
                    if aA > aB:
                        M = aA // g
                    else:
                        M = aB // g
                    hist[M] += 1
                    if ro > 0:
                        ri = R_in[pp + off_in, qq + off_in]
                        histN[ip, M] += ro * ri
                        histQ[ip, M] += ro * ro
            c = 0
            for K in range(Mcap):
                c += hist[K]
                if c > maxdeg[ip, K]:
                    maxdeg[ip, K] = c
                if ro > 0 and c > maxdeg_supp[ip, K]:
                    maxdeg_supp[ip, K] = c
    return maxdeg, maxdeg_supp, histN, histQ


@njit(parallel=True, cache=True)
def max_dir_mult(s_out, s_in):
    """For every outer box point (p,q) (q odd) group the inner box points by the primitive direction of (A,B)
    and return the largest group size (= the largest number of j for one direction), with the K = max(|k|,k')
    of a direction attaining it, per p."""
    n_out = 1 << s_out
    n_in = 1 << s_in
    w = s_out + s_in
    OFF = 1 << w
    BIG = 1 << (w + 2)
    npv = 2 * n_out - 1
    best = np.zeros(npv, dtype=np.int64)
    bestM = np.zeros(npv, dtype=np.int64)
    n_inner = (2 * n_in - 1) * n_in
    for ip in prange(npv):
        p = ip - (n_out - 1)
        keys = np.empty(n_inner, dtype=np.int64)
        for q in range(-(n_out - 1), n_out):
            if (q & 1) == 0:
                continue
            t = 0
            for pp in range(-(n_in - 1), n_in):
                A = pp * n_out + q
                for qq in range(-(n_in - 1), n_in):
                    if (qq & 1) == 0:
                        continue
                    B = p * n_in + qq
                    g = gcd_(abs(A), abs(B))
                    if B < 0:
                        g = -g
                    ka = A // g
                    kb = B // g
                    keys[t] = (ka + OFF) * BIG + kb
                    t += 1
            keys.sort()
            run = 1
            for i in range(1, n_inner):
                if keys[i] == keys[i - 1]:
                    run += 1
                    if run > best[ip]:
                        best[ip] = run
                        ka = keys[i] // BIG - OFF
                        kb = keys[i] % BIG
                        bestM[ip] = max(abs(ka), kb)
                else:
                    run = 1
            if best[ip] == 0:
                best[ip] = 1
    return best, bestM


@njit(parallel=True, cache=True)
def triple_scan(s_out, s_in, Kmax):
    """(D2) mechanism.  For every outer box point (p,q), q odd, and every K <= Kmax:
        Deg_K   = #{inner box points with M <= K}
        T_K     = #{(q', d, e) : |q'| < 2^s_in odd, d >= 1 a divisor of B = p 2^s_in + q', d <= K, j = B/d,
                                |e| <= K, j e = q (mod 2^s_out)}
    Returns  viol = #{(outer, K) : Deg_K > T_K},  maxT[K] = max over outer of T_K,
             maxe = max over (outer, q', d, K) of the number of e,   maxdeg[K].
    """
    n_out = 1 << s_out
    n_in = 1 << s_in
    npv = 2 * n_out - 1
    viol = np.zeros(npv, dtype=np.int64)
    maxe = np.zeros(npv, dtype=np.int64)
    maxT = np.zeros((npv, Kmax + 1), dtype=np.int64)
    maxD = np.zeros((npv, Kmax + 1), dtype=np.int64)
    mask = n_out - 1
    for ip in prange(npv):
        p = ip - (n_out - 1)
        T = np.zeros(Kmax + 1, dtype=np.int64)
        D = np.zeros(Kmax + 1, dtype=np.int64)
        for q in range(-(n_out - 1), n_out):
            if (q & 1) == 0:
                continue
            for t in range(Kmax + 1):
                T[t] = 0
                D[t] = 0
            # exact degree
            for pp in range(-(n_in - 1), n_in):
                A = pp * n_out + q
                aA = abs(A)
                for qq in range(-(n_in - 1), n_in):
                    if (qq & 1) == 0:
                        continue
                    B = p * n_in + qq
                    aB = abs(B)
                    g = gcd_(aA, aB)
                    if aA > aB:
                        M = aA // g
                    else:
                        M = aB // g
                    if M <= Kmax:
                        D[M] += 1
            # triples
            for qq in range(-(n_in - 1), n_in):
                if (qq & 1) == 0:
                    continue
                B = p * n_in + qq
                aB = abs(B)
                for d in range(1, Kmax + 1):
                    if aB % d != 0:
                        continue
                    j = B // d          # exact, signed, odd
                    # c = q * j^{-1} mod 2^s_out  (brute force inverse)
                    jm = j & mask
                    inv = 1
                    for cand in range(1, n_out, 2):
                        if (cand * jm) & mask == 1:
                            inv = cand
                            break
                    if n_out == 1:
                        c = 0
                    else:
                        c = (q * inv) & mask
                    for K in range(d, Kmax + 1):
                        cnt = (K - c) // n_out + (K + c) // n_out + 1
                        if c == 0:
                            cnt = 2 * (K // n_out) + 1
                        T[K] += cnt
                        if cnt > maxe[ip]:
                            maxe[ip] = cnt
            cD = 0
            for K in range(1, Kmax + 1):
                cD += D[K]
                if cD > T[K]:
                    viol[ip] += 1
                if T[K] > maxT[ip, K]:
                    maxT[ip, K] = T[K]
                if cD > maxD[ip, K]:
                    maxD[ip, K] = cD
    return viol, maxe, maxT, maxD


@njit(parallel=True, cache=True)
def height_hist(w, revt):
    """Direct enumeration: histogram of M over the odd pairs (x,x') with top bit of x equal to 1 and top bit of
    x' equal to 0 (the other half is the mirror image (x',x), same M)."""
    n = 1 << w
    half = n >> 1
    nchunk = 64
    hist = np.zeros((nchunk, n), dtype=np.int64)
    for ch in prange(nchunk):
        lo = half + (half * ch) // nchunk
        hi = half + (half * (ch + 1)) // nchunk
        for x in range(lo, hi):
            rx = revt[x]
            x0 = x & 1
            for xp in range(1 - x0, half, 2):
                X = x - xp
                Xi = rx - revt[xp]
                aXi = abs(Xi)
                g = gcd_(X, aXi)
                if X > aXi:
                    M = X // g
                else:
                    M = aXi // g
                hist[ch, M] += 1
    out = np.zeros(n, dtype=np.int64)
    for ch in range(nchunk):
        for m in range(n):
            out[m] += hist[ch, m]
    return out


def check_half_split_identities(w):
    """Brute-force check of rev_w(x) = u 2^s2 + rev_s2(v), X = a' 2^s1 + b, Xi = a 2^s2 + b', oddness,
    multiplicity r_s1 r_s2, on ALL ordered pairs."""
    s1 = (w + 1) // 2
    s2 = w // 2
    n = 1 << w
    R1 = r_table(s1)
    R2 = r_table(s2)
    o1 = (1 << s1) - 1
    o2 = (1 << s2) - 1
    cnt = {}
    bad = 0
    comp = []
    for x in range(n):
        v, y = divmod(x, 1 << s1)
        u = rev(y, s1)
        rx = rev(x, w)
        if rx != u * (1 << s2) + rev(v, s2):
            bad += 1
        comp.append((u, y, v, rev(v, s2), rx))
    if len({(c[0], c[2]) for c in comp}) != n:
        bad += 1
    for x in range(n):
        u, y, v, rv, rx = comp[x]
        for x2 in range(n):
            u2, y2, v2, rv2, rx2 = comp[x2]
            a = u - u2
            b = y - y2
            a_ = v - v2
            b_ = rv - rv2
            X = x - x2
            Xi = rx - rx2
            if X != a_ * (1 << s1) + b or Xi != a * (1 << s2) + b_:
                bad += 1
            if b != rev(u, s1) - rev(u2, s1):
                bad += 1
            odd_pair = (X & 1) and (Xi & 1)
            if bool(odd_pair) != bool((b & 1) and (b_ & 1)):
                bad += 1
            if w >= 2 and bool(odd_pair) != bool(((x ^ x2) & 1) and (((x ^ x2) >> (w - 1)) & 1)):
                bad += 1
            if (b & 1) and a == 0:
                bad += 1
            if (b_ & 1) and a_ == 0:
                bad += 1
            key = (a, b, a_, b_)
            cnt[key] = cnt.get(key, 0) + 1
    for (a, b, a_, b_), c in cnt.items():
        if c != R1[a + o1, b + o1] * R2[a_ + o2, b_ + o2]:
            bad += 1
    # and the product vanishes off the image
    tot = int((R1.astype(object).sum()) * (R2.astype(object).sum()))
    if tot != n * n:
        bad += 1
    return bad


def main():
    wmin = int(sys.argv[1])
    wmax = int(sys.argv[2])
    wide = len(sys.argv) > 3 and sys.argv[3] == "wide"
    print("THM-5141 half split / (1) / (2) / (D1) / (D2) / (B1): independent checks", flush=True)
    print("numba threads:", os.environ.get("NUMBA_NUM_THREADS"), flush=True)
    total_viol = 0
    for w in range(wmin, wmax + 1):
        t0 = time.time()
        s1 = (w + 1) // 2
        s2 = w // 2
        n1 = 1 << s1
        n2 = 1 << s2
        R1 = r_table(s1)
        R2 = r_table(s2)
        e1 = int((R1.astype(object) ** 2).sum())
        e2 = int((R2.astype(object) ** 2).sum())
        print(f"\n=== w = {w}  s1 = {s1}  s2 = {s2} ===", flush=True)
        print(f"energy r_s1: {e1} (6^s1 = {6**s1}) {'OK' if e1 == 6**s1 else 'FAIL'};  "
              f"energy r_s2: {e2} (6^s2 = {6**s2}) {'OK' if e2 == 6**s2 else 'FAIL'};  "
              f"support sizes {int((R1 > 0).sum())} (3^s1 = {3**s1}), {int((R2 > 0).sum())} (3^s2 = {3**s2})")
        if e1 != 6 ** s1 or e2 != 6 ** s2:
            total_viol += 1
        if w <= 8:
            bad = check_half_split_identities(w)
            print(f"half-split identities on all 4^w ordered pairs: {bad} failures")
            total_viol += bad

        # ---------------- box scans ----------------
        Mcap = 1 << w
        mdU, mdUs, hN_u, hQ_u = deg_scan(s1, s2, R1, R2, n1 - 1, n1 - 1, Mcap)
        mdV, mdVs, hN_v, hQ_v = deg_scan(s2, s1, R2, R1, n2 - 1, n2 - 1, Mcap)
        degU = mdU.max(axis=0)      # max over all box P_u (b odd, any a incl. 0)
        degUs = mdUs.max(axis=0)    # max over support of r_s1
        degV = mdV.max(axis=0)
        degVs = mdVs.max(axis=0)
        # N(K) from the split and from the mirrored split
        NK_u = np.cumsum(hN_u.sum(axis=0).astype(object))
        NK_v = np.cumsum(hN_v.sum(axis=0).astype(object))
        # right side of (1)
        RHS1 = (np.cumsum(hQ_u.sum(axis=0).astype(object)) + np.cumsum(hQ_v.sum(axis=0).astype(object)))
        # direct enumeration
        revt = np.array([rev(x, w) for x in range(1 << w)], dtype=np.int64)
        hh = height_hist(w, revt)
        NK = np.cumsum((2 * hh).astype(object))
        n_odd = int(NK[-1])
        print(f"number of odd pairs: {n_odd} (4^(w-1) = {4**(w-1)}) {'OK' if n_odd == 4**(w-1) else 'FAIL'}")
        same = all(int(NK[K]) == int(NK_u[K]) == int(NK_v[K]) for K in range(Mcap))
        print(f"N(K) via the split (both orientations) == N(K) by direct enumeration for every K: "
              f"{'OK' if same else 'FAIL'}")
        if not same or n_odd != 4 ** (w - 1):
            total_viol += 1
        if hh[0] != 0 or any(hh[m] != 0 for m in range(0, Mcap, 2)):
            print("FAIL: even height found")
            total_viol += 1

        # ---- (D1) for every integer K ----
        v_d1 = 0
        v_d1p = 0
        worst_d1 = (0.0, 0)
        worst_d1p = (0.0, 0)
        for K in range(1, Mcap):
            if degU[K] > 4 * K * K:
                v_d1 += 1
            if degV[K] > 8 * K * K:
                v_d1p += 1
            r_ = degU[K] / (4 * K * K)
            if r_ > worst_d1[0]:
                worst_d1 = (r_, K)
            r_ = degV[K] / (8 * K * K)
            if r_ > worst_d1p[0]:
                worst_d1p = (r_, K)
        print(f"(D1) Deg_K <= 4K^2 for all K in [1,2^w): violations {v_d1}; worst ratio {worst_d1[0]:.4f} at K = {worst_d1[1]}")
        print(f"(D1) Deg'_K <= 8K^2 for all K in [1,2^w): violations {v_d1p}; worst ratio {worst_d1p[0]:.4f} at K = {worst_d1p[1]}")
        total_viol += v_d1 + v_d1p

        # ---- (D2) for every integer K <= 2^s2 ----
        v_d2 = 0
        v_d2p = 0
        v_mid = 0
        worst_d2 = (0.0, 0)
        worst_d2p = (0.0, 0)
        for K in range(1, n2 + 1):
            bnd = 3 * n1 * (3 + 2 * math.log(K))
            midU = 3 * (2 * n2 * (1 + math.log(K)) + K)
            midV = 3 * (2 * n1 * (1 + math.log(K)) + K)
            if degU[K] > bnd:
                v_d2 += 1
            if degV[K] > bnd:
                v_d2p += 1
            if degU[K] > midU or degV[K] > midV or midU > bnd * (1 + 1e-12) or midV > bnd * (1 + 1e-12):
                v_mid += 1
            if degU[K] / bnd > worst_d2[0]:
                worst_d2 = (degU[K] / bnd, K)
            if degV[K] / bnd > worst_d2p[0]:
                worst_d2p = (degV[K] / bnd, K)
        print(f"(D2) Deg_K <= 3*2^s1(3+2 ln K) for all K in [1,2^s2]: violations {v_d2}; worst ratio {worst_d2[0]:.4f} at K = {worst_d2[1]}")
        print(f"(D2) Deg'_K <= 3*2^s1(3+2 ln K) for all K in [1,2^s2]: violations {v_d2p}; worst ratio {worst_d2p[0]:.4f} at K = {worst_d2p[1]}")
        print(f"(D2) intermediate bounds 3(2^(s2+1)(1+ln K)+K) resp. 3(2^(s1+1)(1+ln K)+K), and their comparison with the final one: violations {v_mid}")
        total_viol += v_d2 + v_d2p + v_mid
        # is the range K <= 2^s2 needed?  first K (if any) where the (D2) expression fails
        firstfailU = next((K for K in range(1, Mcap) if degU[K] > 3 * n1 * (3 + 2 * math.log(K))), None)
        firstfailV = next((K for K in range(1, Mcap) if degV[K] > 3 * n1 * (3 + 2 * math.log(K))), None)
        print(f"     (for information) first K at which the (D2) expression is exceeded: Deg: {firstfailU}, Deg': {firstfailV}")

        # ---- (1) ----
        v1 = 0
        worst1 = 0.0
        for K in range(1, Mcap):
            if 2 * int(NK[K]) > int(RHS1[K]):
                v1 += 1
            if RHS1[K] > 0:
                worst1 = max(worst1, 2 * int(NK[K]) / int(RHS1[K]))
        print(f"(1) N(K) <= (1/2)[sum r_s1^2 Deg_K + sum r_s2^2 Deg'_K] for all K: violations {v1}; worst ratio N/RHS {worst1:.4f}")
        total_viol += v1

        # ---- N(K) bound of (B1) ----
        vB = 0
        vB_all = 0
        worstB = (0.0, 0)
        for K in range(1, Mcap):
            b_d1 = 6 ** s1 * 8 * K * K
            if int(NK[K]) > b_d1:
                vB_all += 1
            if K <= n2:
                b_d2 = 6 ** s1 * 3 * n1 * (3 + 2 * math.log(K))
                bb = min(b_d1, b_d2)
                if int(NK[K]) > bb:
                    vB += 1
                if int(NK[K]) / bb > worstB[0]:
                    worstB = (int(NK[K]) / bb, K)
        print(f"(B1) N(K) <= 6^s1 min(8K^2, 3*2^s1(3+2 ln K)) for all K in [1,2^s2]: violations {vB}; worst ratio {worstB[0]:.5f} at K = {worstB[1]};  N(K) <= 8K^2 6^s1 for all K: violations {vB_all}")
        total_viol += vB + vB_all

        # ---- table at K = 2^i ----
        print("   i      K   maxDeg  (supp)   4K^2  maxDeg' (supp)   8K^2    D2bound        N(K)      RHS(1)/2   6^s1*min(...)")
        for i in range(0, w + 1):
            K = 1 << i
            Kc = min(K, Mcap - 1)
            d2 = 3 * n1 * (3 + 2 * math.log(K)) if K <= n2 else float('nan')
            bb = min(6 ** s1 * 8 * K * K, 6 ** s1 * d2) if K <= n2 else 6 ** s1 * 8 * K * K
            print(f"  {i:2d} {K:6d} {int(degU[Kc]):8d} {int(degUs[Kc]):7d} {4*K*K:6d} {int(degV[Kc]):8d} {int(degVs[Kc]):6d} {8*K*K:6d} {d2:10.1f} {int(NK[Kc]):11d} {int(RHS1[Kc])/2:13.1f} {bb:14.1f}")

        # ---- (2) and the (B1) chain for S_w ----
        from fractions import Fraction
        S = Fraction(0)
        for m in range(1, Mcap, 2):
            if hh[m]:
                S += Fraction(2 * int(hh[m]), m)
        Sf = float(S)
        rhs2 = Fraction(0)
        for i in range(0, 4 * w + 8):
            K = 1 << (i + 1)
            Nv = int(NK[min(K, Mcap - 1)])
            rhs2 += Fraction(Nv, 1 << i)
        # exact tail: beyond i = 4w+8 the terms are 4^(w-1) 2^-i ; add the geometric remainder
        rhs2 += Fraction(n_odd, 1 << (4 * w + 7))
        i0 = s1 // 2
        mid = 6 ** s1 * (sum(32 * 2 ** i for i in range(0, i0 + 1)) + Fraction(3 * n1 * (3 + 2 * s2), 2 ** i0)) \
            + Fraction(2 * 4 ** w, 2 ** s2)
        fin = (6 * math.sqrt(2)) ** s1 * (82 + 12 * s2) + 2 ** (1.5 * w + 2)
        ok2 = S <= rhs2
        okm = float(rhs2) <= float(mid) or w < 4
        okf = float(mid) <= fin * (1 + 1e-12) or w < 4
        print(f"S_w = {Sf:.6f}   S_w/6^(w/2) = {Sf / 6 ** (w / 2):.4f}   (2) RHS = {float(rhs2):.3f} [S_w <= RHS(2): {'OK' if ok2 else 'FAIL'}]")
        print(f"   (B1) displayed line 1 = {float(mid):.3f}  [RHS(2) <= line 1: {'OK' if okm else 'FAIL'}];  line 2 = (6 sqrt2)^s1 (82+12 s2) + 2^(3w/2+2) = {fin:.3f} [line 1 <= line 2: {'OK' if okf else 'FAIL'}];  S_w/line 2 = {Sf / fin:.5f}")
        print(f"   S_w <= explicit bound: {'OK' if Sf <= fin else 'FAIL'};   line2 / (w (6 sqrt 2)^(w/2)) = {fin / (w * (6 * math.sqrt(2)) ** (w / 2)):.3f}")
        if not ok2 or Sf > fin:
            total_viol += 1
        if w >= 4 and (not okm or not okf):
            total_viol += 1
        # pairs of height 1
        h1 = 2 * int(hh[1])
        print(f"   pairs of height 1: {h1};  h1 / 6^(w/2) = {h1 / 6 ** (w / 2):.4f}"
              f"  (closed form: {4 * 6 ** (w // 2 - 1) if w % 2 == 0 else 2 * 6 ** ((w - 1) // 2)})")

        # ---- multiplicity of j per direction ----
        if w <= 12:
            bj, bjM = max_dir_mult(s1, s2)
            bjp, bjpM = max_dir_mult(s2, s1)
            print(f"max number of j for one direction: Deg: {int(bj.max())} (claimed <= 2), Deg': {int(bjp.max())} (claimed <= 4)")
            if bj.max() > 2 or bjp.max() > 4:
                total_viol += 1

        # ---- (D2) triple mechanism ----
        if w <= 12 and s2 >= 1:
            viol, maxe, maxT, maxD = triple_scan(s1, s2, n2)
            violp, maxep, maxTp, maxDp = triple_scan(s2, s1, n2)
            tb = max(maxT.max(axis=0)[K] / (3 * sum(2 * n2 / d + 1 for d in range(1, K + 1))) for K in range(1, n2 + 1))
            tbp = max(maxTp.max(axis=0)[K] / (3 * sum(2 * n1 / d + 1 for d in range(1, K + 1))) for K in range(1, n2 + 1))
            print(f"(D2) triples: Deg_K > T_K in {int(viol.sum())} cases, Deg'_K > T'_K in {int(violp.sum())} cases; "
                  f"max number of k per (b',k'): {int(maxe.max())} (claimed <= 3); of k' per (b,k): {int(maxep.max())} (claimed <= 3); "
                  f"max T_K / [3 sum_(k'<=K)(2^(s2+1)/k'+1)] = {tb:.4f}; max T'_K / [3 sum_(k<=K)(2^(s1+1)/k+1)] = {tbp:.4f}")
            chk = all(int(maxD.max(axis=0)[K]) == int(degU[K]) for K in range(1, n2 + 1)) and \
                all(int(maxDp.max(axis=0)[K]) == int(degV[K]) for K in range(1, n2 + 1))
            print(f"     consistency of the two degree computations: {'OK' if chk else 'FAIL'}")
            if viol.sum() or violp.sum() or maxe.max() > 3 or maxep.max() > 3 or tb > 1 or tbp > 1 or not chk:
                total_viol += 1

        # ---- wide scan: P_u, P_v outside the box (the proofs of (D1),(D2) do not use the box for the fixed point) ----
        if wide and w <= 10:
            pm = 4 * n1 + 3
            McapW = max((n2 - 1) * n1 + pm, pm * n2 + n2 - 1) + 2
            wU, _, _, _ = deg_scan(s1, s2, R1, R2, pm, pm, McapW)
            pm2 = 4 * n2 + 3
            McapW2 = max((n1 - 1) * n2 + pm2, pm2 * n1 + n1 - 1) + 2
            wV, _, _, _ = deg_scan(s2, s1, R2, R1, pm2, pm2, McapW2)
            dU = wU.max(axis=0)
            dV = wV.max(axis=0)
            vw = 0
            for K in range(1, min(McapW, McapW2)):
                if dU[K] > 4 * K * K or dV[K] > 8 * K * K:
                    vw += 1
                if K <= n2:
                    bnd = 3 * n1 * (3 + 2 * math.log(K))
                    if dU[K] > bnd or dV[K] > bnd:
                        vw += 1
            print(f"wide scan (fixed point up to 4*2^s+3 in both coordinates): (D1)/(D2) violations {vw}")
            total_viol += vw
        print(f"   [{time.time() - t0:.1f} s]", flush=True)
    print(f"\nTOTAL VIOLATIONS: {total_viol}")


if __name__ == "__main__":
    main()
