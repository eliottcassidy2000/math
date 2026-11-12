#!/usr/bin/env python3
"""
b2crt_check.py -- independent check of the paragraph "The CRT form" and of the passage to the degree
bound (D3) in the proof of (B2) of the DRAFT THM-5141 (verifier lane b2crt; written from scratch,
no code shared with the coordinator scripts).

Everything is exact integer / rational arithmetic (Python int, fractions.Fraction).  The only float
is the right-hand side 12 K (1 + ln K)^2, which is compared with an exact rational on the left.

usage:  python b2crt_check.py W [mode] [KFAC]
    W     window length (>= 4); s1 = ceil(W/2), s2 = floor(W/2), ds = s1 - s2
    mode  "diff": P_u runs over the difference set of {(u, rev u)} with b odd (default)
          "box" : P_u runs over ALL (a, b), |a|, |b| < 2^s1, b odd (a = 0 included)   [same for P_v]
    KFAC  directions (k, k') with |k|, k' <= KFAC * 2^s2 are all tested (default 1), and in addition
          EVERY direction that occurs for some point of the box, whatever its height.

What is checked, for Deg (P_u = (a, b) fixed, (a', b') counted):
  [E1] e := 2^ds k' a' - k a satisfies k b' - k' b = 2^s2 e, |e| < |k| + 2^ds k',
       e = -k a (mod 2^ds k'), 2^s2 e = -k' b (mod |k|)                       (every point of the box)
  [E2] the two congruences have exactly one solution e0 modulo m = 2^ds k' |k| (brute force for |k| <= BR)
  [E3] bijection: {(a', b') in the box with direction (k, k')}  =  {(a', b') computed from the integers
       e = e0 (mod m), |e| < |k| + 2^ds k', that give integers with |a'|, |b'| < 2^s2}; e -> (a', b') injective
  [E4] e0/m = -sgn(k) a/(2^ds k') - b inv(2^s1; |k|)/|k|  (mod 1), and
       ||e0/m|| = ||a/(2^ds k') + theta_k||, theta_k = sgn(k) b inv(2^s1; |k|)/|k|
  [E5] count <= 6 M/m + 1 and (count > 0  =>  ||e0/m|| < 3 M/m) with M = max(|k|, k')  (hence for all K >= M)
  [S1] for K = 1..KS:  Deg_K <= sum_{(k,k')} 6K/(|k| k') + #{(k,k') : ||.|| < 3K/(2^ds |k| k')}
  [S2]                 ... <= 12 K (1 + ln K)^2 + sum_{odd k} #{odd k' <= K : ||.|| < 3K/(|k| k')}
  [S3]                 ... <= 12 K (1 + ln K)^2 + sum_{odd k} (1 + sum_{L = 2^i < K} #{odd k' in (L, 2L] : ||.|| <= 3K/(|k| L)})
and the mirror statements for Deg' (P_v = (a', b') fixed, (a, b) counted) with e := k a - 2^ds k' a',
m' = |k| k', condition ||2^ds a'/k + theta'_(k')|| < 3K/(|k| k'), theta'_(k') = b' inv(2^s2; k')/k'.
"""
import sys
import math
import time
from fractions import Fraction
from math import gcd

BR = 70  # brute-force uniqueness of the CRT class when the second modulus is <= BR


def rev(x, w):
    r = 0
    for i in range(w):
        if (x >> i) & 1:
            r |= 1 << (w - 1 - i)
    return r


def diffset(s):
    """multiplicities r_s(a, b) = #{(u, u') : u - u' = a, rev u - rev u' = b}"""
    R = [rev(u, s) for u in range(1 << s)]
    D = {}
    for u in range(1 << s):
        for u2 in range(1 << s):
            key = (u - u2, R[u] - R[u2])
            D[key] = D.get(key, 0) + 1
    return D


def frac1(fr):
    """fractional part in [0, 1) of a Fraction"""
    return fr - (fr.numerator // fr.denominator)


def dist(fr):
    f = frac1(fr)
    return min(f, 1 - f)


def sgn(x):
    return 1 if x > 0 else -1


def inv_mod(x, n):
    """inverse of x modulo n >= 1 (0 when n = 1)"""
    if n == 1:
        return 0
    return pow(x, -1, n)


def direction(X, Xi):
    g = gcd(abs(X), abs(Xi))
    kp = abs(Xi) // g
    k = (X // g) * sgn(Xi)
    return k, kp


class Fail(Exception):
    pass


def check(cond, msg):
    if not cond:
        raise Fail(msg)


def odd_coprime_dirs(Kmax):
    out = []
    for kp in range(1, Kmax + 1, 2):
        for ka in range(1, Kmax + 1, 2):
            if gcd(ka, kp) == 1:
                out.append((ka, kp))
                out.append((-ka, kp))
    return out


def run_deg(W, mode, kfac, out):
    """Deg: P_u = (a, b) fixed, (a', b') in the s2-box counted."""
    s1 = (W + 1) // 2
    s2 = W // 2
    ds = s1 - s2
    B1 = 1 << s1
    B2 = 1 << s2
    Kmax = kfac * B2
    if mode == "diff":
        DU = diffset(s1)
        PU = sorted(p for p in DU if p[1] % 2 != 0)
        check(all(p[0] != 0 for p in PU), "difference point with b odd and a = 0")
    else:
        PU = [(a, b) for a in range(-B1 + 1, B1) for b in range(-B1 + 1, B1) if b % 2 != 0]
    alldirs = odd_coprime_dirs(Kmax)
    KS = list(range(1, B2 + 1)) + [B2 + 1, B2 + 3, 2 * B2 - 1, 2 * B2] if kfac >= 2 else list(range(1, B2 + 1))
    KS = [K for K in KS if K <= Kmax]
    Kd = 2 * Kmax  # theta table needed for odd k' up to the end of the last dyadic block
    stats = dict(points=0, dirs_tested=0, dirs_nonempty=0, max_count=0, e_in_range=0, e_valid=0,
                 e_rejected_range_a=0, e_rejected_range_b=0, brute_unique=0,
                 worst_ratio_6K=Fraction(0), max_ratio_D3={}, max_ratio_S2={}, zero_despite_cond=0,
                 kabs1=0, kneg=0)
    for (a, b) in PU:
        # ---- direct enumeration of the box
        dirs = {}
        for ap in range(-B2 + 1, B2):
            X0 = ap * B1 + b
            for bp in range(-B2 + 1, B2):
                Xi = a * B2 + bp
                if bp % 2 == 0:
                    continue  # Xi even: not an odd pair
                k, kp = direction(X0, Xi)
                # [E1]
                e = (1 << ds) * kp * ap - k * a
                T = abs(k) + (1 << ds) * kp
                check(k * bp - kp * b == B2 * e, f"E1 identity a,b={a,b} a',b'={ap,bp} k,k'={k,kp}")
                check(abs(e) < T, f"E1 |e|<|k|+2^ds k' fails a,b={a,b} a',b'={ap,bp} k,k'={k,kp} e={e}")
                check((e + k * a) % ((1 << ds) * kp) == 0, "E1 congruence 1")
                check((B2 * e + kp * b) % abs(k) == 0, "E1 congruence 2")
                check(k % 2 != 0 and kp % 2 == 1 and gcd(abs(k), kp) == 1, "direction not odd coprime")
                dirs.setdefault((k, kp), set()).add((ap, bp))
                stats["points"] += 1
        # ---- every direction in range, and every occurring direction
        todo = set(alldirs) | set(dirs.keys())
        for (k, kp) in todo:
            m1 = (1 << ds) * kp
            m2 = abs(k)
            m = m1 * m2
            r1 = (-k * a) % m1
            if m2 <= BR:
                sols = [r1 + t * m1 for t in range(m2) if ((r1 + t * m1) * B2 + kp * b) % m2 == 0]
                check(len(sols) == 1, f"E2 uniqueness a,b={a,b} k,k'={k,kp}: {sols}")
                e0 = sols[0]
                stats["brute_unique"] += 1
            else:
                r2 = (-kp * b * inv_mod(B2, m2)) % m2
                t = ((r2 - r1) * inv_mod(m1, m2)) % m2
                e0 = r1 + t * m1
            check((e0 + k * a) % m1 == 0 and (B2 * e0 + kp * b) % m2 == 0 and 0 <= e0 < m, "E2 e0")
            # [E4] the displayed formula
            iv = inv_mod(B1, m2)
            lhs = Fraction(e0, m)
            rhs = -sgn(k) * Fraction(a, m1) - Fraction(b * iv, m2)
            check((lhs - rhs).denominator == 1,
                  f"E4 formula fails a,b={a,b} k,k'={k,kp} e0/m={lhs} rhs={rhs}")
            theta = Fraction(sgn(k) * b * iv, m2)
            d0 = dist(lhs)
            check(d0 == dist(Fraction(a, m1) + theta), f"E4 theta form fails a,b={a,b} k,k'={k,kp}")
            if m2 == 1:
                stats["kabs1"] += 1
                check(frac1(theta) == 0, "theta for |k|=1")
            if k < 0:
                stats["kneg"] += 1
            # [E3] bijection
            T = abs(k) + m1
            e = -T + 1 + ((e0 - (-T + 1)) % m)
            got = set()
            ne = 0
            while e < T:
                ne += 1
                na = e + k * a
                nb = B2 * e + kp * b
                check(na % m1 == 0 and nb % k == 0, "E3 divisibility")
                ap = na // m1
                bp = nb // k
                if abs(ap) < B2 and abs(bp) < B2:
                    check(bp % 2 != 0, "E3 b' odd")
                    check((ap, bp) not in got, "E3 injectivity")
                    got.add((ap, bp))
                else:
                    if abs(ap) >= B2:
                        stats["e_rejected_range_a"] += 1
                    if abs(bp) >= B2:
                        stats["e_rejected_range_b"] += 1
                e += m
            stats["e_in_range"] += ne
            stats["e_valid"] += len(got)
            direct = dirs.get((k, kp), set())
            check(got == direct,
                  f"E3 BIJECTION FAILS W={W} a,b={a,b} k,k'={k,kp} e0={e0} m={m}: from e {sorted(got)} direct {sorted(direct)}")
            # [E5]
            M = max(abs(k), kp)
            cnt = len(direct)
            check(ne * m <= 2 * T + m, "E5 number of e in the open interval > 2T/m + 1")
            check(cnt * m <= 6 * M + m, f"E5 count > 6M/m + 1: a,b={a,b} k,k'={k,kp} cnt={cnt}")
            if cnt > 0:
                check(d0 * m < 3 * M, f"E5 nonempty but ||e0/m|| >= 3M/m: a,b={a,b} k,k'={k,kp}")
                check(d0 * m < T, "E5 sharper: nonempty but ||e0/m|| >= T/m")
                stats["dirs_nonempty"] += 1
                stats["max_count"] = max(stats["max_count"], cnt)
                stats["worst_ratio_6K"] = max(stats["worst_ratio_6K"], Fraction(cnt * m, 6 * M + m))
            else:
                if d0 * m < 3 * M:
                    stats["zero_despite_cond"] += 1
            stats["dirs_tested"] += 1
        # ---- the sums [S1], [S2], [S3]
        # table of ||a/(2^ds k') + theta_k|| for odd |k| <= Kmax, odd k' <= Kd
        dtab = {}
        for ka in range(1, Kmax + 1, 2):
            iv = inv_mod(B1, ka)
            for sg in (1, -1):
                th = Fraction(sg * b * iv, ka)
                for kp in range(1, Kd + 1, 2):
                    dtab[(sg * ka, kp)] = dist(Fraction(a, (1 << ds) * kp) + th)
        for K in KS:
            degK = sum(len(v) for (k, kp), v in dirs.items() if max(abs(k), kp) <= K)
            main = Fraction(0)
            c1 = 0
            c2 = 0
            c3 = 0
            for ka in range(1, K + 1, 2):
                for sg in (1, -1):
                    k = sg * ka
                    for kp in range(1, K + 1, 2):
                        d = dtab[(k, kp)]
                        if gcd(ka, kp) == 1:
                            main += Fraction(6 * K, ka * kp)
                            if d * ((1 << ds) * ka * kp) < 3 * K:
                                c1 += 1
                        if d * (ka * kp) < 3 * K:
                            c2 += 1
                    # blocks
                    c3 += 1  # k' = 1
                    L = 1
                    while L < K:
                        for kp in range(L + 1, 2 * L + 1):
                            if kp % 2 == 1:
                                # every k' of the block has tolerance 3K/(|k| k') < 3K/(|k| L)
                                check(Fraction(3 * K, ka * kp) < Fraction(3 * K, ka * L), "block tolerance")
                                if dtab[(k, kp)] * (ka * L) <= 3 * K:
                                    c3 += 1
                        L *= 2
            check(degK <= main + c1, f"S1 fails W={W} a,b={a,b} K={K}: Deg={degK} main={main} c1={c1}")
            lg = 12 * K * (1 + math.log(K)) ** 2
            check(float(main) <= lg + 1e-9, f"S2 main term W={W} K={K}: {float(main)} > {lg}")
            check(c1 <= c2, "S2 count")
            check(c2 <= c3, f"S3 block count W={W} a,b={a,b} K={K}: c2={c2} c3={c3}")
            if a != 0:
                br = K * (1 + math.log(K)) ** 2 + K * (1 + math.log(K)) * abs(a) ** (1 / 3) + K ** 2.5 * abs(a) ** -0.5
                r = degK / br
                if r > stats["max_ratio_D3"].get(K, (0,))[0]:
                    stats["max_ratio_D3"][K] = (r, (a, b), degK)
            r2 = (float(main) + c2) / max(degK, 1)
            if K not in stats["max_ratio_S2"] or r2 < stats["max_ratio_S2"][K][0]:
                stats["max_ratio_S2"][K] = (r2, (a, b), degK, c2)
    out.append(f"[Deg ] W={W} s1={s1} s2={s2} ds={ds} mode={mode} Kmax={Kmax}: P_u tested {len(PU)}; "
               f"box points {stats['points']}; directions tested {stats['dirs_tested']} "
               f"(nonempty {stats['dirs_nonempty']}, with k=+-1: {stats['kabs1']}, with k<0: {stats['kneg']}); "
               f"CRT class unique by brute force in {stats['brute_unique']} cases; "
               f"e in range {stats['e_in_range']}, giving box points {stats['e_valid']} "
               f"(rejected: |a'|>=2^s2 {stats['e_rejected_range_a']}, |b'|>=2^s2 {stats['e_rejected_range_b']}); "
               f"max points on one direction {stats['max_count']}; max count/(6M/m+1) = {float(stats['worst_ratio_6K']):.4f}; "
               f"directions with ||e0/m|| < 3M/m but no point: {stats['zero_despite_cond']}")
    for K in sorted(stats["max_ratio_D3"]):
        if K & (K - 1) == 0 or K == KS[-1]:
            r, p, d = stats["max_ratio_D3"][K]
            r2, p2, d2, c2 = stats["max_ratio_S2"][K]
            out.append(f"        K={K:5d}: max Deg_K/[K(1+lnK)^2+K(1+lnK)|a|^(1/3)+K^2.5|a|^-.5] = {r:.4f} at P_u={p} (Deg_K={d}); "
                       f"min (S2 bound)/Deg_K = {r2:.2f}")
    return stats


def run_degp(W, mode, kfac, out):
    """Deg': P_v = (a', b') fixed, (a, b) in the s1-box counted."""
    s1 = (W + 1) // 2
    s2 = W // 2
    ds = s1 - s2
    B1 = 1 << s1
    B2 = 1 << s2
    Kmax = kfac * B2
    if mode == "diff":
        DV = diffset(s2)
        PV = sorted(p for p in DV if p[1] % 2 != 0)
        check(all(p[0] != 0 for p in PV), "difference point with b' odd and a' = 0")
    else:
        PV = [(a, b) for a in range(-B2 + 1, B2) for b in range(-B2 + 1, B2) if b % 2 != 0]
    alldirs = odd_coprime_dirs(Kmax)
    KS = list(range(1, B2 + 1)) + [B2 + 1, B2 + 3, 2 * B2 - 1, 2 * B2] if kfac >= 2 else list(range(1, B2 + 1))
    KS = [K for K in KS if K <= Kmax]
    Kd = 2 * Kmax
    stats = dict(points=0, dirs_tested=0, dirs_nonempty=0, max_count=0, e_in_range=0, e_valid=0,
                 e_rejected_range_a=0, e_rejected_range_b=0, brute_unique=0,
                 worst_ratio_6K=Fraction(0), max_ratio_D3={}, zero_despite_cond=0, kabs1=0, kneg=0)
    for (ap, bp) in PV:
        dirs = {}
        for a in range(-B1 + 1, B1):
            Xi = a * B2 + bp
            for b in range(-B1 + 1, B1):
                if b % 2 == 0:
                    continue  # X even
                X = ap * B1 + b
                k, kp = direction(X, Xi)
                e = k * a - (1 << ds) * kp * ap
                T = abs(k) + (1 << ds) * kp
                check(kp * b == k * bp + B2 * e, f"E1' identity a',b'={ap,bp} a,b={a,b}")
                check(abs(e) < T, f"E1' |e| bound fails a',b'={ap,bp} a,b={a,b} k,k'={k,kp} e={e}")
                check((e + (1 << ds) * kp * ap) % abs(k) == 0, "E1' congruence 1")
                check((B2 * e + k * bp) % kp == 0, "E1' congruence 2")
                check(k % 2 != 0 and kp % 2 == 1 and gcd(abs(k), kp) == 1, "direction not odd coprime")
                dirs.setdefault((k, kp), set()).add((a, b))
                stats["points"] += 1
        todo = set(alldirs) | set(dirs.keys())
        for (k, kp) in todo:
            m1 = abs(k)
            m2 = kp
            m = m1 * m2
            r1 = (-(1 << ds) * kp * ap) % m1
            if m2 <= BR:
                sols = [r1 + t * m1 for t in range(m2) if ((r1 + t * m1) * B2 + k * bp) % m2 == 0]
                check(len(sols) == 1, f"E2' uniqueness a',b'={ap,bp} k,k'={k,kp}: {sols}")
                e0 = sols[0]
                stats["brute_unique"] += 1
            else:
                r2 = (-k * bp * inv_mod(B2, m2)) % m2
                t = ((r2 - r1) * inv_mod(m1, m2)) % m2
                e0 = r1 + t * m1
            check((e0 + (1 << ds) * kp * ap) % m1 == 0 and (B2 * e0 + k * bp) % m2 == 0 and 0 <= e0 < m, "E2' e0")
            # [E4'] the condition of the text: ||2^ds a'/k + theta'_(k')||, theta' = b' inv(2^s2; k')/k'
            iv = inv_mod(B2, kp)
            thp = Fraction(bp * iv, kp)
            lhs = Fraction(e0, m)
            rhs = -sgn(k) * (Fraction((1 << ds) * ap, k) + thp)
            check((lhs - rhs).denominator == 1,
                  f"E4' formula fails a',b'={ap,bp} k,k'={k,kp} e0/m={lhs} rhs={rhs}")
            d0 = dist(lhs)
            check(d0 == dist(Fraction((1 << ds) * ap, k) + thp), "E4' theta form")
            # the form used for the lemma: sign * 2^ds a' / |k| + theta'
            check(d0 == dist(Fraction(sgn(k) * (1 << ds) * ap, abs(k)) + thp), "E4' signed form")
            if m1 == 1:
                stats["kabs1"] += 1
            if k < 0:
                stats["kneg"] += 1
            T = abs(k) + (1 << ds) * kp
            e = -T + 1 + ((e0 - (-T + 1)) % m)
            got = set()
            ne = 0
            while e < T:
                ne += 1
                na = e + (1 << ds) * kp * ap
                nb = k * bp + B2 * e
                check(na % k == 0 and nb % kp == 0, "E3' divisibility")
                a = na // k
                b = nb // kp
                if abs(a) < B1 and abs(b) < B1:
                    check(b % 2 != 0, "E3' b odd")
                    check((a, b) not in got, "E3' injectivity")
                    got.add((a, b))
                else:
                    if abs(a) >= B1:
                        stats["e_rejected_range_a"] += 1
                    if abs(b) >= B1:
                        stats["e_rejected_range_b"] += 1
                e += m
            stats["e_in_range"] += ne
            stats["e_valid"] += len(got)
            direct = dirs.get((k, kp), set())
            check(got == direct,
                  f"E3' BIJECTION FAILS W={W} a',b'={ap,bp} k,k'={k,kp} e0={e0} m={m}: from e {sorted(got)} direct {sorted(direct)}")
            M = max(abs(k), kp)
            cnt = len(direct)
            check(ne * m <= 2 * T + m, "E5' number of e")
            check(cnt * m <= 6 * M + m, f"E5' count > 6M/m + 1: a',b'={ap,bp} k,k'={k,kp} cnt={cnt}")
            if cnt > 0:
                check(d0 * m < 3 * M, f"E5' nonempty but ||e0/m|| >= 3M/m: a',b'={ap,bp} k,k'={k,kp}")
                check(d0 * m < T, "E5' sharper")
                stats["dirs_nonempty"] += 1
                stats["max_count"] = max(stats["max_count"], cnt)
                stats["worst_ratio_6K"] = max(stats["worst_ratio_6K"], Fraction(cnt * m, 6 * M + m))
            else:
                if d0 * m < 3 * M:
                    stats["zero_despite_cond"] += 1
            stats["dirs_tested"] += 1
        # sums: table of ||sign 2^ds a'/j + theta'_(k')|| for odd j <= Kd, odd k' <= Kmax
        dtab = {}
        for kp in range(1, Kmax + 1, 2):
            thp = Fraction(bp * inv_mod(B2, kp), kp)
            for j in range(1, Kd + 1, 2):
                for sg in (1, -1):
                    dtab[(sg * j, kp)] = dist(Fraction(sg * (1 << ds) * ap, j) + thp)
        for K in KS:
            degK = sum(len(v) for (k, kp), v in dirs.items() if max(abs(k), kp) <= K)
            main = Fraction(0)
            c1 = 0
            c2 = 0
            c3 = 0
            for kp in range(1, K + 1, 2):
                for sg in (1, -1):
                    for j in range(1, K + 1, 2):
                        d = dtab[(sg * j, kp)]
                        if gcd(j, kp) == 1:
                            main += Fraction(6 * K, j * kp)
                            if d * (j * kp) < 3 * K:
                                c1 += 1
                        if d * (j * kp) < 3 * K:
                            c2 += 1
                    c3 += 1  # |k| = 1
                    L = 1
                    while L < K:
                        for j in range(L + 1, 2 * L + 1):
                            if j % 2 == 1:
                                if dtab[(sg * j, kp)] * (L * kp) <= 3 * K:
                                    c3 += 1
                        L *= 2
            check(degK <= main + c1, f"S1' fails W={W} a',b'={ap,bp} K={K}: Deg'={degK} main={main} c1={c1}")
            lg = 12 * K * (1 + math.log(K)) ** 2
            check(float(main) <= lg + 1e-9, "S2' main term")
            check(c1 <= c2 and c2 <= c3, f"S2'/S3' counts W={W} a',b'={ap,bp} K={K}")
            if ap != 0:
                br = K * (1 + math.log(K)) ** 2 + K * (1 + math.log(K)) * abs(ap) ** (1 / 3) + K ** 2.5 * abs(ap) ** -0.5
                r = degK / br
                if r > stats["max_ratio_D3"].get(K, (0,))[0]:
                    stats["max_ratio_D3"][K] = (r, (ap, bp), degK)
    out.append(f"[Deg'] W={W} s1={s1} s2={s2} ds={ds} mode={mode} Kmax={Kmax}: P_v tested {len(PV)}; "
               f"box points {stats['points']}; directions tested {stats['dirs_tested']} "
               f"(nonempty {stats['dirs_nonempty']}, with k=+-1: {stats['kabs1']}, with k<0: {stats['kneg']}); "
               f"CRT class unique by brute force in {stats['brute_unique']} cases; "
               f"e in range {stats['e_in_range']}, giving box points {stats['e_valid']} "
               f"(rejected: |a|>=2^s1 {stats['e_rejected_range_a']}, |b|>=2^s1 {stats['e_rejected_range_b']}); "
               f"max points on one direction {stats['max_count']}; max count/(6M/m+1) = {float(stats['worst_ratio_6K']):.4f}; "
               f"directions with ||e0/m|| < 3M/m but no point: {stats['zero_despite_cond']}")
    for K in sorted(stats["max_ratio_D3"]):
        if K & (K - 1) == 0 or K == KS[-1]:
            r, p, d = stats["max_ratio_D3"][K]
            out.append(f"        K={K:5d}: max Deg'_K/[K(1+lnK)^2+K(1+lnK)|a'|^(1/3)+K^2.5|a'|^-.5] = {r:.4f} at P_v={p} (Deg'_K={d})")
    return stats


def main():
    W = int(sys.argv[1])
    mode = sys.argv[2] if len(sys.argv) > 2 else "diff"
    kfac = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    which = sys.argv[4] if len(sys.argv) > 4 else "both"
    out = []
    t0 = time.time()
    ok = True
    try:
        if which in ("both", "deg"):
            run_deg(W, mode, kfac, out)
        if which in ("both", "degp"):
            run_degp(W, mode, kfac, out)
    except Fail as ex:
        ok = False
        out.append("FAIL: " + str(ex))
    out.append(f"W={W} mode={mode} kfac={kfac} which={which}: {'ALL CHECKS PASSED' if ok else 'FAILED'} ({time.time() - t0:.1f} s)")
    print("\n".join(out))
    sys.stdout.flush()


if __name__ == "__main__":
    main()
