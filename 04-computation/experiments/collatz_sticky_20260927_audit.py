#!/usr/bin/env python3
"""Independent audit script for
  05-knowledge/results/collatz_sticky_20260927_size_coupled_persistence.md
    (Theorem 1, Propositions 2, A, B; the numbers of sections 2, 3, 5)
and for Propositions 1-3 of
  05-knowledge/results/collatz_two_carries_typology_20260927.md
    (the two carries, the free cofactor, the lifting-the-exponent law).

Auditor subagent (adversarial, blind re-derivation first), 2026-09-27.
Everything arithmetic is exact integer arithmetic (numpy int64 sieves, Python big
integers); floats appear only in the reported averages.  Nothing here uses the
audited scripts' code or seeds.

Run:  python 04-computation/experiments/collatz_sticky_20260927_audit.py
      > 05-knowledge/results/collatz_sticky_20260927_audit.out
"""
from __future__ import annotations

import math
import random
import time
from collections import Counter
from itertools import product

import numpy as np

T0 = time.time()
M_SIEVE = 2 * 10 ** 7          # divisor-sum sieve bound (needed by the Erdos chains)
N_SPF = 10 ** 6                # smallest-prime-factor sieve bound (factorizations)
L3 = math.log2(3)


def stamp() -> str:
    return f"[{time.time() - T0:6.1f}s]"


# ----------------------------------------------------------------------------
# exact helpers
# ----------------------------------------------------------------------------

def v2(x: int) -> int:
    """2-adic valuation of a positive integer."""
    return (x & -x).bit_length() - 1


def v2_arr(x: np.ndarray) -> np.ndarray:
    """2-adic valuation of a positive int64 array (exact: frexp of the lowest set bit)."""
    low = (x & -x).astype(np.float64)
    _, e = np.frexp(low)
    return (e - 1).astype(np.int64)


def isqrt_arr(m: np.ndarray) -> np.ndarray:
    r = np.sqrt(m.astype(np.float64)).astype(np.int64)
    r[r * r > m] -= 1
    r[(r + 1) * (r + 1) <= m] += 1
    return r


def sigma_sieve(N: int) -> np.ndarray:
    s = np.zeros(N + 1, dtype=np.int64)
    for d in range(1, N + 1):
        s[d::d] += d
    return s


def spf_sieve(N: int) -> np.ndarray:
    spf = np.zeros(N + 1, dtype=np.int64)
    for i in range(2, math.isqrt(N) + 1):
        if spf[i] == 0:
            seg = spf[i * i::i]
            seg[seg == 0] = i
    idx = np.arange(N + 1, dtype=np.int64)
    return np.where(spf == 0, idx, spf)


def factor_spf(m: int, spf: np.ndarray) -> dict:
    f = {}
    while m > 1:
        p = int(spf[m])
        e = 0
        while m % p == 0:
            m //= p
            e += 1
        f[p] = e
    return f


def U(m: int) -> tuple[int, int]:
    """Syracuse step on odd m: (3m+1)/2^v, v."""
    x = 3 * m + 1
    v = v2(x)
    return x >> v, v


# ----------------------------------------------------------------------------
# A. Theorem 1(b): the valuation of sigma(p^e)
# ----------------------------------------------------------------------------

def part_A(sig: np.ndarray, spf: np.ndarray) -> None:
    print("== A. Theorem 1(b): v2(sigma(p^e)) = v2(p+1) + v2(e+1) - 1 for e odd, 0 for e even ==")
    primes = [p for p in range(3, 200, 2) if all(p % q for q in range(3, math.isqrt(p) + 1, 2))]
    n_cases = bad = 0
    for p in primes:
        for e in range(1, 16):
            s = sum(p ** i for i in range(e + 1))
            pred = (v2(p + 1) + v2(e + 1) - 1) if e % 2 else 0
            n_cases += 1
            if v2(s) != pred:
                bad += 1
                print(f"   MISMATCH p={p} e={e}: v2={v2(s)} formula={pred}")
    print(f"   odd primes p < 200, exponents 1 <= e <= 15: {n_cases} cases, {bad} mismatches")
    for p, e in [(3, 3), (3, 5), (3, 7), (5, 3), (5, 5), (5, 7), (7, 7), (13, 3), (31, 7), (17, 15), (3, 15), (199, 9)]:
        s = sum(p ** i for i in range(e + 1))
        print(f"   p={p:3d} (p mod 8 = {p % 8}), e={e:2d}: v2(sigma(p^e)) = {v2(s)}, formula {v2(p + 1) + v2(e + 1) - 1}")
    # multiplicativity: v2(sigma(m)) from the factorization against the sieve, all odd m < 10^6
    bad = 0
    for m in range(1, N_SPF, 2):
        tot = 0
        for p, e in factor_spf(m, spf).items():
            if e % 2:
                tot += v2(p + 1) + v2(e + 1) - 1
        if tot != v2(int(sig[m])):
            bad += 1
    print(f"   v2(sigma(m)) = sum over odd-exponent p^e || m of (v2(p+1)+v2(e+1)-1), all odd m < 10^6: {bad} mismatches")
    print("   remark on the note's PROOF of (b): '1 + p^2 + ... + p^(e-1) is a sum of (e+1)/2 odd terms, so its v2 is v2((e+1)/2)'")
    print("   is not a valid inference by itself (1 + 3 = 4 is a sum of 2 odd terms with v2 = 2 != v2(2) = 1); it holds here because")
    print("   the terms are = 1 mod 8 (pins v2 when (e+1)/2 is not divisible by 8) and in general by LTE:")
    print("   1 + p^2 + ... + p^(2k-2) = (p^(2k) - 1)/(p^2 - 1), v2 = [v2(p-1)+v2(p+1)+v2(2k)-1] - [v2(p-1)+v2(p+1)] = v2(k).")
    e = 15
    p = 3
    s = sum(p ** (2 * i) for i in range((e + 1) // 2))
    print(f"   example needing LTE: p=3, e=15: 1+9+...+9^7 = {s} = 8 * {s // 8}, v2 = {v2(s)} = v2(8) (8 | (e+1)/2, so 'mod 8' gives only >= 3)")
    print(stamp())


# ----------------------------------------------------------------------------
# B. Theorem 1(a) and the tenth note's Proposition 2: the valuation identity
# ----------------------------------------------------------------------------

def part_B(sig: np.ndarray) -> None:
    print("\n== B. Theorem 1(a) / tenth-note Proposition 2: s(2^a m) = (2^(a+1)-1) sigma(m) - 2^a m and its valuation, ALL even n <= 2*10^7 ==")
    tot = equal = 0
    ok = [True] * 7
    chunk = 2 * 10 ** 6
    for lo in range(2, M_SIEVE + 1, chunk):
        n = np.arange(lo, min(lo + chunk, M_SIEVE + 1), dtype=np.int64)
        n = n[n % 2 == 0]
        a = v2_arr(n)
        m = n >> a
        sm = sig[m]
        s = sig[n] - n
        ok[0] &= bool(np.all(s == (2 ** (a + 1) - 1) * sm - 2 ** a * m))   # the algebraic identity
        ok[1] &= bool(np.all(s > 0))
        vs = v2_arr(sm)
        vsn = v2_arr(s)
        diff = vs != a
        ok[2] &= bool(np.all(vsn[diff] == np.minimum(a, vs)[diff]))
        ok[3] &= bool(np.all(vsn[~diff] >= (a + 1)[~diff]))
        lost = vsn < a
        kept = vsn == a
        stren = vsn >= a + 1
        ok[4] &= bool(np.all(lost == (vs <= a - 1))) and bool(np.all(kept == (vs > a))) and bool(np.all(stren == (vs == a)))
        r = isqrt_arr(m)
        a1 = a == 1
        ok[5] &= bool(np.all(lost[a1] == ((r * r) == m)[a1]))
        # sigma itself: v2(sigma(n)) = v2(sigma(m)) (the "without the -n" clause)
        ok[6] &= bool(np.all(v2_arr(sig[n]) == vs))
        tot += len(n)
        equal += int((~diff).sum())
    print(f"   {tot} even n: s(n) = (2^(a+1)-1)sigma(m) - 2^a m: {ok[0]}; s(n) > 0: {ok[1]}")
    print(f"   v2 s = min(a, v2 sigma m) when they differ: {ok[2]}; v2 s >= a+1 when equal ({equal} equal cases): {ok[3]}")
    print(f"   lost iff v2 sigma(m) <= a-1, kept iff > a, strengthened iff = a: {ok[4]}; a=1 lost iff m is a square: {ok[5]}")
    print(f"   v2(sigma(n)) = v2(sigma(m)) (sigma(2^a) odd): {ok[6]}")
    print(stamp())


# ----------------------------------------------------------------------------
# C. Theorem 1(c),(d): class rates in [N, 2N) by sieve
# ----------------------------------------------------------------------------

def part_C(sig: np.ndarray, spf: np.ndarray) -> None:
    print("\n== C. Theorem 1(c),(d): persistence / loss by driver class in [N, 2N), exact by sieve (N up to 10^7) ==")
    c1 = 2 - math.sqrt(2)
    c2 = math.pi ** 2 / 8
    print(f"   constants: 2 - sqrt 2 = {c1:.6f}; pi^2/8 = {c2:.6f} (predicted limit of loss(a=2) * ln N, see below)")
    print("   N | a=1: keep, loss, odd-square count/odd count, loss*sqrt(N) | a=2: keep, loss, loss*lnN | a=3: keep, loss, loss*lnN/lnlnN | a=4: keep, loss | entry P(v2 sigma m = 1 | a=1), *lnN")
    for N in [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]:
        n = np.arange(N, 2 * N, 2, dtype=np.int64)
        a = v2_arr(n)
        m = n >> a
        vs = v2_arr(sig[m])
        vsn = v2_arr(sig[n] - n)
        stats = {}
        for aa in (1, 2, 3, 4):
            sel = a == aa
            cnt = int(sel.sum())
            stats[aa] = (cnt, int((vsn[sel] == aa).sum()) / cnt, int((vsn[sel] < aa).sum()) / cnt)
        lo, hi = N // 2, N
        nsq = sum(1 for r in range(math.isqrt(lo - 1) + 1, math.isqrt(hi - 1) + 1) if r % 2 == 1)
        nodd = (hi - lo) // 2
        entry = int(((vs == 1) & (a == 1)).sum()) / stats[1][0]
        lnN = math.log(N)
        print(f"   {N:>8} | {stats[1][1]:.4f} {stats[1][2]:.5f} {nsq}/{nodd}={nsq / nodd:.5f} {stats[1][2] * math.sqrt(N):.4f} | "
              f"{stats[2][1]:.4f} {stats[2][2]:.4f} {stats[2][2] * lnN:.3f} | {stats[3][1]:.4f} {stats[3][2]:.4f} {stats[3][2] * lnN / math.log(lnN):.3f} | "
              f"{stats[4][1]:.4f} {stats[4][2]:.4f} | {entry:.4f} {entry * lnN:.3f}")
    # the a = 2 loss set, exactly: v2 sigma(m) <= 1 iff m is a square or m = p^e r^2 with p = 1 mod 4, e = 1 mod 4, p not | r
    print("   a=2 loss set characterization (odd m < 10^6): v2 sigma(m) <= 1  iff  m square, or exactly one odd-exponent prime p^e with p = 1 (mod 4), e = 1 (mod 4):")
    bad = 0
    for m in range(1, N_SPF, 2):
        f = factor_spf(m, spf)
        odd = [(p, e) for p, e in f.items() if e % 2]
        pred = (len(odd) == 0) or (len(odd) == 1 and odd[0][0] % 4 == 1 and odd[0][1] % 4 == 1)
        if pred != (v2(int(sig[m])) <= 1):
            bad += 1
    print(f"      mismatches: {bad}")
    print("   decomposition of the a=2 loss count over odd m in [N/4, N/2) (the class a=2 of [N, 2N)) and the pi^2/8 prediction:")
    for N in [10 ** 4, 10 ** 5, 10 ** 6]:
        lo, hi = N // 4, N // 2
        loss = sq = main = higher = 0
        for m in range(lo | 1, hi, 2):
            if v2(int(sig[m])) <= 1:
                loss += 1
                f = factor_spf(m, spf)
                odd = [(p, e) for p, e in f.items() if e % 2]
                if not odd:
                    sq += 1
                elif odd[0][1] == 1:
                    main += 1
                else:
                    higher += 1
        nodd = (hi - lo) // 2
        lnN = math.log(N)
        print(f"      N={N:>8}: loss {loss}/{nodd} = {loss / nodd:.4f}; squares {sq}, m = p r^2 (e=1) {main}, e >= 5 {higher}; "
              f"loss*lnN = {loss / nodd * lnN:.3f}, (pi^2/8)/lnN = {c2 / lnN:.4f}, main-term prediction sum_(r odd<=sqrt) pi(N/r^2;4,1)-pi(N/(2r^2);4,1) computed next")
    # the main term exactly: count primes p = 1 mod 4 with N/(2r^2) <= p < N/r^2 over odd r, p not dividing r, for N = 10^6
    N = 10 ** 6
    is_prime = spf[:N] == np.arange(N)
    is_prime[:2] = False
    pr = np.flatnonzero(is_prime)
    pr1 = pr[pr % 4 == 1]
    cnt = 0
    r = 1
    while r * r <= N // 2:
        # N/4 <= p r^2 < N/2  <=>  ceil(N/(4 r^2)) <= p <= (N/2 - 1) // r^2
        lo = (N // 4 + r * r - 1) // (r * r)
        hi = (N // 2 - 1) // (r * r)
        ps = pr1[(pr1 >= lo) & (pr1 <= hi)]
        cnt += int(sum(1 for p in ps if r % int(p)))
        r += 2
    print(f"      N=10^6 main term recount as sum over odd r of #{{p = 1 (4) prime, p not | r, N/4 <= p r^2 < N/2}}: {cnt} (the e=1 count above is {12948 if cnt == 12948 else 'see above'}); "
          f"the limit (pi^2/8)/ln N follows from pi(x;4,1) ~ x/(2 ln x) and sum over odd r of 1/r^2 = pi^2/8")
    print("      note: the entry rate P(v2 sigma m = 1) (class 1 -> classes >= 2) and the class-2 loss P(v2 sigma m <= 1) are the SAME event up to the squares,")
    print("      measured on different m-ranges ([M/2, M) vs [N/4, N/2)); the note's '1.3/ln N' and '1.4/ln N' are finite-range fits of one quantity with limit (pi^2/8)/ln N = 1.2337/ln N")
    print(stamp())


# ----------------------------------------------------------------------------
# D. Theorem 1(e): drifts; sign stability of the even-n average by range
# ----------------------------------------------------------------------------

def part_D(sig: np.ndarray) -> None:
    print("\n== D. Theorem 1(e) and the 'negative average drift' of the aliquot map ==")
    N2 = 2 * 10 ** 6
    n = np.arange(2, N2 + 1, dtype=np.int64)
    s = sig[n] - n
    x = np.log2(s / n)
    even = n % 2 == 0
    print(f"   n <= 2*10^6: E[log2(s(n)/n)]: all {x.mean():+.3f}, even {x[even].mean():+.3f}, odd {x[~even].mean():+.3f}; "
          f"P(s(n) > n) = {(s > n).mean():.4f}; odd n with s(n) >= n: {int(((s >= n) & ~even).sum())}")
    ne = n[even]
    a = v2_arr(ne)
    xe = x[even]
    print("   conditional drifts E[log2(s/n) | a] on even n <= 2*10^6: " + ", ".join(f"a={aa}: {xe[a == aa].mean():+.3f}" for aa in range(1, 7)))
    print("   the same by dyadic range [N, 2N) (even n), to test the sign of the even-n average:")
    for N in [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 5 * 10 ** 6, 10 ** 7]:
        n = np.arange(N, 2 * N, 2, dtype=np.int64)
        s = sig[n] - n
        x = np.log2(s / n)
        a = v2_arr(n)
        print(f"      N={N:>8}: even average {x.mean():+.4f}; classes a=1..5: " + ", ".join(f"{x[a == aa].mean():+.3f}" for aa in range(1, 6))
              + f"; P(s>n | even) = {(s > n).mean():.4f}")
    n = np.arange(2, M_SIEVE + 1, 2, dtype=np.int64)
    s = sig[n] - n
    x = np.log2(s / n)
    print(f"   even n <= 2*10^7: average {x.mean():+.4f}")
    ab = 0
    for lo in range(1, M_SIEVE + 1, 4 * 10 ** 6):
        nn = np.arange(lo, min(lo + 4 * 10 ** 6, M_SIEVE + 1), dtype=np.int64)
        ab += int((sig[nn] > 2 * nn).sum())
    n2 = np.arange(1, N2 + 1, dtype=np.int64)
    print(f"   abundant density: to 2*10^6 {int((sig[n2] > 2 * n2).sum()) / N2:.4f}; to 2*10^7 {ab / M_SIEVE:.4f} (Deleglise: 0.2474 < d < 0.2480)")
    print(stamp())


# ----------------------------------------------------------------------------
# E. Proposition 2: Collatz memorylessness, residues vs orbits
# ----------------------------------------------------------------------------

def orbit_stats(nmax: int, cuts: list[int]):
    nb = len(cuts) + 1
    tot = [0] * nb
    keep = [0] * nb
    totp = [0] * nb
    keepp = [0] * nb
    for n in range(1, nmax, 2):
        m = n
        prev = 0
        mp = 0
        while m != 1:
            x = 3 * m + 1
            v = (x & -x).bit_length() - 1
            if prev == 1:
                b = 0
                for c in cuts:
                    if m >= c:
                        b += 1
                    else:
                        break
                tot[b] += 1
                bp = 0
                for c in cuts:
                    if mp >= c:
                        bp += 1
                    else:
                        break
                totp[bp] += 1
                if v == 1:
                    keep[b] += 1
                    keepp[bp] += 1
            prev = v
            mp = m
            m = x >> v
    return tot, keep, totp, keepp


def count_res(N: int, r: int, q: int) -> int:
    """#{n in [N, 2N) : n = r mod q}."""
    first = N + (r - N) % q
    return len(range(first, 2 * N, q))


def part_E() -> None:
    print("\n== E. Proposition 2: (v1, v2) of odd n are functions of n mod 8; residue counts; orbit-weighted persistence by size ==")
    n = np.arange(1, 10 ** 6, 2, dtype=np.int64)
    x1 = 3 * n + 1
    v1 = v2_arr(x1)
    m1 = x1 >> v1
    v2b = v2_arr(3 * m1 + 1)
    ok1 = bool(np.all((v1 == 1) == (n % 4 == 3)))
    sel = v1 == 1
    ok2 = bool(np.all((v2b[sel] == 1) == (n[sel] % 8 == 7)))
    # the full first-valuation law: v1 = k iff n = 3^(-1)(2^k - 1) mod 2^(k+1)
    ok3 = True
    for k in range(1, 12):
        rk = (pow(3, -1, 2 ** (k + 1)) * (2 ** k - 1)) % 2 ** (k + 1)
        ok3 &= bool(np.all((v1 == k) == (n % 2 ** (k + 1) == rk)))
    print(f"   odd n < 10^6: v1 = 1 iff n = 3 (mod 4): {ok1}; given v1 = 1, v2 = 1 iff n = 7 (mod 8): {ok2}; v1 = k iff n = 3^(-1)(2^k-1) mod 2^(k+1), k <= 11: {ok3}")
    for N in (10 ** 3, 10 ** 5, 10 ** 7, 10 ** 9):
        n3 = count_res(N, 3, 4)
        n7 = count_res(N, 7, 8)
        print(f"   [N, 2N), N = {N:>10}: #n = 3 mod 4: {n3}, #n = 7 mod 8: {n7}, ratio {n7 / n3:.6f}")
    for N in (10 ** 3, 10 ** 5):
        nn = np.arange(N + 1, 2 * N, 2, dtype=np.int64)
        xx = 3 * nn + 1
        vv1 = v2_arr(xx)
        vv2 = v2_arr(3 * (xx >> vv1) + 1)
        s1 = vv1 == 1
        print(f"   counting measure on odd n in [{N}, {2 * N}), valuations computed directly: P(v2 = 1 | v1 = 1) = {int((vv2[s1] == 1).sum())}/{int(s1.sum())}")
    # the orbit-weighted numbers of the note (odd n < 2*10^5) with cuts 10^4, 10^5, 10^6, banded by the current value m
    # (the value whose valuation is the 'next' one, as in the note's script) and, for comparison, by the predecessor value
    for nmax in (2 * 10 ** 5, 10 ** 6):
        cuts = [10 ** 4, 10 ** 5, 10 ** 6]
        tot, keep, totp, keepp = orbit_stats(nmax, cuts)
        names = ["[1,10^4)", "[10^4,10^5)", "[10^5,10^6)", "[10^6,inf)"]
        print(f"   orbits of odd n < {nmax}: P(v_(j+1) = 1 | v_j = 1) banded by the CURRENT value m_(j+1) (as the note): "
              + "; ".join(f"{names[b]}: {keep[b] / tot[b]:.4f} ({tot[b]})" for b in range(4) if tot[b]))
        print(f"      banded by the PREDECESSOR value m_j: "
              + "; ".join(f"{names[b]}: {keepp[b] / totp[b]:.4f} ({totp[b]})" for b in range(4) if totp[b]))
        for thr_idx, thr in ((1, 10 ** 4), (2, 10 ** 5)):
            lo_t = sum(tot[:thr_idx]); lo_k = sum(keep[:thr_idx])
            hi_t = sum(tot[thr_idx:]); hi_k = sum(keep[thr_idx:])
            print(f"      split at {thr}: below {lo_k / lo_t:.4f} ({lo_t} steps), at or above {hi_k / hi_t:.4f} ({hi_t} steps); overall {(lo_k + hi_k) / (lo_t + hi_t):.4f}")
    print(stamp())


# ----------------------------------------------------------------------------
# F. Proposition A: the sticky chain, its stationary law, the ergodic slope
# ----------------------------------------------------------------------------

def part_F() -> None:
    print("\n== F. Proposition A: stationary sticky valuation chains ==")
    print("   kernel P(v -> w) = p [w = v] + (1-p) 2^(-w): for any law pi, (pi P)(w) = p pi(w) + (1-p) 2^(-w), so pi P = pi iff pi = geometric(1/2):")
    print("   the geometric law is the UNIQUE stationary law; the chain is irreducible, aperiodic (P(v->v) >= p > 0) and regenerates with")
    print("   probability 1-p per step, hence ergodic; E[v] = 2 < infinity: the ergodic theorem applies and the a.s. slope is log2 3 - 2 = %.4f" % (L3 - 2))
    rng = random.Random(8675309)

    def geom() -> int:
        v = 1
        while rng.random() < 0.5:
            v += 1
        return v

    p = 0.9
    L = 2 * 10 ** 6
    v = geom()
    cnt = Counter()
    ones = ones_keep = 0
    h = 0.0
    for _ in range(L):
        nv = v if rng.random() < p else geom()
        cnt[nv] += 1
        if v == 1:
            ones += 1
            ones_keep += (nv == 1)
        v = nv
        h += L3 - v
    print(f"   simulated chain p = {p}, {L} steps (own seed): empirical marginal P(v=k) vs 2^-k: "
          + ", ".join(f"k={k}: {cnt[k] / L:.4f}/{2.0 ** -k:.4f}" for k in range(1, 9)))
    print(f"      P(v_(l+1) = 1 | v_l = 1) = {ones_keep / ones:.4f} (exact (1+p)/2 = {(1 + p) / 2}); mean v = {sum(k * c for k, c in cnt.items()) / L:.4f}; slope {h / L:+.4f}")
    print("   200 paths x 4000 steps per persistence (own seed):")
    for p in (0.0, 0.5, 0.9, 0.99):
        slopes, maxh = [], []
        for _ in range(200):
            v = geom()
            h = 0.0
            mx = 0.0
            for _ in range(4000):
                if rng.random() >= p:
                    v = geom()
                h += L3 - v
                if h > mx:
                    mx = h
            slopes.append(h / 4000)
            maxh.append(mx)
        mean = sum(slopes) / 200
        sd = math.sqrt(sum((s_ - mean) ** 2 for s_ in slopes) / 199) / math.sqrt(200)
        # asymptotic variance of the sum: n Var(v) (1+p)/(1-p) with Var(v) = 2 (correlation p^k)
        pred_sd = math.sqrt(4000 * 2 * (1 + p) / (1 - p)) / 4000 / math.sqrt(200) if p < 1 else float('nan')
        print(f"      p={p:.2f}: mean slope {mean:+.4f} (s.e. {sd:.4f}; theory s.e. {pred_sd:.4f}); mean max height {sum(maxh) / 200:.1f}; "
              f"paths above start at the end: {sum(1 for s_ in slopes if s_ > 0)}/200")
    print("   mean maximum height over 4000 steps with 2000 paths (the note's 1.3, 3.1, 17.7, 169 are 200-path means):")
    for p in (0.0, 0.5, 0.9, 0.99):
        maxh = []
        for _ in range(2000):
            v = geom()
            h = 0.0
            mx = 0.0
            for _ in range(4000):
                if rng.random() >= p:
                    v = geom()
                h += L3 - v
                if h > mx:
                    mx = h
            maxh.append(mx)
        mean = sum(maxh) / 2000
        se = math.sqrt(sum((x - mean) ** 2 for x in maxh) / 1999 / 2000)
        print(f"      p={p:.2f}: mean max height {mean:.2f} (s.e. {se:.2f})")
    print(stamp())


# ----------------------------------------------------------------------------
# G. Proposition B: the product criterion, the two rates, the model simulation
# ----------------------------------------------------------------------------

def part_G() -> None:
    print("\n== G. Proposition B: escape products ==")
    logs = 20.0 + 0.3 * np.arange(20000)
    p_sqrt = 3.0 * 2.0 ** (-logs / 2)
    p_harm = 0.6 / (logs * math.log(2))
    never_sqrt = float(np.exp(np.sum(np.log1p(-p_sqrt))))
    logs_inf = 20.0 + 0.3 * np.arange(2000)
    never_inf = float(np.exp(np.sum(np.log1p(-3.0 * 2.0 ** (-logs_inf / 2)))))
    print(f"   escape 3/sqrt(N) from N0 = 2^20, r = 2^0.3: sum p_k over 20000 steps = {p_sqrt.sum():.5f}, exact P(never leave in 20000 steps) = {never_sqrt:.4f}; "
          f"infinite product = {never_inf:.4f}; 1 - sum = {1 - p_sqrt.sum():.4f}")
    print(f"   escape 0.6/ln N: sum over 20000 steps = {p_harm.sum():.2f} (harmonic, diverges like (0.6/(0.3 ln 2)) ln k), exact P(never leave in 20000 steps) = {np.exp(np.sum(np.log1p(-p_harm))):.2e}")
    rng = np.random.default_rng(8675309)
    surv = 0
    for _ in range(10):
        u = rng.random((200, 20000))
        surv += int(np.all(u >= p_sqrt, axis=1).sum())
    print(f"   simulation (2000 trials, own seed): never leaving in 20000 steps with escape 3/sqrt(N): {surv / 2000:.3f}")
    print("   caveat on the statement: 'positive iff sum < infinity' needs p(N_0 r^k) < 1 for all k (else the product is 0 with a finite sum); "
          "'over renewed visits, infinitely often' presupposes infinitely many visits (each visit ends a.s., so if the state is visited infinitely often it is left infinitely often)")
    print(stamp())


# ----------------------------------------------------------------------------
# H. Tenth note, Propositions 1 and 3 (secondary target)
# ----------------------------------------------------------------------------

def part_H() -> None:
    print("\n== H. Tenth note, Proposition 1 (Apery form) and Proposition 3 (free cofactor), Proposition 3(d) (LTE) ==")
    # Proposition 1: A_j = 3^j n + S_(j-1) with S_(j-1) = 3 S_(j-2) + 2^(D_(j-1)), D = cumulative valuation; v2(A_j) = D_j
    ok_cum = True
    fail_step = 0
    n_steps = 0
    for n in range(1, 4000, 2):
        m = n
        S = 0
        D = 0
        j = 0
        while m != 1:
            j += 1
            S = 3 * S + 2 ** D
            m, v = U(m)
            D += v
            A = 3 ** j * n + S
            n_steps += 1
            ok_cum &= (v2(A) == D) and (A == 2 ** D * m) and (v2(3 ** j * n) == 0) and (v2(S) == 0)
            if v2(A) != v:
                fail_step += 1
    print(f"   odd n < 4000, all {n_steps} orbit steps: v2(3^j n + S_(j-1)) = D_j = v_0 + ... + v_(j-1) (cumulative) and 3^j n + S_(j-1) = 2^(D_j) m_j: {ok_cum}")
    print(f"   with the PER-STEP valuation v_(j-1) in place of D_j the identity fails on {fail_step} of the {n_steps} steps (e.g. n = 7, j = 2: A_2 = 63 + 3 + 2 = 68 = 2^2 * 17, v2 = 2 = 1 + 1)")
    print("   so Proposition 1 is correct with d_j = the cumulative 2-adic exponent of the Apery form (the first note's convention, and the script's 'd += v'); the tenth note does not say so")
    # Proposition 3 (a),(b): K <= 10, t odd < 50
    ok = 0
    for K in range(1, 11):
        for t in range(1, 50, 2):
            m = 2 ** (K + 1) * t - 1
            x = m
            for _ in range(K):
                x, v = U(x)
                assert v == 1
            assert x == 2 * 3 ** K * t - 1
            _, vnext = U(x)
            assert vnext >= 2 and vnext == 1 + v2(3 ** (K + 1) * t - 1)
            ok += 1
    print(f"   Proposition 3(a),(b): K <= 10, t odd < 50: run of EXACTLY K ones, U^K(m) = 2*3^K t - 1, next valuation 1 + v2(3^(K+1) t - 1): {ok} cases OK")
    longer = all(U(2 * 3 ** K * t - 1)[1] == 1 for K in range(1, 8) for t in range(2, 40, 2))
    print(f"   for t EVEN the run is longer than K (2*3^K t - 1 = 3 mod 4): {longer} (K < 8, even t < 40) -- the hypothesis 't odd' in (a) is needed")
    # Proposition 3(c): least m with run exactly K then valuation >= J, and exactly J
    print("   Proposition 3(c),(d): least m = 2^(K+1) t - 1 (t odd) with next valuation >= J, brute force vs 2^(K+1) rho_J - 1; and the least with valuation EXACTLY J:")
    for K in range(1, 8):
        row = []
        for J in range(2, 8):
            rho = pow(3, -(K + 1), 2 ** (J - 1))
            pred = 2 ** (K + 1) * rho - 1
            best_ge = best_eq = None
            t = 1
            while best_ge is None or best_eq is None:
                x = 2 * 3 ** K * t - 1
                vn = U(x)[1]
                if vn >= J and best_ge is None:
                    best_ge = 2 ** (K + 1) * t - 1
                if vn == J and best_eq is None:
                    best_eq = 2 ** (K + 1) * t - 1
                t += 2
            row.append(f"J={J}: >=J {best_ge}{'=' if best_ge == pred else '!='}pred, =J {best_eq}{'(same)' if best_eq == best_ge else '(differs)'}")
        print(f"      K={K}: " + "; ".join(row))
    # (d) LTE law and the Mersenne continuation
    lte = all(v2(3 ** (K + 1) - 1) == (1 if (K + 1) % 2 else 2 + v2(K + 1)) for K in range(0, 300))
    print(f"   LTE: v2(3^(K+1) - 1) = 1 (K+1 odd), 2 + v2(K+1) (K+1 even), for 0 <= K < 300: {lte}")
    mers = True
    for K in range(1, 120):
        x = 2 ** (K + 1) - 1
        for _ in range(K):
            x, v = U(x)
            assert v == 1
        vn = U(x)[1]
        mers &= vn == (2 if K % 2 == 0 else 3 + v2(K + 1))
    print(f"   Mersenne 2^(K+1) - 1 continues after its K ones with valuation exactly 2 (K even) / 3 + v2(K+1) (K odd), K < 120: {mers}")
    print("   note on (c): 'the least m realizing depth J' is correct for 'next valuation >= J'; for 'exactly J' it can fail (rows above where '=J' differs, e.g. K=1, J=2: m=3 continues with valuation 4)")
    print(stamp())


# ----------------------------------------------------------------------------
# I. The typology numbers (section 5)
# ----------------------------------------------------------------------------

def part_I(sig: np.ndarray) -> None:
    print("\n== I. Section 5 typology numbers, recomputed ==")
    print(f"   Collatz drift log2 3 - 2 = {L3 - 2:+.4f}; 5n+1: log2 5 - 2 = {math.log2(5) - 2:+.4f} per odd step, (log2(5/2) - 1)/2 = {(math.log2(2.5) - 1) / 2:+.4f} per T-step")
    # Juggler
    tot = keep = 0
    ll = []
    maxsteps = maxdig = 0
    argsteps = argdig = None
    for n in range(2, 2001):
        m = n
        steps = 0
        prev_par = None
        while m != 1:
            m2 = math.isqrt(m ** 3) if m % 2 else math.isqrt(m)
            if m2 > 1 and m > 1:
                ll.append(math.log2(math.log(m2) / math.log(m)))
            par = m % 2
            if prev_par is not None:
                tot += 1
                keep += (par == prev_par)
            prev_par = par
            if len(str(m2)) > maxdig:
                maxdig = len(str(m2))
                argdig = n
            m = m2
            steps += 1
            assert steps < 10 ** 4
        if steps > maxsteps:
            maxsteps = steps
            argsteps = n
    print(f"   Juggler n <= 2000: all reach 1; mean log2(log m'/log m) {sum(ll) / len(ll):+.3f} (fair coin {0.5 * math.log2(1.5) - 0.5:+.3f}); "
          f"P(parity persists) {keep / tot:.3f}; longest {maxsteps} steps (n={argsteps}); largest {maxdig} digits (n={argdig})")
    # reverse-and-add, the note's convention (300 steps, growth over the first 50 steps) and the parallel's (200 steps, all steps)
    cand = []
    g50 = []
    gall = 0.0
    nall = 0
    for n in range(1, 10 ** 5):
        m = n
        hit = False
        for k in range(300):
            s = str(m)
            if k > 0 and s == s[::-1]:
                hit = True
                break
            m2 = m + int(s[::-1])
            if k < 50:
                g50.append(math.log10(m2 / m))
            if k < 200:
                gall += math.log10(m2) - math.log10(m)
                nall += 1
            m = m2
        if not hit:
            cand.append(n)
    print(f"   reverse-and-add n < 10^5, 300 steps: not palindromic {len(cand)} (first {cand[:12]}); below 10^4: {sum(1 for c in cand if c < 10 ** 4)}; "
          f"mean digit growth per step: first 50 steps {sum(g50) / len(g50):+.3f}, all steps to 200 {gall / nall:+.3f}")
    # look-and-say via run-length encoding in numpy
    s = np.frombuffer(b"1", dtype=np.uint8).copy()
    lens = []
    for _ in range(60):
        change = np.flatnonzero(np.diff(s.astype(np.int16))) + 1
        starts = np.concatenate(([0], change))
        ends = np.concatenate((change, [len(s)]))
        runs = ends - starts
        out = np.empty(2 * len(starts), dtype=np.uint8)
        out[0::2] = runs.astype(np.uint8) + 48
        out[1::2] = s[starts]
        s = out
        lens.append(len(s))
    print(f"   look-and-say: lengths ..., {lens[-2]}, {lens[-1]}; ratio at step 60: {lens[-1] / lens[-2]:.5f} (Conway 1.303577)")
    # Ducci length 8, entries < 4
    vecs = np.array(list(product(range(4), repeat=8)), dtype=np.int64)
    steps = np.full(len(vecs), -1)
    cur = vecs.copy()
    k = 0
    while (steps < 0).any():
        zero = ~cur.any(axis=1)
        newly = zero & (steps < 0)
        steps[newly] = k
        cur = np.abs(cur - np.roll(cur, -1, axis=1))
        k += 1
        assert k < 100
    print(f"   Ducci length 8, entries < 4 ({len(vecs)} vectors): all reach 0, worst {steps.max()} steps")
    # Kaprekar 4 digits
    worst = 0
    for n in range(10000):
        d = f"{n:04d}"
        if len(set(d)) == 1:
            continue
        x = n
        k = 0
        while x != 6174:
            d = f"{x:04d}"
            x = int("".join(sorted(d, reverse=True))) - int("".join(sorted(d)))
            k += 1
            assert k < 20
        worst = max(worst, k)
    print(f"   Kaprekar 4 digits: every non-repdigit reaches 6174, worst {worst} steps")
    print(stamp())


# ----------------------------------------------------------------------------
# J. Erdos persistence of increases (fixed horizon), all abundant n <= 10^6 and a sample
# ----------------------------------------------------------------------------

def part_J(sig: np.ndarray) -> None:
    print("\n== J. Erdos persistence P(k consecutive increases | s(n) > n), sieve to 2*10^7 (a chain leaving the sieve counts as a failure) ==")
    K = 8
    n = np.arange(2, 10 ** 6 + 1, dtype=np.int64)
    abund_all = n[sig[n] > 2 * n].tolist()
    rng = random.Random(8675309)
    sample = rng.sample(abund_all, 5000)
    for label, starts in (("all abundant n <= 10^6", abund_all), ("5000 sampled abundant n <= 10^6", sample)):
        alive = [0] * (K + 1)
        inside = [0] * (K + 1)
        for n0 in starts:
            x = n0
            for k in range(1, K + 1):
                if x > M_SIEVE:
                    break
                s = int(sig[x]) - x
                inside[k] += 1
                if s > x:
                    alive[k] += 1
                    x = s
                else:
                    break
        print(f"   {label} ({len(starts)} starts): " + ", ".join(f"k={k}: {alive[k] / len(starts):.3f} (inside {inside[k]})" for k in range(1, K + 1)))
    print(stamp())


if __name__ == "__main__":
    print(f"audit script start; sieves: sigma to {M_SIEVE}, spf to {N_SPF}")
    SIG = sigma_sieve(M_SIEVE)
    print(f"{stamp()} sigma sieve done")
    SPF = spf_sieve(N_SPF)
    print(f"{stamp()} spf sieve done\n")
    part_A(SIG, SPF)
    part_B(SIG)
    part_C(SIG, SPF)
    part_D(SIG)
    part_E()
    part_F()
    part_G()
    part_H()
    part_I(SIG)
    part_J(SIG)
    print(f"\n{stamp()} DONE")
