#!/usr/bin/env python3
"""
lrc14_PC_paircorr_thinning_klein_S215.py

HYP-5790: the (PC) pair-correlation count (THM-677's one ingredient) + the
THINNING ESCAPE.

(PC): #{(j,j') in Good^2, j != j' : ||m(tau_j - tau_j')|| <= delta}
      <= (1+eps) * N^2 * 2delta + c*N,     delta = 1/(2H+1), H = 14.

Studies:
 [1] EXACT (PC) per (killer, bank): the implied (eps, c); sub/super-Poisson?
 [2] the per-difference profile: excess(d) = #near-pairs at difference d
     minus Poisson share -- does the excess concentrate on structured d
     (small d, d ~ multiples of V/m or gcd structure)?
 [3] THE THINNING ESCAPE: greedily select Good' avoiding the top-excess
     differences for ALL killers simultaneously; report density N'/N and the
     resulting D_m on Good' (survival needs D_m <= (1/5-1/7)N' -- margins).
 [4] BERNSTEIN CHECK: is D_m <= sqrt(2H+1) * l2 in practice (the sup step)?
"""

from fractions import Fraction as F
from math import gcd, pi, sin, sqrt
import cmath, itertools as it

ONE14 = F(1, 14)
H = 14
DELTA = 1.0 / (2 * H + 1)


def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))


def is_primitive(S):
    g = 0
    for v in S:
        g = gcd(g, v)
    return g == 1


def good_periods(E, V):
    s = max(E)
    out = {}
    for j in range(V):
        teeth = sorted(set((e * j) % V for e in set(E)))
        n = len(teeth)
        best, ustart = -1, 0
        for i in range(n):
            a = teeth[i]
            b = teeth[(i + 1) % n][0] if False else teeth[(i + 1) % n] + (V if i == n - 1 else 0)
            if b - a > best:
                best, ustart = b - a, a
        if 7 * best <= V:
            continue
        a, g, r = F(ustart, V), F(best, V), F(s, V)
        lo = (a + ONE14) / (1 - r)
        hi = (a + g - ONE14) / (1 + r)
        if lo >= hi:
            continue
        phi = (lo + hi) / 2
        tau = (j + phi) / V
        ok = True
        for e in set(E):
            v = V - e
            x = (v * tau) % 1
            if min(x, 1 - x) < ONE14:
                ok = False
                break
        if ok:
            out[j] = tau
    return out


def near(x, delta=DELTA):
    y = x % 1
    return min(y, 1 - y) <= delta


def study(name, killers, E, V):
    good = good_periods(E, V)
    js = sorted(good)
    N = len(js)
    taus = {j: float(good[j]) for j in js}
    print(f"\n=== {name}: V={V} N={N} ===")
    poisson = N * (N - 1) * 2 * DELTA
    bad_d_all = {}
    for m in killers:
        cnt = 0
        by_d = {}
        for i, j in enumerate(js):
            for j2 in js[i + 1:]:
                if near(m * (taus[j2] - taus[j])):
                    cnt += 2  # ordered pairs
                    d = j2 - j
                    by_d[d] = by_d.get(d, 0) + 1
        eps = cnt / poisson - 1
        # top-excess differences
        rd = {}
        for i, j in enumerate(js):
            for j2 in js[i + 1:]:
                rd[j2 - j] = rd.get(j2 - j, 0) + 1
        excess = {d: by_d.get(d, 0) - rd[d] * 2 * DELTA for d in rd}
        top = sorted(excess.items(), key=lambda kv: -kv[1])[:6]
        print(f"  m={m:5d}: near-pairs={cnt:5d} Poisson={poisson:7.1f} "
              f"ratio={(cnt/poisson):5.2f} (eps={eps:+.2f})")
        print(f"     top-excess d: " + ", ".join(f"d={d}(+{e:.1f})" for d, e in top if e > 0.5))
        for d, e in excess.items():
            if e > 1.0:
                bad_d_all[d] = bad_d_all.get(d, 0) + e
    # [3] thinning: drop j's greedily to avoid the worst differences
    bad_ds = sorted(bad_d_all, key=lambda d: -bad_d_all[d])[:max(3, len(bad_d_all)//3)]
    keep = []
    for j in js:
        if all((j - j2) not in bad_ds and (j2 - j) not in bad_ds for j2 in keep):
            keep.append(j)
    Np = len(keep)
    print(f"  [3] thinning: |bad d set|={len(bad_ds)} -> Good' density {Np}/{N} = {Np/N:.2f}")
    needp = (1 / 5 - 1 / 7) * Np
    ok_all = True
    for m in killers:
        kill = sum(1 for j in keep if near(m * taus[j], 1 / 14))
        D = abs(kill - Np / 7)
        ok = D <= needp
        ok_all &= ok
        print(f"     m={m:5d}: D_m(Good')={D:5.1f} vs needed {needp:.1f}  {'OK' if ok else 'over'}")
    # [4] Bernstein check on full Good
    print(f"  [4] Bernstein sup-vs-L2 (full Good):")
    for m in killers:
        Z = [sum(cmath.exp(2j * pi * h * m * taus[j]) for j in js) for h in range(1, 3 * H + 1)]
        l2 = sqrt(sum((sin(pi * h / 7) / (pi * h)) ** 2 * 2 * abs(Z[h - 1]) ** 2
                      for h in range(1, 3 * H + 1) if h % 7))
        kill = sum(1 for j in js if near(m * taus[j], 1 / 14))
        D = abs(kill - N / 7)
        print(f"     m={m:5d}: D={D:5.1f}  sqrt(2H+1)*l2={sqrt(2*H+1)*l2:6.1f}  "
              f"ratio D/(l2)={D/l2 if l2 else 0:4.2f} (Bernstein factor {sqrt(2*H+1):.1f})")


def main():
    print("=" * 78)
    print("HYP-5790: (PC) exact + per-d profile + thinning escape (klein-S215)")
    print("=" * 78)
    V = 842
    E = [0, 11, 37, 68, 105, 133, 160, 191, 224, 260]
    K = None
    for d in it.product(range(0, 8), repeat=3):
        Kc = [V // 2 - 1 + d[0], V // 3 + 1 + d[1], V // 5 + 2 + d[2]]
        if len(set(Kc)) < 3:
            continue
        S = sorted(Kc + [V - e for e in set(E)])
        if len(set(S)) == 13 and is_covering(S) and is_primitive(S):
            K = Kc
            break
    study("bank2-style 3M", K, E, V)

    V = 1006
    E8 = [0, 17, 55, 96, 141, 190, 243, 300]
    K = None
    for cand in ([12, 13, 251, 402, 503],):
        for d in it.product(range(0, 8), repeat=5):
            Kc = [c + dd for c, dd in zip(cand, d)]
            if len(set(Kc)) < 5:
                continue
            S = sorted(Kc + [V - e for e in set(E8)])
            if len(set(S)) == 13 and is_covering(S) and is_primitive(S):
                K = Kc
                break
        if K:
            break
    study("bank4-style 2P+3M", K, E8, V)
    print("\nDONE.")


if __name__ == '__main__':
    main()
