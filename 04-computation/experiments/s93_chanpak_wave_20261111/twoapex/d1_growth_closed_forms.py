#!/usr/bin/env python3
"""Lane twoapex, deliverables (d) and (e): growth, the limit law of H / (2 * 3^n), and closed forms on sub-families.

usage: python d1_growth_closed_forms.py STATE_DIR
(1) exact first moment: sum over |w| = n of H(w) = 4^n (3^n + 3 * 2^n + 3) / 8   (n <= 13, from all words);
(2) second moment and the constant c2 = lim E[H^2] / 9^n;
(3) H(w) <= 2 * 3^n and the lower bound 2 * 3^n - 4 (3^p 2^q + 2^p 3^q) - 2^n <= H(3^p 0^q), p + q = n <= 24;
(4) limit law: H(w) / (2 * 3^n) for a uniform word of length n against S * E',
        S = sum_{i != j} a_i b_j 2^(|i-j|-1) 3^(-max(i,j)),   E' an independent copy
    (quantiles from the exact histograms of b1_census.py, n = 12..17, and from 2 * 10^6 samples of S * E');
(5) closed forms, checked on ALL words of the family for n <= 11:
      E(a, b)       number of paths [R1] u v [R3]
                    = [b_1 = 1] + [a_n = 0] + sum_k [a_k = 0][b_(k+1) = 1] + sum_{F < L} b_F (1 - a_L) 2^(L-F-1)
      one apex      x(a) = E(a, a) = H(a, 0^n) = H(1^n, a)
      H(0^n, b)     = x(b) + h(b) + b_n,        h(b) = sum_F b_F 2^(n-F)
      H(a, 1^n)     = x(a) + C(a) + 1 - a_1,    C(a) = sum_L (1 - a_L) 2^(L-1)
      H(2^n)        = 2^n + 1
      F1            H(0^p a'', b' 1^q) = x(b') x(a'') + (h(b') + 1)(C(a'') + 1)          (p, q >= 0)
      leading term  | H(w) - 2 * 3^n S_n(w) E_n(w) |  (reported: its largest value relative to n 2^n)
(6) rigidity: for every b and every i < j the second difference d_i d_j H(., b) is <= -1 (n <= 10).
"""
import itertools
import sys

import numpy as np

from ta_fast import FIN24, H24, INIT24, M24, all_alpha


def allH(n):
    """H of all words of length n; index = sum letter_i 4^(n-1-i) (first letter top digit)."""
    return all_alpha(n, M24, INIT24) @ FIN24


def bits(word):
    return [l & 1 for l in word], [(l >> 1) & 1 for l in word]


def E_formula(a, b):
    n = len(a)
    if n == 0:
        return 1
    t = (1 if b[0] == 1 else 0) + (1 if a[-1] == 0 else 0)
    for k in range(n - 1):
        if a[k] == 0 and b[k + 1] == 1:
            t += 1
    for F in range(n):
        if b[F]:
            for L in range(F + 1, n):
                if a[L] == 0:
                    t += 1 << (L - F - 1)
    return t


def main():
    sdir = sys.argv[1] if len(sys.argv) > 1 else None
    ok = True
    print("(1), (2) moments over all words of length n")
    for n in range(1, 14):
        h = allH(n).astype(object) if n <= 10 else allH(n)
        s1 = int(np.sum(allH(n), dtype=np.int64))
        s2 = sum(int(v) * int(v) for v in allH(n)) if n <= 9 else float(np.sum(allH(n).astype(np.float64) ** 2))
        pred = 4 ** n * (3 ** n + 3 * 2 ** n + 3) // 8
        ok &= s1 == pred
        print(f"  n = {n:2d}: sum H = {s1} = 4^n (3^n + 3 2^n + 3)/8: {s1 == pred};  E[H^2] / 9^n = "
              f"{float(s2) / 4 ** n / 9 ** n:.6f};  E[H] / 3^n = {s1 / 4 ** n / 3 ** n:.6f}")
    print("(3) the words 3^p 0^q")
    for n in (8, 12, 16, 20, 24):
        p = n // 2
        q = n - p
        hv = H24([3] * p + [0] * q)
        lb = 2 * 3 ** n - 4 * (3 ** p * 2 ** q + 2 ** p * 3 ** q) - 2 ** n
        ok &= lb <= hv <= 2 * 3 ** n
        print(f"  n = {n}: H(3^{p} 0^{q}) = {hv} = {hv / (2 * 3 ** n):.6f} * 2 * 3^n; lower bound {lb / (2 * 3 ** n):.6f}; "
              f"bound holds: {lb <= hv <= 2 * 3 ** n}")
    # exact ratio sequence for large n: 1 - H / (2 3^n) against (2/3)^(n/2)
    for n in (20, 30, 40, 60, 80):
        p = n // 2
        hv = H24([3] * p + [0] * (n - p))
        print(f"  n = {n}: (1 - H(3^{p} 0^{n - p}) / (2 3^n)) / (2/3)^(n/2) = "
              f"{(1 - hv / (2 * 3 ** n)) / (2 / 3) ** (n / 2):.6f}")
    print("(4) limit law of H / (2 3^n)")
    rng = np.random.default_rng(20261111)
    NS, K = 2_000_000, 48

    def sample_S():
        a = rng.integers(0, 2, size=(NS, K), dtype=np.int8)
        b = rng.integers(0, 2, size=(NS, K), dtype=np.int8)
        S = np.zeros(NS)
        A = np.zeros(NS)          # A_(k-1) / 2^(k-1) scaled: keep A_k = sum_{i<=k} a_i 2^(k-i)
        B = np.zeros(NS)
        w = 1.0
        for k in range(K):
            w /= 3.0
            S += w * (a[:, k] * B + b[:, k] * A)
            A = 2 * A + a[:, k]
            B = 2 * B + b[:, k]
        return S
    S1 = sample_S()
    S2 = sample_S()
    prod = np.sort(S1 * S2)
    qs = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]
    print(f"  S: mean {S1.mean():.5f} (exact 1/4), second moment {np.mean(S1 ** 2):.5f}, max {S1.max():.4f}")
    print(f"  S * E': mean {prod.mean():.5f} (exact 1/16), second moment {np.mean(prod ** 2):.6f}")
    print("  quantiles     " + " ".join(f"{q:8.2f}" for q in qs))
    print("  S * E' (MC)   " + " ".join(f"{prod[int(q * NS)]:8.5f}" for q in qs))
    if sdir:
        for n in (12, 13, 14, 15, 16, 17):
            try:
                bn = np.load(f"{sdir}/valbins_n{n:02d}.npy").astype(np.float64)
            except FileNotFoundError:
                continue
            cdf = np.cumsum(bn) / bn.sum()
            out = []
            for q in qs:
                k = int(np.searchsorted(cdf, q))
                out.append((k + 0.5) / 4096)
            print(f"  H/(2 3^n), n={n} " + " ".join(f"{v:8.5f}" for v in out))
        # small-value exponent: P(H / (2 3^n) <= x) at n = 17
        bn = np.load(f"{sdir}/valbins_n17.npy").astype(np.float64)
        cdf = np.cumsum(bn) / bn.sum()
        print("  n = 17: P(H <= x * 2 * 3^n) for x = 2^-j / 1: " +
              " ".join(f"j={j}:{cdf[(4096 >> j) - 1]:.5f}" for j in range(1, 11)))
        print("     local exponents log2 ratio: " + " ".join(
            f"{np.log2(cdf[(4096 >> j) - 1] / cdf[(4096 >> (j + 1)) - 1]):.3f}" for j in range(1, 10)))
    print("(5) closed forms")
    bad = {k: 0 for k in ("E=x one apex", "H(a,0)=x(a)", "H(1,b)=x(b)", "H(0^n,b)", "H(a,1^n)", "2^n", "F1")}
    cnt = 0
    for n in range(1, 12):
        for a in itertools.product((0, 1), repeat=n):
            x = E_formula(a, a)
            cnt += 1
            if H24([ai for ai in a]) != x:                 # letters (a_i, 0)
                bad["H(a,0)=x(a)"] += 1
            if H24([1 + 2 * ai for ai in a]) != x:         # letters (1, b_i) with b = a
                bad["H(1,b)=x(b)"] += 1
            hb = sum(a[F] << (n - 1 - F) for F in range(n))
            if H24([2 * ai for ai in a]) != x + hb + a[-1]:
                bad["H(0^n,b)"] += 1
            Ca = sum((1 - a[L]) << L for L in range(n))
            if H24([ai + 2 for ai in a]) != x + Ca + 1 - a[0]:
                bad["H(a,1^n)"] += 1
        if H24([2] * n) != 2 ** n + 1:
            bad["2^n"] += 1
    for p in range(0, 7):
        for q in range(0, 7):
            for bp in itertools.product((0, 1), repeat=p):
                for aq in itertools.product((0, 1), repeat=q):
                    word = [2 * t for t in bp] + [t + 2 for t in aq]
                    x1 = E_formula(bp, bp)
                    x2 = E_formula(aq, aq)
                    y1 = sum(bp[F] << (p - 1 - F) for F in range(p)) + 1
                    y2 = sum((1 - aq[L]) << L for L in range(q)) + 1
                    cnt += 1
                    if H24(word) != x1 * x2 + y1 * y2:
                        bad["F1"] += 1
    print(f"  {cnt} instances; violations: {bad}")
    ok &= all(v == 0 for v in bad.values())
    # E(a, b) against the literal count of [R1] u v [R3]
    badE = 0
    for n in range(1, 8):
        for a in itertools.product((0, 1), repeat=n):
            for b in itertools.product((0, 1), repeat=n):
                c = 0
                for cut in itertools.product((0, 1), repeat=n):          # 0: R1, 1: R3
                    R1 = [i for i in range(n) if cut[i] == 0]
                    R3 = [i for i in range(n) if cut[i] == 1]
                    if R1 and a[R1[-1]] != 0:
                        continue
                    if R3 and b[R3[0]] != 1:
                        continue
                    c += 1
                if c != E_formula(a, b):
                    badE += 1
    print(f"  E(a, b) formula against the literal count, all a, b of length <= 7: violations {badE}")
    ok &= badE == 0
    # leading term
    print("  leading term: max over all words of length n of |H - 2 3^n S_n E_n| / (n 2^n)")
    for n in range(2, 11):
        Hn = allH(n)
        worst = 0.0
        idx = np.arange(4 ** n)
        a = np.array([[(i >> (2 * (n - 1 - k))) & 1 for k in range(n)] for i in idx])
        b = np.array([[(i >> (2 * (n - 1 - k) + 1)) & 1 for k in range(n)] for i in idx])
        S = np.zeros(4 ** n)
        Ee = np.zeros(4 ** n)
        for i in range(n):
            for j in range(n):
                if i != j:
                    S += a[:, i] * b[:, j] * 2.0 ** (abs(i - j) - 1) * 3.0 ** (-(max(i, j) + 1))
                    Ee += (1 - a[:, i]) * (1 - b[:, j]) * 2.0 ** (abs(i - j) - 1) * 3.0 ** (-(n - min(i, j)))
        dev = np.abs(Hn - 2 * 3.0 ** n * S * Ee)
        print(f"    n = {n:2d}: {dev.max() / (n * 2 ** n):.4f}   (mean |dev| / 2^n = {dev.mean() / 2 ** n:.4f})")
    print("(6) rigidity of the digits of one apex")
    for n in range(2, 11):
        Hn = allH(n)

        def code(a, b):
            c = 0
            for k in range(n):
                c = (c << 2) | (a[k] + 2 * b[k])
            return c
        mx = -10 ** 18
        for b in itertools.product((0, 1), repeat=n):
            zero = [0] * n
            for i in range(n):
                for j in range(i + 1, n):
                    ai = list(zero)
                    ai[i] = 1
                    aj = list(zero)
                    aj[j] = 1
                    aij = list(zero)
                    aij[i] = aij[j] = 1
                    c = int(Hn[code(aij, b)] - Hn[code(ai, b)] - Hn[code(aj, b)] + Hn[code(zero, b)])
                    mx = max(mx, c)
        ok &= mx <= -1
        print(f"  n = {n:2d}: largest second difference d_i d_j H(., b) over all b and i < j: {mx}")
    print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")


if __name__ == "__main__":
    main()
