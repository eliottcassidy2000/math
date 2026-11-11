#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9665 (A), (C). Is rho(T_k) < 2 for the non-trivial characters k of Z/M, M = 2^(j-1) m ?
# Digit recursion: h -> 2h (q unchanged), h -> 2h+1 (q += ceil(h/2)).  ceil(h/2) mod M is a function of w = h mod 2M.
# T_k[w -> 2w] = 1, T_k[w -> 2w+1] = e(k ceil(w/2)/M), states w mod 2M.
# Also an exact finite check of the statement used in the argument:
#   for M even and k odd, the states w and w+M have the same successors and their phases differ by -1.
# And a direct check of the character sums: max_l |sum_{h<2^n} e(k q(h)/M + l h/(2M))| / 2^n for n = 20.
import numpy as np, math, sys
from f1core import build_psi


def T(k, M):
    n = 2 * M
    A = np.zeros((n, n), dtype=complex)
    for w in range(n):
        A[w, (2 * w) % n] += 1.0
        A[w, (2 * w + 1) % n] += np.exp(2j * np.pi * k * ((w + 1) // 2) / M)
    return A


def main():
    print("# mixed_gap.py (auditor hypotheses; own code)")
    print("# rate(M) = max over k != 0 mod M of rho(T_k)/2;  listed with the worst k, and split by the 2-part of the order of the character")
    for M in (3, 6, 12, 24, 5, 10, 20, 40, 9, 18, 36, 7, 14, 28, 15, 30, 60, 45, 90):
        rates = {}
        worst = (0.0, None)
        for k in range(1, M):
            ev = np.linalg.eigvals(T(k, M))
            r = max(abs(ev)) / 2
            if r > worst[0]:
                worst = (r, k)
            # order of the character k/M
            order = M // math.gcd(k, M)
            two = 0
            o = order
            while o % 2 == 0:
                o //= 2
                two += 1
            key = (two, o)
            rates[key] = max(rates.get(key, 0.0), r)
        desc = "  ".join("ord 2^%d*%d: %.4f" % (a, b, v) for (a, b), v in sorted(rates.items()))
        print("M = %3d  rate %.6f at k = %s   | %s" % (M, worst[0], worst[1], desc))
        # the structural fact used in the gap argument
        if M % 2 == 0:
            ok = True
            for k in range(1, M, 2):
                for w in range(2 * M):
                    p1 = np.exp(2j * np.pi * k * ((w + 1) // 2) / M)
                    p2 = np.exp(2j * np.pi * k * (((w + M) % (2 * M) + 1) // 2) / M)
                    if abs(p1 + p2) > 1e-9:
                        ok = False
            print("         k odd: phase(w + M) = -phase(w) for every w: %s" % ok)
    # direct character sums over h < 2^n
    n = 20
    psi = build_psi(1 << n)
    h = np.arange(1 << n, dtype=np.int64)
    q = h - psi
    print()
    print("# direct sums  max_{k != 0 mod M, l} |sum_{h < 2^n} e(k q/M + l h/M)| / 2^n   at n = 16, 20  (should decay like rate^n)")
    for M in (6, 12, 10, 20, 18, 36, 14, 28):
        out = []
        for nn in (16, 20):
            best = 0.0
            hh = h[: 1 << nn]
            qq = q[: 1 << nn]
            for k in range(1, M):
                for l in range(M):
                    s = np.exp(2j * np.pi * ((k * qq + l * hh) % M) / M).sum()
                    best = max(best, abs(s) / (1 << nn))
            out.append(best)
        print("M = %3d   n=16: %.5f   n=20: %.5f   ratio^(1/4) = %.4f" % (M, out[0], out[1], (out[1] / out[0]) ** 0.25))


if __name__ == "__main__":
    main()
