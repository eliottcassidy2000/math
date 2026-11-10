#!/usr/bin/env python3
"""S91 (THM-5128 / HYP-9660): p-divisible rooted states of one-apex tournaments.

x(h) = 2 psi(h) + 1 (the Hamiltonian-path count of the one-apex tournament with word h), y(h) = 2h + 1,
q(h) = h - psi(h), so x = y - 2q and g(h) := gcd(x(h), y(h)) = gcd(2h + 1, q(h)).

A rooted state (x, y) is p-divisible when p | g(h).  In the F1 family N = x1 x2 + 2^s y1 y2 - [s = 0] a p-divisible
state forces N = 0 (s >= 1) or N = -1 (s = 0) modulo p whatever the other factor is.

Printed (exact integer arithmetic):
  1. for h = 2^j: x = y = 2^(j+1) + 1 (checked), so every prime factor of 2^m + 1 has a divisible state at h = 2^(m-1);
  2. for each odd prime p <= PMAX: the first h with p | g(h), the number of such h below 2^n for several n, and the
     ratio to the equidistribution value 2^n / p^2;
  3. the proportion of h < 2^n with g(h) = 1, against prod_{p odd} (1 - 1/p^2) = 8/pi^2 (the tail of the product is
     NOT proved to be attained: see the note);
  4. the limit law of rho(h) = x(h)/y(h): empirical distribution on [1/2, 1] against the uniform law, and the identity
     rho(h) = 1 - phi(h)/2 + O(log h / h) with phi the binary radical inverse (maximum deviation printed).

Usage: python3 divisible_states.py [NBITS] [PMAX]      (defaults 24, 300)
"""
import math
import sys

import numpy as np


def q_table(n):
    size = 1 << n
    q = np.zeros(size, dtype=np.int64)
    for k in range(n):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0)
    return q


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i:: i] = bytearray(len(s[i * i:: i]))
    return [i for i in range(3, n + 1) if s[i]]


def main():
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    pmax = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    q = q_table(nb)
    h = np.arange(1 << nb, dtype=np.int64)
    y = 2 * h + 1
    x = y - 2 * q
    ok = True
    # psi recursion cross-check on a prefix (independent scalar code)
    psi = [0] * 4096
    for i in range(1, 4096):
        psi[i] = psi[i >> 1] + ((i >> 1) if i % 2 == 0 else (i >> 2) + 1)
    ok &= all(int(x[i]) == 2 * psi[i] + 1 for i in range(4096))
    # 1. powers of two
    for j in range(nb - 1):
        hh = 1 << j
        ok &= int(x[hh]) == int(y[hh]) == (1 << (j + 1)) + 1
    print("1. h = 2^j gives x = y = 2^(j+1) + 1 for j <", nb - 1, ":", ok)
    g = np.gcd(x, y)
    print("   small table  h: (x, y, g)   ", [(int(i), int(x[i]), int(y[i]), int(g[i])) for i in range(1, 17)])
    eq = np.nonzero(x == y)[0]
    ok &= all(int(v) & (int(v) - 1) == 0 for v in eq)
    print("   x = y only at h = 0 and powers of two below 2^%d: %s" % (nb, all(int(v) & (int(v) - 1) == 0 for v in eq)))
    # 2. first divisible states and counts
    print("2. first divisible state and counts  (ratio = count * p^2 / 2^n)")
    ns = [n for n in (10, 14, 18, 22, nb) if n <= nb]
    print(f"   {'p':>5} {'first h':>9} {'bits':>4} {'m: p | 2^m+1':>13} " + " ".join(f"{'n=' + str(n):>9}" for n in ns))
    rows = []
    for p in primes_upto(pmax):
        div = np.nonzero((y % p == 0) & (x % p == 0))[0]
        first = int(div[0]) if len(div) else -1
        m = next((mm for mm in range(1, 2 * p) if pow(2, mm, p) == p - 1), 0)
        ratios = [float(np.searchsorted(div, 1 << n)) * p * p / (1 << n) for n in ns]
        rows.append((p, first, m, ratios))
        ok &= (m == 0) or (first != -1 and first <= (1 << (m - 1))) or ((1 << (m - 1)) >= (1 << nb))
        if p <= 131 or (m and m <= 20) or first == -1:
            print(f"   {p:5d} {first:9d} {first.bit_length() if first > 0 else 0:4d} {m if m else '-':>13} "
                  + " ".join(f"{r:9.3f}" for r in ratios))
    # 3. coprime proportion
    target = 8 / math.pi ** 2
    for n in ns:
        prop = float((g[: 1 << n] == 1).mean())
        print(f"3. proportion of h < 2^{n} with gcd(x, y) = 1: {prop:.5f}   (8/pi^2 = {target:.5f})")
    # 4. limit law of rho = x / y and the radical inverse
    lo = 1 << (nb - 1)
    rho = x[lo:] / y[lo:]
    hist, _ = np.histogram(rho, bins=10, range=(0.5, 1.0))
    print("4. rho = x/y on the top dyadic block: min %.6f max %.6f; decile histogram / uniform:" % (rho.min(), rho.max()),
          [round(float(v) * 10 / len(rho), 4) for v in hist])
    hh = h[lo:]
    phi = np.zeros(len(hh), dtype=np.float64)
    for i in range(nb):
        phi += ((hh >> i) & 1) * 2.0 ** (-(i + 1))
    dev = np.abs(rho - (1 - phi / 2))
    print("   max |rho - (1 - phi/2)| on that block: %.3e   (log h / h = %.3e)" % (dev.max(), nb * math.log(2) / lo))
    ok &= dev.max() < 4 * nb / lo
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
