#!/usr/bin/env python3
"""Coverage of the odd numbers by the type table of order 10.

A root of type (c, r) and exponent j (response c 2^j, offset r, out-set not empty) gives the path counts
c 2^e + r for every e >= j + 1 (tail law, THM-5126 (A)).  This script reads the table of least exponents
(out/blocks_roots_order10.txt: lines "c r j", odd c <= 255, odd |r| < 2047) and reports which odd N in [3, X]
are of the form c 2^e + r with e >= j(c, r) + 1, i.e. are H(L o 0^t), t >= 1, for a rooted tournament (L, v) of
order at most 10 with a non-empty out-set.

usage: coverage.py [X]      (default X = 2^20)
"""
import sys
def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 1 << 20
    types = []
    for line in open('out/blocks_roots_order10.txt'):
        if line.startswith('#'): continue
        c, r, j = map(int, line.split()); types.append((c, r, j))
    cov = bytearray(X + 1)
    for c, r, j in types:
        e = j + 1; v = (c << e) + r
        while v <= X:
            if v >= 1: cov[v] = 1
            e += 1; v = (c << e) + r
    missing = [N for N in range(3, X + 1, 2) if not cov[N]]
    print("types in the table:", len(types))
    print("odd N in [3, %d]: %d ; of the form c 2^e + r with e > j(c, r): %d ; not of this form: %d"
          % (X, (X - 1) // 2, (X - 1) // 2 - len(missing), len(missing)))
    print("first numbers not of this form:", missing[:12])
    for k in range(4, X.bit_length()):
        lo, hi = 1 << k, 1 << (k + 1)
        if hi > X: break
        m = sum(1 for N in range(lo + 1, hi, 2) if not cov[N])
        print("  [2^%d, 2^%d): %d odd numbers, %d not covered" % (k, k + 1, (hi - lo) // 2, m))
    print("1 is covered:", bool(cov[1]), "(a transitive tournament is no L o 0^t with a non-empty out-set)")
if __name__ == "__main__":
    main()
