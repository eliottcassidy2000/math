#!/usr/bin/env python3
"""prime_path_verify.py (klein-2026-11-01) -- independent check of a prime square-sum path file: every prime below N
appears exactly once, nothing else appears, and every consecutive sum is a perfect square.  Standard library only.
usage: python3 prime_path_verify.py data/prime_path_N617059.txt 617059"""
import math, sys
path = [int(x) for x in open(sys.argv[1]).read().split()]; N = int(sys.argv[2])
s = bytearray([1]) * N; s[0:2] = b'\x00\x00'
for i in range(2, int(N ** 0.5) + 1):
    if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
primes = [i for i in range(N) if s[i]]
ok_set = sorted(path) == primes
sums = [a + b for a, b in zip(path, path[1:])]
ok_sq = all(math.isqrt(t) ** 2 == t for t in sums)
roots = sorted({math.isqrt(t) for t in sums})
print('file %s: %d terms; set = all %d primes below %d: %s; all %d consecutive sums are squares: %s' % (
    sys.argv[1].split('/')[-1], len(path), len(primes), N, ok_set, len(sums), ok_sq))
print('largest prime %d (the %dth prime); first terms %s; last terms %s' % (primes[-1], len(primes), path[:8], path[-5:]))
print('square roots used: %d distinct, from %d to %d; most frequent root %d (%d times)' % (
    len(roots), roots[0], roots[-1], max(roots, key=lambda r: sums.count(r * r)), max(sums.count(r * r) for r in roots)))
print('VERIFIED' if ok_set and ok_sq else 'FAILED'); sys.exit(0 if ok_set and ok_sq else 1)
