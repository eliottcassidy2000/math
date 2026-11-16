#!/usr/bin/env python3
"""Wave 42: parity of the number of tournament classes with H = h at h's least order (orders 3..9), and least orders of
seed counts.  Needs gentourng (nauty) and a compiled w38_hq.c (H and Q from out-masks).
Usage: python3 w42_class_parity.py HQ_BINARY"""
import os, subprocess, sys
from collections import Counter, defaultdict
HQ = sys.argv[1]
cnt = defaultdict(Counter)
for n in range(3, 10):
    g = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, check=True).stdout.decode().split()
    lines = []
    for s in g:
        out = [0] * n; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if s[k] == '1': out[i] |= 1 << j
                else: out[j] |= 1 << i
                k += 1
        lines.append(f"{n} " + " ".join(map(str, out)))
    res = subprocess.run([HQ], input=("\n".join(lines) + "\n").encode(), capture_output=True, check=True).stdout.decode().split("\n")
    for l in res:
        if l.strip():
            cnt[int(l.split()[0])][n] += 1
least = {h: min(c) for h, c in cnt.items()}
odd = sum(1 for h in cnt if cnt[h][least[h]] % 2)
primes = [h for h in cnt if h > 1 and all(h % d for d in range(2, int(h ** .5) + 1))]
print('distinct H values (orders 3-9):', len(cnt), '; odd class count at least order:', odd)
print('primes:', len(primes), '; odd at least order:', sum(1 for p in primes if cnt[p][least[p]] % 2))
print('least orders:', {h: least.get(h) for h in (3, 5, 17, 49, 63, 67, 69, 343, 7, 21)})
