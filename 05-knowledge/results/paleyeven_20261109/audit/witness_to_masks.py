#!/usr/bin/env python3
"""Audit P1: convert witness lines "W h code" (code = upper triangle row-major, '1' iff i -> j) to the
"name n out_0 ... out_{n-1}" format, deducing n from len(code) = n(n-1)/2."""
import sys
for k, line in enumerate(open(sys.argv[1])):
    parts = line.split()
    if len(parts) < 3 or parts[0] != 'W':
        continue
    h, code = int(parts[1]), parts[2]
    L = len(code)
    n = next(m for m in range(2, 40) if m * (m - 1) // 2 == L)
    out = [0] * n
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code[idx] == '1':
                out[i] |= 1 << j
            elif code[idx] == '0':
                out[j] |= 1 << i
            else:
                raise SystemExit('bad char')
            idx += 1
    print(f'wit{k}_claimed_{h}', n, *out)
