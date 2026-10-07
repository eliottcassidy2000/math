# Independent dense DP for Z_k = #(zeroless k-digit x with 2^k | x), k <= 24, from the least significant digit:
# state r_i = (sum_{l<i} a_l 10^l)/2^i mod 2^(k-i); next digit a must have a == r_i (mod 2); r_{i+1} = (r_i + a 5^i)/2.
import numpy as np, json, sys
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 24
def Z(k):
    cnt = np.zeros(1 << k, dtype=np.int64); cnt[0] = 1
    for i in range(k):
        size = 1 << (k - i); nsize = size >> 1; mask = nsize - 1
        new = np.zeros(nsize, dtype=np.int64)
        r = np.arange(size, dtype=np.int64)
        p5 = pow(5, i, size)
        for a in range(1, 10):
            sel = r[(r & 1) == (a & 1)]
            idx = ((sel + a * p5) >> 1) & mask     # bijection on this parity class -> no index collisions
            new[idx] += cnt[sel]
        cnt = new
    return int(cnt[0])
vals = [Z(k) for k in range(1, KMAX + 1)]
print("Z_1..Z_%d =" % KMAX, vals)
oeis = [int(l.split()[1]) for l in open('b181610.txt') if l.strip() and not l.startswith('#')]
print("matches OEIS A181610 b-file:", vals == oeis[:KMAX], " (b-file has %d terms)" % len(oeis))
repo = json.load(open('/Users/e/Documents/GitHub/math-wt-chessboard-20261006/04-computation/experiments/oai3_20261007_readers/zeroless/zk_values.json'))
rz = {int(k): int(v) for k, v in repo['Z'].items()}
print("matches repo zk_values.json:", all(rz[k] == vals[k-1] for k in range(1, KMAX + 1)))
print("repo Z_25, Z_26 vs OEIS:", rz[25] == oeis[24], rz[26] == oeis[25])
